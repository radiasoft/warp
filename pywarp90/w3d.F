#include "top.h"
c=============================================================================
c@(#) File W3D.F, version $Revision: 3.196 $, $Date: 2005/03/24 00:32:48 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package W3D of code WARP
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      subroutine w3dinit
      use Subtimers3d
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Called at first reference to package (not nec. a "run" etc.).

      call w3dvers (STDOUT)

!$OMP MASTER
      if (lw3dtimesubs) timew3dinit = timew3dinit + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine w3dvers (iout)
      use Subtimers3d
      use W3Dversion
      integer(ISZ):: iout
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
c  Echoes code version, etc. to output files as they're created
      call printpkgversion(iout,"Particle package W3D",versw3d)
!$OMP MASTER
      if (lw3dtimesubs) timew3dvers = timew3dvers + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine w3dgen()
      use GlobalVars
      use Subtimers3d
      use Ch_var
      use Constant
      use InGen
      use InGen3d
      use InDiag
      use InPart
      use InPart3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles
      use Efields3d
      use GridBoundary3d
      use Multipole
      use Io
      use Lattice
      use LatticeInternal
      use Particles
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
      use InjectVars
      use InjectVars3d
      use DKInterp
      use InDiag3d
      use ExtPart

c  Invoked by the GENERATE command, it sets up the problem
c  This routine allots all of the neccesary dynamic arrays, calls the
c  particle loader and does the initial load onto the charge density
c  mesh, initializes arrays for the field solver and sets the mesh arrays,
c  does the initial field solve, and sets up other arrays that are needed.


      integer(ISZ):: i,j,k,is,ipmin,ip,iwin,nl
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Announce that we're starting up

      call remark(" ***  particle simulation package W3D generating")

c  Put ny to zero for RZ geometry and setup a few flags

      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        ny = 0
        fstype = 10
        lsezax3d = .false.
        lsrhoax3d = .false.
        lsphiax3d = .false.
        lrhodia3d = .false.
      elseif(solvergeom==Zgeom) then
        nx = 0
        ny = 0
        fstype = 10
      elseif(solvergeom==Rgeom) then
        ny = 0
        nz = 0
        fstype = 10
      elseif(solvergeom==AMRgeom) then
        fstype = 11
      endif

c Set wpid if needed and not set by user
      if((l_inj_rz.or.l_inj_regular) .and. wpid==0) then
        wpid = npid+1
        npid = npid+1
      endif

c If delay on setting temperature of injected particles, then
c time of particle creation is recorded initially with a negative value
c switched to positive when temperature has been added.
      if(l_inj_delay_temp) l_inj_rec_inittime = .true.

c Set tpid if needed and not set by user
      if(l_inj_rec_inittime .and. tpid==0) then
        tpid = npid+1
        npid = npid+1
      end if

c Set rpid if needed and not set by user
      if(l_inj_rec_initradius .and. rpid==0) then
        rpid = npid+1
        npid = npid+1
      end if

c Set epflagpid if needed
      if (lepsaveonce) then
        epflagpid = npid + 1
        npid = npid + 1
      endif

c initializes boundary conditions

c --- if using old variables
      if(.not. periinz .or. stickyz .or. .not. stickyxy) then
        call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        call remark("Notice: The variables periinz, stickyz and stickyxy are still ")
        call remark("        valid but will be obsolete soon. They are replaced ")
        call remark("        by the variables pbound0, pboundnz and pboundxy ")
        call remark("        which can be set to absorb, reflect or periodic.")
        call remark("        Note that the periodicity is now set independently")
        call remark("        for fields and for particles.")
        call remark("        Please fix your input files accordingly.")
        call remark("        Your simulation will run ok now, but in future")
        call remark("        versions, it will stop with an error.")
        call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
       if(periinz) then
          bound0  = periodic
          boundnz = periodic
          if(.not. stickyz) then
            pbound0  = periodic
            pboundnz = periodic
          end if
        end if
        if(stickyz) then
          pbound0  = absorb
          pboundnz = absorb
        end if
        if(stickyxy) then
          pboundxy  = absorb
        end if
      end if

c --- using new variables
      if((bound0==periodic .and. boundnz/=periodic) .or.
     &   (bound0/=periodic .and. boundnz==periodic)) then
         write(0,*) 'Error in setup of bound0 and boundnz: either both or none should be periodic, stopping.'
         stop
      end if
      if((pbound0==periodic .and. pboundnz/=periodic) .or.
     &   (pbound0/=periodic .and. pboundnz==periodic)) then
         write(0,*) 'Error in setup of pbound0 and pboundnz: either both or none should be periodic, stopping.'
         stop
      end if
      if(bound0==periodic .or. boundnz==periodic) then
        bound0  = periodic
        boundnz = periodic
      end if
      if(pbound0==periodic .or. pboundnz==periodic) then
        pbound0  = periodic
        pboundnz = periodic
      end if

c  Set global values of zmmin and zmmax before the decomposition is done
c  for the parallel code.
      zmminglobal = zmmin
      zmmaxglobal = zmmax

c------- Parallel only -----------------------------------------------------
#ifdef MPIPARALLEL
      call init_w3d_parallel()
#endif
c---------------------------------------------------------------------------

c  Estimate wall radius, needed for g-factor calc
c  (rwallfac = 1 is probably NOT a good guess)

      rwall = rwallfac * sqrt( xmmax**2 + xmmin**2 )

c  Calculate derived quantities and species related arrays  (species are 
c  set from an internal call to "species" in derivqty).   

      call derivqty

c  Load text lines which are printed at the bottom of all plots.

      call stepid (it, time, zbeam)

c  Print values of input variables, other interesting things to text file
      if (warpout > -1) then
        call edit (warpout, "runid")
        call edit (warpout, "it")
        call edit (warpout, "time")
        call edit (warpout, "InGen")
        call edit (warpout, "InDiag")
        call edit (warpout, "InPart3d")
        call edit (warpout, "InMesh3d")
      endif

c  Create the dynamic arrays for fields, contour plot workspace
c  Also set vbeamfrm, adjust x and ymmin for symmetries and set
c  so that by default izfsmax = nz (for full field solve).

      if (nzfull == 0) nzfull = nz
      nmxy  = max(nx,ny)
      nmxyz = max(nx,ny,nzfull)
      call gchange("Fields3d", 0)
      vbeamfrm = vbeam
      if(solvergeom==XYZgeom .or. solvergeom == XYZgeomMR) then
        if (l2symtry) then
          ymmin = 0.
        elseif (l4symtry) then
          xmmin = 0.
          ymmin = 0.
        endif
      else if(solvergeom==XZgeom) then
        ymmin = 0.
        if (l2symtry .or. l4symtry) xmmin = 0.
      else if(solvergeom==RZgeom .or. solvergeom==Zgeom .or. solvergeom==Rgeom) then
        l2symtry=.false.
        l4symtry=.false.
        xmmin = 0.
        ymmin = 0.
      endif
      if (izfsmax == 0) izfsmax = nz

c  Calculate mesh dimensioning quantities

      if (nx .ne. 0) dx = (xmmax - xmmin) / nx
      if (ny .ne. 0) dy = (ymmax - ymmin) / ny
      if (nzfull .ne. 0) dz = (zmmaxglobal - zmminglobal) / nzfull
      if(solvergeom==RZgeom) then
        dy = dx
      elseif(solvergeom==XZgeom) then
        dy = 1.
      elseif(solvergeom==Zgeom) then
        dx = 1.
        dy = 1.
      elseif(solvergeom==Rgeom) then
        dy = 1.
        dz = 1.
      end if
      do i = 0, nx
         xmesh(i) = i*dx + xmmin
      enddo
      do j = 0, ny
         ymesh(j) = j*dy + ymmin
      enddo
      do k = 0, nz
         zmesh(k) = k*dz + zmmin
      enddo

c  Calculate location of axis in mesh, the term dx*1.e-5 acts as fuzp
      ix_axis = nint(-xmmin/dx)
      if(solvergeom==XYZgeom .or. solvergeom == XYZgeomMR .or. solvergeom==AMRgeom) then
        iy_axis = nint(-ymmin/dy)
      else ! RZgeom, XZgeom, Zgeom, Rgeom
        iy_axis = 0
      endif
      iz_axis = nint(-zmmin/dz)

c  Setup arrays for potential and charge density for the particles.

      nxp = nx
      nyp = ny
#ifndef MPIPARALLEL
      nzp = nz
      zmminp = zmmin
      zmmaxp = zmmax
#endif

c  Initialize base grid for RZ solver
      if(solvergeom==RZgeom .or.
     &   solvergeom==XZgeom .or.
     &   solvergeom==Rgeom  .or.
     &   solvergeom==Zgeom) then
        call init_base(nx,nz,dx,dz,xmmin,zmmin)
      end if

c  Initial call to fieldsolver in order to initialize attx, kxsq, etc.
      call vp3d (1)

c  Create the dynamic arrays in Z_arrays; set the z mesh

      if (nzzarr == 0) nzzarr = nzfull
      call gchange("Z_arrays", 0)
      if (zzmin == 0.) zzmin = zmminglobal
      if (zzmax == 0.) zzmax = zmmaxglobal
      dzz = (zzmax - zzmin)/nzzarr
      dzzi = 1./dzz
      do k = 0, nzzarr
         zplmesh(k) = zzmin + k*dzz
      enddo

c  Set the value of prwall, radius at which particles are lost
      if (prwall == 0.) prwall = xmmax - dx
      do k=0,nzzarr
        prwallz(k) = prwall
        if (prwallxz(k) == 0.) prwallxz(k) = prwallx
        if (prwallyz(k) == 0.) prwallyz(k) = prwally
        if (prwelipz(k) == 1.) prwelipz(k) = prwelip
      enddo

c  Re-size the dynamic arrays for the lattice (scan for true length, first).
c  Also allocate the dynamic internal lattice arrays.

      call remark(" ---  Resetting lattice array sizes")
      call resetlat
      if (nzlmax == 0) nzlmax = nzfull
      if (nzl    == 0 .or. nzl > nzlmax) nzl = nzlmax
      call gchange("LatticeInternal", 0)
      if (zlmin == 0.) zlmin = zmminglobal
      if (zlmax == 0.) zlmax = zmmaxglobal
      dzl = (zlmax - zlmin)/nzl
      dzli = 1./dzl
      do k = 0, nzl
         zlmesh(k) = zlmin + k*dzl
      enddo

c Allocate space for the selfe array, which contains the finite differenced
c phi (Ex, Ey, Ez) on the mesh.
      if (efetch == 3 .and.
     &    (solvergeom == XYZgeom .or.
     &     solvergeom == RZgeom .or.
     &     solvergeom == XZgeom .or.
     &     solvergeom == Rgeom  .or.
     &     solvergeom == Zgeom)) then
        nx_selfe = nxp
        ny_selfe = nyp
        nz_selfe = nzp
        if(solvergeom==RZgeom .or.
     &     solvergeom==XZgeom .or.
     &     solvergeom==Rgeom  .or.
     &     solvergeom==Zgeom) ny_selfe = 0
        if(solvergeom==Zgeom) nx_selfe = 0
        call gchange("Efields3d",0)
      endif

c  Create dynamic arrays for 3D multipole moments of the electrostatic
c  potential
      call gchange("Multipole", 0)

c  Create the dynamic arrays for particles (set npmax to an estimated
c  length for now, for those loading schemes that don't actually
c  use a user-set npmax directly)

      call remark(" ---  Allocating space for particles")
      if (xrandom == "grid") npmax = nxstripe*nystripe*nzstripe
      if (xrandom == "fibonacc") npmax = nfibgrps*fibg1

      call alotpart
      call alotlostpart

c  Load the particles, calculate the charge density

      call remark(" ---  Loading particles")
      call stptcl3d
      call setgamma(lrelativ)
      call remark(" ---  Setting charge density")
      call loadrho3d(-1,-1,-1,.true.)
      call remark(" ---  done")

c  Set up for injection
      if (inject > 0) call injctint

c  Create the dynamic arrays for the partcle qtys needed for the
c  damped mover.  They are always allocated, but with length 1 if not used.

      npdamp = 1
      if (eomdamp /= 0.) npdamp = npmax
      call gchange("Damped_eom", 0)

c  Create the dynamic arrays for "window" moments

      call remark(" ---  Allocating Win_Moments")
      zwindows(1,0) = zmminglobal
      zwindows(2,0) = zmmaxglobal
      nzwind = 0
      do iwin = 1, NWINDOWS
         if (zwindows(1,iwin) /= zwindows(2,iwin)) nzwind = nzwind + 1
      enddo
      call gchange("Win_Moments", 0)

c  Create the dynamic arrays for z moments

      call remark(" ---  Allocating Z_Moments")
      if (nzmmnt == 0) nzmmnt = nzfull
      call gchange("Z_Moments", 0)
      if (zmmntmax == 0.) zmmntmax = zmmaxglobal
      if (zmmntmin == 0.) zmmntmin = zmminglobal
      dzm = (zmmntmax - zmmntmin)/nzmmnt
      dzmi = 1./dzm
      do k = 0, nzmmnt
         zmntmesh(k) = zmmntmin + k*dzm
      enddo

c  Create the dynamic arrays for lab frame moments

      call remark(" ---  Allocating Lab_Moments")
      if (vbeam*dt /= 0.) then
        nl = (zmmax-zmmin)/(vbeam*dt) +.5
      else
        nl = 1
      endif
      call initlabwn(nl)

c  Create the scratch arrays for phase space plots (permanent, for now)
c  and set limits for plots

c     call remark(" ---  Allocating scratch space for plots")
c     if (npsplt == 0) npsplt = nparpgrp-1
c     call gchange("Pspwork", 0)
      if (xplmin == 0.) xplmin = -xmmax
      if (xplmax == 0.) xplmax = xmmax
      if (yplmin == 0.) yplmin = -ymmax
      if (yplmax == 0.) yplmax = ymmax
      if (zplmin == 0.) zplmin = zmmin
      if (zplmax == 0.) zplmax = zmmax

c  Choose which particles will be plotted in "subset" plots

      call psubsets

c  Setup history mechanism

      if (nhist > 0) then
c       --- create the dynamic arrays for history data; set pointer into them
        call remark(" ---  Allocating history arrays")
        if (lenhist == 0) lenhist = min ( nt/nhist + 1, 100)
        call gchange("Hist", 0)
        jhist = -1
      elseif (nhist < 0) then
c       --- call interpreter routine to setup hst package
c       --- setup_hst is in bas.wrp
        call execuser("setup_hst")
      endif

c  Print interesting things to plot file and teletype

      call prntpara(dx,dy,dz,lprntpara)
      call prntpa3d(lprntpara)

c  Initial fieldsolve, diagnostics

      call step3d ("w3dgen")

c  Set up Eears of z

      call seteears()

!$OMP MASTER
      if (lw3dtimesubs) timew3dgen = timew3dgen + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine w3dexe()
      use Subtimers3d
      use Picglb
      use InGen
      use InPart
      use Picglb3d
      use Ctl_to_pic
      use Particles

c  Takes a time step.
c  This routine advances the mesh in the lab frame, sets the logicals
c  which control how this next step is to be done, and then calls
c  the routine STEP3D to do the step.


      real(kind=8):: zcorrection
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c --- Announce that we're running

      if (it == 0) call remark(" ***  particle simulation package W3D running")

c --- Accelerate grid frame.
      call acclbfrm(zcorrection)

c --- Set timestep counter, time, and advance grid frame. The grid
c --- frame is advanced here to be with the particles after the
c --- position advance.
c --- The zcorrection from the accleration of the beam frame is added on
c --- by adding it to zgrid.
c --- zgridprv is set here (as well as in padvnc3d) in case the
c --- user has changed zbeam.

      it = it + 1
c     time = time + dt ! Moved to end of padvnc3d
      if (lbeamcom) then
c       --- In this case, it is still up to the user to set zgrid if zbeam is
c       --- set. Here, zbeam is not necessarily equal to zgrid.
        zgridprv = zgrid
      else
c       --- Set zgridprv to zbeam so that the user only has to set zbeam.
c       --- Otherwise, zbeam is the same as zgrid, as set at the
c       --- end of padvnc3d.
        zgridprv = zbeam
      endif
      if (lgridqnt) then
c       --- zgrid is integer number of dz's
        zgrid = int((zbeam + dt*vbeamfrm + zcorrection)/dz + .5)*dz
      else
c       --- zgrid is the same as the beam frame
        zgrid = zbeam + dt*vbeamfrm + zcorrection
      endif
      call stepid (it, time+dt, zgrid)

c  set logicals

      lfirst = .false.
      if (ncall == 1) lfirst = .true.
      llast = .false.
      if (ncall == maxcalls) llast = .true.

c  call the routine that does the actual work

      call step3d ("w3dexe")

c  Have we reached the end of the run or run out of particles?

      if ( lfinishd ) then
         call remark("w3dexe: problem completed.")
!$OMP MASTER
         if (lw3dtimesubs) timew3dexe = timew3dexe + wtime() - substarttime
!$OMP END MASTER
         return
      elseif (nplive <= 0) then
c Check if there are any live particles left.
c If on first time step, set nplive to one and continue running since if using
c injection, there may not be any particles yet but still want to run.
         if (it == 1) then
           nplive = 1
         else
!           print*," *** W3DEXE: stopping, nplive =",nplive
!$OMP MASTER
           if (lw3dtimesubs) timew3dexe = timew3dexe + wtime() - substarttime
!$OMP END MASTER
           return
         endif
      endif

c  Continue to run.
!$OMP MASTER
      if (lw3dtimesubs) timew3dexe = timew3dexe + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine w3dfin()
      use Subtimers3d
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
      use Picglb
      use Picglb3d
      use Win_Moments
      use Z_Moments
      use Z_arrays
      use Hist
c     use Pspwork
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Finish up at end of RUN, or on receipt of FIN
c  This routine is never called, at present; history plots are
c  made using a Python interpreter script (histplot), and we just
c  end the run.  If we wanted to chain runs so that an output qty
c  might be plotted vs a parameter, this routine might be useful.
c  It would be needed for a non-Python version of WARP.
c  For now it serves as a place-holder.


c  print final edits

c  perform diagnostics (unless we just did)

c  Create history plots

c  Make a restart dump (unless we just did, or the user inhibits it)

c  create final printouts

c  release storage

      call gfree ("Fields3d")
      call gfree ("Fields3dParticles")
      call gfree ("Hist")
      call gfree ("Win_Moments")
      call gfree ("Z_Moments")
      call gfree ("Lab_Moments")
      call gfree ("Moments")
      call gfree ("Lattice")
      call gfree ("LatticeInternal")
      call gfree ("Particles")
c     call gfree ("Pspwork")
      call gfree ("Z_arrays")

!$OMP MASTER
      if (lw3dtimesubs) timew3dfin = timew3dfin + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine step3d (caller)
      use Subtimers3d
      use Constant
      use InGen
      use InGen3d
      use InDiag
      use InDiag3d
      use InPart
      use InMesh3d
      use Fields3d
      use Fields3dParticles
      use Efields3d
      use Io
      use Particles
      use Picglb
      use Picglb3d
      use InjectVars
      use LatticeInternal
      use GridBoundary3d
      use Timers
      character(*):: caller
      logical(ISZ):: thisstep,thiszbeam,dolabwn

c  When called by W3DEXE, step3d advances the system forward in time one 
c  timestep and gathers diagnostics.  When called by W3DGEN, step3d takes 
c  a step of zero size, to compute fields, and gather diagnostics at start 
c  of run.


      real(kind=8):: zbeaml,zbeamr
      integer(ISZ):: k
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c --- Set the internal lattice variables. This is not generally necessary at
c --- this point (it is redundant most of the time, the next call to
c --- setlatt in this subroutine is sufficient). There are cases where
c --- this is required for consistency. Since it is cheap (time wise),
c --- it is better to make sure the data is consistent than to save a
c --- little bit of time. The value of nzl etc must be checked since other
c --- packages (like WXY or ENV) may have reset it. For example, if the
c --- ENV package is generated after the W3D package, nzl will be set to
c --- zero. Switching back to W3D and running step, the internal lattice
c --- would still be setup for the ENV package and so the step would produce
c --- erroneaous results.
      if (nzl == 0) then
        nzl = nzlmax
        zlmin = zmmin
        zlmax = zmmax
        dzl = (zlmax - zlmin)/nzl
        dzli = 1./dzl
        do k = 0, nzl
          zlmesh(k) = zlmin + k*dzl
        enddo
      endif
      call setlatt()
c  Main particle advance: x to t.l. it; v to t.l. it-1/2
c  Half-step in v from t.l. it-1   if last step was "special"
c  Full-step in v from t.l. it-3/2 if last step not "special"
c  No step at all if generating.

      if (caller == "w3dexe") then
        if (lspecial) then
            call padvnc3d ("halfv")
         else
            call padvnc3d ("fullv")
         endif
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

      lfinishd = (it >= nt) .or. (time >= tstop*(1.-MACHEPS)) .or.
     &                             (zbeam >= zstop)
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
      if (nhist /= 0) then
        lhist  = mod(it,nhist) == 0
      else
        lhist  = .false.
      endif
      ldump    = mod(it, itdump) == 0
      llabwn   = dolabwn()
      lspecial = (lfinishd .or. lalways .or. lseldom .or. ldump .or.
     &            lmoments .or. lhist .or. llabwn .or. llast .or.
     &            (it == 0) .or. allspecl)

c  Set the "gap" electric field.

      call setegap

c  Gather moments used in diagnostics at "special" timesteps only. 
c  Compute line charge density (gtlchg3d) and the axial line charge 
c  (srhoax3d) on 1-d meshes.  Note -- these moments accumulations are
c  done at this phase of the particle advance to allow for the eventual 
c  use of a single array for rho and phi.   

      if (lspecial) then
        call gtlchg3d
        call srhoax3d
      endif

c  Charge density contour plot diagnostics.  Note -- these diagnostics 
c  are done at this phase of the particle advance to allow for the eventual 
c  use of a single array for rho and phi.   

      if (lalways .or. lseldom) call pltfld3d("rho",ALWAYS)
      if (lseldom)            call pltfld3d("rho",SELDOM)

c  Set lattice; this is done just before field solve, and so is
c  relative to ZBEAM in the same way that self-fields are.

      call setlatt
      if (linbend) call setrstar(rstar(-1),nz,dz,zmmin,zgrid)

c  Field-solve for potential.

      if (lbeforefs) call execuser("beforefs")
      call fieldsol3d(-1)
      if (lafterfs) call execuser("afterfs")

c  Pre-calculate the self-E if it is needed for sete3d. This is done
c  after the call to afterfs in case some manipulation is done to phi.
      if (efetch == 3 .and.
     &    (solvergeom == XYZgeom .or.
     &     solvergeom == RZgeom .or.
     &     solvergeom == XZgeom .or.
     &     solvergeom == Rgeom  .or.
     &     solvergeom == Zgeom)) then
        call getselfe3d(phip(0,0,-1),nxp,nyp,nzp,selfe(1,0,0,0),
     &                  nx_selfe,ny_selfe,nz_selfe,dx,dy,dz,
     &                  boundxy,boundxy,boundxy,boundxy)
      endif

c  Set the potential near the emitting surface.
      if (inject > 0 .and. inject < 4) then
        call getinj_phi()
      endif

c  Complete constant current and axially directed space-charge limited
c  injection with new fields including injected particles.
      if (inject > 0 .and. caller == "w3dexe")
     &  call inject3d(2)

c  If a flag was set making this a "special" step,
c  do a half-advance to bring v to t.l. it 

      if (caller == "w3dexe" .and. lspecial) then
         call padvnc3d ("synchv")
      elseif (caller == "w3dgen") then
         call padvnc3d ("gen")
      endif

c  Finish transverse space-charge limited injection: change ins and nps.
c  This is done here so the newly injected particles don't go through
c  the 'synchv' advance.
      if (inject > 0 .and. caller == "w3dexe") then
         call inject3d(3)
      end if

c  Gather moments used in diagnostics at "special" timesteps only. 
c  Compute mean beam z velocity from current and line charge density 
c  on a 1-d mesh.  Also, calculate the electrostatic energy (getese3d), 
c  electrostatic potential on axis (sphiax3d), and the axial electric 
c  field (sezax3d).  

      if (lspecial) then
        if (lgetvzofz) call getvzofz
        call getese3d 
        call sphiax3d 
        call  sezax3d
      endif 

c  Electrostatic potential contour plot diagnostics 

      if (lalways .or. lseldom) call pltfld3d("phi",ALWAYS)
      if (lseldom)            call pltfld3d("phi",SELDOM)

c  1d array plot diagnostics.

      if (lalways .or. lseldom) call onedplts(ALWAYS)
      if (lseldom)            call onedplts(SELDOM)

c  Phase space diagnostics

      if (lalways .or. lseldom) call psplots (ALWAYS)
      if (lseldom)            call psplots (SELDOM)

c  Finally, moment diagnostic printout and history storage

      if (caller == "w3dgen" .or. lspecial)
     &  call minidiag (it,time,lspecial)
!$OMP MASTER
      if (lw3dtimesubs) timestep3d = timestep3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine exteb3d(np,xp,yp,zp,uzp,gaminv,dtl,dtr,
     &                   bx,by,bz,ex,ey,ez,m,q,bendres,bendradi,gammabar,dt)
      use Subtimers3d
      use Timers
      integer(ISZ):: np
      real(kind=8):: dtl,dtr,m,q,gammabar,dt
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uzp(np), gaminv(np)
      real(kind=8):: bx(np), by(np), bz(np), bendres(np), bendradi(np)
      real(kind=8):: ex(np), ey(np), ez(np)

c  Calculates "external" E, B fields
c  Calculates electric or magnetic AG focusing fields, bending and dipole
c  fields, and accelerating fields.
c  Includes back-rotation associated with coordinate transformation into By. 

c  NOTE: When we (someday) set B_self from a Lorentz transformation
c  on E_self, we'll have to carefully work out a sequence of calls,
c  since this routine is called more than once in the PADVNC3D loop on
c  a single step at present.  Perhaps we will also have to compute
c  B_self more than once.

      integer(ISZ):: ip
      real(kind=8):: qoverm
      real(kind=8):: timetemp,wtime
      real(kind=8):: substarttime
      if (lw3dtimesubs) substarttime = wtime()

      timetemp = wtime()

      qoverm = q/m

c     --- handle uniform fields
      call applyuniformfields(np,ez,bz)

c     --- handle quads
      call applyquad(np,xp,yp,np,zp,uzp,gaminv,dtl,dtr,dt,.false.,ex,ey,bx,by)

c     --- handle dipos 
      call applydipo(np,np,zp,uzp,gaminv,dtl,dtr,dt,.false.,ex,ey,bx,by)

c     --- handle sexts
      call applysext(np,xp,yp,np,zp,uzp,gaminv,dtl,dtr,dt,.false.,ex,ey,bx,by)

c     --- handle hard-edge electric and magnetic multipoles
      call applyhele(np,xp,yp,np,zp,uzp,gaminv,dtl,dtr,dt,.false.,
     &               ex,ey,ez,bx,by,bz)

c     --- fold in coordinate transformation associated with bends
      call applybend(np,xp,uzp,np,bendres,bendradi,m,q,.false.,by)

c     --- handle acceleration (and calculate position correction)
      call applyaccl(np,xp,zp,uzp,gaminv,dtl,dtr,dt,qoverm,.false.,ez)

c     --- handle electrostatic multipole components
      call applyemlt(np,xp,yp,np,zp,dtl,dtr,dt,.false.,ex,ey,ez)

c     --- handle magnetostatic multipole components
      call applymmlt(np,xp,yp,np,zp,dtl,dtr,dt,.false.,bx,by,bz)

c     --- handle magnetic fields from 3-D grid
      call applybgrd(np,xp,yp,np,zp,.false.,bx,by,bz)

c     --- handle electrostatic potential from 3-D grid
      call applypgrd(np,xp,yp,np,zp,.false.,ex,ey,ez)

c     --- Set E to zero for lost particles
      do ip=1,np
        if (uzp(ip) == 0) then
          ex(ip) = 0.
          ey(ip) = 0.
          ez(ip) = 0.
        endif
      enddo

c     --- Accumulate time for applying the fields of the lattice
!$OMP MASTER
      latticetime = latticetime + (wtime() - timetemp)
!$OMP END MASTER

!$OMP MASTER
      if (lw3dtimesubs) timeexteb3d = timeexteb3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine othere3d (np,xp,yp,zp,zbeam,zimax,zimin,straight,
     &                     ifeears,eears,eearsofz,dzzi,nzzarr,zzmin,
     &                     dedr,dexdx,deydy,dbdr,ex,ey,ez,bx,by,bz)
      use Subtimers3d
      integer(ISZ):: np,ifeears,nzzarr
      real(kind=8):: zbeam,zimax,zimin,straight,eears,dzzi,zzmin
      real(kind=8):: dedr,dexdx,deydy,dbdr
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: ex(np), ey(np), ez(np)
      real(kind=8):: bx(np), by(np), bz(np)
      real(kind=8):: eearsofz(0:nzzarr)

      integer(ISZ):: ip,iz
      real(kind=8):: zs,wz
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Set the electric fields from external sources,
c  inculding axial confining fields, and uniform focusing fields.

c  Apply axial confining ears for finite beam.

      if (ifeears == 1) then
        zs = (zimax - zimin)*straight*0.5
        do ip=1,np
          if ((zp(ip)-zbeam) > zs) then
            ez(ip) = ez(ip) + eears*(zp(ip) - zbeam - zs)
          elseif ((zp(ip)-zbeam) < -zs) then
            ez(ip) = ez(ip) + eears*(zp(ip) - zbeam + zs)
          endif
        enddo
      endif

c  axial confining ears for finite beam as a function of z

      if (ifeears == 2) then
        do ip=1,np
          iz = (zp(ip) - zbeam - zzmin)*dzzi
          if (0 <= iz .and. iz <= nzzarr) then
            wz = (zp(ip) - zbeam - zzmin)*dzzi - iz
            ez(ip) = ez(ip) + eearsofz(iz)*(1.-wz) + eearsofz(iz+1)*wz
          endif
        enddo
      endif

c  uniform focusing forces

c     --- radial electric field
      if (dedr /= 0.) then
        do ip=1,np
          ex(ip) = ex(ip) + dedr*xp(ip)
          ey(ip) = ey(ip) + dedr*yp(ip)
        enddo
      endif
c     --- x- and y-electric fields
      if ((dexdx /= 0.) .or. (deydy /= 0.)) then
        do ip=1,np
          ex(ip) = ex(ip) + dexdx*xp(ip)
          ey(ip) = ey(ip) + deydy*yp(ip)
        enddo
      endif
c     --- azimuthal magnetic field
      if (dbdr /= 0.) then
        do ip=1,np
          bx(ip) = bx(ip) - dbdr*yp(ip)
          by(ip) = by(ip) + dbdr*xp(ip)
        enddo
      endif

!$OMP MASTER
      if (lw3dtimesubs) timeothere3d = timeothere3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine getese3d
      use Subtimers3d
      use Beam_acc
      use Picglb
      use Picglb3d
      use InMesh3d
      use Fields3d
      use Moments
      use InPart
      use InGen
      use InGen3d
      use Z_arrays
      use InDiag3d

c Calculate the electrostatic energy, sum of rho*phi*dx*dy*dz/2
c Includes external fields: Uniform focusing and eears for cigars
c This does it in a partially vectorized manner.

      integer(ISZ):: ix,iy,iz,izb
      real(kind=8):: zm,zz,phiextun,wzb
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (.not. lgetese3d) then
!$OMP MASTER
        if (lw3dtimesubs) timegetese3d = timegetese3d + wtime() - substarttime
!$OMP END MASTER
        return
      endif

      ese = 0.
      zm = (zimax - zimin)*(1. - straight)*0.5

c     conversion factor to go from grid frame to beam frame
      zz = zmmin + zgrid - zbeam

      do iy=0,ny
        do ix=0,nx
          scrtch(ix,iy) = 0.0
        enddo
      enddo

c     --- Loop order switched for the OpenMP parallel version - this
c     --- avoids requiring the accumulation into scrtch to be put
c     --- into a critical section.
!$OMP PARALLEL
!$OMP DO PRIVATE(phiextun,izb,wzb,ix,iy,iz)
#ifdef _OPENMP
      do iy=0,ny
#endif
      do iz=izfsmin,izfsmax-1
        phiextun = 0.
        if (ifeears == 1) then
          if (iz*dz+zz < -zm) phiextun=eears*.5*(iz*dz+zz+zm)**2
          if (iz*dz+zz >  zm) phiextun=eears*.5*(iz*dz+zz-zm)**2
        elseif (ifeears == 2) then
c         --- This is not quite right for the parallel version since
c         --- the longitudinal extent of eearsofz (the same as the extent
c         --- of the particles) will not necessarily overlap the extent of
c         --- the field solver.
          izb = (iz*dz+zz-zzmin)*dzzi
          wzb = (iz*dz+zz-zzmin)*dzzi - izb
          if (0 <= iz .and. iz < nzzarr)
     &      phiextun = eearsofz(izb)*(1. - wzb) + eearsofz(izb+1)*wzb
        endif
#ifndef _OPENMP
        do iy=0,ny
#endif
          do ix=0,nx
            scrtch(ix,iy) = scrtch(ix,iy) + rho(ix,iy,iz)*(phi(ix,iy,iz)
     &                       - dedr*.5*((ix*dx+xmmin)**2 + (iy*dy+ymmin)**2)
     &                       - dexdx*.5*(ix*dx+xmmin)**2 
     &                       - deydy*.5*(iy*dy+ymmin)**2
     &                       + phiextun)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

c     --- Apply appropriate multiplication factors for symmetries
      if (l4symtry) then
c       --- Center point is unchanged
c       --- Edges are doubled
c       --- Bulk is quadrupoled
        scrtch(0,1:ny) = 2.*scrtch(0,1:ny)
        scrtch(1:nx,0) = 2.*scrtch(1:nx,0)
        scrtch(1:nx,1:ny) = 4.*scrtch(1:nx,1:ny)
      else if (l2symtry) then
c       --- Edge is unchanged
c       --- Bulk is doubled
        scrtch(1:nx,0:ny) = 2.*scrtch(1:nx,0:ny)
      endif

      do iy=1,ny
        do ix=0,nx
          scrtch(ix,0) = scrtch(ix,0) + scrtch(ix,iy)
        enddo
      enddo

      do ix=0,nx
        ese = ese + scrtch(ix,0)
      enddo

      ese = dx*dy*dz*0.5*ese

!$OMP MASTER
      if (lw3dtimesubs) timegetese3d = timegetese3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine gtlchg3d
      use Subtimers3d
      use Picglb
      use InGen3d
      use Picglb3d
      use InMesh3d
      use Fields3d
      use Z_arrays
      use InDiag3d

c  Calculates the line charge density from rho.

      real(kind=8):: zz,zg,wzg,dzi
      integer(ISZ):: ix,iy,iz,izg
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      if (.not. lgtlchg3d) return

      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        call gtlchgrz
        return
      end if

      dzi = 1./dz

c     conversion factor to go from grid frame to beam frame
      zz = zgrid + zmmin - zzmin - zbeam

c First, sum up the charge density in each axial plane.

!$OMP PARALLEL DO
      do iz=0,nz

        if (l2symtry) then
c         --- 2-fold: account for only one side of rho being stored
          scrtch(iz,0) = sum(rho(:,0,iz)) + 2.*sum(rho(:,1:,iz))

        elseif (l4symtry) then
c         --- 4-fold: account for only one quadrant of rho being stored
          scrtch(iz,0) = rho(0,0,iz) + 4.*sum(rho(1:,1:,iz)) +
     &                   2.*sum(rho(1:,0,iz)) + 2.*sum(rho(0,1:,iz))

        else
c         --- no symmetries...
          scrtch(iz,0) = sum(rho(:,:,iz))

        endif

      enddo
!$OMP END PARALLEL DO

c Then, linearly interpolate the data into the beam frame array, linechg.

      do iz=0,nzzarr

        zg = (iz*dzz - zz)*dzi
        izg = int(zg + nz) - nz
        wzg = zg - izg

        if (0 <= izg .and. izg <= nz) then
          linechg(iz) = scrtch(izg,0)*dx*dy*(1. - wzg)
        else
          linechg(iz) = 0.
        endif
        if (0 <= izg+1 .and. izg+1 <= nz) then
          linechg(iz) = linechg(iz) + scrtch(izg+1,0)*dx*dy*wzg
        endif

      enddo

#ifdef MPIPARALLEL
      call parallelnonzerorealarray(linechg,nzzarr+1)
#endif

!$OMP MASTER
      if (lw3dtimesubs) timegtlchg3d = timegtlchg3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine seteears()
      use Subtimers3d
      use Beam_acc
      use InGen
      use InPart
      use InMesh3d
      use Picglb3d
      use Constant
      use Particles
      use Z_arrays
      use Fields3d

c Calculate ear fields to confine the beam axially.
c Controlled by ifears.  If zero, no ear fields.
c If one, strictly linear ear fields, calculated from expression from Nueffer.
c If two, ear fields obtained from axial E fields on axis with linear pressure
c         term added on.
c Calculated for first species.


      real(kind=8):: zm
      integer(ISZ):: nzla,nzlb,nzma,nzmb,iz
      real(kind=8):: zs,zl,eearsprs
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      zm = (zimax - zimin)*(1. - straight)*0.5
      emitlong = zm*2.*vthz/dvnz(vbeam)

c     set to be linear in z (eears is actually slope of Eears)
      if (ifeears == 1) then
        eears = eears*(- 2*ibeam/(vbeam*zm**2)*gfactor/(2.*Pi*eps0)
     &                  - sm(1)*vbeam**2*emitlong**2/(sq(1)*zm**4))

      elseif (ifeears == 2) then

c     set to initial field on axis (with the center part zero)
c     the '5' is just a guestimate, it should be nz dependent
c     The ends (beyond the beam) are set with the E at the end of the beam
        zl = (zimax - zimin)*0.5
        nzla = nzzarr/2-int(zl*dzzi)
        nzlb = nzzarr/2+int(zl*dzzi)
        zs = (zimax - zimin)*straight*0.5
        nzma = nzzarr/2-int(zs*dzzi)
        nzmb = nzzarr/2+int(zs*dzzi)
        do iz=0,nzla
          eearsofz(iz) = - ezax(nzla)
        enddo
        do iz=nzla+1,nzma+5
          eearsofz(iz) = - ezax(iz)
        enddo
        do iz=nzmb-5,nzlb
          eearsofz(iz) = - ezax(iz)
        enddo
        do iz=nzlb+1,nzzarr
          eearsofz(iz) = - ezax(nzlb)
        enddo
c       --- Add on linear pressure term
        if (emitlong /= 0.) then
          eearsprs = - sm(1)*vbeam**2*emitlong**2/(sq(1)*zm**4)
          do iz=0,nzma
            eearsofz(iz) = eearsofz(iz) + eearsprs*(zplmesh(iz) + zs)
          enddo
          do iz=nzmb,nzzarr
            eearsofz(iz) = eearsofz(iz) + eearsprs*(zplmesh(iz) - zs)
          enddo
        endif
      endif

!$OMP MASTER
      if (lw3dtimesubs) timeseteears = timeseteears + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine padvnc3d(center_in)
      use GlobalVars
      use Subtimers3d
      use InMesh3d
      use InGen
      use InGen3d
      use InPart
      use InGaps
      use InDiag
      use Lattice
      use LatticeInternal
      use Particles
      use Fields3d
      use Fields3dParticles
      use Efields3d
      use Picglb
      use Picglb3d
      use Beam_acc
      use Z_arrays
      use Z_Moments, only: nzmmnt,nszmmnt,tempmaxp,tempminp,tempzmmnts0,tempzmmnts
      use Damped_eom
      use InjectVars
      use InjectVars3d
      use Apertures
      use GridBoundary3d
      use DKInterp
      use DKInterptmp
      character(*):: center_in
      character(len=len(center_in)):: center

c  Advances the particles position and velocity according to CENTER,
c  and also loads RHO at the new time level.

      integer(ISZ):: is,ip,ipmin,i
      real(kind=8):: uxpadv,uypadv,uzpadv,halfdt
      real(kind=8):: fulldt_s,halfdt_s
      real(kind=8):: zpo(nparpgrp)
c RHC: moved ex, ey, ez, bx, by, bz to DKInterptmp in w3d.v  3/11/04
c      real(kind=8):: ex(nparpgrp), ey(nparpgrp), ez(nparpgrp)
c      real(kind=8):: bx(nparpgrp), by(nparpgrp), bz(nparpgrp)
      real(kind=8):: uxpo(nparpgrp), uypo(nparpgrp), uzpo(nparpgrp)
      real(kind=8):: bendres(nparpgrp), bendradi(nparpgrp)
#ifdef _OPENMP
      real(kind=8),allocatable:: threadmaxp(:,:,:),threadminp(:,:,:)
      real(kind=8),allocatable:: threadzmmnts0(:,:,:),threadzmmnts(:,:,:,:)
      integer(ISZ):: ithread,omp_get_thread_num
      integer(ISZ):: nthread,omp_get_num_threads
      integer(ISZ):: allocerror
#endif
      real(kind=8),pointer:: maxp(:,:),minp(:,:)
      real(kind=8),pointer:: zmmnts0(:,:),zmmnts(:,:,:)
      real(kind=8):: getbeamcom,zbeamcom
      real(kind=8):: substarttime,wtime

      if (lw3dtimesubs) substarttime = wtime()

      halfdt = 0.5*dt
      center = center_in

c  Zero the bend radius and residence arrays
c     call zeroarry (bendres,nparpgrp)
c     call zeroarry (bendradi,nparpgrp)
      bendres = 0.
      bendradi = 0.

c  Zero the moments if center is synchv
      if ((center == "synchv" .or. center == "gen" ).and. npmax>0) then
        call getzmmnt(1,xp(1),yp(1),zp(1),uxp(1),uyp(1),uzp(1),gaminv(1),
     &                sq(1),sm(1),sw(1),dt*0.5,1,
     &                nplive,uxpo,uypo,uzpo,1,ns,
     &                tempmaxp,tempminp,tempzmmnts0,tempzmmnts)
#ifdef _OPENMP
        allocate(maxp(6,0:nszmmnt),minp(6,0:nszmmnt),
     &           zmmnts0(NUMZMMNT,0:nszmmnt),
     &           zmmnts(0:nzmmnt,NUMZMMNT,0:nszmmnt),stat=allocerror)
        if (allocerror /= 0) then
          print*,"padvnc3d: allocation error ",allocerror,
     &           ": could not allocate temp arrays to shape ",nszmmnt
          stop
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
c  Set the transverse E fields near any defined apertures.
      if (napertures > 0) call set_aperture_e()

!$OMP PARALLEL
!$OMP&PRIVATE(ip,zpo,ex,ey,ez,bx,by,bz,uxpo,uypo,uzpo,bendres,bendradi,
!$OMP&        uxpadv,uypadv,uzpadv)
!$OMP&FIRSTPRIVATE(maxp,minp,zmmnts0,zmmnts)

c  Set storage for field arrays if not adequate already
      if (npfield < nparpgrp) then
         npfield = nparpgrp
         call gchange("DKInterptmp",0)
      endif

c  Loop over species
      do is=1,ns

c  Advance species according to subcycling rule governed by ndts.
        if (mod(it,ndts(is))/=0) then
          cycle
        else
          if (ndts(is)==1) then
            center = center_in
          else
            if (center_in=="fullv" .and. ldts(is)) center="halfv"
            if (center_in=="synchv") then
              ldts(is)=.true.
            else
              ldts(is)=.false.
            end if
          end if
        end if

c Include the scaling factors in the time step size, both the one for
c subcycling and step size scaling for steady-state and slice modes.
        fulldt_s    = dt*ndts(is)*dtscale(is)
        halfdt_s    = 0.5*fulldt_s

c If this species involves drift-kinetic interpolation, allocate
c  space for temporary arrays if sufficient space not already allocated
c Also set m_over_q, q_over_m**2
         if (interpdk(is) .ne. 0) then
            if (npint < nparpgrp) then
               npint = nparpgrp
               call gchange("DKInterptmp",0)
            endif
            notusealphcalc(is) = (1.-usealphacalc(is))*alpha0(is)
            if (sq(is) .ne. 0.) m_over_q(is)=sm(is)/sq(is)
            if (sm(is) .ne. 0.) qovermsq(is)= (sq(is)/sm(is))**2
         endif

c  Loop over particle blocks; move each block separately
!$OMP DO
         do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)

c           --- Obtain the self-field from the electrostatic potential
            call fetche3d(ipmin,ip,is)

c           --- Get E field for particles near the  injection surface.
            if (inject > 0) then
             call inj_sete3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                       pid(ipmin,injpid),inj_dz,ex,ey,ez)
             if(l_inj_delay_temp) call inj_addtemp3d(ip,ipmin,dz)
            endif
c           --- Get transverse E field for particles near any defined apertures.
            if (napertures > 0) then
              call sete3d_aperture(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                             zgridprv,xmmin,ymmin,zmminp,zmmax,dx,dy,dz,
     &                             nxp,nyp,nzp,ex,ey,l2symtry,l4symtry)
            end if
c           --- Scale the self E-field to get the lowest order relativistic
c           --- correction.
            if (relativity == 1) then
              call sete3d_relativity(ip,ex,ey,vbeam)
            end if

c           --- Zero out B field arrays, but only if B is going
c             to be recomputed; B doesn't get recomputed if 
c             center is "synchv" and interpdk is nonzero and allspecl
c             is 0
            if (.not. (center .eq. "synchv" .and. interpdk(is) .ne. 0 
     &          .and. .not. allspecl)) then
               do i = 1,ip
                  bx(i) = 0.
                  by(i) = 0.
                  bz(i) = 0.
               enddo
            endif
c           --- Compute lag average for experimental damping algorithm
            if (eomdamp /= 0.) then
              call edamp (eomdamp,it,itdamp,center,ip,ex,ey,ez,
     &                    exold(ipmin),eyold(ipmin),ezold(ipmin),
     &                    exlag(ipmin),eylag(ipmin),ezlag(ipmin))
            end if
c           --- FULLV
c           Note if this is a species with interpolation between
c           drift kinetics and full PIC, and allspecl
c           is not 1, then ALWAYS do full-v steps.
c           That is, v is always interpreted as at a half step.
c           In this case fullv and halfv do the same thing and
c           synchv does nothing
            if (center == "fullv" .or. (center == "halfv" .and.
     &             interpdk(is) == 1 .and. .not. allspecl)) then
               do i = 1,ip
                 zpo(i) = zp(ipmin+i-1) 
               enddo 
c              --- Obtain bend radii and residence factors
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),bendres,
     &                       bendradi,-halfdt_s,halfdt_s,.false.)
c              --- Correct Ez_self for warped mesh effect
               call bendez3d (ip,xp(ipmin),zp(ipmin),ez,
     &                        bendres,bendradi,bends,bnezflag,linbend)
c              --- Add in Ez from axially-smoothed gaps 
               call gapfield (ip,zp(ipmin),ez,zbeam,zzmin,egap(0),dzz)
c
c              --- Add in ears and uniform focusing E field pieces
               call othere3d (ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                        zbeam,zimax,zimin,straight,ifeears,eears,
     &                        eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,
     &                        ex,ey,ez,bx,by,bz)
c              --- Set quad, dipole E and B; All: Bz
               call exteb3d (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),-halfdt_s,halfdt_s,
     &                       bx,by,bz,ex,ey,ez,sm(is),sq(is),
     &                       bendres,bendradi,gammabar,fulldt_s)
c              --- Correction to z on entry/exit to accelerator gap
               call zgapcorr(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       -halfdt_s, halfdt_s, fulldt_s, sm(1), sq(1), time)
c              --- First half-electric field increment to momenta
               call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                       ex, ey, ez, sq(is), sm(is), halfdt_s)
c              --- Magnetic field increment to momenta
c              --- If interpolating with drift kinetics, must do half
c              ---  a mu grad B correction before and half after Bpush
               if (interpdk(is) .ne. 0) then
c                  call setptrs(bx,by,bz,ex,ey,ez)
                  call mugrdbpush(ip,is,ipmin,halfdt_s,1)
               endif
               call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bx, by, bz, sq(is), sm(is), fulldt_s, ibpush)
               if (interpdk(is) .ne. 0) call mugrdbpush(ip,is,ipmin,halfdt_s,0)
c              --- Final half-electric field increment to momenta
               call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                       ex, ey, ez, sq(is), sm(is), halfdt_s)
c              --- Advance relativistic Gamma factor
               call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                       gamadv,lrelativ)
c              --- Position advance
               if (interpdk(is) .ne. 0) then
                 call xpush3dintrp(ip,is,ipmin)
               else
                 call xpush3d (ip, xp(ipmin), yp(ipmin), zp(ipmin), uxp(ipmin),
     &                       uyp(ipmin),uzp(ipmin),gaminv(ipmin),fulldt_s)
               endif
c              --- Bend residence factor over POSITION advance step
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bendres,bendradi,-fulldt_s,0.,.false.)
c              --- Correct position advance for warped mesh effect
               call zbendcor(ip,xp(ipmin),zp(ipmin),uxp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),fulldt_s,bendres,bendradi,bends,linbend)
c              --- Correct position advance for slanted dipole entry/exit
               call sledgcor(ip,xp(ipmin),yp(ipmin),ip,zp(ipmin),zpo,
     &                       uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       sm(is),sq(is),.false.)

c           --- HALFV
            elseif ((center == "halfv" .and. (interpdk(is) .ne. 1 .or.
     &               allspecl)) .or. (center == "hbor2" .and.
     &               interpdk(is) .ne. 0)) then
c           --- drift-kinetic-interpolated species are always treated
c               with full-v steps except when initalizing, which is
c               optionally taken care of by "hbor2" option
c               The hbor1 option does 1st half of a Boris step,
c               hbor2 does 2nd half.
               do i = 1,ip
                 zpo(i) = zp(ipmin+i-1) 
               enddo 
c              --- Obtain bend radii and residence factors
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bendres,bendradi,0.,halfdt_s,.false.)
c              --- Correct Ez_self for warped mesh effect
               call bendez3d (ip,xp(ipmin),zp(ipmin),ez,
     &                        bendres,bendradi,bends,bnezflag,linbend)
c              --- Add in Ez from axially-smoothed gaps 
               call gapfield (ip,zp(ipmin),ez,zbeam,zzmin,egap(0),dzz)
c              --- Add in ears and uniform focusing E field pieces
               call othere3d (ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                        zbeam,zimax,zimin,straight,ifeears,eears,
     &                        eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,
     &                        ex,ey,ez,bx,by,bz)
c              --- Set quad, dipole E and B;  All: Bz
               call exteb3d (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),0.,halfdt_s,bx,by,bz,
     &                       ex,ey,ez,sm(is),sq(is),bendres,bendradi,
     &                       gammabar,fulldt_s)
c              --- Correction to z on entry/exit to accelerator gap
               call zgapcorr(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       0., halfdt_s, fulldt_s, sm(1), sq(1), time)
c              --- Magnetic field increment to momenta
               if (interpdk(is) .ne. 0) then
                  call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bx, by, bz, sq(is), sm(is), fulldt_s, -1)
c              --- If interpolating with drift kinetics, must do 
c              ---  a mu grad B correction.  Since may be calling
c              ---  this on first timestep, must make sure grad B's are
c              ---  calculated, so last arg = 1.
                  call mugrdbpush(ip,is,ipmin,halfdt_s,1)
               else
                  call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bx, by, bz, sq(is), sm(is), halfdt_s, ibpush)
c              --- Final half-electric field increment to momenta
               endif
               call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                       ex, ey, ez, sq(is), sm(is), halfdt_s)
c              --- Advance relativistic Gamma factor
               call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                       gamadv,lrelativ)
c              --- Position advance
               if (center .ne. "hbor2") then
c              hbor2 just does a v advance, no x push
                 if (interpdk(is) .ne. 0) then
                   call xpush3dintrp(ip,is,ipmin)
                 else
                   call xpush3d (ip, xp(ipmin), yp(ipmin), zp(ipmin), uxp(ipmin),
     &                       uyp(ipmin),uzp(ipmin),gaminv(ipmin),fulldt_s)
                 endif
c                --- Bend residence factor over POSITION advance step
                 call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                         bendres,bendradi,-fulldt_s,0.,.false.)
c                --- Correct position advance for warped mesh effect
                 call zbendcor(ip,xp(ipmin),zp(ipmin),uxp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),fulldt_s,bendres,bendradi,bends,linbend)
c                --- Correct position advance for slanted dipole entry/exit
                 call sledgcor(ip,xp(ipmin),yp(ipmin),ip,zp(ipmin),zpo,
     &                       uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       sm(is),sq(is),.false.)
               endif
c           --- SYNCHV or GEN
c           Double "if" follows so that moments are calculated
c           for center = "synchv" or "gen" regardless of interpdk,
c           but the vpush is only done for interpdk=1 if allspecl=1.
            elseif (center == "synchv" .or. center == "gen"
     &        .or. center == "hbor1") then
              if ((center == "synchv" .and. (interpdk(is) == 0
     &          .or. allspecl)) 
     &         .or. center == "gen" .or. (center == "hbor1"
     &          .and. interpdk(is) .ne. 0)) then
c           --- drift-kinetic-interpolated species are always treated
c               with full-v steps except in initial setup, when we want
c               a half-v step packwards.  center = "hbor2" provides
c               the latter capability with dt set to be negative.
c               Exception, if allspecl = 1, then we synchronize every
c               step, using half Boris if interpdk = 1.
c               More generally hbor1 provides 1st half of a Boris
c               step, hbor2 provides 2nd half.
c              --- Copy 'old' velocity into uxpo, uypo, and uzpo
               do i=1,ip
                  uxpo(i) = uxp(ipmin+i-1)
                  uypo(i) = uyp(ipmin+i-1)
                  uzpo(i) = uzp(ipmin+i-1)
               enddo
c              --- Obtain bend radii and residence factors
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),bendres,
     &                       bendradi,-halfdt_s,0.,.false.)
c              --- Correct Ez_self for warped mesh effect
               call bendez3d (ip,xp(ipmin),zp(ipmin),ez,
     &                        bendres,bendradi,bends,bnezflag,linbend)
c              --- Add in Ez from axially-smoothed gaps 
               call gapfield (ip,zp(ipmin),ez,zbeam,zzmin,egap(0),dzz)
c              --- Add in ears and uniform focusing E field pieces
               call othere3d (ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                        zbeam,zimax,zimin,straight,ifeears,eears,
     &                        eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,
     &                        ex,ey,ez,bx,by,bz)
c              --- Set quad, dipole E and B; All: Bz
               call exteb3d (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),-halfdt_s,0.,
     &                       bx,by,bz,ex,ey,ez,sm(is),sq(is),
     &                       bendres,bendradi,gammabar,fulldt_s)
c              --- Half electric field increment to momenta
               call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                       ex, ey, ez, sq(is), sm(is), halfdt_s)
c              --- If interpolating with drift kinetics, must do half
c              ---  a mu grad B correction
c              --- Half magnetic field increment to momenta
               if (interpdk(is) .ne. 0) then
                  call mugrdbpush(ip,is,ipmin,halfdt_s,1)
                  call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bx, by, bz, sq(is), sm(is), fulldt_s, -1)
c                 Half a Boris push
               else
                  call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bx, by, bz, sq(is), sm(is), halfdt_s, ibpush)
               endif
               if (center == "gen") then
c                --- Reset uxp to uxpo, set uxpo to half step backward
c                --- for interpolation in moments calculation
c                --- Don't do this if we're doing drift-kinetic interp
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
              endif
c             --- Calculate moments over particles, now that we're sync'd
c             --- THIS IS NOT DONE FOR center = hbor1.
              if (center == "gen" .or. center == "synchv") then
                if(wpid==0) then
                 call getzmmnt(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                         uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                         gaminv(ipmin),sq(is),sm(is),sw(is),
     &                         halfdt_s,2,nplive,uxpo,uypo,uzpo,is,ns,
     &                         maxp,minp,zmmnts0,zmmnts)
                else
                 call getzmmnt_weights(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                         uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                         gaminv(ipmin),pid(ipmin,wpid),
     &                         sq(is),sm(is),sw(is),halfdt_s,2,nplive,
     &                         uxpo,uypo,uzpo,is,ns,
     &                         maxp,minp,zmmnts0,zmmnts)
                endif
                if (npid > 0) then
                  call getextrapolatedparticles(ip,xp(ipmin),yp(ipmin),
     &                           zp(ipmin),
     &                           uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                           npidmax,pid(ipmin:ipmin+ip-1,:),
     &                           dt/2.,uxpo,uypo,uzpo,is,time,zbeam)
                else
                  call getextrapolatedparticles(ip,xp(ipmin),yp(ipmin),
     &                           zp(ipmin),
     &                           uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                           0,0.,
     &                           dt/2.,uxpo,uypo,uzpo,is,time,zbeam)
                endif
              endif
            endif
         enddo
!$OMP END DO

      enddo

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
          print*,"padvnc3d: allocation error ",allocerror,
     &           ": could not allocate temp arrays to shape ",nszmmnt,nthread
          stop
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

c----------------------------------------------------------------------------
c     ---  Do final stuff for moments calculation
      if ((center == "synchv" .or. center == "gen") .and. npmax>0) then
#ifdef _OPENMP
        tempmaxp = max(tempmaxp,maxval(threadmaxp(:,1:nthread),3))
        tempminp = min(tempminp,minval(threadminp(:,1:nthread),3))
        tempzmmnts0 = tempzmmnts0 + sum(threadzmmnts0(:,1:nthread),3)
        tempzmmnts = tempzmmnts + sum(threadzmmnts(:,:,1:nthread),4)
        deallocate(threadmaxp,threadminp,threadzmmnts0,threadzmmnts)
#endif
        call getzmmnt(1,xp(1),yp(1),zp(1),uxp(1),uyp(1),uzp(1),gaminv(1),
     &                sq(1),sm(1),sw(1),halfdt,3,nplive,uxpo,uypo,uzpo,
     &                1,ns,tempmaxp,tempminp,tempzmmnts0,tempzmmnts)
        call getlabwn()
c       --- load rho diagnostic quantities
        call rhodia3d
      endif

      if (center == "synchv" .or. center == "gen" .or. ifgap) then

c       --- Calculate current

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

        if (.not. laccumulate_rho) call zeroarry(curr,(nzzarr+1)*(1+nszarr))

!$OMP PARALLEL PRIVATE(ip)
        do is=1,ns
!$OMP DO
          do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)
            if(wpid==0) then
              call setcurr(nzzarr,nszarr,curr,ip,zp(ipmin),uzp(ipmin),
     &                     gaminv(ipmin),sq(is),sw(is)*dtscale(is),is,
     &                     zbeam,dzz,zzmin,ns,lspeciesmoments)
            else
              call setcurr_weights(nzzarr,nszarr,curr,ip,zp(ipmin),uzp(ipmin),
     &                             pid(ipmin,wpid),gaminv(ipmin),
     &                             sq(is),sw(is)*dtscale(is),is,zbeam,dzz,zzmin,
     &                             ns,lspeciesmoments)
            end if
          enddo
!$OMP END DO
        enddo
!$OMP END PARALLEL
        call percurr(curr,nzzarr,nszarr,bound0)
      endif
c----------------------------------------------------------------------------
      if (center == "fullv" .or. center == "halfv") then
c       --- Inject more particles. This adds particles and alters ins and nps.

        if (inject > 0) then
           call inject3d(1)
        end if

c       --- Treat particles at boundaries
c        if (any((mod(it,ndts)==0) .and. (nps>0))) then
        call particleboundaries3d()
 
c       --- Collect charge density (3d)
        call loadrho3d(-1,-1,-1,.not. laccumulate_rho)

        if (lbeamcom) then
c         --- Set zbeam so that it follows the center of mass of the beam.
          zbeam = getbeamcom() - zbeamcomoffset

        else
c         --- Advance beam frame location using the nominal beam frame velocity.
          zbeam = zgrid
        endif

c       --- zgridprv needs to be updated for the "synchv" step
c       --- Note that zgridprv is also set at the beginning of w3dexe
        zgridprv = zgrid

        time = time + dt
      endif
c----------------------------------------------------------------------------

!$OMP MASTER
      if (lw3dtimesubs) timepadvnc3d = timepadvnc3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine getextpart()
      use GlobalVars
      use InGen
      use InPart
      use Particles
      use Picglb

      integer(ISZ):: is,ip,ipmin
      real(kind=8):: fulldt_s,halfdt_s

c  Loop over species
      do is=1,ns

c  Advance species according to subcycling rule governed by ndts.
c  Also include time step size scaling for steady-state and slice modes.
        if (mod(it,ndts(is))/=0) cycle
        fulldt_s    = dt*ndts(is)*dtscale(is)
        halfdt_s    = 0.5*fulldt_s

         do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)
                if (npid > 0) then
                  call getextrapolatedparticles(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                           uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                           npidmax,pid(ipmin:ipmin+ip-1,:),
     &                           halfdt_s,uxp(ipmin),uyp(ipmin),uzp(ipmin),is,time,zbeam)
                else
                  call getextrapolatedparticles(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                           uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                           0,0.,
     &                           halfdt_s,uxp(ipmin),uyp(ipmin),uzp(ipmin),is,time,zbeam)
                endif
         enddo
       enddo

      return
      end
c=============================================================================
      subroutine perphi3d(phi,nx,ny,nz)
      use Subtimers3d
      use InGen3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      integer(ISZ):: ix,iy
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Sets the slices on the exterior of phi for periodicity
c  sets slice at -1 equal to the slice at nz-1
c  sets slice at nz+1 equal to the slice at 1

      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        call perphirz()
        return
      end if

#ifdef MPIPARALLEL
      call perphi3d_slave(phi(0,0,-1),nx,ny,nz)
#else
      do iy=0,ny
        do ix=0,nx
          phi(ix,iy,-1)   = phi(ix,iy,nz-1)
          phi(ix,iy,nz)   = phi(ix,iy,0)
          phi(ix,iy,nz+1) = phi(ix,iy,1)
        enddo
      enddo
#endif

!$OMP MASTER
      if (lw3dtimesubs) timeperphi3d = timeperphi3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine perrho3d(rho,nx,ny,nz,bound0,boundxy)
      use GlobalVars
      use Subtimers3d
      use InGen3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      integer(ISZ):: bound0,boundxy
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Sums the first and last slices of rho for periodicity
c  and puts the result into both slices.

      if (solvergeom==AMRgeom) return

c  Enforce transverse periodicity if requested.
      if (solvergeom==XYZgeom) then
        if (boundxy==periodic) then
          rho(:,0,:) = rho(:,0,:) + rho(:,ny,:)
          rho(0,:,:) = rho(0,:,:) + rho(nx,:,:)
          rho(:,ny,:) = rho(:,0,:)
          rho(nx,:,:) = rho(0,:,:)
        endif
        if (boundxy==neumann) then
          if (.not. (l2symtry .or. l4symtry)) rho(:,0,:) = 2.*rho(:,0,:)
          if (.not. l4symtry) rho(0,:,:) = 2.*rho(0,:,:)
          rho(:,ny,:) = 2.*rho(:,ny,:)
          rho(nx,:,:) = 2.*rho(nx,:,:)
        endif
      endif

c  Distribute rho for 2d solver
      if(solvergeom==RZgeom .or. solvergeom==XZgeom .or.
     &   solvergeom==XYgeom .or. solvergeom==Zgeom) then
        call distribute_rho_rz()
      end if

c  Enforce longitudinal periodicity.
      if(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==XYgeom) then
        call rhobndrz()
#ifdef MPIPARALLEL
        if(bound0==periodic) call perrhorz()
#endif
      end if

c  Copy charge density from frz.basegrid to w3d.rho
      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        call get_rho_rz(rho(0,0,0),nx,nz,1,0)
      else if(solvergeom==XYgeom) then
        call get_rho_rz(rho(0,0,0),nx,ny,1,0)
      else if(solvergeom==Zgeom) then
        call get_rho_z(rho(0,0,0),nz,1,0)
      else if(solvergeom==Rgeom) then
        call get_rho_r(rho(0,0,0),nx,1,0)
      end if

      if (solvergeom==XYZgeom .and. bound0==periodic) then
#ifdef MPIPARALLEL
        call perrho3d_slave(rho,nx,ny,nz)
#else
        rho(:,:,0)  = rho(:,:,0) + rho(:,:,nz)
        rho(:,:,nz) = rho(:,:,0)
#endif
      endif

!$OMP MASTER
      if (lw3dtimesubs) timeperrho3d = timeperrho3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine epush3d(np,uxp,uyp,uzp,ex,ey,ez,q,m,ddt)
      use Subtimers3d
      integer(ISZ):: np
      real(kind=8):: uxp(np),uyp(np),uzp(np)
      real(kind=8):: ex(np),ey(np),ez(np)
      real(kind=8):: q,m,ddt

c  Push the particle velocity with E field

      integer(ISZ):: ip
      real(kind=8):: const
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      const = q*ddt/m

      do ip=1,np
        uxp(ip) = uxp(ip) + ex(ip)*const
        uyp(ip) = uyp(ip) + ey(ip)*const
        uzp(ip) = uzp(ip) + ez(ip)*const
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timeepush3d = timeepush3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine epusht3d(np,uxp,uyp,uzp,ex,ey,ez,q,m,dtp,fdt)
      use Subtimers3d
      integer(ISZ):: np
      real(kind=8):: uxp(np),uyp(np),uzp(np)
      real(kind=8):: ex(np),ey(np),ez(np)
      real(kind=8):: dtp(np)
      real(kind=8):: q,m,fdt

c Push the particle velocity with E field using a different time step for
c each particle.

      integer(ISZ):: ip
      real(kind=8):: const
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      const = q*fdt/m

      do ip=1,np
        uxp(ip) = uxp(ip) + ex(ip)*const*dtp(ip)
        uyp(ip) = uyp(ip) + ey(ip)*const*dtp(ip)
        uzp(ip) = uzp(ip) + ez(ip)*const*dtp(ip)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timeepusht3d = timeepusht3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c============================================================================= 
      subroutine bpush3d(np,uxp,uyp,uzp,gaminv,bx,by,bz,q,m,ddt,ibpush)
      use Subtimers3d
      integer(ISZ):: np,ibpush
      real(kind=8):: uxp(np),uyp(np),uzp(np),gaminv(np)
      real(kind=8):: bx(np),by(np),bz(np)
      real(kind=8):: q,m,ddt

c  Push the particle velocity with B field

      integer(ISZ):: ip
      real(kind=8):: btot,btotinv,tanalphab
      real(kind=8):: tx,ty,tz,tsqi,sx,sy,sz,uxppr,uyppr,uzppr
      real(kind=8):: const,t,tiny,const1
      real(kind=8):: substarttime,wtime

      if (lw3dtimesubs) substarttime = wtime()
      const = q*ddt*0.5/m

      if (ibpush == 1) then
c        --- fast b-field rotation algorithm
         do ip=1,np
            tx = gaminv(ip)*bx(ip)*const
            ty = gaminv(ip)*by(ip)*const
            tz = gaminv(ip)*bz(ip)*const
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
            tanalphab = tan(gaminv(ip)*btot*const)/btot
            tx = bx(ip)*tanalphab
            ty = by(ip)*tanalphab
            tz = bz(ip)*tanalphab
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
      elseif (ibpush == -1) then
c        --- fast b-field rotation algorithm with half Boris angle
         tiny = 1.e-20
         do ip=1,np
            btot = sqrt(bx(ip)**2 + by(ip)**2 + bz(ip)**2)
            t=gaminv(ip)*btot*const
            if (btot == 0 .or. t == 0.) cycle
            btot=max(btot,tiny)
            const1=(sqrt(1.+t*t)-1.)/(t*btot)
            tx = bx(ip)*const1
            ty = by(ip)*const1
            tz = bz(ip)*const1
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

!$OMP MASTER
      if (lw3dtimesubs) timebpush3d = timebpush3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c============================================================================= 
      subroutine bpusht3d(np,uxp,uyp,uzp,gaminv,bx,by,bz,q,m,dtp,fdt,ibpush)
      use Subtimers3d
      integer(ISZ):: np,ibpush
      real(kind=8):: uxp(np),uyp(np),uzp(np),gaminv(np)
      real(kind=8):: bx(np),by(np),bz(np)
      real(kind=8):: dtp(np)
      real(kind=8):: q,m,fdt

c Push the particle velocity with B field using a different time step
c for each particle.

      integer(ISZ):: ip
      real(kind=8):: btot,btotinv,tanalpha
      real(kind=8):: tx,ty,tz,tsqi,sx,sy,sz,uxppr,uyppr,uzppr
      real(kind=8):: const
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
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

!$OMP MASTER
      if (lw3dtimesubs) timebpusht3d = timebpusht3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=========================================================================== 
      subroutine xpush3d(np,xp,yp,zp,uxp,uyp,uzp,gaminv,ddt)
      use Subtimers3d
      integer(ISZ):: np
      real(kind=8):: ddt
      real(kind=8):: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np),gaminv(np)

c  Advance particle positions

      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      do ip=1,np
        xp(ip) = xp(ip) + uxp(ip)*ddt*gaminv(ip)
        yp(ip) = yp(ip) + uyp(ip)*ddt*gaminv(ip)
        zp(ip) = zp(ip) + uzp(ip)*ddt*gaminv(ip)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timexpush3d = timexpush3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c===========================================================================
      subroutine xpusht3d(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dtp)
      use Subtimers3d
      integer(ISZ):: np
      real(kind=8):: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np),gaminv(np)
      real(kind=8):: dtp(np)

c  Advance particle positions

      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      do ip=1,np
        xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dtp(ip)
        yp(ip) = yp(ip) + uyp(ip)*gaminv(ip)*dtp(ip)
        zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dtp(ip)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timexpusht3d = timexpusht3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c============================================================================
      subroutine sete3d_relativity(np,ex,ey,vbeam)
      use Subtimers3d
      use Constant
      integer(ISZ):: np
      real(kind=8):: ex(np),ey(np)
      real(kind=8):: vbeam
c Applies the first order relativistic correction to the self E-field.

      real(kind=8):: gammabarisq
      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      gammabarisq = 1. -  (vbeam/clight)**2

      do ip=1,np
        ex(ip) = ex(ip)*gammabarisq
        ey(ip) = ey(ip)*gammabarisq
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesete3d_relativity = timesete3d_relativity + wtime() - substarttime
!$OMP END MASTER
      return
      end
c============================================================================ 
      subroutine edamp(eomdamp,it,itdamp,center,np,
     &                 ex,ey,ez,exold,eyold,ezold,exlag,eylag,ezlag)
      use Subtimers3d
      real(kind=8):: eomdamp
      integer(ISZ):: it,itdamp,np
      real(kind=8):: ex(np), ey(np), ez(np), exold(np), eyold(np), ezold(np)
      real(kind=8):: exlag(np), eylag(np), ezlag(np)
      character(*):: center

c  Computes modified electric field to apply "adjustably damped" mover.
c  Updates "old" and "lag" qtys only after "fullv" or "halfv" advance,
c    never after "synchv" since they are needed by the "halfv" to come.

      integer(ISZ):: ip
      real(kind=8):: exip,eyip,ezip
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (eomdamp == 0.) then
!$OMP MASTER
        if (lw3dtimesubs) timeedamp = timeedamp + wtime() - substarttime
!$OMP END MASTER
        return
      endif

c  Protect against doing algebra on garbage; give lag qtys a good start.

      if (it == 1) then
         do ip = 1, np
            exold(ip) = ex(ip)
            eyold(ip) = ey(ip)
            ezold(ip) = ez(ip)
            exlag(ip) = ex(ip)
            eylag(ip) = ey(ip)
            ezlag(ip) = ez(ip)
         enddo
      endif

c  Main loop to set e, and perhaps update old and lag qtys

      if (center == "synchv") then
c        --- only modify e if late enough that lags are well established
         if (it >= itdamp) then
            if (eomdamp > 0.) then
c              --- second order scheme
               do ip = 1, np
                  ex(ip) = 0.5 * ( (2.+0.5*eomdamp)*ex(ip) - exold(ip)
     &                     + (1. - 0.5*eomdamp)*exlag(ip) )
                  ey(ip) = 0.5 * ( (2.+0.5*eomdamp)*ey(ip) - eyold(ip)
     &                     + (1. - 0.5*eomdamp)*eylag(ip) )
                  ez(ip) = 0.5 * ( (2.+0.5*eomdamp)*ez(ip) - ezold(ip)
     &                     + (1. - 0.5*eomdamp)*ezlag(ip) )
               enddo
            else
c              --- first order backward biased scheme
               do ip = 1, np
                  ex(ip) = (1.-eomdamp)*ex(ip) + eomdamp*exold(ip)
                  ey(ip) = (1.-eomdamp)*ey(ip) + eomdamp*eyold(ip)
                  ez(ip) = (1.-eomdamp)*ez(ip) + eomdamp*ezold(ip)
               enddo
            endif
         endif
      else
c        --- modify e and update old and lag qtys
         if (eomdamp > 0.) then
c           --- second order scheme
            do ip = 1, np
               exip = ex(ip)
               eyip = ey(ip)
               ezip = ez(ip)
               ex(ip) = 0.5 * ( (2.+0.5*eomdamp)*ex(ip) - exold(ip)
     &                  + (1. - 0.5*eomdamp)*exlag(ip) )
               ey(ip) = 0.5 * ( (2.+0.5*eomdamp)*ey(ip) - eyold(ip)
     &                  + (1. - 0.5*eomdamp)*eylag(ip) )
               ez(ip) = 0.5 * ( (2.+0.5*eomdamp)*ez(ip) - ezold(ip)
     &                  + (1. - 0.5*eomdamp)*ezlag(ip) )
               exold(ip) = exip
               eyold(ip) = eyip
               ezold(ip) = ezip
               exlag(ip) = (1.-0.5*eomdamp)*exip + 0.5*eomdamp*exlag(ip)
               eylag(ip) = (1.-0.5*eomdamp)*eyip + 0.5*eomdamp*eylag(ip)
               ezlag(ip) = (1.-0.5*eomdamp)*ezip + 0.5*eomdamp*ezlag(ip)
            enddo
         else
c           --- first order backward biased scheme
            do ip = 1, np
               exip = ex(ip)
               eyip = ey(ip)
               ezip = ez(ip)
               ex(ip) = (1.-eomdamp)*ex(ip) + eomdamp*exold(ip)
               ey(ip) = (1.-eomdamp)*ey(ip) + eomdamp*eyold(ip)
               ez(ip) = (1.-eomdamp)*ez(ip) + eomdamp*ezold(ip)
               exold(ip) = exip
               eyold(ip) = eyip
               ezold(ip) = ezip
            enddo
         endif
c        --- don't use modified field if too early in run
         if (it < itdamp+1) then
            do ip = 1, np
               ex(ip) = exold(ip)
               ey(ip) = eyold(ip)
               ez(ip) = ezold(ip)
            enddo
         endif
      endif

!$OMP MASTER
      if (lw3dtimesubs) timeedamp = timeedamp + wtime() - substarttime
!$OMP END MASTER
      return
      end
c============================================================================ 
      subroutine bendez3d(np,xp,zp,ez,bendres,bendradi,bends,bnezflag,linbend)
      use Subtimers3d
      integer(ISZ):: np
      real(kind=8):: xp(np), zp(np), ez(np), bendres(np), bendradi(np)
      logical(ISZ):: bends, bnezflag, linbend

c  Corrects axial electric field at particle position for warped geometry
c  via multiplying by r_star/r = 1 - x/r, in a residence-corrected way;
c  at smaller radii, zones are closer together, so field is larger.

      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (.not. (bends.and.bnezflag) .or. .not. linbend) then
!$OMP MASTER
        if (lw3dtimesubs) timebendez3d = timebendez3d + wtime() - substarttime
!$OMP END MASTER
        return
      endif

      do ip=1,np
         ez(ip) = ez(ip)*(1. - bendres(ip)*xp(ip)/(bendradi(ip) + xp(ip)))
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timebendez3d = timebendez3d + wtime() - substarttime
!$OMP END MASTER
      return
      end      
c=========================================================================== 
      subroutine zbendcor(np,xp,zp,uxp,uzp,gaminv,ddt,bendres,bendradi,bends,
     &                    linbend)
      use Subtimers3d
      integer(ISZ):: np
      real(kind=8):: ddt
      real(kind=8):: xp(np), zp(np), uxp(np), uzp(np), gaminv(np)
      real(kind=8):: bendres(np), bendradi(np)
      logical(ISZ):: bends,linbend

c  Corrects particle axial positions for bending; zp is position
c  "along the centerline."

      integer(ISZ):: ip
      real(kind=8):: xprv,xc
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (.not. bends .or. .not. linbend) then
!$OMP MASTER
        if (lw3dtimesubs) timezbendcor = timezbendcor + wtime() - substarttime
!$OMP END MASTER
        return
      endif

      do ip=1,np
        xprv = xp(ip) - ddt * uxp(ip) * gaminv(ip)
        xc = 0.5 * (xp(ip) + xprv)
        zp(ip) = zp(ip) + ddt*uzp(ip)*gaminv(ip)*bendres(ip)
     &                    *(bendradi(ip)/(bendradi(ip) + xc) - 1.)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timezbendcor = timezbendcor + wtime() - substarttime
!$OMP END MASTER
      return
      end      
c=========================================================================== 
      subroutine sete3d(phi1d,selfe,np,xp,yp,zp,zgrid,xmmin,ymmin,zmmin,
     &                  dx,dy,dz,nx,ny,nz,efetch,ex,ey,ez,l2symtry,l4symtry)
      use Subtimers3d
      integer(ISZ):: np,nx,ny,nz
      real(kind=8):: zgrid,xmmin,ymmin,zmmin,dx,dy,dz
      real(kind=8):: phi1d(0:*),selfe(3,0:nx,0:ny,0:nz),xp(np),yp(np),zp(np)
      real(kind=8):: ex(np),ey(np),ez(np)
      integer(ISZ):: efetch
      logical(ISZ):: l2symtry,l4symtry

c  Sets self electric field for particles
c  Note that the phi1d passed in is assumed to start at phi(0,0,-1).
c  This is a change is behavior from what is was earlier. This change
c  was made to avoid problems with bounds checking. phi needs to be accessed
c  at the plane iz=-1, so when phi was passed in starting at phi(0,0,0), with
c  bounds checking turned on, this would be caught (even though the code was
c  technically correct).

c  Algorithm notes: phi array is dimensioned (0:nx,0:ny,-1:nz+1) outside,
c  but is made one dimensional in this routine
c  so cell index into 1d phi array for vectorized deposition is:
c     i + j*(nx+1) + k*(nx+1)*(ny+1)
c  The field is:
c     Ex = u0 * v0 * w0 * ex(i  ,j  ,k  )
c        + u1 * v0 * w0 * ex(i+1,j  ,k  )
c        + u0 * v1 * w0 * ex(i  ,j+1,k  )
c        + ...

      integer(ISZ):: nnx,nnxy,ip,i,j,k,ind0,m
      real(kind=8):: dxi,dyi,dzi,tdxi,tdyi,tdzi,u0,u1,v0,v1,w0,w1,ysign,xsign
      integer(ISZ):: ox,oy,sox,soy
      real(kind=8):: sx,sy
      integer(ISZ):: noff(32)
      save noff

c     --- For second version
c     --- Use a first dimension of 33 not 32 to avoid bank conflicts
      integer(ISZ),allocatable:: indx(:,:)
      real(kind=8),allocatable:: p(:,:), u1p(:), v1p(:), w1p(:)
      integer(ISZ):: allocerror

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      nnx  = nx + 1
      nnxy = (nx + 1) * (ny + 1)

c  Evaluation of offsets occurs on the first call only provided
c  local arrays are SAVE'd. nnxy is added to all offsets to account for
c  the fact that the phi1d passed begins at phi(0,0,-1), so the location
c  of phi(0,0,0) is equivalent to phi1d(nnxy).

c     if (noff(1) /= -nnxy) then
         noff(1)  = nnxy - nnxy
         noff(2)  = nnxy - nnxy   + 1
         noff(3)  = nnxy - nnxy   + nnx
         noff(4)  = nnxy - nnxy   + nnx   + 1
         noff(5)  = nnxy          - nnx
         noff(6)  = nnxy          - nnx   + 1
         noff(7)  = nnxy                  - 1
         noff(8)  = nnxy                  + 0
         noff(9)  = nnxy                  + 1
         noff(10) = nnxy                  + 2
         noff(11) = nnxy          + nnx   - 1
         noff(12) = nnxy          + nnx
         noff(13) = nnxy          + nnx   + 1
         noff(14) = nnxy          + nnx   + 2
         noff(15) = nnxy          + 2*nnx
         noff(16) = nnxy          + 2*nnx + 1
         noff(17) = nnxy + nnxy   - nnx
         noff(18) = nnxy + nnxy   - nnx   + 1
         noff(19) = nnxy + nnxy           - 1
         noff(20) = nnxy + nnxy
         noff(21) = nnxy + nnxy           + 1
         noff(22) = nnxy + nnxy           + 2
         noff(23) = nnxy + nnxy   + nnx   - 1
         noff(24) = nnxy + nnxy   + nnx
         noff(25) = nnxy + nnxy   + nnx   + 1
         noff(26) = nnxy + nnxy   + nnx   + 2
         noff(27) = nnxy + nnxy   + 2*nnx
         noff(28) = nnxy + nnxy   + 2*nnx + 1
         noff(29) = nnxy + 2*nnxy
         noff(30) = nnxy + 2*nnxy         + 1
         noff(31) = nnxy + 2*nnxy + nnx
         noff(32) = nnxy + 2*nnxy + nnx   + 1
c     endif

c  Evaluation of E, vectorized over particles
      tdxi = 1. / (2. * dx)
      tdyi = 1. / (2. * dy)
      tdzi = 1. / (2. * dz)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      if (efetch == 1) then

      if (.not. (l2symtry .or. l4symtry)) then
        do ip = 1, np

           i =  (xp(ip) - xmmin) * dxi
           j =  (yp(ip) - ymmin) * dyi
           k =  (zp(ip) - zgrid - zmmin) * dzi
           ind0 = i + j*nnx + k*nnxy

           u1 = (xp(ip) - xmmin) * dxi - i
           v1 = (yp(ip) - ymmin) * dyi - j
           w1 = (zp(ip) - zgrid - zmmin) * dzi - k

           u0 = 1. - u1
           v0 = 1. - v1
           w0 = 1. - w1

           ex(ip)=tdxi*(u0*v0*w0*(phi1d(noff( 7)+ind0) - phi1d(noff( 9)+ind0))
     &                + u1*v0*w0*(phi1d(noff( 8)+ind0) - phi1d(noff(10)+ind0))
     &                + u0*v1*w0*(phi1d(noff(11)+ind0) - phi1d(noff(13)+ind0))
     &                + u1*v1*w0*(phi1d(noff(12)+ind0) - phi1d(noff(14)+ind0))
     &                + u0*v0*w1*(phi1d(noff(19)+ind0) - phi1d(noff(21)+ind0))
     &                + u1*v0*w1*(phi1d(noff(20)+ind0) - phi1d(noff(22)+ind0))
     &                + u0*v1*w1*(phi1d(noff(23)+ind0) - phi1d(noff(25)+ind0))
     &                + u1*v1*w1*(phi1d(noff(24)+ind0) - phi1d(noff(26)+ind0)))

           ey(ip)=tdyi*(u0*v0*w0*(phi1d(noff( 5)+ind0) - phi1d(noff(12)+ind0))
     &                + u1*v0*w0*(phi1d(noff( 6)+ind0) - phi1d(noff(13)+ind0))
     &                + u0*v1*w0*(phi1d(noff( 8)+ind0) - phi1d(noff(15)+ind0))
     &                + u1*v1*w0*(phi1d(noff( 9)+ind0) - phi1d(noff(16)+ind0))
     &                + u0*v0*w1*(phi1d(noff(17)+ind0) - phi1d(noff(24)+ind0))
     &                + u1*v0*w1*(phi1d(noff(18)+ind0) - phi1d(noff(25)+ind0))
     &                + u0*v1*w1*(phi1d(noff(20)+ind0) - phi1d(noff(27)+ind0))
     &                + u1*v1*w1*(phi1d(noff(21)+ind0) - phi1d(noff(28)+ind0)))

           ez(ip)=tdzi*(u0*v0*w0*(phi1d(noff( 1)+ind0) - phi1d(noff(20)+ind0))
     &                + u1*v0*w0*(phi1d(noff( 2)+ind0) - phi1d(noff(21)+ind0))
     &                + u0*v1*w0*(phi1d(noff( 3)+ind0) - phi1d(noff(24)+ind0))
     &                + u1*v1*w0*(phi1d(noff( 4)+ind0) - phi1d(noff(25)+ind0))
     &                + u0*v0*w1*(phi1d(noff( 8)+ind0) - phi1d(noff(29)+ind0))
     &                + u1*v0*w1*(phi1d(noff( 9)+ind0) - phi1d(noff(30)+ind0))
     &                + u0*v1*w1*(phi1d(noff(12)+ind0) - phi1d(noff(31)+ind0))
     &                + u1*v1*w1*(phi1d(noff(13)+ind0) - phi1d(noff(32)+ind0)))
        enddo

      else

c       --- Set offsets for indices on axis of symmetry.  The offsets change
c       --- the sign of the grid cells which are on the negative side
c       --- of the axis of symmetry.
        soy = 2*nnx
        sox = 0
        if (l4symtry) sox = 2

c       --- Set the signs of the E field for particles on negative side of
c       --- the axis of symmetry.
        sy = -1.
        sx = 1.
        if (l4symtry) sx = -1.

c       --- special loop symmetry is used
        do ip = 1, np

           i =  (abs(xp(ip)) - xmmin)*dxi
           j =  (abs(yp(ip)) - ymmin)*dyi
           k =  (zp(ip) - zgrid - zmmin)*dzi
           ind0 = i + j*nnx + k*nnxy

           u1 = (abs(xp(ip)) - xmmin)*dxi - i
           v1 = (abs(yp(ip)) - ymmin)*dyi - j
           w1 = (zp(ip) - zgrid - zmmin)*dzi - k

           u0 = 1. - u1
           v0 = 1. - v1
           w0 = 1. - w1

c          --- Set offsets for points on symmetry axis of grid.  The offset
c          --- for points off the axis is zero.
           ox = 0
           oy = 0
           if (i == 0 .and. xmmin == 0.) ox = sox
           if (j == 0 .and. ymmin == 0.) oy = soy

c          --- Adjust sign of E field for appropriate quadrant.
           xsign = tdxi
           ysign = tdyi
           if (xp(ip) < 0.) xsign = sx*tdxi
           if (yp(ip) < 0.) ysign = sy*tdyi
          ex(ip)=xsign*(u0*v0*w0*(phi1d(noff( 7)+ind0+ox)-phi1d(noff( 9)+ind0))
     &                + u1*v0*w0*(phi1d(noff( 8)+ind0   )-phi1d(noff(10)+ind0))
     &                + u0*v1*w0*(phi1d(noff(11)+ind0+ox)-phi1d(noff(13)+ind0))
     &                + u1*v1*w0*(phi1d(noff(12)+ind0   )-phi1d(noff(14)+ind0))
     &                + u0*v0*w1*(phi1d(noff(19)+ind0+ox)-phi1d(noff(21)+ind0))
     &                + u1*v0*w1*(phi1d(noff(20)+ind0   )-phi1d(noff(22)+ind0))
     &                + u0*v1*w1*(phi1d(noff(23)+ind0+ox)-phi1d(noff(25)+ind0))
     &                + u1*v1*w1*(phi1d(noff(24)+ind0   )-phi1d(noff(26)+ind0)))

          ey(ip)=ysign*(u0*v0*w0*(phi1d(noff( 5)+ind0+oy)-phi1d(noff(12)+ind0))
     &                + u1*v0*w0*(phi1d(noff( 6)+ind0+oy)-phi1d(noff(13)+ind0))
     &                + u0*v1*w0*(phi1d(noff( 8)+ind0   )-phi1d(noff(15)+ind0))
     &                + u1*v1*w0*(phi1d(noff( 9)+ind0   )-phi1d(noff(16)+ind0))
     &                + u0*v0*w1*(phi1d(noff(17)+ind0+oy)-phi1d(noff(24)+ind0))
     &                + u1*v0*w1*(phi1d(noff(18)+ind0+oy)-phi1d(noff(25)+ind0))
     &                + u0*v1*w1*(phi1d(noff(20)+ind0   )-phi1d(noff(27)+ind0))
     &                + u1*v1*w1*(phi1d(noff(21)+ind0   )-phi1d(noff(28)+ind0)))

           ez(ip)=tdzi*(u0*v0*w0*(phi1d(noff( 1)+ind0) - phi1d(noff(20)+ind0))
     &                + u1*v0*w0*(phi1d(noff( 2)+ind0) - phi1d(noff(21)+ind0))
     &                + u0*v1*w0*(phi1d(noff( 3)+ind0) - phi1d(noff(24)+ind0))
     &                + u1*v1*w0*(phi1d(noff( 4)+ind0) - phi1d(noff(25)+ind0))
     &                + u0*v0*w1*(phi1d(noff( 8)+ind0) - phi1d(noff(29)+ind0))
     &                + u1*v0*w1*(phi1d(noff( 9)+ind0) - phi1d(noff(30)+ind0))
     &                + u0*v1*w1*(phi1d(noff(12)+ind0) - phi1d(noff(31)+ind0))
     &                + u1*v1*w1*(phi1d(noff(13)+ind0) - phi1d(noff(32)+ind0)))
        enddo

      endif


      elseif (efetch == 2) then

c     --- Allocate the temporary arrays
c     --- The arrays are made allocatable so that memory is only created for them when
c     --- efetch is 2.
        allocate(indx(33,np),p(33,np),u1p(np),v1p(np),w1p(np),stat=allocerror)
        if (allocerror /= 0) then
          print*,"sete3d: allocation error ",allocerror,
     &           ": could not allocate temp arrays to shape ",np
          stop
      endif


C  WARNING: This coding does not include symmetries!
c  NOTE: This section can probably be deleted since it is not really neccesary.
c  Vectorized loop over particles setting up indices of cells from which
c  to fetch phi, and weights for evaluating E at particle positions 
      if (l2symtry .or. l4symtry) then
        call remark("WARNING:sete3d - the algorithm does not take into account")
        call remark("               symmetries. Set efetch=1 or 3 and restart.")
        call kaboom(0)
      endif

      do ip = 1, np

         i =  (xp(ip) - xmmin) * dxi
         j =  (yp(ip) - ymmin) * dyi
         k =  (zp(ip) - zgrid - zmmin) * dzi
         ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)

         indx( 1,ip) = ind0 + noff(1)
         indx( 2,ip) = ind0 + noff(2)
         indx( 3,ip) = ind0 + noff(3)
         indx( 4,ip) = ind0 + noff(4)
         indx( 5,ip) = ind0 + noff(5)
         indx( 6,ip) = ind0 + noff(6)
         indx( 7,ip) = ind0 + noff(7)
         indx( 8,ip) = ind0 + noff(8)
         indx( 9,ip) = ind0 + noff(9)
         indx(10,ip) = ind0 + noff(10)
         indx(11,ip) = ind0 + noff(11)
         indx(12,ip) = ind0 + noff(12)
         indx(13,ip) = ind0 + noff(13)
         indx(14,ip) = ind0 + noff(14)
         indx(15,ip) = ind0 + noff(15)
         indx(16,ip) = ind0 + noff(16)
         indx(17,ip) = ind0 + noff(17)
         indx(18,ip) = ind0 + noff(18)
         indx(19,ip) = ind0 + noff(19)
         indx(20,ip) = ind0 + noff(20)
         indx(21,ip) = ind0 + noff(21)
         indx(22,ip) = ind0 + noff(22)
         indx(23,ip) = ind0 + noff(23)
         indx(24,ip) = ind0 + noff(24)
         indx(25,ip) = ind0 + noff(25)
         indx(26,ip) = ind0 + noff(26)
         indx(27,ip) = ind0 + noff(27)
         indx(28,ip) = ind0 + noff(28)
         indx(29,ip) = ind0 + noff(29)
         indx(30,ip) = ind0 + noff(30)
         indx(31,ip) = ind0 + noff(31)
         indx(32,ip) = ind0 + noff(32)

         u1p(ip) = (xp(ip) - xmmin) * dxi - i
         v1p(ip) = (yp(ip) - ymmin) * dyi - j
         w1p(ip) = (zp(ip) - zgrid - zmmin) * dzi - k

      enddo

c  Vectorized fetch of phi over the 32 cells touched;

      do ip = 1, np
         do m = 1, 32
            p(m,ip) = phi1d(indx(m,ip))
         enddo
      enddo

c  Evaluation of E, vectorized over particles

      do ip = 1, np

         u1 = u1p(ip)
         v1 = v1p(ip)
         w1 = w1p(ip)
         u0 = 1. - u1
         v0 = 1. - v1
         w0 = 1. - w1

         ex(ip) = tdxi * (u0 * v0 * w0 * (p( 7,ip) - p( 9,ip))
     &                  + u1 * v0 * w0 * (p( 8,ip) - p(10,ip))
     &                  + u0 * v1 * w0 * (p(11,ip) - p(13,ip)) 
     &                  + u1 * v1 * w0 * (p(12,ip) - p(14,ip)) 
     &                  + u0 * v0 * w1 * (p(19,ip) - p(21,ip))
     &                  + u1 * v0 * w1 * (p(20,ip) - p(22,ip))
     &                  + u0 * v1 * w1 * (p(23,ip) - p(25,ip)) 
     &                  + u1 * v1 * w1 * (p(24,ip) - p(26,ip)))

         ey(ip) = tdyi * (u0 * v0 * w0 * (p( 5,ip) - p(12,ip)) 
     &                  + u1 * v0 * w0 * (p( 6,ip) - p(13,ip))
     &                  + u0 * v1 * w0 * (p( 8,ip) - p(15,ip)) 
     &                  + u1 * v1 * w0 * (p( 9,ip) - p(16,ip)) 
     &                  + u0 * v0 * w1 * (p(17,ip) - p(24,ip))
     &                  + u1 * v0 * w1 * (p(18,ip) - p(25,ip))
     &                  + u0 * v1 * w1 * (p(20,ip) - p(27,ip)) 
     &                  + u1 * v1 * w1 * (p(21,ip) - p(28,ip)))

         ez(ip) = tdzi * (u0 * v0 * w0 * (p( 1,ip) - p(20,ip)) 
     &                  + u1 * v0 * w0 * (p( 2,ip) - p(21,ip))
     &                  + u0 * v1 * w0 * (p( 3,ip) - p(24,ip)) 
     &                  + u1 * v1 * w0 * (p( 4,ip) - p(25,ip)) 
     &                  + u0 * v0 * w1 * (p( 8,ip) - p(29,ip))
     &                  + u1 * v0 * w1 * (p( 9,ip) - p(30,ip))
     &                  + u0 * v1 * w1 * (p(12,ip) - p(31,ip)) 
     &                  + u1 * v1 * w1 * (p(13,ip) - p(32,ip)))

      enddo

c       --- Free the memory
        deallocate(indx)
        deallocate(p, u1p, v1p, w1p)

      elseif (efetch == 3) then

        if (.not. (l2symtry .or. l4symtry)) then
          do ip = 1, np

            i = (xp(ip) - xmmin) * dxi
            j = (yp(ip) - ymmin) * dyi
            k = (zp(ip) - zgrid - zmmin) * dzi

            u1 = (xp(ip) - xmmin) * dxi - i
            v1 = (yp(ip) - ymmin) * dyi - j
            w1 = (zp(ip) - zgrid - zmmin) * dzi - k

            u0 = 1. - u1
            v0 = 1. - v1
            w0 = 1. - w1

            ex(ip) = u0*v0*w0*selfe(1,i  ,j  ,k  )
     &             + u1*v0*w0*selfe(1,i+1,j  ,k  )
     &             + u0*v1*w0*selfe(1,i  ,j+1,k  )
     &             + u1*v1*w0*selfe(1,i+1,j+1,k  )
     &             + u0*v0*w1*selfe(1,i  ,j  ,k+1)
     &             + u1*v0*w1*selfe(1,i+1,j  ,k+1)
     &             + u0*v1*w1*selfe(1,i  ,j+1,k+1)
     &             + u1*v1*w1*selfe(1,i+1,j+1,k+1)

            ey(ip) = u0*v0*w0*selfe(2,i  ,j  ,k  )
     &             + u1*v0*w0*selfe(2,i+1,j  ,k  )
     &             + u0*v1*w0*selfe(2,i  ,j+1,k  )
     &             + u1*v1*w0*selfe(2,i+1,j+1,k  )
     &             + u0*v0*w1*selfe(2,i  ,j  ,k+1)
     &             + u1*v0*w1*selfe(2,i+1,j  ,k+1)
     &             + u0*v1*w1*selfe(2,i  ,j+1,k+1)
     &             + u1*v1*w1*selfe(2,i+1,j+1,k+1)

            ez(ip) = u0*v0*w0*selfe(3,i  ,j  ,k  )
     &             + u1*v0*w0*selfe(3,i+1,j  ,k  )
     &             + u0*v1*w0*selfe(3,i  ,j+1,k  )
     &             + u1*v1*w0*selfe(3,i+1,j+1,k  )
     &             + u0*v0*w1*selfe(3,i  ,j  ,k+1)
     &             + u1*v0*w1*selfe(3,i+1,j  ,k+1)
     &             + u0*v1*w1*selfe(3,i  ,j+1,k+1)
     &             + u1*v1*w1*selfe(3,i+1,j+1,k+1)

          enddo

        else

c         --- Set the signs of the E field for particles on negative side of
c         --- the axis of symmetry.
          sy = -1.
          sx = 1.
          if (l4symtry) sx = -1.

c         --- special loop symmetry is used
          do ip = 1, np

            i = (abs(xp(ip)) - xmmin)*dxi
            j = (abs(yp(ip)) - ymmin)*dyi
            k = (zp(ip) - zgrid - zmmin)*dzi

            u1 = (abs(xp(ip)) - xmmin)*dxi - i
            v1 = (abs(yp(ip)) - ymmin)*dyi - j
            w1 = (zp(ip) - zgrid - zmmin)*dzi - k

            u0 = 1. - u1
            v0 = 1. - v1
            w0 = 1. - w1

c           --- Adjust sign of E field for appropriate quadrant.
            xsign = +1.
            ysign = +1.
            if (xp(ip) < 0.) xsign = sx
            if (yp(ip) < 0.) ysign = sy

            ex(ip) = xsign*(u0*v0*w0*selfe(1,i  ,j  ,k  )
     &                    + u1*v0*w0*selfe(1,i+1,j  ,k  )
     &                    + u0*v1*w0*selfe(1,i  ,j+1,k  )
     &                    + u1*v1*w0*selfe(1,i+1,j+1,k  )
     &                    + u0*v0*w1*selfe(1,i  ,j  ,k+1)
     &                    + u1*v0*w1*selfe(1,i+1,j  ,k+1)
     &                    + u0*v1*w1*selfe(1,i  ,j+1,k+1)
     &                    + u1*v1*w1*selfe(1,i+1,j+1,k+1))

            ey(ip) = ysign*(u0*v0*w0*selfe(2,i  ,j  ,k  )
     &                    + u1*v0*w0*selfe(2,i+1,j  ,k  )
     &                    + u0*v1*w0*selfe(2,i  ,j+1,k  )
     &                    + u1*v1*w0*selfe(2,i+1,j+1,k  )
     &                    + u0*v0*w1*selfe(2,i  ,j  ,k+1)
     &                    + u1*v0*w1*selfe(2,i+1,j  ,k+1)
     &                    + u0*v1*w1*selfe(2,i  ,j+1,k+1)
     &                    + u1*v1*w1*selfe(2,i+1,j+1,k+1))

            ez(ip) =        u0*v0*w0*selfe(3,i  ,j  ,k  )
     &                    + u1*v0*w0*selfe(3,i+1,j  ,k  )
     &                    + u0*v1*w0*selfe(3,i  ,j+1,k  )
     &                    + u1*v1*w0*selfe(3,i+1,j+1,k  )
     &                    + u0*v0*w1*selfe(3,i  ,j  ,k+1)
     &                    + u1*v0*w1*selfe(3,i+1,j  ,k+1)
     &                    + u0*v1*w1*selfe(3,i  ,j+1,k+1)
     &                    + u1*v1*w1*selfe(3,i+1,j+1,k+1)

          enddo

        endif

      elseif (efetch == 4) then
c      --- Energy conserving

        if (.not. (l2symtry .or. l4symtry)) then
          do ip = 1, np

             i =  (xp(ip) - xmmin) * dxi
             j =  (yp(ip) - ymmin) * dyi
             k =  (zp(ip) - zgrid - zmmin) * dzi
             ind0 = i + j*nnx + k*nnxy

             u1 = (xp(ip) - xmmin) * dxi - i
             v1 = (yp(ip) - ymmin) * dyi - j
             w1 = (zp(ip) - zgrid - zmmin) * dzi - k

             u0 = 1. - u1
             v0 = 1. - v1
             w0 = 1. - w1

             ex(ip)=dxi*(+v0*w0*phi1d(noff( 8)+ind0)
     &                   -v0*w0*phi1d(noff( 9)+ind0)
     &                   +v1*w0*phi1d(noff(12)+ind0)
     &                   -v1*w0*phi1d(noff(13)+ind0)
     &                   +v0*w1*phi1d(noff(20)+ind0)
     &                   -v0*w1*phi1d(noff(21)+ind0)
     &                   +v1*w1*phi1d(noff(24)+ind0)
     &                   -v1*w1*phi1d(noff(25)+ind0))

             ey(ip)=dyi*(+u0*w0*phi1d(noff( 8)+ind0)
     &                   +u1*w0*phi1d(noff( 9)+ind0)
     &                   -u0*w0*phi1d(noff(12)+ind0)
     &                   -u1*w0*phi1d(noff(13)+ind0)
     &                   +u0*w1*phi1d(noff(20)+ind0)
     &                   +u1*w1*phi1d(noff(21)+ind0)
     &                   -u0*w1*phi1d(noff(24)+ind0)
     &                   -u1*w1*phi1d(noff(25)+ind0))

             ez(ip)=dzi*(+u0*v0*phi1d(noff( 8)+ind0)
     &                   +u1*v0*phi1d(noff( 9)+ind0)
     &                   +u0*v1*phi1d(noff(12)+ind0)
     &                   +u1*v1*phi1d(noff(13)+ind0)
     &                   -u0*v0*phi1d(noff(20)+ind0)
     &                   -u1*v0*phi1d(noff(21)+ind0)
     &                   -u0*v1*phi1d(noff(24)+ind0)
     &                   -u1*v1*phi1d(noff(25)+ind0))

          enddo

        else

c         --- Set the signs of the E field for particles on negative side of
c         --- the axis of symmetry.
          sy = -1.
          sx = 1.
          if (l4symtry) sx = -1.

c         --- special loop symmetry is used
          do ip = 1, np
  
             i =  (abs(xp(ip)) - xmmin)*dxi
             j =  (abs(yp(ip)) - ymmin)*dyi
             k =  (zp(ip) - zgrid - zmmin)*dzi
             ind0 = i + j*nnx + k*nnxy

             u1 = (abs(xp(ip)) - xmmin)*dxi - i
             v1 = (abs(yp(ip)) - ymmin)*dyi - j
             w1 = (zp(ip) - zgrid - zmmin)*dzi - k

             u0 = 1. - u1
             v0 = 1. - v1
             w0 = 1. - w1

c            --- Adjust sign of E field for appropriate quadrant.
             xsign = dxi
             ysign = dyi
             if (xp(ip) < 0.) xsign = sx*dxi
             if (yp(ip) < 0.) ysign = sy*dyi

             ex(ip)=xsign*(+v0*w0*phi1d(noff( 8)+ind0)
     &                     -v0*w0*phi1d(noff( 9)+ind0)
     &                     +v1*w0*phi1d(noff(12)+ind0)
     &                     -v1*w0*phi1d(noff(13)+ind0)
     &                     +v0*w1*phi1d(noff(20)+ind0)
     &                     -v0*w1*phi1d(noff(21)+ind0)
     &                     +v1*w1*phi1d(noff(24)+ind0)
     &                     -v1*w1*phi1d(noff(25)+ind0))

             ey(ip)=ysign*(+u0*w0*phi1d(noff( 8)+ind0)
     &                     +u1*w0*phi1d(noff( 9)+ind0)
     &                     -u0*w0*phi1d(noff(12)+ind0)
     &                     -u1*w0*phi1d(noff(13)+ind0)
     &                     +u0*w1*phi1d(noff(20)+ind0)
     &                     +u1*w1*phi1d(noff(21)+ind0)
     &                     -u0*w1*phi1d(noff(24)+ind0)
     &                     -u1*w1*phi1d(noff(25)+ind0))

             ez(ip) = dzi*(+u0*v0*phi1d(noff( 8)+ind0)
     &                     +u1*v0*phi1d(noff( 9)+ind0)
     &                     +u0*v1*phi1d(noff(12)+ind0)
     &                     +u1*v1*phi1d(noff(13)+ind0)
     &                     -u0*v0*phi1d(noff(20)+ind0)
     &                     -u1*v0*phi1d(noff(21)+ind0)
     &                     -u0*v1*phi1d(noff(24)+ind0)
     &                     -u1*v1*phi1d(noff(25)+ind0))

          enddo

        endif

      endif

!$OMP MASTER
      if (lw3dtimesubs) timesete3d = timesete3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=========================================================================== 
      subroutine getselfe3d(phi,nx,ny,nz,selfe,nx_selfe,ny_selfe,nz_selfe,
     &                      dx,dy,dz,boundx0,boundxn,boundy0,boundyn)
      use GlobalVars
      use Subtimers3d
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: nx_selfe,ny_selfe,nz_selfe
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),selfe(3,0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz
      integer(ISZ):: boundx0,boundxn,boundy0,boundyn

c Calculate the self-E via finite differences of phi.

      integer(ISZ):: ix,iy,iz
      real(kind=8):: tdxi,tdyi,tdzi
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Do some checking first.
      if (nx /= nx_selfe .or. ny /= ny_selfe .or. nz /= nz_selfe) then
        call remark("ERROR: the dimensions of phi and selfe must be the same")
        call kaboom(1)
      endif

      tdxi = 0.5/dx
      tdyi = 0.5/dy
      tdzi = 0.5/dz

c     --- Do the calculation for the bulk of the array
      do iz=0,nz
        do iy=1,ny-1
          do ix=1,nx-1
            selfe(1,ix,iy,iz) = tdxi*(phi(ix-1,iy  ,iz  ) - phi(ix+1,iy  ,iz  ))
            selfe(2,ix,iy,iz) = tdyi*(phi(ix  ,iy-1,iz  ) - phi(ix  ,iy+1,iz  ))
            selfe(3,ix,iy,iz) = tdzi*(phi(ix  ,iy  ,iz-1) - phi(ix  ,iy  ,iz+1))
          enddo
        enddo
      enddo

c     --- E on boundaries. These are only really needed when symmetry is used.
c     --- It is done always anyway since it is little work and also forces the
c     --- value at the boundaries.
c     --- Note that the z boundaries at taken care of above, taking advantage
c     --- of the guard cells.
      selfe(1,1:nx-1, 0,:) = tdxi*(phi(0:nx-2, 0,:) - phi(2:nx, 0,:))
      selfe(1,1:nx-1,ny,:) = tdxi*(phi(0:nx-2,ny,:) - phi(2:nx,ny,:))
      selfe(2, 0,1:ny-1,:) = tdyi*(phi( 0,0:ny-2,:) - phi( 0,2:ny,:))
      selfe(2,nx,1:ny-1,:) = tdyi*(phi(nx,0:ny-2,:) - phi(nx,2:ny,:))

      if (boundx0 == dirichlet .or. boundx0 == neumann) then
        selfe(1, 0, :,:) = 0.
      else if (boundx0 == periodic) then
        selfe(1, 0, :,:) = tdxi*(phi(nx-1,:,:) - phi(1,:,:))
      endif

      if (boundxn == dirichlet .or. boundxn == neumann) then
        selfe(1,nx, :,:) = 0.
      else if (boundxn == periodic) then
        selfe(1,nx, :,:) = tdxi*(phi(nx-1,:,:) - phi(1,:,:))
      endif

      if (boundy0 == dirichlet .or. boundy0 == neumann) then
        selfe(2, :, 0,:) = 0.
      else if (boundy0 == periodic) then
        selfe(2, :, 0,:) = tdyi*(phi(:,ny-1,:) - phi(:,1,:))
      endif

      if (boundyn == dirichlet .or. boundyn == neumann) then
        selfe(2, :,ny,:) = 0.
      else if (boundyn == periodic) then
        selfe(2, :,ny,:) = tdyi*(phi(:,ny-1,:) - phi(:,1,:))
      endif

c     --- Axial self field does not depend on the transverse boundary condition
      selfe(3, 0,:,0:nz) = tdzi*(phi( 0,:,-1:nz-1) - phi( 0,:,1:nz+1))
      selfe(3,nx,:,0:nz) = tdzi*(phi(nx,:,-1:nz-1) - phi(nx,:,1:nz+1))
      selfe(3,1:nx-1, 0,0:nz) = tdzi*(phi(1:nx-1, 0,-1:nz-1) -
     &                                phi(1:nx-1, 0,1:nz+1))
      selfe(3,1:nx-1,ny,0:nz) = tdzi*(phi(1:nx-1,ny,-1:nz-1) -
     &                                phi(1:nx-1,ny,1:nz+1))

!$OMP MASTER
      if (lw3dtimesubs) timegetselfe3d = timegetselfe3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine stptcl3d
      use GlobalVars
      use Subtimers3d
      use Beam_acc
      use Constant
      use InGen
      use InGen3d
      use InPart
      use InPart3d
      use InMesh3d
      use Particles
      use Picglb
      use Picglb3d
#ifdef MPIPARALLEL
      use Parallel
#endif

c  Loads particles.
c  At end, particle boundary conditions are enforced.
c
c  For the parallel version, each processor loads the fraction of the
c  longitudinally uniform beam of length (zimax-zimin) which is within its
c  region. The routine sw_globalsum recalculates the particle weights (sw)
c  appropriately.
c
c  Picglb is included only to use zbeam


      integer(ISZ):: npm,ip,is,i,j,k,i1,nxstri_s,nystri_s,nzstri_s,izstripe,icheck
      integer(ISZ):: nfib2,nfib3,nfib4,ii,ioff
      real(kind=8):: rnpmi,zmid,zlen,zleni,vtx,vty,vtz,rr,r0,rpp,phi1,phi2,rp,wz,wr
      real(kind=8):: vthr,vrbar
      real(kind=8):: zoff
      real(kind=8):: sumzdist,zdistmax
      real(kind=8):: sumrdist
      integer(ISZ):: iz,ipp,izz
      integer(ISZ):: ir,irr
      integer(ISZ):: envxport
      character(72):: errline
      real(kind=8):: wranf,wrandom,wrandomgauss

      logical(ISZ):: firstk
      integer(ISZ):: npgrp,iptotal
      integer(ISZ),allocatable,dimension(:):: indx
      real(kind=8),allocatable,dimension(:):: xt,yt,zt,rt,tt,uxt,uyt,uzt
      real(kind=8),allocatable,dimension(:):: perpscal,at,apt,bt,bpt
      real(kind=8),allocatable,dimension(:):: epsxt,epsyt,vzt,vtzt,ibeamt
      real(kind=8),allocatable,dimension(:):: xct,xpct,yct,ypct
      integer(ISZ):: allocerror

#ifdef MPIPARALLEL
      real(kind=8):: zmin,zmax
#endif
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  If distrbtn is not set by the user, exit
      if (distrbtn == "none") then
!$OMP MASTER
        if (lw3dtimesubs) timestptcl3d = timestptcl3d + wtime() - substarttime
!$OMP END MASTER
        return
      endif

c  Set npgrp so arrays are correct size
      if (xrandom == "fibonacc" .or.
     &    xrandom == "digitrev" .or.
     &    xrandom == "pseudo") then
        npgrp = nparpgrp
      elseif (xrandom == "grid") then
        npgrp = nxstripe*nystripe
      else
        write (errline,'("ERROR: xrandom has an improper value = ",a8)')
     &         xrandom
        call remark (errline)
        call kaboom (0)
      endif

c  Set npm to loop through all possible particles
      if (xrandom == "fibonacc") then
        npm = nfibgrps*fibg1
      elseif (xrandom == "digitrev") then
        npm = npmax
      elseif (xrandom == "pseudo") then
        npm = npmax
      elseif (xrandom == "grid") then
        npm = nxstripe*nystripe*nzstripe
      endif

c  If not loading any particles, skip the next big section, but still calculate
c  values for sq and sm. For parallel version, also update the value for sw.
      if (npm > 0) then

c  Prepare for arbitrary particle distribution in z.  This coding will only
c  load a beam approximately matching the distribution.  A gathering of
c  the particles on a z mesh will not exactly reproduce the z distribution.
c  This problem is being studied.
      if (nzdist > 0) then
c       --- Find sum of zdist so it can be normalized
        sumzdist = 0.5*zdist(0)
        do iz=1,nzdist-1
          sumzdist = sumzdist + zdist(iz)
        enddo
        sumzdist = sumzdist + 0.5*zdist(nzdist)
c       --- Calculate integral of zdist from 0 to iz
c       --- Integral = 0.5*zdist(0) + sum(zdist(1:iz-1)) + 0.5*zdist(iz)
c       --- Also, normalize zdist.
        call gchange("InPart3d",0)
        nrmzdist(0) = zdist(0)/sumzdist
        intzdist(0) = 0.
        do iz=1,nzdist
          nrmzdist(iz) = zdist(iz)/sumzdist
          intzdist(iz) = intzdist(iz-1) + 0.5*(nrmzdist(iz-1) + nrmzdist(iz))
        enddo
      endif

c  Prepare for arbitrary particle distribution in r.  This coding will only
c  load a beam approximately matching the distribution.  A gathering of
c  the particles on a r mesh will not exactly reproduce the r distribution.
c  This problem is being studied.
      if (nrdist > 0) then
        call gchange("InPart3d",0)
c       --- Multiply rdist by r.
        do ir=0,nrdist
          nrmrdist(ir) = rdist(ir)*ir/nrdist
        enddo
c       --- Integrate rdist so it can be normalized. Note that nrmrdist(0)
c       --- is ignored since it is always zero.
        sumrdist = 0.
        do ir=1,nrdist-1
          sumrdist = sumrdist + nrmrdist(ir)
        enddo
        sumrdist = sumrdist + 0.5*nrmrdist(nrdist)
c       --- Calculate integral of rdist from 0 to ir
c       --- Integral = 0.5*rdist(0) + sum(rdist(1:ir-1)) + 0.5*rdist(ir)
c       --- Also, normalize rdist.
        intrdist(0) = 0.
        do ir=1,nrdist
          nrmrdist(ir) = nrmrdist(ir)/sumrdist
          intrdist(ir) = intrdist(ir-1) + 0.5*(nrmrdist(ir-1) + nrmrdist(ir))
        enddo
      endif

c     Calculate rms equivalent matched beam equilibria in a continuous 
c     focusing channel for use in generalized psudo-equilibrium distribution 
c     loads.  
c     
c     As presently implemented, this only calculates the equilibrium 
c     once for a single beam slice, making it useful only for xy slice code 
c     simulations.  To get this to work in 3D, the following should be done:
c       1) Modify the main particle loading loop to first load the particle 
c          z coordinates of the full beam distribution before any transverse 
c          loading.  
c       2) Using the loaded z coordinates and the desired transverse loading 
c          option, calculate the range max -> min of all continuous focusing 
c          eqivalent beam parameters for all particles 
c       3) Based on the ranges of parameters in 2) calculate needed equilibria 
c          for the transverse load by discretizing the parameter range and 
c          then calculating the equilibria in scaled coordinates over the 
c          discretized range.  
c       4) Go through the main particle loading loop again to load the 
c          transverse distribution using the equilibria calculated in 3).  
c          Interpolation can be used for equilibria not saved at exactly the 
c          parameter values needed.


      if (distrbtn == "WB" .or. distrbtn == "Waterbag"  .or.  
     &    distrbtn == "PA" .or. distrbtn == "Parabolic" .or. 
     &    distrbtn == "TE" .or. distrbtn == "ThermalEquilibrium") then 

        call  perp_cfe_den()

      endif 

c----------------------------------
c  Begin main loop over species
c----------------------------------

!$OMP PARALLEL
!$OMP+PRIVATE(k,i,j,ipp,i1,izz,nzstri_s,izstripe,vtz,iz,wz,irr,ii,rp,rpp,r0)
!$OMP+PRIVATE(firstk)
!$OMP+PRIVATE(ir,wr,vthr,vrbar,vtx,vty,ip,rr,phi1,phi2)
!$OMP+PRIVATE(indx,xt,yt,zt,rt,tt,uxt,uyt,uzt,perpscal,at,apt,bt,bpt)
!$OMP+PRIVATE(xct,xpct,yct,ypct,icheck,nfib2,nfib3,nfib4,nxstri_s,nystri_s)
!$OMP+PRIVATE(epsxt,epsyt,vzt,vtzt,ibeamt)
!$OMP+SHARED(iptotal)

c  Allocate scratch space for load. Must be done in the parallel block
      allocate(indx(npgrp),xt(npgrp),yt(npgrp),zt(npgrp),rt(npgrp),tt(npgrp),
     &         uxt(npgrp),uyt(npgrp),uzt(npgrp),perpscal(npgrp),
     &         at(npgrp),apt(npgrp),bt(npgrp),bpt(npgrp),
     &         epsxt(npgrp),epsyt(npgrp),vzt(npgrp),vtzt(npgrp),ibeamt(npgrp),
     &         xct(npgrp),xpct(npgrp),yct(npgrp),ypct(npgrp),stat=allocerror)
      if (allocerror /= 0) then
        print*,"stptcl3d: allocation error ",allocerror,
     &         ": could not allocate temp arrays to shape ",npgrp
        stop
      endif

      do is = 1, ns

c       --- Get number of particles to be loaded for this species.
        if (xrandom == "fibonacc") then
          npm = nfibgrps*fibg1*sp_fract(is)
        elseif (xrandom == "digitrev") then
          npm = np_s(is)
        elseif (xrandom == "pseudo") then
          npm = np_s(is)
        elseif (xrandom == "grid") then
          npm = nxstripe*nystripe*nzstripe*sp_fract(is)
        endif

        zmid = .5 * (zimax_s(is) + zimin_s(is)) + zbeam
        zlen = zimax_s(is) - zimin_s(is)
        zleni = 1./zlen
        ip = ins(is) - 1
!$OMP SINGLE
#ifdef _OPENMP
        iptotal = ins(is) - 1
#endif
!$OMP END SINGLE

c       --- Inverse of the number of particles. When multiplied by the particle
c       --- index and added to the normalized value of zmin, the normalized 
c       --- axial location of the particle in the beam, between 0 and 1, 
c       --- results.
        rnpmi = 1./dvnz(real(npm,kind=8))

c       --- Reset counter for arbitrary particle distribution in z.
c       --- It is reset for each species since each species has its own
c       --- initial uniform distribution.
c       --- Resetting izz here and starting each particle with its previous
c       --- value only works since the initial uniform z positions are
c       --- monotonically increasing.  Using this feature greatly increases
c       --- the speed of the code.
        if (nzdist > 0) izz = 0

c       --- Calculate zoff, the axial offset of particles for processors, and
c       --- ioff, the offset for the random number generators for processors.
c       --- This ensures that each processor has the correct z distribution and
c       --- seeds for the generators.
        zoff = 0.
        ioff = 0

#ifdef MPIPARALLEL
c       --- first, calculate z extent of beam within this processor
        zmin = max(zimin_s(is),zpslmin(my_index))
        zmax = min(zimax_s(is),zpslmax(my_index))
c       --- calculate number of particles that are loaded by processors
c       --- to the left
        ioff = int(npm*(zmin - zimin_s(is))*zleni + 1.e-5)
c       --- calculate axial offset so the first particle for this processor
c       --- is the correct distance from the last particle of the processor
c       --- to the left.
        zoff = ioff*rnpmi
c       --- number of particles to be loaded in this processor: the number
c       --- of particles between zmax and zimin minus the number between
c       --- zmin and zimin.
        npm = int(npm*(zmax - zimin_s(is))*zleni + 1.e-5) - ioff
c       --- Make sure that there is enough room to load the particles.
c       --- This is only needed here in the parallel version since it is
c       --- difficult to gaurantee that there is room when an arbitrary
c       --- domain decomposition is allowed.
!$OMP MASTER
        call chckpart(is,0,npm,.false.)
!$OMP END MASTER
#endif

c       --- Shift ioff by randoffset so that the seed for the quiet start
c       --- random number generators can be controlled.
        ioff = ioff + randoffset

c       --- Set firstk = true. This is needed so that the code knows that it
c       --- is in the first iteration of the loop over k. It is for the
c       --- OpenMP code since there, each processor won't execute the
c       --- iteration where k=1 and so some things would otherwise not be
c       --- initialized properly, for example the code for a cylinder beam.
        firstk = .true.

c  MAIN INITIALIZATION LOOP

!$OMP DO
      do k = 1, npm, npgrp
         i1 = min(npm, k+npgrp-1)

c  Load longitudinal stuff
c        --- load normalized z (0 < z < 1)
         if (xrandom == "fibonacc" .or.
     &       xrandom == "digitrev" .or.
     &       xrandom == "pseudo") then
           do i=k,i1
              zt(i-k+1) = (i-.5)*rnpmi + zoff
           enddo
c          --- Transform the longitudinal positions to match the specified
c          --- longitudinal distribution.
           if (nzdist > 0) then
             do ipp=k,i1
c              --- Find what grid cell the particle is in.
c              --- Note that the previous value of izz is used (see note above)
               do while (zt(ipp-k+1) > intzdist(izz+1))
                 izz = izz + 1
               enddo
c              --- Do the transformation.  The transformation assumes a
c              --- linearly varying distribution between grid locations.
c              --- The transformation is derived by equating an integral
c              --- over the uniform distribution to an integral over the
c              --- specified distribution.  The particle remains within
c              --- the range [0,1].
               if (nrmzdist(izz+1) == nrmzdist(izz)) then
                 zt(ipp-k+1) = (zt(ipp-k+1) - intzdist(izz))/nrmzdist(izz)
               else
                 zt(ipp-k+1) = ((sqrt(nrmzdist(izz)**2 -
     &                     2.*(nrmzdist(izz+1)-nrmzdist(izz))*(intzdist(izz) -
     &                     zt(ipp-k+1)))-nrmzdist(izz))/
     &                     (nrmzdist(izz+1)-nrmzdist(izz)) + izz)/nzdist
               endif
             enddo
           endif

         elseif (xrandom == "grid") then
           nzstri_s = nzstripe*sp_fract(is)**(1./3.)
           izstripe = (k-1)/npgrp + 1
           do i=k,i1
              zt(i-k+1) = (izstripe-.5)/nzstri_s + zjig*(wranf()-.5)/nzstri_s +
     &                    zoff
           enddo
         endif

c        --- zero vz if vthz is zero only the first time through
         if (vthz_s(is) == 0. .and. firstk) then
           do i=1,npgrp
              uzt(i) = 0.
           enddo
         endif

         if (cigarld) then
c          --- set vtz to vthz*2 since vtz is vz_max - vz_bar for cigar
           vtz = 2.*vthz_s(is)
c          --- load uzt with linear distribution (-.5 < uzt < .5)
           if (vtz /= 0.) then
             do i=1,i1-k+1
               uzt(i) = wrandom(vzrandom,i+k-1+ioff,dig7,fibg1,npm) - .5
             enddo
           endif
c          --- use cigar to adjust the z's
c          --- xt,yt,uxt,uyt passed as scratch arrays
           call cigar(i1-k+1,zt(1),uzt(1),zt(1),uzt(1),perpscal(1),
     &                straight_s(is),xt(1),yt(1),uxt(1),uyt(1))

c          --- Loading gaussian longitudinal velocity distribution.
           if (distr_l == "gaussian") then
             do i=1,i1-k+1
               uzt(i) = wrandomgauss(vzrandom,i+k-1+ioff,dig7,dig8,1,1,.true.)
               uzt(i) = 0.5*perpscal(i)*uzt(i)
             enddo
           endif

         else
c          --- non-cigar load
c          --- set vtz to vthz
           vtz = vthz_s(is)
c          --- load uzt with gaussian distribution (mean 0, variance 1)
           if (vtz /= 0.) then
             do i=1,i1-k+1
               uzt(i) = wrandomgauss(vzrandom,i+k-1+ioff,dig7,dig8,1,1,.true.)
             enddo
           endif
c          --- set perpscal to one first time through
           if (firstk) perpscal = 1.0
         endif

c        --- unormalize z to fetch envelope
         zt = zbeam + zimin_s(is) + zlen*zt

c        --- fetch envelope, centroid, and emittance.  Reset values 
c        --- after initial set if axial varying emittance, axial velocity, 
c        --- and axial thermal velocity are input
         if (nenvofz > 0) then
           do i=1,i1-k+1
c            --- envelope interpolations 
             iz = (zt(i) - zimin_s(is) - zbeam)*zleni*nenvofz
             wz = (zt(i) - zimin_s(is) - zbeam)*zleni*nenvofz - iz
c            --- envelope 
             at(i)  = aofz(iz)*(1.-wz)  + aofz(iz+1)*wz
             bt(i)  = bofz(iz)*(1.-wz)  + bofz(iz+1)*wz
             apt(i) = apofz(iz)*(1.-wz) + apofz(iz+1)*wz
             bpt(i) = bpofz(iz)*(1.-wz) + bpofz(iz+1)*wz
c            --- centroid 
             xct(i)  = xofz(iz)*(1.-wz)  + xofz(iz+1)*wz
             yct(i)  = yofz(iz)*(1.-wz)  + yofz(iz+1)*wz
             xpct(i) = xpofz(iz)*(1.-wz) + xpofz(iz+1)*wz
             ypct(i) = ypofz(iz)*(1.-wz) + ypofz(iz+1)*wz
           enddo
         elseif (.not. cylinder) then
           icheck = envxport(i1-k+1,zt,at,apt,bt,bpt,xct,xpct,yct,ypct,
     &                       vzt,epsxt,epsyt,ibeamt)
           if (icheck == 1) then
              call remark ("ERROR: out-of-range z sent to ENVXPORT")
              call kaboom (0)
           endif
c          --- set axial thermal velocity based on species .... note that 
c          --- at present this is inconsistent with envxport not knowing about 
c          --- individual species.  Eventually this might be fixed with a 
c          --- generalized envelope solver 
           do i=1,npgrp
             vtzt(i) = vthz_s(is)
           enddo 
         else
c          --- set only first time through
           if (firstk) then
             do i=1,npgrp
               at(i)   = a0_s(is)
               bt(i)   = b0_s(is)
               apt(i)  = ap0_s(is)
               bpt(i)  = bp0_s(is)
c
               xct(i)  = x0_s(is)  + xcent_s(is)
               yct(i)  = y0_s(is)  + ycent_s(is)
               xpct(i) = xp0_s(is) + xpcent_s(is)
               ypct(i) = yp0_s(is) + ypcent_s(is)
c
               epsxt(i) = emitx_s(is) 
               epsyt(i) = emity_s(is) 
c
               vzt(i)  = vbeam_s(is) 
               vtzt(i) = vthz_s(is) 
c
               ibeamt(i) = ibeam_s(is)
             enddo
           endif
         endif
c        --- reset emittance interpolations               
         if (nemitofz > 0) then 
           do i=1,i1-k+1
             iz = (zt(i) - zimin_s(is) - zbeam)*zleni*nemitofz
             wz = (zt(i) - zimin_s(is) - zbeam)*zleni*nemitofz - iz
c
             epsxt(i) = emitxofz(iz)*(1.-wz) + emitxofz(iz+1)*wz
             epsyt(i) = emityofz(iz)*(1.-wz) + emityofz(iz+1)*wz
           enddo
         endif  
c        --- reset axial velocity interpolations 
         if (nvbeamofz > 0) then 
           do i=1,i1-k+1
c
             iz = (zt(i) - zimin_s(is) - zbeam)*zleni*nvbeamofz
             wz = (zt(i) - zimin_s(is) - zbeam)*zleni*nvbeamofz - iz
c
             vzt(i) = vbeamofz(iz)*(1.-wz) + vbeamofz(iz+1)*wz
           enddo 
         endif 
c        --- reset axial thermal velocity interpolations 
         if (nvthzofz > 0) then 
           do i=1,i1-k+1
c
             iz = (zt(i) - zimin_s(is) - zbeam)*zleni*nvthzofz
             wz = (zt(i) - zimin_s(is) - zbeam)*zleni*nvthzofz - iz
c
             vtzt(i) = vthzofz(iz)*(1.-wz) + vthzofz(iz+1)*wz
           enddo
         endif 
c        --- XXX beam current interpolations (placeholder)  
c         if () then 
c           do i=1,i1-k+1 
c             ibeamt(i) = .... 
c           enddo 
c         endif  


c  Load tranverse stuff
c  semi-gaussian distribution
         if (distrbtn == "semigaus" .or. distrbtn == "SG" .or. 
     &       distrbtn == "SemiGaussian") then

c          --- load normalized tranverse space variables
           if (xrandom .ne. "grid") then
             nfib2 = nfibgrps*fibg2
             nfib3 = nfibgrps*fibg3
             do i=1,i1-k+1
               xt(i) = 2.*wrandom(xrandom,i+k-1+ioff,dig1,nfib2,npm) - 1.
               yt(i) = 2.*wrandom(xrandom,i+k-1+ioff,dig2,nfib3,npm) - 1.
             enddo
           elseif (xrandom == "grid") then
             nxstri_s = nxstripe*sp_fract(is)**(1./3.)
             nystri_s = nystripe*sp_fract(is)**(1./3.)
             do j=1,nystri_s
               do i=1,nxstri_s
                 xt(k+i-1+(j-1)*nxstri_s) = 2.*(i-.5)/nxstri_s - 1.
                 yt(k+i-1+(j-1)*nystri_s) = 2.*(j-.5)/nystri_s - 1.
               enddo
             enddo
           endif

           if (ldprfile == "polar") then
c          --- use random numbers to load particles in polar coordinates
             do i=1,i1-k+1
               rt(i) = 0.5*(xt(i) + 1.)
               tt(i) = yt(i)*pi

c              --- Transform the radial positions to match the specified
c              --- distribution.
               if (nrdist > 0) then
c                --- Find what grid cell the particle is in.
                 irr = 0
                 do while (rt(i) > intrdist(irr+1))
                   irr = irr + 1
                 enddo
c                --- Do the transformation.  The transformation assumes a
c                --- linearly varying distribution between grid locations.
c                --- The transformation is derived by equating an integral
c                --- over the uniform distribution to an integral over the
c                --- specified distribution.  The particle remains within
c                --- the range [0,1].
                 if (nrmrdist(irr+1) == nrmrdist(irr)) then
                   rt(i) = (rt(i) - intrdist(irr))/nrmrdist(irr)
                 else
                   rt(i) = ((sqrt(nrmrdist(irr)**2 -
     &              2.*(nrmrdist(irr+1)-nrmrdist(irr))*(intrdist(irr) - rt(i)))-
     &              nrmrdist(irr))/(nrmrdist(irr+1)-nrmrdist(irr)) + irr)/nrdist
                 endif
c                --- The above gives the radius, but below, rt(i) is assumed to
c                --- hold the r squared.
                 rt(i) = rt(i)*rt(i)
               endif

c              --- Hollow beam using an analytic function.
               if (hollow == 2) then
c                --- f(r)~(h+(1-h)r^2)
                 rt(i) = (1 + hollow_h)*rt(i)/
     &                (hollow_h + sqrt(hollow_h**2 + (1. - hollow_h**2)*rt(i)))
               endif

c              --- Now convert radius and theta to x and y.
               rt(i) = sqrt(rt(i))
               xt(i) = rt(i)*cos(tt(i))
               yt(i) = rt(i)*sin(tt(i))
               indx(i) = i
             enddo
             j = i1-k+1

           elseif (ldprfile == "streamls") then
c            --- carve into cylinder
             j=0
             do i=1,i1-k+1
               if (xt(i)**2 + yt(i)**2 < 1.) then
                 rt(i) = sqrt(xt(i)**2 + yt(i)**2)
                 tt(i) = atan2(yt(i),xt(i))
                 j=j+1
                 indx(j) = i
               endif
             enddo

c            --- transform to hollow beam of type one
c            --- f(r^2)  =    f0 * (r/rmax)^2                r < rmax/2
c            ---              f0 * (1 - (r/rmax)^2)/3        r > rmax/2
c            --- x and y are multiplied by 1.08 to keep rbar and rrms roughly
c            --- the same as in the uniform beam
             if (hollow == 1)  then
               do i=1,i1-k+1
                 ii = indx(i)
                 rp = sqrt(xt(ii)**2 + yt(ii)**2)
                 rpp = sqrt(0.5*sqrt(rp**2))
                 if (rp > 0.5) rpp = sqrt(1. - 0.5*sqrt(3. - 3.*rp**2))
                 xt(ii) = rpp*xt(ii)/rp*1.08
                 yt(ii) = rpp*yt(ii)/rp*1.08
               enddo
             endif

           elseif (ldprfile == "stripes") then
c            --- carve into normalized envelope
             j=0
             r0 = max(a0_s(is),b0_s(is))
             do i=1,i1-k+1
               if ((xt(i)*bt(i))**2 + (yt(i)*at(i))**2
     &                                < (at(i)*bt(i)/r0)**2) then
                 j=j+1
                 indx(j) = i
                 xt(i) = xt(i)*r0/at(i)
                 yt(i) = yt(i)*r0/bt(i)
                 rt(i) = sqrt(xt(i)**2 + yt(i)**2)
                 tt(i) = atan2(yt(i),xt(i))
               endif
             enddo

           else
             write (errline,'("ERROR: ldprfile has an improper value = ",a8)')
     &           ldprfile
             call remark (errline)
             call kaboom (0)
           endif

c          --- load normalized transverse velocity variables
           if (distr_t == "gaussian") then
c             --- Note that the loops must be broken up in this way since the
c             --- routine wrandomgauss must be called in the correct order.
             do i=1,i1-k+1
               uxt(i) = wrandomgauss(vtrandom,i+k-1+ioff,dig3,dig4,1,1,.true.)
             enddo
             do i=1,i1-k+1
               uyt(i) = wrandomgauss(vtrandom,i+k-1+ioff,dig5,dig6,1,1,.true.)
             enddo
           else if (distr_t == "uniform") then
             do i=1,i1-k+1
               uxt(i) = 2.*wrandom(vtrandom,i+k-1+ioff,dig3,fibg1,npm) - 1.
               uyt(i) = 2.*wrandom(vtrandom,i+k-1+ioff,dig5,fibg1,npm) - 1.
             enddo
           endif

c          --- Modify radial velocity distribution
           if (nvrdist > 0) then
             do i=1,i1-k+1
               ir = rt(i)*nvrdist
               wr = rt(i)*nvrdist - ir
               vthr = (vthrofr(ir)*(1.-wr) + vthrofr(ir+1)*wr)
               vrbar = (vrbarofr(ir)*(1.-wr) + vrbarofr(ir+1)*wr)
               uxt(i) = vthr*uxt(i) + vrbar*cos(tt(i))
               uyt(i) = vthr*uyt(i) + vrbar*sin(tt(i))
             enddo
           endif

!$OMP CRITICAL (STPTCL3D1)
#ifdef _OPENMP
           ip = iptotal
           iptotal = iptotal + j
#endif
!$OMP END CRITICAL (STPTCL3D1)

c          --- unnormalize everything and load into particle arrays
           do i=1,j
             ii = indx(i)
c            --- put position into particle arrays
             xp(ip+i) = at(ii)*xt(ii)*perpscal(ii) + xct(ii)
             yp(ip+i) = bt(ii)*yt(ii)*perpscal(ii) + yct(ii)
             zp(ip+i) = zt(ii)
c            --- put velocity into particle arrays
             vtx = .5*vzt(ii)*epsxt(ii)/at(ii)*perpscal(ii) + vthperp_s(is)
             vty = .5*vzt(ii)*epsyt(ii)/bt(ii)*perpscal(ii) + vthperp_s(is)
             uxp(ip+i) = vzt(ii)*apt(ii)*xt(ii)*perpscal(ii) + vtx*uxt(ii) +
     &                   vzt(ii)*xpct(ii)
             uyp(ip+i) = vzt(ii)*bpt(ii)*yt(ii)*perpscal(ii) + vty*uyt(ii) +
     &                   vzt(ii)*ypct(ii)
             uzp(ip+i) = vzt(ii)*(1. + vtilt_s(is)*(zmid - zt(ii))*zleni) +
     &                   vtzt(ii)*uzt(ii)
             if (npid > 0) pid(ip+i,:) = 0.
           enddo

c          --- add sinusoidal perturbation to uzp
           if (vzperamp /= 0.) then
             do i=ip+1,ip+j
               uzp(i) = uzp(i) + vzperamp*sin(2.*pi*zp(i)/vzperlam + vzperphs)
             enddo
           endif

c          --- increment number of particles by size of current group
           ip = ip + j

c        --- K-V distribution
         elseif (distrbtn == "K-V" .or. distrbtn == "KV") then

c          --- fetch random numbers and put into xt,yt, & uxt  temporarily
           nfib2 = nfibgrps*fibg2
           nfib3 = nfibgrps*fibg3
           nfib4 = nfibgrps*fibg4
           do i=1,i1-k+1
             xt(i) = wrandom(xrandom,i+k-1+ioff,dig1,nfib2,npm)
             yt(i) = wrandom(xrandom,i+k-1+ioff,dig2,nfib3,npm)
             uxt(i) = wrandom(xrandom,i+k-1+ioff,dig3,nfib4,npm)
           enddo

!$OMP CRITICAL (STPTCL3D2)
#ifdef _OPENMP
           ip = iptotal
           iptotal = iptotal + i1-k+1
#endif
!$OMP END CRITICAL (STPTCL3D2)

c          --- load x,y,ux, and uy evenly onto a 4-D ellipsoid
           do i=1,i1-k+1
             rr = sqrt(xt(i))
             phi1 = 2.*Pi*yt(i)
             phi2 = 2.*Pi*uxt(i)
c            --- put position into particle arrays
c            --- Centroids are added in below, since positions are used
c            --- to set the coherent velocity - the positions used for that
c            --- are relative to the centroid.
             xp(ip+i) = rr*cos(phi1)*at(i)*perpscal(i)
             yp(ip+i) = rr*sin(phi1)*bt(i)*perpscal(i)
             zp(ip+i) = zt(i)
c            --- put velocity into particle arrays
             rr = sqrt(1.-rr*rr)
             vtx = vzt(i)*epsxt(i)/at(i)*perpscal(i) + vthperp_s(is)
             vty = vzt(i)*epsyt(i)/bt(i)*perpscal(i) + vthperp_s(is)
             uxp(ip+i) = vzt(i)*xp(ip+i)*apt(i)/at(i) + vtx*rr*cos(phi2) +
     &                   vzt(i)*xpct(i)
             uyp(ip+i) = vzt(i)*yp(ip+i)*bpt(i)/bt(i) + vty*rr*sin(phi2) +
     &                   vzt(i)*ypct(i)
             uzp(ip+i) = vzt(i)*(1. + vtilt_s(is)*(zmid - zt(i))*zleni) +
     &                   vtzt(i)*uzt(i)
c            --- The centroids are added to the transverse positions
c            --- after the variables are used above to calculate the
c            --- coherent velocity component.
             xp(ip+i) = xp(ip+i) + xct(i)
             yp(ip+i) = yp(ip+i) + yct(i)
           enddo

c          --- add sinusoidal perturbation to uzp
           if (vzperamp /= 0.) then
             do i=ip+1,ip+i1-k+1
               uzp(i) = uzp(i) + vzperamp*sin(2.*pi*zp(i)/vzperlam + vzperphs)
             enddo
           endif

c          --- increment number of particles by size of current group
           ip = ip + i1-k+1

c        --- Load various transverse distributions 
c        --- using local Courant-Snyder zero-current invariants (0 suffixes) 
c        --- appropriate for periodic focusing channels or 

         elseif (distrbtn == "WB"  .or. distrbtn == "Waterbag"  .or.  
     &           distrbtn == "PA"  .or. distrbtn == "Parabolic" .or. 
     &           distrbtn == "TE"  .or. distrbtn == "ThermalEquilibrium" .or. 
     &           distrbtn == "KV0" .or. 
     &           distrbtn == "WB0" .or. distrbtn == "Waterbag0"  .or.  
     &           distrbtn == "PA0" .or. distrbtn == "Parabolic0" .or. 
     &           distrbtn == "GA0" .or. distrbtn == "Gaussian0"  .or. 
     &           distrbtn == "TE"  .or. distrbtn == "ThermalEquilibrium") 
     &           then 

c          --- Load transverse particles in x-x', y-y' phase space using 
c          --- the local beam envelope and emittance in the slice of each 
c          --- particle 

           if (distrbtn == "WB" .or. distrbtn == "Waterbag"  .or.  
     &         distrbtn == "PA" .or. distrbtn == "Parabolic" .or. 
     &         distrbtn == "TE" .or. distrbtn == "ThermalEquilibrium") then 
c            --- use transforms of continuous focusing equilibria 
             call loadperpdist(k+ioff,i1-k+1,xt,yt,rt,tt,uxt,uyt,
     &                         at,bt,apt,bpt,epsxt,epsyt)
           elseif (distrbtn == "KV0" .or. 
     &             distrbtn == "WB0" .or. distrbtn == "Waterbag0"  .or.  
     &             distrbtn == "PA0" .or. distrbtn == "Parabolic0" .or. 
     &             distrbtn == "GA0" .or. distrbtn == "Gaussian0") then 
c            --- use zero applied field Courant-Snyder invariants 
             call loadperpdist0(k+ioff,i1-k+1,xt,yt,uxt,uyt,
     &                          at,bt,apt,bpt,epsxt,epsyt)
           endif 

!$OMP CRITICAL (STPTCL3D2)
#ifdef _OPENMP
           ip = iptotal
           iptotal = iptotal + i1-k+1
#endif
!$OMP END CRITICAL (STPTCL3D2)

c          --- load coordinates 
           do i=1,i1-k+1

c            --- set longitudinal particle coordinates 
c                  1st term: distribution about centroid 
c                  2nd term: centroid component  
             xp(ip+i) = xt(i) + xct(i)
             yp(ip+i) = yt(i) + yct(i)
c            --- set the longitudinal particle coordinate 
             zp(ip+i) = zt(i)
c            --- set transverse particle gamma*velocity   
c                  1st term: rescale uxp returned in x' units to uxp  
c                  2nd term: centroid component 
             uxp(ip+i) = vzt(i)*uxt(i) + vzt(i)*xpct(i)
             uyp(ip+i) = vzt(i)*uyt(i) + vzt(i)*ypct(i)
c            --- set longitudinal particle gamma*velocity 
             uzp(ip+i) = vzt(i)*(1. + vtilt_s(is)*(zmid - zt(i))*zleni) +
     &                   vtzt(i)*uzt(i)
             if (npid > 0) pid(ip+i,:) = 0.
           enddo

c          --- add sinusoidal perturbation to uzp
           if (vzperamp /= 0.) then
             do i=ip+1,ip+i1-k+1
               uzp(i) = uzp(i) + vzperamp*sin(2.*pi*zp(i)/vzperlam + vzperphs)
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

c        --- Set flag to indicate that the loop has been passed through once
        firstk = .false.
      enddo
!$OMP END DO

c  Set particle number, indices, etc.

!$OMP MASTER
#ifdef _OPENMP
      ip = iptotal
#endif
      nps(is) = ip - ins(is) + 1
!$OMP END MASTER

      enddo

c     deallocate scratch space for load
c      --- XXX These were commented out - I'm not sure why.
       deallocate(indx,xt,yt,zt,rt,tt)
       deallocate(uxt,uyt,uzt,perpscal)
       deallocate(at,apt,bt,bpt)
       deallocate(epsxt,epsyt,vzt,vtzt,ibeamt)
       deallocate(xct,xpct,yct,ypct)

!$OMP END PARALLEL

c--------------------------------
c  End main loop over species
c--------------------------------

c     --- end of if checking for positive npm
      endif

c     --- Find maximum of zdist, which is used to scale sw.  See below.
      if (nzdist > 0) then
        zdistmax = nrmzdist(0)
        do iz=1,nzdist
          if (nrmzdist(iz) > zdistmax) zdistmax = nrmzdist(iz)
        enddo
      endif

c     --- calculate values for sm, sq, and sw for each species.
      do is = 1, ns
        sq(is) = zion_s(is) * echarge
        sm(is) = aion_s(is) * amu
        zlen = zimax_s(is) - zimin_s(is)
        if (vbeam_s(is) /= 0. .and. nps(is) /= 0) then
          sw(is)=ibeam_s(is)*zlen*sp_fract(is)/
     &                   (vbeam_s(is)*echarge*zion_s(is)*nps(is))
        else
          sw(is) = 0.
        endif

c       --- Adjust weighting of particles since cigar() makes beam more dense
        if (cigarld) then
          sw(is) = sw(is)*(straight_s(is) + (1. - straight_s(is))*2./3.)
        endif

c       --- Scale sw when using zdist so that ibeam/vbeam gives the maximum
c       --- linecharge of the distribution.  sw is scaled by the ratio of the
c       --- sum of zdist and the the product of the max of zdist and the
c       --- number of grid points in zdist.  That is the ratio of the loaded
c       --- linecharge and what the linecharge would have been with a uniform
c       --- distribution.
        if (nzdist > 0) then
          sw(is) = sw(is)/(zdistmax*nzdist)
        endif

      enddo

c     --- if slave, recalculate sw
#ifdef MPIPARALLEL
      call sw_globalsum(ns,sw(1))
#endif

c     --- Apply particle boundary conditions
      call zpartbnd(zmmax,zmmin,dz,zgrid)
      do is=1,ns
        if (nps(is) > 0) then
           call stckxy3d(nps(is),xp(ins(is)),xmmax,xmmin,dx,yp(ins(is)),ymmax,
     &                   ymmin,dy,zp(ins(is)),zmmin,dz,uxp(ins(is)),
     &                   uyp(ins(is)),uzp(ins(is)),gaminv(ins(is)),
     &                   zgrid,zbeam,l2symtry,l4symtry,pboundxy,.true.)
        endif
        call processlostpart(is,clearlostpart,time,zbeam)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timestptcl3d = timestptcl3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setrho3d(rho,rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,depos,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry)
      use GlobalVars
      use Subtimers3d
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: rho1d(0:(1+nx)*(1+ny)*(1+nz)-1)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      character(8):: depos
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c  Sets charge density

c  Algorithm notes: rho array is dimensioned (0:nx,0:ny,0:nz) outside,
c  but is made one dimensional in this routine
c  so cell index into 1d rho array for vectorized deposition is:
c     i + j*(nx+1) + k*(nx+1)*(ny+1)
c  In each case,
c     rho(i  ,j  ,k  ) = rho(i  ,j  ,k  ) + u0 * v0 * w0 * g
c     rho(i+1,j  ,k  ) = rho(i+1,j  ,k  ) + u1 * v0 * w0 * g
c  Note that many changes are possible; for example, we might define
c  ind0(ir) and not use indx; this saves some store operations but
c  leads to a more complicated indirect address for the vectorized
c  gather-add-scatter loop.  It seems about 3% slower than the present way.
c  RHO must be zeroed in ZERORHO since it is not zeroed here (to allow
c  handling of blocks of particles at a time)

c     --- For vectorized algorithm
      integer(ISZ):: moff(0:7), indx(0:7,0:nparpgrp-1)
      save moff
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ):: ii(0:nparpgrp-1), jj(0:nparpgrp-1), kk(0:nparpgrp-1)
c     --- For both algorithms
      real(kind=8):: s(0:7,0:nparpgrp-1)

      integer(ISZ):: ipmin,nptmp,ip,i,j,k,ind0,m,ir
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Set up offset array for vectorized deposition:

c      if (moff(4) /= (nx+1)*(ny+1)) then
        moff(0) = 0
        moff(1) = 1
        moff(2) = nx+1
        moff(3) = nx+2
        moff(4) = (nx+1)*(ny+1)
        moff(5) = (nx+1)*(ny+1)+1
        moff(6) = (nx+1)*(ny+1)+nx+1
        moff(7) = (nx+1)*(ny+1)+nx+2
c      endif

      g = wght * q / (dx * dy * dz)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      if (l2symtry) then
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        g = g*0.5
      elseif (l4symtry) then
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        g = g*0.25
      endif

! np was made FIRSTPRIVATE to get around a bug when the expression
! np+1-ipmin was evaluating to 1-ipmin (as if np was zero).
! I don't know why it works, but it does.
!$OMP PARALLEL PRIVATE(ipmin,nptmp,i,j,k,u1,u0,v1,v0,w1,w0,ir,ip,ind0,indx,
!$OMP&                 gyfact,gxfact,s,m,ii,jj,kk)
!$OMP&FIRSTPRIVATE(np)

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
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        do ip = ipmin,ipmin+nptmp-1
           i  = int((xp(ip) - xmmin) * dxi)
           u1 =     (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = int((abs(yp(ip)) - ymmin)*dyi)
           v1 =     (abs(yp(ip)) - ymmin)*dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 =     (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ir = ip - ipmin
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           gyfact = 1.
           if (j == 0 .and. ymmin == 0.) gyfact = 2.
           s(0,ir) = u0 * v0 * w0 * g*gyfact
           s(1,ir) = u1 * v0 * w0 * g*gyfact
           s(2,ir) = u0 * v1 * w0 * g
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g*gyfact
           s(5,ir) = u1 * v0 * w1 * g*gyfact
           s(6,ir) = u0 * v1 * w1 * g
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      elseif (l4symtry) then
c       --- special loop for 4-fold symmetry
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        do ip = ipmin,ipmin+nptmp-1
           i  = int((abs(xp(ip)) - xmmin)*dxi)
           u1 =     (abs(xp(ip)) - xmmin)*dxi - i
           u0 = 1. - u1
           j  = int((abs(yp(ip)) - ymmin)*dyi)
           v1 =     (abs(yp(ip)) - ymmin)*dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 =     (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ir = ip - ipmin
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           gxfact = 1.
           gyfact = 1.
           if (i == 0 .and. xmmin == 0.) gxfact = 2.
           if (j == 0 .and. ymmin == 0.) gyfact = 2.
           s(0,ir) = u0 * v0 * w0 * g*gxfact*gyfact
           s(1,ir) = u1 * v0 * w0 * g*gyfact
           s(2,ir) = u0 * v1 * w0 * g*gxfact
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g*gxfact*gyfact
           s(5,ir) = u1 * v0 * w1 * g*gyfact
           s(6,ir) = u0 * v1 * w1 * g*gxfact
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
           i  = int((xp(ip) - xmmin) * dxi)
           u1 =     (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = int((yp(ip) - ymmin) * dyi)
           v1 =     (yp(ip) - ymmin) * dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 =     (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ir = ip - ipmin
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           s(0,ir) = u0 * v0 * w0 * g
           s(1,ir) = u1 * v0 * w0 * g
           s(2,ir) = u0 * v1 * w0 * g
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g
           s(5,ir) = u1 * v0 * w1 * g
           s(6,ir) = u0 * v1 * w1 * g
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      endif
c     --- vectorized deposition over the 8 cells touched;
c     --- there'd be a hazard if we interchanged the loops.
!$OMP CRITICAL (CRITICAL_SETRHO3D1)
      do ir = 0,nptmp-1
        if (uzp(ir+ipmin) /= 0.) then
          do m = 0, 7
             rho1d(indx(m,ir)) = rho1d(indx(m,ir)) + s(m,ir)
          enddo
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETRHO3D1)

c--------------------------------------
c  Begin scalar deposition loop
c--------------------------------------
      elseif (depos == "scalar") then

c  Begin main loop over species

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        do ip = ipmin,ipmin+nptmp-1
          ir = ip - ipmin
          ii(ir) = int((xp(ip) - xmmin) * dxi)
          u1     =     (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(yp(ip)) - ymmin)*dyi)
          v1     =     (abs(yp(ip)) - ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     =     (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gyfact = 1.
          if (jj(ir) == 0 .and. ymmin == 0.) gyfact = 2.
          s(0,ir) = u0 * v0 * w0 * g*gyfact
          s(1,ir) = u1 * v0 * w0 * g*gyfact
          s(2,ir) = u0 * v1 * w0 * g
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g*gyfact
          s(5,ir) = u1 * v0 * w1 * g*gyfact
          s(6,ir) = u0 * v1 * w1 * g
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      elseif (l4symtry) then
        do ip = ipmin,ipmin+nptmp-1
          ir = ip - ipmin
          ii(ir) = int((abs(xp(ip)) - xmmin)*dxi)
          u1     =     (abs(xp(ip)) - xmmin)*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(yp(ip)) - ymmin)*dyi)
          v1     =     (abs(yp(ip)) - ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     =     (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii(ir) == 0 .and. xmmin == 0.) gxfact = 2.
          if (jj(ir) == 0 .and. ymmin == 0.) gyfact = 2.
          s(0,ir) = u0 * v0 * w0 * g*gxfact*gyfact
          s(1,ir) = u1 * v0 * w0 * g*gyfact
          s(2,ir) = u0 * v1 * w0 * g*gxfact
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g*gxfact*gyfact
          s(5,ir) = u1 * v0 * w1 * g*gyfact
          s(6,ir) = u0 * v1 * w1 * g*gxfact
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
          ir = ip - ipmin
          ii(ir) = int((xp(ip) - xmmin) * dxi)
          u1     = (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((yp(ip) - ymmin) * dyi)
          v1     = (yp(ip) - ymmin) * dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     = (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          s(0,ir) = u0 * v0 * w0 * g
          s(1,ir) = u1 * v0 * w0 * g
          s(2,ir) = u0 * v1 * w0 * g
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g
          s(5,ir) = u1 * v0 * w1 * g
          s(6,ir) = u0 * v1 * w1 * g
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      endif
c     --- scalar loop does the actual deposition
!$OMP CRITICAL (CRITICAL_SETRHO3D2)
      do ir = 0, nptmp-1
         if (uzp(ir+ipmin) /= 0) then
         rho(ii(ir)  ,jj(ir)  ,kk(ir)  )=rho(ii(ir)  ,jj(ir)  ,kk(ir)  )+s(0,ir)
         rho(ii(ir)+1,jj(ir)  ,kk(ir)  )=rho(ii(ir)+1,jj(ir)  ,kk(ir)  )+s(1,ir)
         rho(ii(ir)  ,jj(ir)+1,kk(ir)  )=rho(ii(ir)  ,jj(ir)+1,kk(ir)  )+s(2,ir)
         rho(ii(ir)+1,jj(ir)+1,kk(ir)  )=rho(ii(ir)+1,jj(ir)+1,kk(ir)  )+s(3,ir)
         rho(ii(ir)  ,jj(ir)  ,kk(ir)+1)=rho(ii(ir)  ,jj(ir)  ,kk(ir)+1)+s(4,ir)
         rho(ii(ir)+1,jj(ir)  ,kk(ir)+1)=rho(ii(ir)+1,jj(ir)  ,kk(ir)+1)+s(5,ir)
         rho(ii(ir)  ,jj(ir)+1,kk(ir)+1)=rho(ii(ir)  ,jj(ir)+1,kk(ir)+1)+s(6,ir)
         rho(ii(ir)+1,jj(ir)+1,kk(ir)+1)=rho(ii(ir)+1,jj(ir)+1,kk(ir)+1)+s(7,ir)
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETRHO3D2)

      endif

      enddo
!$OMP END DO

!$OMP END PARALLEL

      if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine setrho3dselect(rho,rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,depos,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry)
      use GlobalVars
      use Subtimers3d
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: rho1d(0:(1+nx)*(1+ny)*(1+nz)-1)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      character(8):: depos
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c  Sets charge density

c  Algorithm notes: rho array is dimensioned (0:nx,0:ny,0:nz) outside,
c  but is made one dimensional in this routine
c  so cell index into 1d rho array for vectorized deposition is:
c     i + j*(nx+1) + k*(nx+1)*(ny+1)
c  In each case,
c     rho(i  ,j  ,k  ) = rho(i  ,j  ,k  ) + u0 * v0 * w0 * g
c     rho(i+1,j  ,k  ) = rho(i+1,j  ,k  ) + u1 * v0 * w0 * g
c  Note that many changes are possible; for example, we might define
c  ind0(ir) and not use indx; this saves some store operations but
c  leads to a more complicated indirect address for the vectorized
c  gather-add-scatter loop.  It seems about 3% slower than the present way.
c  RHO must be zeroed in ZERORHO since it is not zeroed here (to allow
c  handling of blocks of particles at a time)

c     --- For vectorized algorithm
      integer(ISZ):: moff(0:7), indx(0:7,0:nparpgrp-1)
      save moff
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ):: ii(0:nparpgrp-1), jj(0:nparpgrp-1), kk(0:nparpgrp-1)
c     --- For both algorithms
      real(kind=8):: s(0:7,0:nparpgrp-1)

      integer(ISZ):: ipmin,nptmp,ip,i,j,k,ind0,m,ir
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Set up offset array for vectorized deposition:

c      if (moff(4) /= (nx+1)*(ny+1)) then
        moff(0) = 0
        moff(1) = 1
        moff(2) = nx+1
        moff(3) = nx+2
        moff(4) = (nx+1)*(ny+1)
        moff(5) = (nx+1)*(ny+1)+1
        moff(6) = (nx+1)*(ny+1)+nx+1
        moff(7) = (nx+1)*(ny+1)+nx+2
c      endif

      g = wght * q / (dx * dy * dz)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      if (l2symtry) then
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        g = g*0.5
      elseif (l4symtry) then
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        g = g*0.25
      endif

! np was made FIRSTPRIVATE to get around a bug when the expression
! np+1-ipmin was evaluating to 1-ipmin (as if np was zero).
! I don't know why it works, but it does.
!$OMP PARALLEL PRIVATE(ipmin,nptmp,i,j,k,u1,u0,v1,v0,w1,w0,ir,ip,ind0,indx,
!$OMP&                 gyfact,gxfact,s,m,ii,jj,kk)
!$OMP&FIRSTPRIVATE(np)

!$OMP DO
      do ipmin = 1,np,nparpgrp
        nptmp = min(nparpgrp, np+1-ipmin)
        s = 0.

c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
      if (depos == "vector") then

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
c       --- special loop for 2-fold symmetry
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        do ip = ipmin,ipmin+nptmp-1
           if (uzp(ip) == 0.) cycle
           i  = int((xp(ip) - xmmin) * dxi)
           u1 =     (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = int((abs(yp(ip)) - ymmin)*dyi)
           v1 =     (abs(yp(ip)) - ymmin)*dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 =     (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ir = ip - ipmin
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           gyfact = 1.
           if (j == 0 .and. ymmin == 0.) gyfact = 2.
           s(0,ir) = u0 * v0 * w0 * g*gyfact
           s(1,ir) = u1 * v0 * w0 * g*gyfact
           s(2,ir) = u0 * v1 * w0 * g
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g*gyfact
           s(5,ir) = u1 * v0 * w1 * g*gyfact
           s(6,ir) = u0 * v1 * w1 * g
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      elseif (l4symtry) then
c       --- special loop for 4-fold symmetry
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        do ip = ipmin,ipmin+nptmp-1
           if (uzp(ip) == 0.) cycle
           i  = int((abs(xp(ip)) - xmmin)*dxi)
           u1 =     (abs(xp(ip)) - xmmin)*dxi - i
           u0 = 1. - u1
           j  = int((abs(yp(ip)) - ymmin)*dyi)
           v1 =     (abs(yp(ip)) - ymmin)*dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 =     (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ir = ip - ipmin
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           gxfact = 1.
           gyfact = 1.
           if (i == 0 .and. xmmin == 0.) gxfact = 2.
           if (j == 0 .and. ymmin == 0.) gyfact = 2.
           s(0,ir) = u0 * v0 * w0 * g*gxfact*gyfact
           s(1,ir) = u1 * v0 * w0 * g*gyfact
           s(2,ir) = u0 * v1 * w0 * g*gxfact
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g*gxfact*gyfact
           s(5,ir) = u1 * v0 * w1 * g*gyfact
           s(6,ir) = u0 * v1 * w1 * g*gxfact
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
           if (uzp(ip) == 0.) cycle
           i  = int((xp(ip) - xmmin) * dxi)
           u1 =     (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = int((yp(ip) - ymmin) * dyi)
           v1 =     (yp(ip) - ymmin) * dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 =     (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ir = ip - ipmin
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           s(0,ir) = u0 * v0 * w0 * g
           s(1,ir) = u1 * v0 * w0 * g
           s(2,ir) = u0 * v1 * w0 * g
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g
           s(5,ir) = u1 * v0 * w1 * g
           s(6,ir) = u0 * v1 * w1 * g
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      endif
c     --- vectorized deposition over the 8 cells touched;
c     --- there'd be a hazard if we interchanged the loops.
!$OMP CRITICAL (CRITICAL_SETRHO3D1)
      do ir = 0,nptmp-1
        if (uzp(ir+ipmin) /= 0.) then
          do m = 0, 7
             rho1d(indx(m,ir)) = rho1d(indx(m,ir)) + s(m,ir)
          enddo
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETRHO3D1)

c--------------------------------------
c  Begin scalar deposition loop
c--------------------------------------
      elseif (depos == "scalar") then

c  Begin main loop over species

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        do ip = ipmin,ipmin+nptmp-1
          if (uzp(ip) == 0.) cycle
          ir = ip - ipmin
          ii(ir) = int((xp(ip) - xmmin) * dxi)
          u1     =     (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(yp(ip)) - ymmin)*dyi)
          v1     =     (abs(yp(ip)) - ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     =     (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gyfact = 1.
          if (jj(ir) == 0 .and. ymmin == 0.) gyfact = 2.
          s(0,ir) = u0 * v0 * w0 * g*gyfact
          s(1,ir) = u1 * v0 * w0 * g*gyfact
          s(2,ir) = u0 * v1 * w0 * g
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g*gyfact
          s(5,ir) = u1 * v0 * w1 * g*gyfact
          s(6,ir) = u0 * v1 * w1 * g
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      elseif (l4symtry) then
        do ip = ipmin,ipmin+nptmp-1
          if (uzp(ip) == 0.) cycle
          ir = ip - ipmin
          ii(ir) = int((abs(xp(ip)) - xmmin)*dxi)
          u1     =     (abs(xp(ip)) - xmmin)*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(yp(ip)) - ymmin)*dyi)
          v1     =     (abs(yp(ip)) - ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     =     (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii(ir) == 0 .and. xmmin == 0.) gxfact = 2.
          if (jj(ir) == 0 .and. ymmin == 0.) gyfact = 2.
          s(0,ir) = u0 * v0 * w0 * g*gxfact*gyfact
          s(1,ir) = u1 * v0 * w0 * g*gyfact
          s(2,ir) = u0 * v1 * w0 * g*gxfact
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g*gxfact*gyfact
          s(5,ir) = u1 * v0 * w1 * g*gyfact
          s(6,ir) = u0 * v1 * w1 * g*gxfact
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
          if (uzp(ip) == 0.) cycle
          ir = ip - ipmin
          ii(ir) = int((xp(ip) - xmmin) * dxi)
          u1     = (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((yp(ip) - ymmin) * dyi)
          v1     = (yp(ip) - ymmin) * dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     = (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          s(0,ir) = u0 * v0 * w0 * g
          s(1,ir) = u1 * v0 * w0 * g
          s(2,ir) = u0 * v1 * w0 * g
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g
          s(5,ir) = u1 * v0 * w1 * g
          s(6,ir) = u0 * v1 * w1 * g
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      endif
c     --- scalar loop does the actual deposition
!$OMP CRITICAL (CRITICAL_SETRHO3D2)
      do ir = 0, nptmp-1
         if (uzp(ir+ipmin) /= 0) then
         rho(ii(ir)  ,jj(ir)  ,kk(ir)  )=rho(ii(ir)  ,jj(ir)  ,kk(ir)  )+s(0,ir)
         rho(ii(ir)+1,jj(ir)  ,kk(ir)  )=rho(ii(ir)+1,jj(ir)  ,kk(ir)  )+s(1,ir)
         rho(ii(ir)  ,jj(ir)+1,kk(ir)  )=rho(ii(ir)  ,jj(ir)+1,kk(ir)  )+s(2,ir)
         rho(ii(ir)+1,jj(ir)+1,kk(ir)  )=rho(ii(ir)+1,jj(ir)+1,kk(ir)  )+s(3,ir)
         rho(ii(ir)  ,jj(ir)  ,kk(ir)+1)=rho(ii(ir)  ,jj(ir)  ,kk(ir)+1)+s(4,ir)
         rho(ii(ir)+1,jj(ir)  ,kk(ir)+1)=rho(ii(ir)+1,jj(ir)  ,kk(ir)+1)+s(5,ir)
         rho(ii(ir)  ,jj(ir)+1,kk(ir)+1)=rho(ii(ir)  ,jj(ir)+1,kk(ir)+1)+s(6,ir)
         rho(ii(ir)+1,jj(ir)+1,kk(ir)+1)=rho(ii(ir)+1,jj(ir)+1,kk(ir)+1)+s(7,ir)
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETRHO3D2)

      endif

      enddo
!$OMP END DO

!$OMP END PARALLEL

      if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine loadrho3d(ins_i,nps_i,is_i,lzero)
      use GlobalVars
      use Subtimers3d
      use InGen
      use InGen3d
      use InPart
      use InMesh3d
      use Picglb
      use Picglb3d
      use Particles
      use Fields3d
      use Fields3dParticles
      use GridBoundary3d
      integer(ISZ):: ins_i,nps_i,is_i
      logical(ISZ):: lzero

c --- This routine provides a simple call from the interpreter to load the
c --- rhop array.  The value '-1' is used as a flag in the input to use
c --- all of the particles, otherwise the specified particles are loaded.

      integer(ISZ):: ins_u,nps_u
      integer(ISZ):: is1,is2
      integer(ISZ):: ip,ipmin,is
      real(kind=8):: swtmp
      real(kind=8):: substarttime,wtime
      real(kind=8),dimension(:,:,:),pointer::rhopcopy
      integer(ISZ):: allocerror

      if (lw3dtimesubs) substarttime = wtime()

c     --- Ensure that the rhop array is setup properly
      if(solvergeom==XYZgeom) then
        call setupFields3dParticles(ns,ndts,it)
      endif

c     --- zero rhop if requested
      if (lzero) then
        if (solvergeom == XYZgeom) then
          call zeroarry(rhop(0,0,0),(nxp+1)*(nyp+1)*(nzp+1))
        endif
        if(solvergeom==RZgeom .or. solvergeom==XZgeom .or.
     &     solvergeom==Rgeom  .or. solvergeom==Zgeom) call reset_rzmgrid_rho()
      end if

c     --- change AMR grid if necessary
      if(solvergeom==RZgeom .or. solvergeom==XZgeom .or.
     &   solvergeom==Rgeom  .or. solvergeom==Zgeom) call change_loc_part()

c     --- set limits on loop over species
      if (is_i == -1) then
        is1 = 1
        is2 = ns
      else
        is1 = is_i
        is2 = is_i
      endif

      if(depos /= 'none') then

c     --- set initial limits from input
c     --- (will be changed if necessary in the loop)
      ins_u = ins_i
      nps_u = nps_i

c     --- loop over species
      do is=is1,is2

c        --- For the AMR/Chombo version, all of the lost particles must be
c        --- cleared out.
         if (solvergeom == AMRgeom) call clearpart(is,1)

c        --- get loop limits for particles if needed
         if (ins_i == -1) ins_u = ins(is)
         if (nps_i == -1) nps_u = nps(is)

c        --- Scale the weight, sw, by the time step scale size. This only
c        --- makes sense for steady-state and slice modes. In time-dependent
c        --- mode, it is assumed that dtscale has not been changed from 1.
         swtmp = sw(is)*dtscale(is)

         if(solvergeom==XYZgeom) then
           ipmin = ins_u
           if (nps_u > 0) then
             if (ndts(is)>1) then
c              --- This code would be better if it were in setrho3d, but since
c              --- rhop is passed into it, this must be done here.
               if (mod(it+1,ndts(is))==0) then
c                 --- If this species is being advanced this step, then
c                 --- deposit its rho into rhoptmp.
                  rhopcopy => rhop
                  rhop => rhoptmp(:,:,:,jsrhoptmp(is-1))
                  rhop = 0.
               else
c                 --- If this species is not being advanced this step, then
c                 --- just add its saved rhoptmp into rhop.
                  rhop = rhop + rhoptmp(:,:,:,jsrhoptmp(is-1))
                  cycle
                end if
             end if
             call setrho3d(rhop,rhop,nps_u,xp(ipmin),yp(ipmin),zp(ipmin),
     &                     zgrid,uzp(ipmin),sq(is),swtmp,depos,
     &                     nxp,nyp,nzp,dx,dy,dz,xmmin,ymmin,zmminp,
     &                     l2symtry,l4symtry)
             if (ndts(is)>1 .and. mod(it+1,ndts(is))==0) then
c               --- If this species is being advanced this step, then restore
c               --- rhop and copy its rho into rhop.
                rhop => rhopcopy
                rhop = rhop + rhoptmp(:,:,:,jsrhoptmp(is-1))
             end if
           endif
         elseif(solvergeom==RZgeom .or. solvergeom==XZgeom) then
c          --- loop over particle blocks
           do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
             ip = min(nparpgrp, ins_u+nps_u-ipmin)
             if(wpid==0) then
               call rhoweightrz(xp(ipmin),yp(ipmin),zp(ipmin),ip,
     &                          sq(is)*swtmp,nx,nz,dx,dz,xmmin,zgrid)
             else
               call rhoweightrz_weights(xp(ipmin),yp(ipmin),zp(ipmin),pid(ipmin,wpid),ip,
     &                          sq(is)*swtmp,nx,nz,dx,dz,xmmin,zgrid)
             end if
           enddo
         elseif(solvergeom==Zgeom) then
c          --- loop over particle blocks
           do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
              ip = min(nparpgrp, ins_u+nps_u-ipmin)
              call rhoweightz(zp(ipmin),ip,sq(is)*swtmp,nz,dz,zgrid)
           enddo
         elseif(solvergeom==Rgeom) then
c          --- loop over particle blocks
           do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
              ip = min(nparpgrp, ins_u+nps_u-ipmin)
              if(wpid==0) then
                call rhoweightr(xp(ipmin),yp(ipmin),ip,sq(is)*swtmp,nx,dx,xmmin)
              else
                call rhoweightr_weights(xp(ipmin),yp(ipmin),pid(ipmin,wpid),ip,
     &                                  sq(is)*swtmp,nx,dx,xmmin)
              end if
           enddo
         elseif(solvergeom==AMRgeom) then
c          --- loop over particle blocks
           do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
              ip = min(nparpgrp, ins_u+nps_u-ipmin)
              call cho_setrho3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),zgrid,
     &                          sq(is),swtmp,is,(ipmin-ins(is))/nparpgrp)
           enddo
         end if
      enddo

      if (solvergeom == XYZgeomMR) then
        call execuser("loadrhoMR")
      endif

      endif

c  For parallel version, each processor sends rho to neighboring processors
c  whose field solve region overlap its particle region.
      if (lzero) call getrhoforfieldsolve()

c     --- enforce periodicity
c     --- This is done on the rho array, rather than rhop, since that will
c     --- be used for the field solve. Also, in some cases, the rhop array
c     --- will not cover the full axial extent of the system, so 
c     --- periodicity can over be enforced on rho which does cover the full
c     --- system.
      if (lzero) call perrho3d(rho(0,0,0),nx,ny,nz,bound0,boundxy)

!$OMP MASTER
      if (lw3dtimesubs) timeloadrho3d = timeloadrho3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine getrhoforfieldsolve()
      use InGen
      use InGen3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles

#ifndef MPIPARALLEL

c     --- If rhop is not associated with rho, then copy the data.
      if (.not. associated(rho,rhop)) then
        if (nxp == nx .and. nyp == ny .and. nzp == nz) then
          rho = rhop
        else
          call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
          call remark("ERROR!! rhop and rho are not the same shape!       ")
          call remark("        The space charge will not be properly used!")
          call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        endif
      endif

#else

      use Parallel

c  For parallel version, each processor sends rho to neighboring processors
c  whose field solve region overlap its particle region.
      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        call getrhoforfieldsolverz(nx,nz,rho(0,0,0))
      else if(solvergeom==Zgeom) then
        call getrhoforfieldsolvez(nz,rho(0,0,0))
      elseif (solvergeom==XYZgeom) then
        call getrhoforfieldsolve3d(nx,ny,nz,rho,nxp,nyp,nzp,rhop,
     &                             my_index,nslaves,
     &                             izfsslave,nzfsslave,izpslave,nzpslave)
      end if

#endif

      return
      end
c=============================================================================
      subroutine getrhoforfieldsolve3d(nx,ny,nz,rho,nxp,nyp,nzp,rhop,
     &                                 my_index,nslaves,
     &                                 izfsslave,nzfsslave,izpslave,nzpslave)
      integer(ISZ):: nx,ny,nz,nxp,nyp,nzp
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)

#ifdef MPIPARALLEL

      call sumrhoondomainboundaries(rhop,nxp,nyp,nzp,
     &                              my_index,nslaves,izpslave,nzpslave)
      call getrhoforfieldsolve3d_parallel(nx,ny,nz,rho,nxp,nyp,nzp,rhop,
     &                                    my_index,nslaves,
     &                                    izfsslave,nzfsslave,izpslave,nzpslave)

#endif

      return
      end
c=============================================================================
      subroutine setupFields3dParticles(ns,ndts,it)
      use Picglb3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles
#ifdef MPIPARALLEL
      use Parallel
#endif

      integer(ISZ):: ns,ndts(0:ns-1),it

c Ensures that the phip and rhop arrays are setup properly.
c Also, setup rhoptmp for the case when species have different time step sizes.

      integer(ISZ):: js,i,oldnsndts
      integer(ISZ),allocatable:: lirhoptmp(:)

      nxp = nx
      nyp = ny
#ifndef MPIPARALLEL
      nzp = nz
      zmminp = zmmin
      zmmaxp = zmmax
#else
      nzp = nzpslave(my_index)
      zmminp = zmminglobal + izpslave(my_index)*dz
      zmmaxp = zmminp + nzp*dz
#endif

c     --- If either rhop of phip are not associated or either does not have
c     --- the correct dimensions as given by nxp, etc, something 
c     --- needs to be done for fix them.
      if (.not. associated(rhop) .or. .not. associated(phip) .or.
     &    ANY(ubound(phip) .ne. (/ nxp,nyp,nzp+1 /) ) .or.
     &    ANY(ubound(rhop) .ne. (/ nxp,nyp,nzp   /) )) then

c       --- Free whatever the rhop and phip arrays are now refering to.
        call gfree('Fields3dParticles')

        if (ALL( (/ nxp,nyp,nzp /) == (/ nx,ny,nzfull /) )) then
c         --- If the sizes of rhop and rho, and phip and phi, are the same,
c         --- then associate rhop with rho and phip with phi
          rhop => rho
          phip => phi
        else
c         --- If the sizes are not the same, then allocate the group
          call gallot('Fields3dParticles',0)
        endif
      endif

c     --- Setup rhoptmp

c     --- First, count how many species take larger time steps and increase
c     --- the size of the arrays appropriately.
      oldnsndts = nsndts
      nsndts = 0
      do js = 0, ns-1
        if (ndts(js) > 1) nsndts = nsndts + 1
      enddo

      if (nsndts > 0) then

c       --- Set local copy of the number of species, used to size jsrhoptmp
        nsrhoptmp = ns

c       --- Only change size if it has increased. If it is decreased, do the
c       --- gchange afterward so that data that may still be needed at the end
c       --- of the array is not thrown away.
        if (nsndts > oldnsndts) call gchange("Fields3dParticles",0)

c       --- Now set the species index into the rhoptmp array's last dimension.
c       --- Clear out any that have have ndts set back to 1 and find which
c       --- indices are already taken. If ndts is changed and a value is
c       --- changed from or to 1, the species that already have ndts > 1
c       --- must keep the existing data in rhoptmp. This code fills in any
c       --- empty spaces with new species or ones at the end of the list.
c       --- The lirhoptmp array flags which indices have already been claimed.
c       --- Also, jsrhoptmp defaults to -1 for species which have ndts == 1.
        allocate(lirhoptmp(0:nsndts-1))
        lirhoptmp = 0
        do js = 0, ns-1
          if (ndts(js) == 1) jsrhoptmp(js) = -1
          if (ndts(js) > 0 .and. jsrhoptmp(js) >= 0)
     &       lirhoptmp(jsrhoptmp(js)) = 1
        enddo

        i = 0
        do js = 0, ns-1
          if (ndts(js) > 1) then
c           --- If ndts is now set, or if this species index falls beyond the
c           --- new size of the array, find a spot for it.
            if (jsrhoptmp(js) == -1 .or. jsrhoptmp(js) >= nsndts) then
c             --- Find the next available index
              do while (lirhoptmp(i) == 1)
                i = i + 1
              enddo
c             --- If the species was beyond the end of the array, copy the
c             --- data into the new place, unless it will be recalculated
c             --- anyway this step.
              if (jsrhoptmp(js) >= nsndts .and. mod(it+1,ndts(js)) /= 0) then
                rhoptmp(:,:,:,i) = rhoptmp(:,:,:,jsrhoptmp(js))
              endif
c             --- Finally, do the assignment
              jsrhoptmp(js) = i
              lirhoptmp(i) = 1
            endif
          endif
        enddo
        deallocate(lirhoptmp)
      endif

c     --- If the new size is smaller, the change can now be done since the
c     --- data that needed to be saved is copied to its new place.
c     --- Note that this is done outside the if-block above to catch the case
c     --- when nsndts==0 and oldnsndts>0.
      if (nsndts < oldnsndts) call gchange("Fields3dParticles",0)

      return
      end
c=============================================================================
      subroutine fetche3d(ipmin,ip,is)
      use GlobalVars
      use Subtimers3d
      use Picglb
      use Picglb3d
      use InGen
      use InGen3d
      use InMesh3d
      use Particles
      use Fields3d
      use Fields3dParticles
      use Efields3d
      use FieldSolveAPI
      use DKInterptmp
      integer(ISZ):: ipmin,ip,is

c     --- Obtain the self-field from the electrostatic potential

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if(solvergeom==XYZgeom) then
        call sete3d(phip(0,0,-1),selfe(1,0,0,0),ip,
     &              xp(ipmin),yp(ipmin),zp(ipmin),
     &              zgridprv,xmmin,ymmin,zmminp,dx,dy,dz,nxp,nyp,nzp,
     &              efetch,ex,ey,ez,l2symtry,l4symtry)
      elseif(solvergeom==XYZgeomMR) then
        xfsapi => xp(ipmin:ipmin+ip-1)
        yfsapi => yp(ipmin:ipmin+ip-1)
        zfsapi => zp(ipmin:ipmin+ip-1)
        exfsapi => ex
        eyfsapi => ey
        ezfsapi => ez
        call execuser("fetcheMR")
        nullify(xfsapi)
        nullify(yfsapi)
        nullify(zfsapi)
        nullify(exfsapi)
        nullify(eyfsapi)
        nullify(ezfsapi)
      elseif(solvergeom==RZgeom) then
        call setemgridrz(ipmin,ip,is,ex,ey,ez)
      elseif(solvergeom==XZgeom) then
         call fieldweightxz(xp(ipmin),zp(ipmin),ex,ez,ip,zgridprv)
      elseif(solvergeom==Zgeom) then
         call fieldweightz(zp(ipmin),ez,ip,zgridprv)
      elseif(solvergeom==Rgeom) then
         call fieldweightr(xp(ipmin),yp(ipmin),ex,ey,ip)
      elseif(solvergeom==AMRgeom) then
        call cho_gete3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),zgridprv,
     &                  ex,ey,ez,is,(ipmin-ins(is))/nparpgrp)
      else
        ex = 0.
        ey = 0.
        ez = 0.
      endif

!$OMP MASTER
      if (lw3dtimesubs) timefetche3d = timefetche3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine particleboundaries3d()
      use Subtimers3d
      use InGen
      use InGen3d
      use InPart
      use InMesh3d
      use Picglb
      use Picglb3d
      use Particles

c Does everything needed to apply particle boundary conditions.

      integer(ISZ):: is
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Apply particle boundary conditions
      call zpartbnd(zmmax,zmmin,dz,zgrid)

c     --- absorb particles passing through semitransparent disc
      call semitransparent_disc(dz)

      do is=1,ns
        if (mod(it,ndts(is))==0 .and. nps(is) > 0)
     &     call stckxy3d(nps(is),xp(ins(is)),xmmax,xmmin,dx,yp(ins(is)),ymmax,
     &                   ymmin,dy,zp(ins(is)),zmmin,dz,uxp(ins(is)),
     &                   uyp(ins(is)),uzp(ins(is)),gaminv(ins(is)),
     &                   zgrid,zbeam,l2symtry,l4symtry,pboundxy,.true.)
      enddo

      if (lcallscraper) call execuser("callscraper")

      do is=1,ns
c       --- do it only if species was advanced on this time step
        if (mod(it,ndts(is))/=0) cycle
        call processlostpart(is,clearlostpart,time+dt*ndts(is),zbeam)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timeparticleboundaries3d = timeparticleboundaries3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine stckxy3d(np,xp,xmmax,xmmin,dx,yp,ymmax,ymmin,dy,
     &                    zp,zmmin,dz,uxp,uyp,uzp,gaminv,zgrid,zbeam,
     &                    l2symtry,l4symtry,pboundxy,lcountaslost)
      use Subtimers3d
      use GlobalVars
      use Z_arrays
      use InGen3d, only: solvergeom, RZgeom, XZgeom, Zgeom
      integer(ISZ):: np,pboundxy
      real(kind=8):: xmmax,xmmin,dx,ymmax,ymmin,dy,zmmin,dz,zgrid,zbeam
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np)
      logical(ISZ):: l2symtry,l4symtry,lcountaslost

c Enforces transverse sticky boundary conditions.
c Particles are considered stuck when outside a z-dependent aperture of
c radius prwallz, centered at (prwallxz, prwallyz). Particles are also
c lost if outside of the computational domain.
c Lost particles are moved to the edge transversely and to the left end of
c the grid longitudinally, as well as having the velocity set to zero. Zero
c longitudinal velocity is the flag for lost particles. Note that particles
c previously lost are continually moved to the left end of the grid (which
c is important with a moving grid).
c If pboundxy == periodic applies periodic transverse bc.
c Rectangular boundary capability added, where the values of the xmaxz,xminz,
c ymaxz,yminz arrays are used as a rectangular sticky boundary when
c their values(z) are within the grid.
c If lcountaslost is not true, then the uzp of the lost particles is
c set to zero - those particles then won't be added to the count of
c lost particles.


      integer(ISZ):: ip,iz,izz
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: xper,yper
      logical(ISZ):: lstckxy
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      call getparticleextant(1,xmin,xmax,ymin,ymax,
     &                         xmmin,xmmax,ymmin,ymmax)
      call getparticleextant(nzzarr+1,xminz,xmaxz,yminz,ymaxz,
     &                                xminz,xmaxz,yminz,ymaxz)

      xper = xmax - xmin
      yper = ymax - ymin

c     --- Apply periodic boundary condition
      if(pboundxy == periodic) then
!$OMP PARALLEL DO
        do ip=1,np
          if(xp(ip) > xmax) then
            xp(ip) = xp(ip) - xper
          elseif(xp(ip) < xmin) then
            xp(ip) = xp(ip) + xper
          endif
          if(yp(ip) > ymax) then
            yp(ip) = yp(ip) - yper
          elseif(yp(ip) < ymin) then
            yp(ip) = yp(ip) + yper
          endif
        enddo
!$OMP END PARALLEL DO
      endif

c     --- Apply reflection boundary condition
      if(pboundxy == reflect) then
!$OMP PARALLEL DO
        do ip=1,np
          if(xp(ip) > xmax) then
            xp(ip) = xmax - (xp(ip) - xmax)
            uxp(ip) = -uxp(ip)
          elseif(xp(ip) < xmin) then
            xp(ip) = xmin + (xmin - xp(ip))
            uxp(ip) = -uxp(ip)
          endif
          if(yp(ip) > ymax) then
            yp(ip) = ymax - (yp(ip) - ymax)
            uyp(ip) = -uyp(ip)
          elseif(yp(ip) < ymin) then
            yp(ip) = ymin + (ymin - yp(ip))
            uyp(ip) = -uyp(ip)
          endif
        enddo
!$OMP END PARALLEL DO
      endif

c     --- Check values of rectangular boundary arrays to see if they will
c     --- be used.
      lstckxy = .false.
      do izz = 0, nzzarr
        if(xmaxz(izz) <= xmax) lstckxy = .true.
        if(xminz(izz) >= xmin) lstckxy = .true.
        if(ymaxz(izz) <= ymax) lstckxy = .true.
        if(yminz(izz) >= ymin) lstckxy = .true.
      enddo

      if(lstckxy) then
!$OMP PARALLEL DO PRIVATE(iz)
        do ip=1,np
          iz = max(0,min(nzzarr,int((zp(ip) - zzmin - zbeam)*dzzi)))
          if (  xp(ip) >= xmaxz(iz)
     &        .or. yp(ip) >= ymaxz(iz)
     &        .or. xp(ip) <= xminz(iz)
     &        .or. yp(ip) <= yminz(iz) ) then
            gaminv(ip) = 0.
            if (.not. lcountaslost) uzp(ip) = 0.
          endif
        enddo
!$OMP END PARALLEL DO
      endif
      
!$OMP PARALLEL DO PRIVATE(iz)
      do ip=1,np
        iz = max(0,min(nzzarr,int((zp(ip) - zzmin - zbeam)*dzzi)))
        if ((((xp(ip)-prwallxz(iz))*prwelipz(iz))**2 +
     &        (yp(ip)-prwallyz(iz))**2) >= (prwallz(iz)*prwelipz(iz))**2 .or.
     &      (xmin > xp(ip) .or. xp(ip) > xmax .or.
     &       ymin > yp(ip) .or. yp(ip) > ymax)) then
          gaminv(ip) = 0.
          if (.not. lcountaslost) uzp(ip) = 0.
        endif
      enddo
!$OMP END PARALLEL DO

!$OMP MASTER
      if (lw3dtimesubs) timestckxy3d = timestckxy3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine getparticleextant(n,xmin,xmax,ymin,ymax,
     &                               xmmin,xmmax,ymmin,ymmax)
      use InGen3d
      integer(ISZ):: n
      real(kind=8):: xmin(n),xmax(n),ymin(n),ymax(n)
      real(kind=8):: xmmin(n),xmmax(n),ymmin(n),ymmax(n)

      xmin = xmmin
      xmax = xmmax
      ymin = ymmin
      ymax = ymmax
      if (l2symtry) then
        ymin = -ymmax
      elseif (l4symtry) then
        xmin = -xmmax
        ymin = -ymmax
      endif
      if(solvergeom==RZgeom) then
        xmin = -xmmax
        ymin = -xmmax
        ymax =  xmmax
      else if(solvergeom==XZgeom) then
        if(l2symtry .or. l4symtry) then
          xmin = -xmmax
        endif
      endif

      return
      end
c=============================================================================
      subroutine fixgridextent()
      use InMesh3d
      use InGen3d

c      --- Put ny to zero for RZ geometry
      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        ny = 0
      elseif(solvergeom==Zgeom) then
        nx = 0
        ny = 0
      elseif(solvergeom==Rgeom) then
        ny = 0
        nz = 0
      endif

c     --- Put grid min's to zero as appropriate for the symmetry
      if(solvergeom==XYZgeom .or. solvergeom==XYZgeomMR) then
        if (l2symtry) then
          ymmin = 0.
        elseif (l4symtry) then
          xmmin = 0.
          ymmin = 0.
        endif
      else if(solvergeom==XZgeom) then
        ymmin = 0.
        if (l2symtry .or. l4symtry) xmmin = 0.
      else if(solvergeom==RZgeom .or. solvergeom==Zgeom .or. solvergeom==Rgeom) then
        l2symtry=.false.
        l4symtry=.false.
        xmmin = 0.
        ymmin = 0.
      endif

      return
      end
c=============================================================================
      subroutine setrstar(rstar,nz,dz,zmmin,zgrid)
      use Subtimers3d
      integer(ISZ):: nz
      real(kind=8):: rstar(-1:nz+1)
      real(kind=8):: dz,zmmin,zgrid

c  Loads radius of curvature of reference orbit into rstar array, so that
c  r_star is known on mesh points for purposes of field solution.

      integer(ISZ):: iz
      real(kind=8):: zz(-1:nz+1),vz(-1:nz+1),gi(-1:nz+1)
      real(kind=8):: bendres(-1:nz+1)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Fill temporary arrays. Note that getbend expects the pair vz and dt,
c     --- but this passes in the pair 1 and dz which is equivalent.
      do iz = -1,nz+1
        zz(iz) = zmmin + iz*dz + zgrid
        vz(iz) = 1.
        gi(iz) = 1.
      enddo

c     --- Kludgy fix for the end point. The getbend routine relies on the
c     --- internal lattice arrays which only range from 0 to nz. The points
c     --- 0 and nz+1 are outside that range. The fix is to replace then with
c     --- the neighboring point inside the internal lattice range. Note
c     --- that this may give incorrect results, but less so than using
c     --- garbage out of the array bounds.
      zz(-1) = zz(0)
      zz(nz+1) = zz(nz)

c     --- Get the bend radius of curvature and residence fraction at the
c     --- list of z points.
      call getbend(nz+3,nz+3,zz,vz,gi,bendres,rstar,-0.5*dz,0.5*dz,.false.)

c     --- Scale the radius by the reciprical of the residence fraction.
c     --- Also, make sure the rstar is non-zero.
      do iz = -1,nz+1
        if (bendres(iz) > 0.) rstar(iz) = rstar(iz)/bendres(iz)
        if (rstar(iz) == 0.) rstar(iz) = LARGEPOS
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesetrstar = timesetrstar + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine fetchphi(n,x,y,z,p)
      use InGen3d
      use Picglb
      use FieldSolveAPI
      integer(ISZ):: n
      real(kind=8),target:: x(n),y(n),z(n),p(n)

c     --- Call the appropriate routine to get the phi
      if(solvergeom==XYZgeom) then
        call fetchphi3d(n,x,y,z,p)
      elseif(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        call setphirz(n,x,y,z,p,zgrid)
      elseif(solvergeom==Zgeom) then
        call setphiz(n,z,p,zgrid)
      elseif(solvergeom==AMRgeom) then
        call cho_getphi3d(n,x,y,z,0.,p,-1,-1)
      elseif(solvergeom==XYZgeomMR) then
        xfsapi => x
        yfsapi => y
        zfsapi => z
        phifsapi => p
        call execuser("fetchphiMR")
        nullify(xfsapi)
        nullify(yfsapi)
        nullify(zfsapi)
        nullify(phifsapi)
      endif

      call fetchphi_from_pgrd(n,x,y,z,p)

      return
      end
c=============================================================================
      subroutine fetchphi3d(n,x,y,z,p)
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Fields3d
      use Fields3dParticles
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),p(n)

c Fetch phi at points.

      integer(ISZ):: ix,iy,iz,ii
      real(kind=8):: wx,wy,wz
      real(kind=8):: dxi,dyi,dzi

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      if (l4symtry) then
        do ii=1,n
          ix = (abs(x(ii)) - xmmin)*dxi
          iy = (abs(y(ii)) - ymmin)*dyi
          iz =     (z(ii)  - zmminp)*dzi
          wx = (abs(x(ii)) - xmmin)*dxi - ix
          wy = (abs(y(ii)) - ymmin)*dyi - iy
          wz =     (z(ii)  - zmminp)*dzi - iz
          if (xmmin <= abs(x(ii)) .and. abs(x(ii)) < xmmax .and.
     &        ymmin <= abs(y(ii)) .and. abs(y(ii)) < ymmax .and.
     &        zmminp <= z(ii) .and. z(ii) < zmmaxp) then
            p(ii) = phip(ix  ,iy,  iz  )*(1.-wx)*(1.-wy)*(1.-wz) +
     &              phip(ix+1,iy,  iz  )*    wx *(1.-wy)*(1.-wz) +
     &              phip(ix  ,iy+1,iz  )*(1.-wx)*    wy *(1.-wz) +
     &              phip(ix+1,iy+1,iz  )*    wx *    wy *(1.-wz) +
     &              phip(ix  ,iy,  iz+1)*(1.-wx)*(1.-wy)*    wz  +
     &              phip(ix+1,iy,  iz+1)*    wx *(1.-wy)*    wz  +
     &              phip(ix  ,iy+1,iz+1)*(1.-wx)*    wy *    wz  +
     &              phip(ix+1,iy+1,iz+1)*    wx *    wy *    wz
          endif
        enddo
      else if (l2symtry) then
        do ii=1,n
          ix =     (x(ii)  - xmmin)*dxi
          iy = (abs(y(ii)) - ymmin)*dyi
          iz =     (z(ii)  - zmminp)*dzi
          wx =     (x(ii)  - xmmin)*dxi - ix
          wy = (abs(y(ii)) - ymmin)*dyi - iy
          wz =     (z(ii)  - zmminp)*dzi - iz
          if (xmmin <= x(ii) .and. x(ii) < xmmax .and.
     &        ymmin <= abs(y(ii)) .and. abs(y(ii)) < ymmax .and.
     &        zmminp <= z(ii) .and. z(ii) < zmmaxp) then
            p(ii) = phip(ix  ,iy,  iz  )*(1.-wx)*(1.-wy)*(1.-wz) +
     &              phip(ix+1,iy,  iz  )*    wx *(1.-wy)*(1.-wz) +
     &              phip(ix  ,iy+1,iz  )*(1.-wx)*    wy *(1.-wz) +
     &              phip(ix+1,iy+1,iz  )*    wx *    wy *(1.-wz) +
     &              phip(ix  ,iy,  iz+1)*(1.-wx)*(1.-wy)*    wz  +
     &              phip(ix+1,iy,  iz+1)*    wx *(1.-wy)*    wz  +
     &              phip(ix  ,iy+1,iz+1)*(1.-wx)*    wy *    wz  +
     &              phip(ix+1,iy+1,iz+1)*    wx *    wy *    wz
          endif
        enddo
      else
        do ii=1,n
          ix = (x(ii) - xmmin)*dxi
          iy = (y(ii) - ymmin)*dyi
          iz = (z(ii) - zmminp)*dzi
          wx = (x(ii) - xmmin)*dxi - ix
          wy = (y(ii) - ymmin)*dyi - iy
          wz = (z(ii) - zmminp)*dzi - iz
          if (xmmin  <= x(ii) .and. x(ii) < xmmax .and.
     &        ymmin  <= y(ii) .and. y(ii) < ymmax .and.
     &        zmminp <= z(ii) .and. z(ii) < zmmaxp) then
            p(ii) = phip(ix  ,iy,  iz  )*(1.-wx)*(1.-wy)*(1.-wz) +
     &              phip(ix+1,iy,  iz  )*    wx *(1.-wy)*(1.-wz) +
     &              phip(ix  ,iy+1,iz  )*(1.-wx)*    wy *(1.-wz) +
     &              phip(ix+1,iy+1,iz  )*    wx *    wy *(1.-wz) +
     &              phip(ix  ,iy,  iz+1)*(1.-wx)*(1.-wy)*    wz  +
     &              phip(ix+1,iy,  iz+1)*    wx *(1.-wy)*    wz  +
     &              phip(ix  ,iy+1,iz+1)*(1.-wx)*    wy *    wz  +
     &              phip(ix+1,iy+1,iz+1)*    wx *    wy *    wz
          endif
        enddo
      endif

      return
      end
c===========================================================================
      subroutine sete3d_aperture(np,xp,yp,zp,zgrid,xmmin,ymmin,zmmin,
     &                           zmmax,dx,dy,dz,nx,ny,nz,ex,ey,
     &                           l2symtry,l4symtry)
      use Subtimers3d
      use Apertures
      integer(ISZ):: np,nx,ny,nz
      real(kind=8):: zgrid,xmmin,ymmin,zmmin,zmmax,dx,dy,dz
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: ex(np), ey(np)
      logical(ISZ):: l2symtry,l4symtry

c Replaces the transverse E field for particles near an aperture described by
c the variables in the group Apertures. This relys on the routine
c set_aperture_e being called to fill the aper_ex and aper_ey arrays.


      real(kind=8):: zs,ze
      integer(ISZ):: ip,i,j,k,ia,izs,ize
      real(kind=8):: dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,ysign,xsign
      real(kind=8):: sx,sy
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Evaluation of E, vectorized over particles

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c     --- Loop over apertures
      do ia=1,napertures
        zs = aper_zs(ia)
        ze = aper_ze(ia)
        izs = nz - int((zmmax - zs + zgrid) * dzi)
        ize = int((ze - zgrid - zmmin) * dzi)

c       --- Only do calculation if aperture is within the grid.
        if (0 < ize .and. izs < nz .and. izs <= ize) then

          if (.not. (l2symtry .or. l4symtry)) then
            do ip = 1, np

              k =  (zp(ip) - zgrid - zmmin) * dzi
              if (izs-1 <= k .and. k < ize+1) then

                i =  (xp(ip) - xmmin) * dxi
                j =  (yp(ip) - ymmin) * dyi

                u1 = (xp(ip) - xmmin) * dxi - i
                v1 = (yp(ip) - ymmin) * dyi - j
                w1 = (zp(ip) - zgrid - zmmin) * dzi - k

                u0 = 1. - u1
                v0 = 1. - v1
                w0 = 1. - w1

                ex(ip) = u0*v0*w0*aper_ex(i  ,j  ,k-izs  ,ia)
     &                 + u1*v0*w0*aper_ex(i+1,j  ,k-izs  ,ia)
     &                 + u0*v1*w0*aper_ex(i  ,j+1,k-izs  ,ia)
     &                 + u1*v1*w0*aper_ex(i+1,j+1,k-izs  ,ia)
     &                 + u0*v0*w1*aper_ex(i  ,j  ,k-izs+1,ia)
     &                 + u1*v0*w1*aper_ex(i+1,j  ,k-izs+1,ia)
     &                 + u0*v1*w1*aper_ex(i  ,j+1,k-izs+1,ia)
     &                 + u1*v1*w1*aper_ex(i+1,j+1,k-izs+1,ia)

                ey(ip) = u0*v0*w0*aper_ey(i  ,j  ,k-izs  ,ia)
     &                 + u1*v0*w0*aper_ey(i+1,j  ,k-izs  ,ia)
     &                 + u0*v1*w0*aper_ey(i  ,j+1,k-izs  ,ia)
     &                 + u1*v1*w0*aper_ey(i+1,j+1,k-izs  ,ia)
     &                 + u0*v0*w1*aper_ey(i  ,j  ,k-izs+1,ia)
     &                 + u1*v0*w1*aper_ey(i+1,j  ,k-izs+1,ia)
     &                 + u0*v1*w1*aper_ey(i  ,j+1,k-izs+1,ia)
     &                 + u1*v1*w1*aper_ey(i+1,j+1,k-izs+1,ia)

              endif

            enddo

          else

c           --- Set the signs of the E field for particles on negative side of
c           --- the axis of symmetry.
            sy = -1.
            sx = 1.
            if (l4symtry) sx = -1.

c           --- special loop symmetry is used
            do ip = 1, np

              k =  (zp(ip) - zgrid - zmmin)*dzi
              if (izs-1 <= k .and. k < ize+1) then
                i =  (abs(xp(ip)) - xmmin)*dxi
                j =  (abs(yp(ip)) - ymmin)*dyi

                u1 = (abs(xp(ip)) - xmmin)*dxi - i
                v1 = (abs(yp(ip)) - ymmin)*dyi - j
                w1 = (zp(ip) - zgrid - zmmin)*dzi - k

                u0 = 1. - u1
                v0 = 1. - v1
                w0 = 1. - w1

c               --- Adjust sign of E field for approiate quadrant.
                xsign = 1.
                ysign = 1.
                if (xp(ip) < 0.) xsign = sx
                if (yp(ip) < 0.) ysign = sy

                ex(ip) = xsign*(u0*v0*w0*aper_ex(i  ,j  ,k-izs  ,ia)
     &                        + u1*v0*w0*aper_ex(i+1,j  ,k-izs  ,ia)
     &                        + u0*v1*w0*aper_ex(i  ,j+1,k-izs  ,ia)
     &                        + u1*v1*w0*aper_ex(i+1,j+1,k-izs  ,ia)
     &                        + u0*v0*w1*aper_ex(i  ,j  ,k-izs+1,ia)
     &                        + u1*v0*w1*aper_ex(i+1,j  ,k-izs+1,ia)
     &                        + u0*v1*w1*aper_ex(i  ,j+1,k-izs+1,ia)
     &                        + u1*v1*w1*aper_ex(i+1,j+1,k-izs+1,ia))

                ey(ip) = ysign*(u0*v0*w0*aper_ey(i  ,j  ,k-izs  ,ia)
     &                        + u1*v0*w0*aper_ey(i+1,j  ,k-izs  ,ia)
     &                        + u0*v1*w0*aper_ey(i  ,j+1,k-izs  ,ia)
     &                        + u1*v1*w0*aper_ey(i+1,j+1,k-izs  ,ia)
     &                        + u0*v0*w1*aper_ey(i  ,j  ,k-izs+1,ia)
     &                        + u1*v0*w1*aper_ey(i+1,j  ,k-izs+1,ia)
     &                        + u0*v1*w1*aper_ey(i  ,j+1,k-izs+1,ia)
     &                        + u1*v1*w1*aper_ey(i+1,j+1,k-izs+1,ia))

              endif
            enddo

c         --- End of if for types of symmetry
          endif

c       --- End of if checking if aperture is within the grid
        endif

c     --- End of loop over apertures
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesete3d_aperture = timesete3d_aperture + wtime() - substarttime
!$OMP END MASTER
      return
      end
c===========================================================================
      subroutine set_aperture_e()
      use Subtimers3d
      use InGen3d
      use Picglb
      use Picglb3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles
      use Apertures

c Explicitly calculate the E field in the plane of the apertures which
c scrape particles. This method was found empirically to give the
c smoothest variation of E near the edge of the aperture (and therefore
c the lowest emittance). The resulting phase space agrees fairly well with
c the phase space from a much higher resolution simulation, except that in
c the high resolution simulation, there are particles right at the edge of
c the beam that are pushed out the phase space ellipse from the fields
c right at the aperture. Which version more closely models the actual
c behavior of particles near the edge of the aperture I don't know.
c
c There are four different regions which are treated differently.
c   - for point far inside the aperture, the standard finite difference is used
c   - for points right next to the aperture, the standard finite difference is
c     also used, but assuming that the potential at the nearby point within
c     the conductor is obtained from exptrapolation from the current point to
c     the conductor surface.
c   - for points in the conductor that are near the perpendicular axis
c     (y-axis for Ex and x-axis for Ey), the E field is extrapolated from the
c     two points off to one side along the axis of the E-field which are not
c     in the conductor
c   - for the rest of the points in the conductor, the extrapolation is done
c     from the two points off to one side along the axis perpendicular to the
c     E-field.


      integer(ISZ):: midx,midy,minx,maxx,miny,maxy
      integer(ISZ):: minxm1,maxxp1,minym1,maxyp1
      integer(ISZ):: ixmin,ixmax,iymin,iymax
      integer(ISZ):: ia,ix,iy,iz,ixs,iys,izs,ize,iza
      real(kind=8):: zs,ze,xx,yy,wx,wy
      real(kind=8):: tdxi,tdyi,dzi
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      tdxi = 0.5/dx
      tdyi = 0.5/dy
      dzi = 1./dz

c     --- Zero out E field arrays
c     call zeroarry(aper_ex,(nx+1)*(ny+1)*(2+aper_zmax)*napertures)
c     call zeroarry(aper_ey,(nx+1)*(ny+1)*(2+aper_zmax)*napertures)
      aper_ex = 0.
      aper_ey = 0.

c     --- Loop over apertures.
      do ia=1,napertures
        zs = aper_zs(ia)
        ze = aper_ze(ia)
        izs = nzp - int((zmmaxp - zs + zgrid) * dzi)
        ize = int((ze - zgrid - zmminp) * dzi)

c       --- Only do calculation if aperture is within the grid.
        if (0 < ize .and. izs < nz .and. izs <= ize) then

        midx = max(0,min(nxp,int((aper_x(ia)                - xmmin)/dx + 0.5)))
        midy = max(0,min(nyp,int((aper_y(ia)                - ymmin)/dy + 0.5)))
        minx = max(0,min(nxp,int((aper_x(ia) - aper_rad(ia) - xmmin)/dx + 1  )))
        miny = max(0,min(nyp,int((aper_y(ia) - aper_rad(ia) - ymmin)/dy + 1  )))
        maxx = max(0,min(nxp,int((aper_x(ia) + aper_rad(ia) - xmmin)/dx      )))
        maxy = max(0,min(nyp,int((aper_y(ia) + aper_rad(ia) - ymmin)/dy      )))
        minxm1 = max(0,minx-1)
        maxxp1 = min(nxp,maxx+1)
        minym1 = max(0,miny-1)
        maxyp1 = min(nyp,maxy+1)

c       -----------------------------------------------------------------------
c       --- Calculate Ex and Ey for all points using standard finite
c       --- difference.  For the planes before and after the aperture,
c       --- do the calculation over the full transverse extent of the grid.
c       --- Those points near the conductor will be recalculated below.
        do iz=izs-1,ize+1
          iza = iz - izs
          if (iz == izs-1 .or. iz == ize+1) then
            ixmin = 0
            ixmax = nxp-1
            iymin = 0
            iymax = nyp-1
          else
            ixmin = minxm1
            ixmax = maxxp1
            iymin = minym1
            iymax = maxyp1
          endif
          do iy=iymin,iymax
            do ix=ixmin,ixmax
             aper_ex(ix,iy,iza,ia)=tdxi*(phip(abs(ix-1),iy,iz)-phip(ix+1,iy,iz))
             aper_ey(ix,iy,iza,ia)=tdyi*(phip(ix,abs(iy-1),iz)-phip(ix,iy+1,iz))
            enddo
          enddo
        enddo

        do iz=izs,ize
          iza = iz - izs

c         -------------------------------------------------------------------
c         --- Calculate Ex, first looping over the range of y which is within
c         --- the aperture.

c         --- First quadrant (x > midx, y > midy)
          do iy=midy,maxy
            xx = sqrt(abs(aper_rad(ia)**2 - (iy*dy - aper_y(ia) + ymmin)**2)) +
     &           aper_x(ia) - xmmin
            ixs = int(xx/dx)
            ixmax = min(ixs+2,maxxp1)
            do ix=ixs,ixmax
              if (ix*dx < xx) then
                wx = xx/dx - ix
                aper_ex(ix,iy,iza,ia) = tdxi*(phip(ix-1,iy,iz) - phip(ix,iy,iz)+
     &                                  (phip(ix  ,iy,iz) - aper_volt(ia))/wx)
              elseif (ix*dx == xx) then
                aper_ex(ix,iy,iza,ia) = tdxi*2.*(phip(ix-1,iy,iz)-aper_volt(ia))
              else
                if (iy-midy < 2) then
                  aper_ex(ix,iy,iza,ia) = 2.*aper_ex(ix-1,iy,iza,ia) -
     &                                       aper_ex(ix-2,iy,iza,ia)
                else
                  aper_ex(ix,iy,iza,ia) = 2.*aper_ex(ix,iy-1,iza,ia) -
     &                                       aper_ex(ix,iy-2,iza,ia)
                endif
              endif
            enddo
          enddo

c         --- Second quadrant (x < midx, y > midy)
          do iy=midy,maxy
            xx = -sqrt(abs(aper_rad(ia)**2 - (iy*dy - aper_y(ia) + ymmin)**2))+
     &           aper_x(ia) - xmmin
            ixs = nx - int(nx-xx/dx)
            ixmin = max(ixs - 2,minxm1)
            do ix=ixs,ixmin,-1
              if (ix*dx > xx) then
                wx = ix - xx/dx
                aper_ex(ix,iy,iza,ia) = tdxi*(phip(ix,iy,iz) - phip(ix+1,iy,iz)-
     &                                  (phip(ix,iy,iz) - aper_volt(ia))/wx)
              elseif (ix*dx == xx) then
                aper_ex(ix,iy,iza,ia) = tdxi*2.*(aper_volt(ia)-phip(ix+1,iy,iz))
              else
                if (iy-midy < 2) then
                  aper_ex(ix,iy,iza,ia) = 2.*aper_ex(ix+1,iy,iza,ia) -
     &                                       aper_ex(ix+2,iy,iza,ia)
                else
                  aper_ex(ix,iy,iza,ia) = 2.*aper_ex(ix,iy-1,iza,ia) -
     &                                       aper_ex(ix,iy-2,iza,ia)
                endif
              endif
            enddo
          enddo

c         --- Third quadrant (x > midx, y < midy)
          do iy=midy-1,miny,-1
            xx = sqrt(abs(aper_rad(ia)**2 - (iy*dy - aper_y(ia) + ymmin)**2)) +
     &           aper_x(ia) - xmmin
            ixs = int(xx/dx)
            ixmax = min(ixs+2,maxxp1)
            do ix=ixs,ixmax
              if (ix*dx < xx) then
                wx = xx/dx - ix
                aper_ex(ix,iy,iza,ia) = tdxi*(phip(ix-1,iy,iz) - phip(ix,iy,iz)+
     &                                  (phip(ix  ,iy,iz) - aper_volt(ia))/wx)
              elseif (ix*dx == xx) then
                aper_ex(ix,iy,iza,ia) = tdxi*2.*(phip(ix-1,iy,iz)-aper_volt(ia))
              else
                if (iy-midy > -2) then
                  aper_ex(ix,iy,iza,ia) = 2.*aper_ex(ix-1,iy,iza,ia) -
     &                                       aper_ex(ix-2,iy,iza,ia)
                else
                  aper_ex(ix,iy,iza,ia) = 2.*aper_ex(ix,iy+1,iza,ia) -
     &                                       aper_ex(ix,iy+2,iza,ia)
                endif
              endif
            enddo
          enddo

c         --- Fourth quadrant (x < midx, y < midy)
          do iy=midy-1,miny,-1
            xx = -sqrt(abs(aper_rad(ia)**2 - (iy*dy - aper_y(ia) + ymmin)**2))+
     &           aper_x(ia) - xmmin
            ixs = nx - int(nx-xx/dx)
            ixmin = max(ixs - 2,minxm1)
            do ix=ixs,ixmin,-1
              if (ix*dx > xx) then
                wx = ix - xx/dx
                aper_ex(ix,iy,iza,ia) = tdxi*(phip(ix,iy,iz) - phip(ix+1,iy,iz)-
     &                                  (phip(ix,iy,iz) - aper_volt(ia))/wx)
              elseif (ix*dx == xx) then
                aper_ex(ix,iy,iza,ia) = tdxi*2.*(aper_volt(ia)-phip(ix+1,iy,iz))
              else
                if (iy-midy < 2) then
                  aper_ex(ix,iy,iza,ia) = 2.*aper_ex(ix+1,iy,iza,ia) -
     &                                       aper_ex(ix+2,iy,iza,ia)
                else
                  aper_ex(ix,iy,iza,ia) = 2.*aper_ex(ix,iy+1,iza,ia) -
     &                                       aper_ex(ix,iy+2,iza,ia)
                endif
              endif
            enddo
          enddo

c         --- Now set the lines just above and below the aperture, only going
c         --- as far out in x as needed.
          if (maxyp1 == maxy + 1) then
c           --- Find min x
            if (minxm1 == minx - 1) then
              xx = -sqrt(abs(aper_rad(ia)**2-(maxy*dy-aper_y(ia)+ymmin)**2))+
     &             aper_x(ia) - xmmin
              ixmin = int(xx/dx)
            else
              ixmin = minx
            endif
c           --- Find max x
            if (maxxp1 == maxx + 1) then
              xx = sqrt(abs(aper_rad(ia)**2-(maxy*dy-aper_y(ia)+ymmin)**2))+
     &             aper_x(ia) - xmmin
              ixmax = int(xx/dx) + 1
            else
              ixmax = maxx
            endif
c           --- Do the work
            do ix=ixmin,ixmax
              aper_ex(ix,maxyp1,iza,ia) = 2.*aper_ex(ix,maxyp1-1,iza,ia) -
     &                                       aper_ex(ix,maxyp1-2,iza,ia)
            enddo
          endif

          if (minym1 == miny - 1) then
c           --- Find min x
            if (minxm1 == minx - 1) then
              xx = -sqrt(abs(aper_rad(ia)**2-(miny*dy-aper_y(ia)+ymmin)**2))+
     &             aper_x(ia) - xmmin
              ixmin = int(xx/dx)
            else
              ixmin = minx
            endif
c           --- Find max x
            if (maxxp1 == maxx + 1) then
              xx = sqrt(abs(aper_rad(ia)**2-(miny*dy-aper_y(ia)+ymmin)**2))+
     &             aper_x(ia) - xmmin
              ixmax = int(xx/dx) + 1
            else
              ixmax = maxx
            endif
c           --- Do the work
            do ix=ixmin,ixmax
              aper_ex(ix,minym1,iza,ia) = 2.*aper_ex(ix,minym1+1,iza,ia) -
     &                                       aper_ex(ix,minym1+2,iza,ia)
            enddo
          endif

c         -------------------------------------------------------------------
c         --- Calculate Ey, first looping over the range of x which is within
c         --- the aperture.

c         --- First quadrant (x > midx, y > midy)
          do ix=midx,maxx
            yy = sqrt(abs(aper_rad(ia)**2 - (ix*dx - aper_x(ia) + xmmin)**2)) +
     &           aper_y(ia) - ymmin
            iys = int(yy/dy)
            iymax = min(iys+2,maxyp1)
            do iy=iys,iymax
              if (iy*dy < yy) then
                wy = yy/dy - iy
                aper_ey(ix,iy,iza,ia) = tdyi*(phip(ix,iy-1,iz) - phip(ix,iy,iz)+
     &                                  (phip(ix,iy  ,iz) - aper_volt(ia))/wy)
              elseif (iy*dy == yy) then
                aper_ey(ix,iy,iza,ia) = tdyi*2.*(phip(ix,iy-1,iz)-aper_volt(ia))
              else
                if (ix-midx < 2) then
                  aper_ey(ix,iy,iza,ia) = 2.*aper_ey(ix,iy-1,iza,ia) -
     &                                       aper_ey(ix,iy-2,iza,ia)
                else
                  aper_ey(ix,iy,iza,ia) = 2.*aper_ey(ix-1,iy,iza,ia) -
     &                                       aper_ey(ix-2,iy,iza,ia)
                endif
              endif
            enddo
          enddo

c         --- Second quadrant (x < midx, y > midy)
          do ix=midx-1,minx,-1
            yy = sqrt(abs(aper_rad(ia)**2 - (ix*dx - aper_x(ia) + xmmin)**2)) +
     &           aper_y(ia) - ymmin
            iys = int(yy/dy)
            iymax = min(iys+2,maxyp1)
            do iy=iys,iymax
              if (iy*dy < yy) then
                wy = yy/dy - iy
                aper_ey(ix,iy,iza,ia) = tdyi*(phip(ix,iy-1,iz) - phip(ix,iy,iz)+
     &                                  (phip(ix,iy  ,iz) - aper_volt(ia))/wy)
              elseif (iy*dy == yy) then
                aper_ey(ix,iy,iza,ia) = tdyi*2.*(phip(ix,iy-1,iz)-aper_volt(ia))
              else
                if (ix-midx > -2) then
                  aper_ey(ix,iy,iza,ia) = 2.*aper_ey(ix,iy-1,iza,ia) -
     &                                       aper_ey(ix,iy-2,iza,ia)
                else
                  aper_ey(ix,iy,iza,ia) = 2.*aper_ey(ix+1,iy,iza,ia) -
     &                                       aper_ey(ix+2,iy,iza,ia)
                endif
              endif
            enddo
          enddo

c         --- Third quadrant (x > midx, y < midy)
          do ix=midx,maxx
            yy = -sqrt(abs(aper_rad(ia)**2 - (ix*dx - aper_x(ia) + xmmin)**2)) +
     &           aper_y(ia) - ymmin
            iys = ny - int(ny-yy/dy)
            iymin = max(iys-2,minym1)
            do iy=iys,iymin,-1
              if (iy*dy > yy) then
                wy = iy - yy/dy
                aper_ey(ix,iy,iza,ia) = tdyi*(phip(ix,iy,iz) - phip(ix,iy+1,iz)-
     &                                  (phip(ix,iy,iz) - aper_volt(ia))/wy)
              elseif (iy*dy == yy) then
                aper_ey(ix,iy,iza,ia) = tdyi*2.*(aper_volt(ia)-phip(ix,iy+1,iz))
              else
                if (ix-midx < 2) then
                  aper_ey(ix,iy,iza,ia) = 2.*aper_ey(ix,iy+1,iza,ia) -
     &                                       aper_ey(ix,iy+2,iza,ia)
                else
                  aper_ey(ix,iy,iza,ia) = 2.*aper_ey(ix-1,iy,iza,ia) -
     &                                       aper_ey(ix-2,iy,iza,ia)
                endif
              endif
            enddo
          enddo

c         --- Fourth quadrant (x < midx, y < midy)
          do ix=midx-1,minx,-1
            yy = -sqrt(abs(aper_rad(ia)**2 - (ix*dx - aper_x(ia) + xmmin)**2)) +
     &           aper_y(ia) - ymmin
            iys = ny - int(ny-yy/dy)
            iymin = max(iys-2,minym1)
            do iy=iys,iymin,-1
              if (iy*dy > yy) then
                wy = iy - yy/dy
                aper_ey(ix,iy,iza,ia) = tdyi*(phip(ix,iy,iz) - phip(ix,iy+1,iz)-
     &                                  (phip(ix,iy,iz) - aper_volt(ia))/wy)
              elseif (iy*dy == yy) then
                aper_ey(ix,iy,iza,ia) = tdyi*2.*(aper_volt(ia)-phip(ix,iy+1,iz))
              else
                if (ix-midx < 2) then
                  aper_ey(ix,iy,iza,ia) = 2.*aper_ey(ix,iy+1,iza,ia) -
     &                                       aper_ey(ix,iy+2,iza,ia)
                else
                  aper_ey(ix,iy,iza,ia) = 2.*aper_ey(ix+1,iy,iza,ia) -
     &                                       aper_ey(ix+2,iy,iza,ia)
                endif
              endif
            enddo
          enddo

c         --- Now set the lines just to the left and right of the aperture
          if (maxxp1 == maxx + 1) then
c           --- Find min y
            if (minym1 == miny - 1) then
              yy = -sqrt(abs(aper_rad(ia)**2-(maxx*dx-aper_x(ia)+xmmin)**2))+
     &             aper_y(ia) - ymmin
              iymin = int(yy/dy)
            else
              iymin = miny
            endif
c           --- Find max y
            if (maxyp1 == maxy + 1) then
              yy = sqrt(abs(aper_rad(ia)**2-(maxx*dx-aper_x(ia)+xmmin)**2))+
     &             aper_y(ia) - ymmin
              iymax = int(yy/dy) + 1
            else
              iymax = maxy
            endif
c           --- Do the work
            do iy=iymin,iymax
              aper_ey(maxxp1,iy,iza,ia) = 2.*aper_ey(maxxp1-1,iy,iza,ia) -
     &                                       aper_ey(maxxp1-2,iy,iza,ia)
            enddo
          endif
          if (minxm1 == minx - 1) then
c           --- Find min y
            if (minym1 == miny - 1) then
              yy = -sqrt(abs(aper_rad(ia)**2-(minx*dx-aper_x(ia)+xmmin)**2))+
     &             aper_y(ia) - ymmin
              iymin = int(yy/dy)
            else
              iymin = miny
            endif
c           --- Find max y
            if (maxyp1 == maxy + 1) then
              yy = sqrt(abs(aper_rad(ia)**2-(minx*dx-aper_x(ia)+xmmin)**2))+
     &             aper_y(ia) - ymmin
              iymax = int(yy/dy) + 1
            else
              iymax = maxy
            endif
c           --- Do the work
            do iy=iymin,iymax
              aper_ey(minxm1,iy,iza,ia) = 2.*aper_ey(minxm1+1,iy,iza,ia) -
     &                                       aper_ey(minxm1+2,iy,iza,ia)
            enddo
          endif

c       --- End of loop over z planes
        enddo

c       --- End of if checking if aperture is within the grid
        endif

c     --- End of loop over apertures
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timeset_aperture_e = timeset_aperture_e + wtime() - substarttime
!$OMP END MASTER
      return
      end
c===========================================================================
      subroutine bendfieldsol3d()
      use Constant
      use InGen
      use InGen3d
      use InMesh3d
      use Picglb
      use Fields3d
      use Picglb3d

c Does the FFT field solve in bends.

      integer(ISZ):: i,j,k
      real(kind=8):: ccmult,rs,dhdz,cutoff,x,r,phiref,rskm1,rskp1,dxi,dzi
      character(80):: outstr
      dxi = 1./dx
      dzi = 1./dz

      bndfit = 0
 1000 continue
      bndfit = bndfit + 1

c     --- save phi on midplane for error measure
      do k = 0, nz
         do i = 0, nx
            phiprv(i,k) = phi(i,ny/2,k)
         enddo
      enddo

c     --- set multiplier for "jump term" (curvature change term)
      ccmult = 0.
      if (bnjtflag) ccmult = 1.
c     --- loop over slices in y, computing source with bend correction
      do j = 0, ny
c        --- save phi this slice, for future algebra
         do k = 0, nz
            do i = 0, nx
               phisav(i,k) = phi(i,j,k)
            enddo
         enddo
         do i = 0, nx
            phisav(i,-1) = phisav(i,nz-1) 
         enddo
c        --- loop over mesh points in z
         do k = 0, nz-1
            rs = rstar(k)
            rskm1 = rstar(k-1)
            rskp1 = rstar(k+1)
            dhdz = (1./rskp1 - 1./rskm1) *0.5*dzi
c           --- only modify source if actually in a bend, or at entrance/exit
            cutoff = LARGEPOS*1.e-6
            if ((abs(rskm1) < cutoff).or.(abs(rskp1) < cutoff))then
               do i = 1, nx-1
                  x = xmmin + i*dx
                  r = 1./(rs + x)
                  phi(i,j,k) = rho(i,j,k)*rs*r 
     &             + eps0 * ( (phisav(i+1,k) - phisav(i-1,k)) *0.5*dxi*r
     &                  + (phisav(i,k+1) - 2.*phisav(i,k) + phisav(i,k-1))
     &                  * (-2.*x*r + x**2*r**2) * dzi**2
     &                  - ccmult * (phisav(i,k+1) - phisav(i,k-1)) * 0.5*dzi
     &                  * (rs*r)**3 * x * dhdz 
     &                      )
               enddo
            else
               do i = 1, nx-1
                  phi(i,j,k) = rho(i,j,k)
               enddo
            endif
         enddo
      enddo

c     --- call Cartesian field solver
      call vp3d (-1)

c     --- compute error
      bndferr = 0.
      phiref = 0.
      do k = 0, nz
         do i = 0, nx
            bndferr = max( bndferr, abs(phi(i,ny/2,k)-phiprv(i,k)) )
            phiref = max( phiref, phi(i,ny/2,k) )
         enddo
      enddo

c     --- For the parallel version, find global max of errors.
#ifdef MPIPARALLEL
      call parallelmaxrealarray(phiref,1)
#endif
      bndferr = bndferr / dvnz(phiref)
#ifdef MPIPARALLEL
      call parallelmaxrealarray(bndferr,1)
#endif

c     --- for debug, print out the error
      if (bnprflag) then
        write (outstr,9985) it, bndfit, bndferr
 9985   format ("It =",i7," Bent field iteration",i3," Rel Change = ", 1pe12.4)
        call remark(outstr)
      endif

c     --- if error too big, repeat main loop - provided iters remain
      if ( (bndferr > bndftol) .and. (bndfit < bndfitmx) ) go to 1000

c     --- if failure to converge, report the bad news to user
      if (bndferr > bndftol) then
        print*,"*** NONCONVERGENCE in bent field iteration"
        print*,"Relative change = ",bndferr," after ",bndfit," iterations."
      endif

      return
      end
c===========================================================================
      subroutine getphiforparticles()
      use InGen3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles

c     --- Ensure that the phip array is setup properly
      if(solvergeom==XYZgeom .or. solvergeom==XYZgeomMR) then
c       --- zeros are passed in since the arguments are only needed for
c       --- setting up rhoptmp.
        call setupFields3dParticles(0,0,0)
      endif

#ifndef MPIPARALLEL

c     --- If phip is not associated with phi, then copy the data.
      if (.not. associated(phi,phip)) then
        if (nxp == nx .and. nyp == ny .and. nzp == nz) then
          phip = phi
        else
          call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
          call remark("ERROR!! phip and phi are not the same shape!        ")
          call remark("        The self field will not be properly applied!")
          call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        endif
      endif

#else

c     --- Distribute phi among the processes so each has phi in its
c     --- particle domain.
      if(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==Zgeom) then
        call getphiforparticlesrz()
c      else if(solvergeom==Zgeom) then
c        call getphiforparticlesz()
      elseif (solvergeom==XYZgeom) then
        call getphiforparticles3d(nx,ny,nz,phi(0,0,-1),nxp,nyp,nzp,phip)
      end if
#endif

      return
      end
c===========================================================================
      subroutine fieldsol3d(iwhich)
      use GlobalVars
      use Subtimers3d
      use Timers
      use InGen
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Fields3d
      use Fields3dParticles
      use Efields3d
      use LatticeInternal
      use InjectVars
      use GridBoundary3d
      integer(ISZ):: iwhich

c  Field solver for 3d warped Cartesian geometry. 
c  Enter with charge density in rho array, old potential in phi array.
c  Exit with new potential in phi array, and rho unchanged.
c  
c  For field solve in a bend there are two possible cases:
c      1) Call to an SOR field solver.  
c      2) Iteration loop for a perturbative calculation of fields. 
c         The loop is in bendfieldsol3d

      real(kind=8):: timetemp
      real(kind=8):: substarttime,wtime
      integer(ISZ):: izmin,izmax,tnz
      if (lw3dtimesubs) substarttime = wtime()
      timetemp = wtime()

c     --- If no field solve, return 
      if ( fstype == -1) return 

      if (fstype == 3 .or. fstype == 7 .or. fstype == 11 .or. fstype == 12) then
c       --- SOR and multigrid field solvers - deal with bends directly and set
c       --- axial boundary conditions (call to perphi3d is not needed)
c       --- Chombo knows nothing of bends and does its own b.c.'s
        call vp3d (iwhich)

      else if (.not. linbend) then
c       --- If not in bend, call VP3D
        if (fstype == 0 .or. fstype == 1 .or. fstype == 2 .or.
     &      fstype == 4 .or. fstype == 5 .or. fstype == 6 .or.
     &      fstype == 8 .or. fstype == 9) then
          if (bound0  == dirichlet) izmin = 1
          if (bound0  == neumann)   izmin = 0
          if (bound0  == periodic)  izmin = 0
          if (boundnz == dirichlet) izmax = nz-1
          if (boundnz == neumann)   izmax = nz
          if (boundnz == periodic)  izmax = nz
          tnz = izmax - izmin + 1
          call copyarry (rho(0,0,izmin), phi(0,0,izmin), (nx+1)*(ny+1)*(tnz))
        endif
        call vp3d (iwhich)
        if (bound0==periodic) call perphi3d(phi(0,0,-1),nx,ny,nz)

      else
c       --- Call perturbative bent beam field solver
        call bendfieldsol3d
        if (bound0==periodic) call perphi3d(phi(0,0,-1),nx,ny,nz)

      endif

c     --- Distribute phi among the processes so each has phi in its
c     --- particle domain.
      call getphiforparticles()

!$OMP MASTER
      if (lw3dtimesubs) timefieldsol3d = timefieldsol3d + wtime() - substarttime
!$OMP END MASTER
      fstime = fstime + (wtime() - timetemp)
      return
      end
c=============================================================================
      subroutine vp3d (iwhich)
      use Subtimers3d
      use Constant
      use InGen
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Picglb
      use GridBoundary3d
      use Fields3d
      use LatticeInternal
      use Parallel
      integer(ISZ):: iwhich

      real(kind=8):: xlen,ylen,zlen
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Interface to VPOIS3D using variables from database of package 3D

c LatticeInternal is included to have access to the variable linbend.

c fstype
c   0  only 3d fft on rectangular mesh, square pipe
c   1  3d fft and Kz capacity matrix for variable pipe shape
c   2  3d fft and 3d capacity matrix for internal quadrupole conductors
c   3  3d SOR field solve 
c   4  2d fft in transverse plane, tridiag solve in z

      if(solvergeom==RZgeom .or.
     &   solvergeom==XZgeom .or.
     &   solvergeom==XYgeom .or.
     &   solvergeom==Rgeom  .or.
     &   solvergeom==Zgeom) fstype=10
      if(solvergeom==AMRgeom) fstype = 11
      if (solvergeom == XYZgeomMR) fstype = 12

      if (iwhich==0 .or. iwhich==1) call execuser('initfieldsolver')

c Note that the select case statement was giving the f90 compiler on HPUX
c fits and so was replaced with if statements.

      xlen = xmmax - xmmin
      ylen = ymmax - ymmin
      zlen = zmmaxglobal - zmminglobal

      if (fstype == 0) then
          call vpois3d(iwhich,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq(0),kysq(0),kzsq(0),
     &                 attx(0),atty(0),attz(0),filt(1,1),xlen,ylen,zlen,
     &                 nx,ny,izfsmax-izfsmin,nzfull,
     &                 scrtch(0,0),xywork,zwork,0,l2symtry,l4symtry,
     &                 bound0,boundnz,boundxy)
      elseif (fstype == 1) then
          call pipe3df(iwhich,pipeshpe,rho(0,0,izfsmin),phi(0,0,izfsmin),
     &                kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,
     &                nx,ny,izfsmax-izfsmin,nzfull,scrtch(0,0),xywork,
     &                zwork(1,0,0),l2symtry,l4symtry,bound0,boundnz,boundxy)
      elseif (fstype == 2) then
          call vcap3d(iwhich,rho(0,0,izfsmin),phi(0,0,izfsmin),
     &                kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xlen,ylen,zlen,nx,ny,izfsmax-izfsmin,nzfull,
     &                scrtch(0,0),xywork,zwork,
     &                xmmax,zmmin+izfsmin*dz,zgrid,pipeshpe,
     &                l2symtry,l4symtry,bound0,boundnz,boundxy)
      elseif (fstype == 3) then
c          --- Note that by default, izfsmax=nz, and izfsmin=0.
c          --- izfsmin and izfsmax allow field solves on a partial grid.
           call psor3df(iwhich,nx,ny,nz,phi(0,0,-1),rho(0,0,0),
     &                 phi(0,0,-1),rho(0,0,0),rstar, 
     &                 dx,dy,dz,xmmin,ymmin,zmmin,zbeam,zgrid,linbend,
     &                 bound0,boundnz,boundxy,l2symtry,l4symtry,lzerophiedge,
     &                 scrtch(0,0),izfsmin,izfsmax)
      elseif (fstype == 4) then
          if (iwhich == 1 .or. iwhich == 0) then
            call vpois3d(1,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq,kysq,kzsq,
     &                   attx,atty,attz,filt,xlen,ylen,zlen,
     &                   nx,ny,izfsmax-izfsmin,nzfull,scrtch(0,0),
     &                   xywork,zwork,0,l2symtry,l4symtry,
     &                   bound0,boundnz,boundxy)
          endif
          if (iwhich == -1 .or. iwhich == 0) then
            call vpois3d(12,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq,kysq,kzsq,
     &                   attx,atty,attz,filt,xlen,ylen,zlen,
     &                   nx,ny,izfsmax-izfsmin,nzfull,
     &                   scrtch(0,0),xywork,zwork,0,l2symtry,l4symtry,
     &                   bound0,boundnz,boundxy)
            call vpois3d(14,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq,kysq,kzsq,
     &                   attx,atty,attz,filt,xlen,ylen,zlen,
     &                   nx,ny,izfsmax-izfsmin,nzfull,
     &                   scrtch(0,0),xywork,zwork,0,l2symtry,l4symtry,
     &                   bound0,boundnz,boundxy)
            call vpois3d(13,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq,kysq,kzsq,
     &                   attx,atty,attz,filt,xlen,ylen,zlen,
     &                   nx,ny,izfsmax-izfsmin,nzfull,
     &                   scrtch(0,0),xywork,zwork,0,l2symtry,l4symtry,
     &                   bound0,boundnz,boundxy)
          endif

      elseif (fstype == 5) then
c         --- General capacity matrix solver in kz space
          call capmatkz3d(iwhich,phi(0,0,0),rho,kxsq,kysq,kzsq,attx,atty,attz,
     &                    filt,xlen,ylen,zlen,
     &                    nx,ny,izfsmax-izfsmin,nzfull,dx,dy,dz,
     &                    xmmin,ymmin,zmmin+izfsmin*dz,scrtch(0,0),
     &                    xywork,zwork(1,0,0),l2symtry,l4symtry,
     &                    bound0,boundnz,boundxy)
      elseif (fstype == 6) then
c         --- General capacity matrix solver
          call capmat3df(iwhich,phi(0,0,0),rho,kxsq,kysq,kzsq,attx,atty,attz,
     &                  filt,xlen,ylen,zlen,
     &                  nx,ny,izfsmax-izfsmin,nzfull,dx,dy,dz,
     &                  xmmin,ymmin,zmmin+izfsmin*dz,scrtch(0,0),
     &                  xywork,zwork(1,0,0),l2symtry,l4symtry,
     &                  bound0,boundnz,boundxy)

      elseif (fstype == 7) then
          call multigrid3df(iwhich,nx,ny,izfsmax-izfsmin,nzfull,dx,dy,dz,
     &            phi(0,0,izfsmin-1),rho(0,0,izfsmin),rstar(izfsmin-1),
     &            linbend,bound0,boundnz,boundxy,l2symtry,l4symtry,
     &            xmmin,ymmin,zmmin,zbeam,zgrid)

#ifdef MPIPARALLEL
      elseif (fstype == 8) then
c         --- Experimental parallel solver
c         --- Combines local tridiag solves with iteration to exchange
c         --- boundary data.
          call paralleltridiag(phi(0,0,izfsmin-1),kxsq,kysq,kzsq,
     &                         attx,atty,attz,filt,xlen,ylen,zlen,
     &                         nx,ny,izfsmax-izfsmin,
     &                         scrtch(0,0),xywork,zwork,l2symtry,l4symtry,
     &                         bound0,boundnz,boundxy)
      elseif (fstype == 9) then
c         --- Experimental parallel solver based on the method
c         --- put forward by Lantz
          call lantzsolver(iwhich,phi(0,0,izfsmin-1),kxsq,kysq,kzsq,
     &                     attx,atty,attz,filt,xlen,ylen,zlen,
     &                     nx,ny,izfsmax-izfsmin,
     &                     scrtch(0,0),xywork,zwork,l2symtry,l4symtry,
     &                     bound0,boundnz,boundxy)

#endif
      elseif (fstype == 10) then
c         --- RZ full-multigrid solver
       call multigridrzf(iwhich,phi(0,0,-1),rho(0,0,0),nx,nz)
      elseif (fstype == 11) then
c         --- Chombo AMR  full-multigrid solver
        call cho_solve3d(iwhich,nx,ny,nz,nzfull,dx,dy,dz,
     &                   l2symtry,l4symtry,xmmin,ymmin,zmmin,zmminglobal)
      elseif (fstype == 12) then
        if (iwhich <= 0) then
          call execuser('fieldsolMR')
        endif

      endif

!$OMP MASTER
      if (lw3dtimesubs) timevp3d = timevp3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
