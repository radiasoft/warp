#include "top.h"
c=============================================================================
c@(#) File W3D.M, version $Revision: 3.43 $, $Date: 2002/03/18 21:36:19 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package W3D of code WARP
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      subroutine w3dinit
      use Subtimers
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Called at first reference to package (not nec. a "run" etc.).

      call w3dvers (STDOUT)

      if (lw3dtimesubs) timew3dinit = timew3dinit + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine w3dvers (iout)
      use Subtimers
      use W3Dversion
      integer(ISZ):: iout
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
c  Echoes code version, etc. to output files as they're created
      call printpkgversion(iout,"Particle package W3D",versw3d)
      if (lw3dtimesubs) timew3dvers = timew3dvers + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine w3dgen()
      use Subtimers
      use Ch_var
      use Constant
      use InGen
      use InGen3d
      use InDiag
      use InPart
      use InPart3d
      use InMesh3d
      use Fields3d
      use Efields3d
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
      use Lab_Moments
      use Moments
      use Damped_eom
      use Hist
      use Pspwork
      use InjectVars
      use FRZmgrid

c  Invoked by the GENERATE command, it sets up the problem
c  This routine allots all of the neccesary dynamic arrays, calls the
c  particle loader and does the initial load onto the charge density
c  mesh, initializes arrays for the field solver and sets the mesh arrays,
c  does the initial field solve, and sets up other arrays that are needed.


      integer(ISZ):: i,j,k,is,ipmin,ip,iwin
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Announce that we're starting up

      call remark(" ***  particle simulation package W3D generating")

c  Put ny to zero for RZ geometry

      if(solvergeom==RZgeom) then
        ny=0
        fstype=10
      endif

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

c  Create the dynamic arrays in Z_arrays; set the z mesh

      if (nzzarr == 0) nzzarr = nz
      call gchange("Z_arrays", 0)
      if (zzmin == 0.) zzmin = zmmin
      if (zzmax == 0.) zzmax = zmmax
      dzz = (zzmax - zzmin)/nzzarr
      dzzi = 1./dzz
      do k = 0, nzzarr
         zplmesh(k) = zzmin + k*dzz
      enddo

c  Re-size the dynamic arrays for the lattice (scan for true length, first).
c  Also allocate the dynamic internal lattice arrays.

      call remark(" ---  Resetting lattice array sizes")
      call resetlat
      if (nzlmax == 0) nzlmax = nz
      if (nzl    == 0 .or. nzl > nzlmax) nzl = nzlmax
      call gchange("LatticeInternal", 0)
      if (zlmin == 0.) zlmin = zmmin
      if (zlmax == 0.) zlmax = zmmax
      dzl = (zlmax - zlmin)/nzl
      dzli = 1./dzl
      do k = 0, nzl
         zlmesh(k) = zlmin + k*dzl
      enddo

c  Create the dynamic arrays for fields, contour plot workspace
c  Also set vbeamfrm, adjust x and ymmin for symmetries and set
c  so that by default izfsmax = nz (for full field solve).

      if (nzfull == 0) nzfull = nz
      nmxy  = max(nx,ny)
      nmxyz = max(nx,ny,nzfull)
      call gallot("Fields3d", 0)
      vbeamfrm = vbeam
      if (l2symtry) then
        ymmin = 0.
      elseif (l4symtry) then
        xmmin = 0.
        ymmin = 0.
      endif
      if(solvergeom==RZgeom) then
        l2symtry=.false.
        l4symtry=.false.
        xmmin = 0.
        ymmin = 0.
      endif
      if (izfsmax == 0) izfsmax = nz

c Allocate space for the selfe array, which contains the finite differenced
c phi (Ex, Ey, Ez) on the mesh.
      if (efetch == 3) then
        nx_selfe = nx
        ny_selfe = ny
        nz_selfe = nz
        call gallot("Efields3d",0)
      endif

      if(solvergeom==RZgeom) then
        nx_selfe = nx
        ny_selfe = 0
        nz_selfe = nz
        call gallot("Efields3d",0)
      endif

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
      if(solvergeom==XYZgeom) then
        dy = (ymmax - ymmin) / ny
      elseif(solvergeom==RZgeom) then
        dy = dx
      endif
      dz = (zmmax - zmmin) / nz
      do i = 0, nx
         xmesh(i) = i*dx + xmmin
      enddo
      do j = 0, ny
         ymesh(j) = j*dy + ymmin
      enddo
      do k = 0, nz
         zmesh(k) = k*dz + zmmin
      enddo

c Calculate location of axis in mesh, the term dx*1.e-5 acts as fuzz
      ix_axis = nint(-xmmin/dx)
      if(solvergeom==XYZgeom) then
        iy_axis = nint(-ymmin/dy)
      elseif(solvergeom==RZgeom) then
        iy_axis = 0
      endif
      iz_axis = nint(-zmmin/dz)

c  Create dynamic arrays for 3D multipole moments of the electrostatic
c  potential
      call gallot("Multipole", 0)

c  Create the dynamic arrays for particles (set npmax to an estimated
c  length for now, for those loading schemes that don't actually
c  use a user-set npmax directly)

      call remark(" ---  Allocating space for particles")
      if (xrandom == "grid") npmax = nxstripe*nystripe*nzstripe
      if (xrandom == "fibonacc") npmax = nfibgrps*fibg1

      call alotpart

c  Load the particles, calculate the charge density

      call remark(" ---  Loading particles")
      call stptcl3d
      call setgamma(lrelativ)
      call remark(" ---  Setting charge density")
      call zeroarry(rho(0,0,0),(nx+1)*(ny+1)*(nz+1))
      call zeroarry(curr(0),nzzarr+1)
c     rho = 0.
c     curr = 0.
      do is=1,ns
         do ipmin = ins(is), ins(is) + nps(is) - 1, NPARPGRP
           ip = min(NPARPGRP, ins(is)+nps(is)-ipmin)
           if(solvergeom==XYZgeom) then
              call setrho3d (rho(0,0,0), ip, xp(ipmin), yp(ipmin), zp(ipmin),
     &                       zgrid, uzp(ipmin), sq(is), sw(is), depos)
            elseif(solvergeom==RZgeom) then
              call rhoweightrz(xp(ipmin),yp(ipmin),zp(ipmin),ip,
     &                         sq(is)*sw(is),rho(0,0,0),nx,nz,dx,dz,zmmin)
            end if
         enddo
      enddo
      call perrho3d(rho(0,0,0),nx,ny,nz,periinz)

c  Set up for injection
      if (inject .gt. 0) call injctint

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

      call remark(" ---  Allocating Z_Moments")
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
      if (iflabwn .ne. 0) then
c       --- Note that the original code used abs(nhist), but it was giving
c       --- the HPUX f90 compiler fits. The negative nhist feature was
c       --- made obsolete anyway.
        if (itlabwn == 0) itlabwn = max(1,nhist)
        if (ntlabwn == 0) then
          if (vbeam*dt .ne. 0.) then
            ntlabwn = (zmmax-zmmin)/(vbeam*dt*itlabwn) +.5
          else
            call remark("WARNING: Unable to calculate lab windows array sizes")
            call remark("since either vbeam and/or dt are zero.")
            ntlabwn = 1
          endif
        endif
        if (iflabwn .ne. 0) call gallot("Lab_Moments", 0)
      endif

c  Create the scratch arrays for phase space plots (permanent, for now)
c  and set limits for plots

      call remark(" ---  Allocating scratch space for plots")
      if (npsplt == 0) npsplt = NPARPGRP-1
      call gallot("Pspwork", 0)
      if (xplmin == 0.) xplmin = -xmmax
      if (xplmax == 0.) xplmax = xmmax
      if (yplmin == 0.) yplmin = -ymmax
      if (yplmax == 0.) yplmax = ymmax
      if (zplmin == 0.) zplmin = zmmin
      if (zplmax == 0.) zplmax = zmmax

c  Choose which particles will be plotted in "subset" plots

      call psubsets

c  Setup history mechanism

      if (nhist .gt. 0) then
c       --- create the dynamic arrays for history data; set pointer into them
        call remark(" ---  Allocating history arrays")
        if (lenhist == 0) lenhist = min ( nt/nhist + 1, 100)
        call gallot("Hist", 0)
        jhist = -1
      elseif (nhist .lt. 0) then
c       --- call interpreter routine to setup hst package
c       --- setup_hst is in bas.wrp
        call execuser("setup_hst")
      endif

c  Print interesting things to plot file and teletype

      call prntpara(dx,dy,dz,lprntpara)
      call prntpa3d(lprntpara)

c  Initial call to fieldsolver in order to initialize attx, kxsq, etc.

      call vp3d (1)

c  Initial fieldsolve, diagnostics

      call step3d ("w3dgen")

c  Set up Eears of z

      call seteears()

      if (lw3dtimesubs) timew3dgen = timew3dgen + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine w3dexe()
      use Subtimers
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
      time = time + dt
      if (lgridqnt) then
c       --- zgrid is integer number of dz's
        zgridprv = int((zbeam)/dz + .5)*dz
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

c  call the routine that does the actual work

      call step3d ("w3dexe")

c  Have we reached the end of the run or run out of particles?

      if ( lfinishd ) then
         call remark("w3dexe: problem completed.")
         if (lw3dtimesubs) timew3dexe = timew3dexe + wtime() - substarttime
         return
      elseif (nplive .le. 0) then
c Check if there are any live particles left.
c If on first time step, set nplive to one and continue running since if using
c injection, there may not be any particles yet but still want to run.
         if (it == 1) then
           nplive = 1
         else
           print*," *** W3DEXE: stopping, nplive =",nplive
           if (lw3dtimesubs) timew3dexe = timew3dexe + wtime() - substarttime
           return
         endif
      endif

c  Continue to run.
      if (lw3dtimesubs) timew3dexe = timew3dexe + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine w3dfin()
      use Subtimers
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
      use Pspwork
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

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
      call gfree ("Pspwork")
      call gfree ("Z_arrays")

      if (lw3dtimesubs) timew3dfin = timew3dfin + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine step3d (caller)
      use Subtimers
      use Constant
      use InGen
      use InGen3d
      use InDiag
      use InPart
      use InMesh3d
      use Fields3d
      use Efields3d
      use Io
      use Particles
      use Picglb
      use Picglb3d
      use InjectVars
      use LatticeInternal
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
      llabwn   = dolabwn()
      lspecial = (lfinishd .or. lalways .or. lseldom .or. ldump .or. lmoments
     &            .or. llabwn .or. llast .or. (it .eq. 0) .or. allspecl)

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
      if (efetch == 3)
     &  call getselfe3d(phi(0,0,-1),nx,ny,nz,selfe(1,0,0,0),
     &                  nx_selfe,ny_selfe,nz_selfe,dx,dy,dz)
c  Set the potential near the emitting surface.
      if (inject > 0) then
        if(solvergeom==XYZgeom) then
          call getinj_phi(nx,ny,nz,phi(0,0,0),dx,dy,dz,xmmin,ymmin)
        elseif(solvergeom==RZgeom) then
          call getinj_phi_rz(nx,nz,phi(0,0,0),dx,dz,xmmin,ymmin)
        endif
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
      if (inject > 0 .and. caller == "w3dexe")
     &  call inject3d(3)

c  Gather moments used in diagnostics at "special" timesteps only. 
c  Compute mean beam z velocity from current and line charge density 
c  on a 1-d mesh.  Also, calculate the electrostatic energy (getese3d), 
c  electrostatic potential on axis (sphiax3d), and the axial electric 
c  field (sezax3d).  

      if (lspecial) then 
        call getvzofz
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
     &  call minidiag (it,time,lmoments)

      if (lw3dtimesubs) timestep3d = timestep3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine exteb3d(np,xp,yp,zp,uzp,gaminv,dtl,dtr,
     &                   bx,by,bz,ex,ey,ez,m,q,bendres,bendradi,gammabar,dt)
      use Subtimers
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

c     --- Set E to zero for dead particles
      do ip=1,np
        if (uzp(ip) == 0) then
          ex(ip) = 0.
          ey(ip) = 0.
          ez(ip) = 0.
        endif
      enddo

c     --- Accumulate time for applying the fields of the lattice
      latticetime = latticetime + (wtime() - timetemp)

      if (lw3dtimesubs) timeexteb3d = timeexteb3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine othere3d (np,xp,yp,zp,zbeam,zimax,zimin,straight,
     &                     ifeears,eears,eearsofz,dzzi,nzzarr,zzmin,
     &                     dedr,dexdx,deydy,dbdr,ex,ey,ez,bx,by,bz)
      use Subtimers
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
          if ((zp(ip)-zbeam) .gt. zs) then
            ez(ip) = ez(ip) + eears*(zp(ip) - zbeam - zs)
          elseif ((zp(ip)-zbeam) .lt. -zs) then
            ez(ip) = ez(ip) + eears*(zp(ip) - zbeam + zs)
          endif
        enddo
      endif

c  axial confining ears for finite beam as a function of z

      if (ifeears == 2) then
        do ip=1,np
          iz = (zp(ip) - zbeam - zzmin)*dzzi
          if (0 .le. iz .and. iz .le. nzzarr) then
            wz = (zp(ip) - zbeam - zzmin)*dzzi - iz
            ez(ip) = ez(ip) + eearsofz(iz)*(1.-wz) + eearsofz(iz+1)*wz
          endif
        enddo
      endif

c  uniform focusing forces

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

      if (lw3dtimesubs) timeothere3d = timeothere3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine getese3d
      use Subtimers
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
        if (lw3dtimesubs) timegetese3d = timegetese3d + wtime() - substarttime
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

      do iz=izfsmin,izfsmax-1
        phiextun = 0.
        if (ifeears == 1) then
          if (iz*dz+zz .lt. -zm) phiextun=eears*.5*(iz*dz+zz+zm)**2
          if (iz*dz+zz .gt.  zm) phiextun=eears*.5*(iz*dz+zz-zm)**2
        elseif (ifeears == 2) then
c         --- This is not quite right for the parallel version since
c         --- the longitudinal extent of eearsofz (the same as the extent
c         --- of the particles) will not necessarily overlap the extent of
c         --- the field solver.
          izb = (iz*dz+zz-zzmin)*dzzi
          wzb = (iz*dz+zz-zzmin)*dzzi - izb
          if (0 .le. iz .and. iz .lt. nzzarr)
     &      phiextun = eearsofz(izb)*(1. - wzb) + eearsofz(izb+1)*wzb
        endif
        do iy=0,ny
          do ix=0,nx
            scrtch(ix,iy) = scrtch(ix,iy) + rho(ix,iy,iz)*(phi(ix,iy,iz)
     &                       - dedr*.5*((ix*dx+xmmin)**2 + (iy*dy+ymmin)**2)
     &                       - dexdx*.5*(ix*dx+xmmin)**2 
     &                       - deydy*.5*(iy*dy+ymmin)**2
     &                       + phiextun)
          enddo
        enddo
      enddo

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

      if (lw3dtimesubs) timegetese3d = timegetese3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine gtlchg3d
      use Subtimers
      use Picglb
      use InGen3d
      use Picglb3d
      use InMesh3d
      use Fields3d
      use Z_arrays
      use InDiag3d

c  Calculates the line charge density from rho.

      real(kind=8):: zz,wzg,dzi
      integer(ISZ):: ix,iy,iz,izg
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      if (.not. lgtlchg3d) return

      dzi = 1./dz

c     conversion factor to go from grid frame to beam frame
      zz = zgrid + zmmin - zzmin - zbeam

c First, sum up the charge density in each axial plane.

c If no symmetries...
      if (.not. (l2symtry .or. l4symtry)) then

c       --- loop over z planes
        do iz=0,nz

c         --- start with rho along iy=0 line
          do ix=0,nx
            scrtch(iz,ix) = rho(ix,0,iz)
          enddo

c         --- add rho along lines in y
          do iy=1,ny
            do ix=0,nx
              scrtch(iz,ix) = scrtch(iz,ix) + rho(ix,iy,iz)
            enddo
          enddo

c         --- add up sums from each line in y
          do ix=1,nx
            scrtch(iz,0) = scrtch(iz,0) + scrtch(iz,ix)
          enddo


        enddo

      elseif (l2symtry) then

c  special loop for 2-fold to account for only one side of rho being stored

c       --- loop over z planes
        do iz=0,nz

c         --- start with rho along iy=0 line
          do ix=0,nx
            scrtch(iz,ix) = rho(ix,0,iz)
          enddo

c         --- add rho along lines in y
          do iy=1,ny
            do ix=0,nx
              scrtch(iz,ix) = scrtch(iz,ix) + 2.*rho(ix,iy,iz)
            enddo
          enddo

c         --- add up sums from each line in y
          do ix=1,nx
            scrtch(iz,0) = scrtch(iz,0) + scrtch(iz,ix)
          enddo

        enddo

      elseif (l4symtry) then

c  special loop for 4-fold to account for only one quadrant of rho being
c  stored

c       --- loop over z planes
        do iz=0,nz

c         --- start with rho along iy=0 line
          scrtch(iz,0) = rho(0,0,iz)
          do ix=1,nx
            scrtch(iz,ix) = 2.*rho(ix,0,iz)
          enddo

c         --- add rho along lines in y
          do iy=1,ny
            scrtch(iz,0) = scrtch(iz,0) + 2.*rho(0,iy,iz)
            do ix=1,nx
              scrtch(iz,ix) = scrtch(iz,ix) + 4.*rho(ix,iy,iz)
            enddo
          enddo

c         --- add up sums from each line in y
          do ix=1,nx
            scrtch(iz,0) = scrtch(iz,0) + scrtch(iz,ix)
          enddo

        enddo

      endif

c Then, linearly interpolate the data into the beam frame array, linechg.

      do iz=0,nzzarr

        izg = (iz*dzz - zz)*dzi
        wzg = (iz*dzz - zz)*dzi - izg

        if (0 .le. izg .and. izg .le. nz) then
          linechg(iz) = scrtch(izg,0)*dx*dy*(1. - wzg)
        else
          linechg(iz) = 0.
        endif
        if (0 .le. izg+1 .and. izg+1 .le. nz) then
          linechg(iz) = linechg(iz) + scrtch(izg+1,0)*dx*dy*wzg
        endif

      enddo

      if (lw3dtimesubs) timegtlchg3d = timegtlchg3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine seteears()
      use Subtimers
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
#ifndef MPIPARALLEL
      integer(ISZ):: nzla,nzlb,nzma,nzmb,iz
      real(kind=8):: zs,zl,eearsprs
#endif
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      zm = (zimax - zimin)*(1. - straight)*0.5
      emitlong = zm*2.*vthz/vbeam

c     set to be linear in z (eears is actually slope of Eears)
      if (ifeears == 1) then
        eears = eears*(- 2*ibeam/(vbeam*zm**2)*gfactor/(2.*Pi*eps0)
     &                  - sm(1)*vbeam**2*emitlong**2/(sq(1)*zm**4))

      elseif (ifeears == 2) then
c     --- For Parallel slave, message passing must be done so seperate
c     --- routine is used.  Note that second half of ifelse macro
c     --- includes the next 30 lines.
#ifdef MPIPARALLEL
        call seteears_slave(zwork(1,0,0))
#else
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
        if (emitlong .ne. 0.) then
          eearsprs = - sm(1)*vbeam**2*emitlong**2/(sq(1)*zm**4)
          do iz=0,nzma
            eearsofz(iz) = eearsofz(iz) + eearsprs*(zplmesh(iz) + zs)
          enddo
          do iz=nzmb,nzzarr
            eearsofz(iz) = eearsofz(iz) + eearsprs*(zplmesh(iz) - zs)
          enddo
        endif
#endif
      endif

      if (lw3dtimesubs) timeseteears = timeseteears + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine padvnc3d(center)
      use Subtimers
      use InMesh3d
      use InGen
      use InGen3d
      use InPart
      use InGaps
      use Lattice
      use LatticeInternal
      use Particles
      use Fields3d
      use Efields3d
      use Picglb
      use Picglb3d
      use Beam_acc
      use Z_arrays
      use Damped_eom
      use InjectVars
      use Apertures
      character(*):: center

c  Advances the particles position and velocity according to CENTER,
c  and also loads RHO at the new time level.


      integer(ISZ):: is,ip,ipmin,i
      real(kind=8):: uxpadv,uypadv,uzpadv
      real(kind=8):: zpo(NPARPGRP)
      real(kind=8):: ex(NPARPGRP), ey(NPARPGRP), ez(NPARPGRP)
      real(kind=8):: bx(NPARPGRP), by(NPARPGRP), bz(NPARPGRP)
      real(kind=8):: uxpo(NPARPGRP), uypo(NPARPGRP), uzpo(NPARPGRP)
      real(kind=8):: bendres(NPARPGRP), bendradi(NPARPGRP)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Zero the bend radius and residence arrays
c     call zeroarry (bendres,NPARPGRP)
c     call zeroarry (bendradi,NPARPGRP)
      bendres = 0.
      bendradi = 0.

c  Zero rho only if particles are advanced
      if (center == "fullv" .or. center == "halfv" .and. .not. laccumulate_rho)
     &  call zeroarry(rho(0,0,0),(nx+1)*(ny+1)*(nz+1))
c    &  rho = 0.

c  Zero curr if center is synchv or if using gaps
      if (center == "synchv" .or. ifgap .and. .not. laccumulate_rho)
     &  call zeroarry(curr(0),nzzarr+1)
c    &  curr = 0.

c  Zero the moments if center is synchv
      if (center == "synchv" .or. center == "gen") then
        call getzmmnt(1,xp(1),yp(1),zp(1),uxp(1),uyp(1),uzp(1),gaminv(1),
     &                sq(1),sm(1),sw(1),dt*0.5,1,
     &                nplive,uxpo,uypo,uzpo,1,ns)
      endif

c  Set the transverse E fields near any defined apertures.
      if (napertures > 0) call set_aperture_e()

c  Loop over species
      do is=1,ns

c  Loop over particle blocks; move each block separately
         do ipmin = ins(is), ins(is) + nps(is) - 1, NPARPGRP
            ip = min(NPARPGRP, ins(is)+nps(is)-ipmin)

c           --- Obtain the self-field from the electrostatic potential
            if(solvergeom==XYZgeom) then
              call sete3d(phi(0,0,-1),selfe(1,0,0,0),ip,
     &                    xp(ipmin),yp(ipmin),zp(ipmin),
     &                    zgridprv,xmmin,ymmin,zmmin,dx,dy,dz,nx,ny,nz,
     &                    efetch,ex,ey,ez,l2symtry,l4symtry)
            elseif(solvergeom==RZgeom) then
             if(is==1.and.ipmin==ins(is)) then
               call fieldweightrz(xp(ipmin),yp(ipmin),zp(ipmin),ex,ey,ez,
     &                            ip,phi(0,0,-1),selfe(1,0,0,0),nx,nz,dx,dz,
     &                            zmmin,.true.)
             else
               call fieldweightrz(xp(ipmin),yp(ipmin),zp(ipmin),ex,ey,ez,
     &                            ip,phi(0,0,-1),selfe(1,0,0,0),nx,nz,dx,dz,
     &                            zmmin,.false.)
             end if
            endif
c           --- Get E field for particles near the  injection surface.
            if (inject > 0) then
              call inj_sete3d(ip,ipmin,xmmin,ymmin,zmmin,dx,dy,dz,ex,ey,ez)
            endif
c           --- Get transverse E field for particles near any defined apertures.
            if (napertures .gt. 0)
     &        call sete3d_aperture(phi(0,0,0),ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                             zgridprv,xmmin,ymmin,zmmin,zmmax,dx,dy,dz,
     &                             nx,ny,nz,ex,ey,l2symtry,l4symtry)
c           --- Scale the self E-field to get the lowest order relativistic
c           --- correction.
            if (relativity == 1)
     &        call sete3d_relativity(ip,ex,ey,vbeam)

c           --- Zero out B field arrays
            do i = 1,ip
               bx(i) = 0.
               by(i) = 0.
               bz(i) = 0.
            enddo
c           --- Compute lag average for experimental damping algorithm
            if (eomdamp .ne. 0.)
     &        call edamp (eomdamp,it,itdamp,center,ip,ex,ey,ez,
     &                    exold(ipmin),eyold(ipmin),ezold(ipmin),
     &                    exlag(ipmin),eylag(ipmin),ezlag(ipmin))

c           --- FULLV
            if (center == "fullv") then
               do i = 1,ip
                 zpo(i) = zp(ipmin+i-1) 
               enddo 
c              --- Obtain bend radii and residence factors
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),bendres,
     &                       bendradi,-dt*0.5,dt*0.5,.false.)
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
     &                       gaminv(ipmin),-dt*0.5,dt*0.5,
     &                       bx,by,bz,ex,ey,ez,sm(is),sq(is),
     &                       bendres,bendradi,gammabar,dt)
c              --- Correction to z on entry/exit to accelerator gap
               call zgapcorr(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       -dt*0.5, dt*0.5, dt, sm(1), sq(1), time)
c              --- First half-electric field increment to momenta
               call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                       ex, ey, ez, sq(is), sm(is), dt*0.5)
c              --- Magnetic field increment to momenta
               call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bx, by, bz, sq(is), sm(is), dt, ibpush)
c              --- Final half-electric field increment to momenta
               call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                       ex, ey, ez, sq(is), sm(is), dt*0.5)
c              --- Advance relativistic Gamma factor
               call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                       gamadv,lrelativ)
c              --- Position advance
               call xpush3d (ip, xp(ipmin), yp(ipmin), zp(ipmin), uxp(ipmin),
     &                       uyp(ipmin),uzp(ipmin),gaminv(ipmin),dt)
c              --- Bend residence factor over POSITION advance step
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bendres,bendradi,-dt,0.,.false.)
c              --- Correct position advance for warped mesh effect
               call zbendcor(ip,xp(ipmin),zp(ipmin),uxp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),dt,bendres,bendradi,bends,linbend)
c              --- Correct position advance for slanted dipole entry/exit
               call sledgcor(ip,xp(ipmin),yp(ipmin),ip,zp(ipmin),zpo,
     &                       uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       sm(is),sq(is),.false.)

c           --- HALFV
            elseif (center == "halfv") then
               do i = 1,ip
                 zpo(i) = zp(ipmin+i-1) 
               enddo 
c              --- Obtain bend radii and residence factors
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bendres,bendradi,0.,dt*0.5,.false.)
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
     &                       gaminv(ipmin),0.,dt*0.5,bx,by,bz,
     &                       ex,ey,ez,sm(is),sq(is),bendres,bendradi,
     &                       gammabar,dt)
c              --- Correction to z on entry/exit to accelerator gap
               call zgapcorr(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       0., dt*0.5, dt, sm(1), sq(1), time)
c              --- Magnetic field increment to momenta
               call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bx, by, bz, sq(is), sm(is), dt*0.5, ibpush)
c              --- Final half-electric field increment to momenta
               call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                       ex, ey, ez, sq(is), sm(is), dt*0.5)
c              --- Advance relativistic Gamma factor
               call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                       gamadv,lrelativ)
c              --- Position advance
               call xpush3d (ip, xp(ipmin), yp(ipmin), zp(ipmin), uxp(ipmin),
     &                       uyp(ipmin),uzp(ipmin),gaminv(ipmin),dt)
c              --- Bend residence factor over POSITION advance step
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bendres,bendradi,-dt,0.,.false.)
c              --- Correct position advance for warped mesh effect
               call zbendcor(ip,xp(ipmin),zp(ipmin),uxp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),dt,bendres,bendradi,bends,linbend)
c              --- Correct position advance for slanted dipole entry/exit
               call sledgcor(ip,xp(ipmin),yp(ipmin),ip,zp(ipmin),zpo,
     &                       uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       sm(is),sq(is),.false.)

c           --- SYNCHV or GEN
            elseif (center == "synchv" .or. center == "gen") then
c              --- Copy 'old' velocity into uxpo, uypo, and uzpo
               do i=1,ip
                  uxpo(i) = uxp(ipmin+i-1)
                  uypo(i) = uyp(ipmin+i-1)
                  uzpo(i) = uzp(ipmin+i-1)
               enddo
c              --- Obtain bend radii and residence factors
               call getbend (ip,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),bendres,
     &                       bendradi,-dt*0.5,0.,.false.)
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
     &                       gaminv(ipmin),-dt*0.5,0.,
     &                       bx,by,bz,ex,ey,ez,sm(is),sq(is),
     &                       bendres,bendradi,gammabar,dt)
c              --- Half electric field increment to momenta
               call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                       ex, ey, ez, sq(is), sm(is), dt*0.5)
c              --- Half magnetic field increment to momenta
               call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                       bx, by, bz, sq(is), sm(is), dt*0.5, ibpush)
               if (center == "gen") then
c                --- Reset uxp to uxpo, set uxpo to half step backward
c                --- for interpolation in moments calculation
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
c              --- Calculate moments over particles, now that we're sync'd
               call getzmmnt(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                       uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),sq(is),sm(is),sw(is),dt*0.5,2,nplive,
     &                       uxpo,uypo,uzpo,is,ns)
            endif
         enddo

      enddo

c----------------------------------------------------------------------------
c     ---  Do final stuff for moments calculation
      if (center == "synchv" .or. center == "gen") then
        call getzmmnt(1,xp(1),yp(1),zp(1),uxp(1),uyp(1),uzp(1),gaminv(1),
     &                sq(1),sm(1),sw(1),
     &                dt*0.5,3,nplive,uxpo,uypo,uzpo,1,ns)
        call getlabwn()
c       --- load rho diagnostic quantities
        call rhodia3d
      endif

      if (center == "synchv" .or. center == "gen" .or. ifgap) then
c       --- Calculate current
        do is=1,ns
          do ipmin = ins(is), ins(is) + nps(is) - 1, NPARPGRP
            ip = min(NPARPGRP, ins(is)+nps(is)-ipmin)
            call setcurr(curr(0),ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                   sq(is),sw(is),zbeam,dzz,zzmin)
          enddo
        enddo
        call percurr(curr(0),nzzarr,periinz)
      endif

c----------------------------------------------------------------------------
      if (center == "fullv" .or. center == "halfv") then
c       --- Inject more particles. This addes particles and alters ins and nps.
        if (inject .gt. 0) call inject3d(1)

c       --- Apply particle boundary conditions
        call zpartbnd(zmmax,zmmin,dz,zgrid)
        do is=1,ns
          if (stickyxy .and. nps(is) > 0)
     &       call stckxy3d(nps(is),xp(ins(is)),xmmax,xmmin,dx,yp(ins(is)),ymmax,
     &                     ymmin,dy,zp(ins(is)),zmmin,dz,uxp(ins(is)),
     &                     uyp(ins(is)),uzp(ins(is)),zgrid,zbeam,nzzarr,
     &                     prwallz(0),prwallxz(0),prwallyz(0),prwelips(0),
     &                     l2symtry,l4symtry,lostpars(0),zzmin,dzzi)
        enddo

c       --- Collect charge density (3d) and current (1d)
        do is=1,ns
          do ipmin = ins(is), ins(is) + nps(is) - 1, NPARPGRP
            ip = min(NPARPGRP, ins(is)+nps(is)-ipmin)
            if(solvergeom==XYZgeom) then
              call setrho3d (rho(0,0,0), ip, xp(ipmin), yp(ipmin), zp(ipmin),
     &                       zgrid, uzp(ipmin), sq(is), sw(is), depos)
            elseif(solvergeom==RZgeom) then
              call rhoweightrz(xp(ipmin),yp(ipmin),zp(ipmin),ip,
     &                         sq(is)*sw(is),rho(0,0,0),nx,nz,dx,dz,zmmin)
            endif
          enddo
        enddo
        if(solvergeom==RZgeom) call rhobndrz(rho(0,0,0),nx,nz)

c       --- Make rho periodic if periodic in z (if not accumulating rho)
        if (.not. laccumulate_rho)
     &    call perrho3d(rho(0,0,0),nx,ny,nz,periinz)

        if (lbeamcom) then
c         --- Set zbeam so that it follows the center of mass of the beam.
          call zbeamcom(zbeam)
        else
c         --- Advance beam frame location using the nominal beam frame velocity.
          zbeam = zgrid
        endif

c       --- zgridprv needs to be updated for the "synchv" step
c       --- Note that zgridprv is also set at the beginning of w3dexe
        zgridprv = zgrid
      endif


      if (lw3dtimesubs) timepadvnc3d = timepadvnc3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine perphi3d(phi,nx,ny,nz)
      use Subtimers
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      integer(ISZ):: ix,iy
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Sets the slices on the exterior of phi for periodicity
c  sets slice at -1 equal to the slice at nz-1
c  sets slice at nz+1 equal to the slice at 1

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

      if (lw3dtimesubs) timeperphi3d = timeperphi3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine perrho3d(rho,nx,ny,nz,periinz)
      use Subtimers
      integer(ISZ):: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      logical(ISZ):: periinz
      integer(ISZ):: ix,iy
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Sums the first and last slices of rho for periodicity
c  and puts the result into both slices.


c Parallel version makes call to special routine perrho3d_slave
#ifdef MPIPARALLEL

c This sums rho in the overlapping region of neighboring processes.
      call perrho3d_slave(rho(0,0,0),nx,ny,nz,periinz)

c  For parallel version, each processor sends rho to neighboring processors
c  whose field solve region overlap its particle region.
      call getrhoforfieldsolve(nx,ny,nz,rho(0,0,0))

#else

      if (periinz) then

        do iy=0,ny
          do ix=0,nx
            rho(ix,iy,0)  = rho(ix,iy,0) + rho(ix,iy,nz)
            rho(ix,iy,nz) = rho(ix,iy,0)
          enddo
        enddo

      endif

#endif

      if (lw3dtimesubs) timeperrho3d = timeperrho3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine epush3d(np,uxp,uyp,uzp,ex,ey,ez,q,m,ddt)
      use Subtimers
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

      if (lw3dtimesubs) timeepush3d = timeepush3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine epusht3d(np,uxp,uyp,uzp,ex,ey,ez,q,m,dtp,fdt)
      use Subtimers
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

      if (lw3dtimesubs) timeepusht3d = timeepusht3d + wtime() - substarttime
      return
      end
c============================================================================= 
      subroutine bpush3d(np,uxp,uyp,uzp,gaminv,bx,by,bz,q,m,ddt,ibpush)
      use Subtimers
      integer(ISZ):: np,ibpush
      real(kind=8):: uxp(np),uyp(np),uzp(np),gaminv(np)
      real(kind=8):: bx(np),by(np),bz(np)
      real(kind=8):: q,m,ddt

c  Push the particle velocity with B field

      integer(ISZ):: ip
      real(kind=8):: tx,ty,tz,tsqi,sx,sy,sz,uxppr,uyppr,uzppr
      real(kind=8):: const
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
            tx = tan (gaminv(ip)*bx(ip)*const)
            ty = tan (gaminv(ip)*by(ip)*const)
            tz = tan (gaminv(ip)*bz(ip)*const)
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

      if (lw3dtimesubs) timebpush3d = timebpush3d + wtime() - substarttime
      return
      end
c============================================================================= 
      subroutine bpusht3d(np,uxp,uyp,uzp,gaminv,bx,by,bz,q,m,dtp,fdt,ibpush)
      use Subtimers
      integer(ISZ):: np,ibpush
      real(kind=8):: uxp(np),uyp(np),uzp(np),gaminv(np)
      real(kind=8):: bx(np),by(np),bz(np)
      real(kind=8):: dtp(np)
      real(kind=8):: q,m,fdt

c Push the particle velocity with B field using a different time step
c for each particle.

      integer(ISZ):: ip
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
            tx = tan(gaminv(ip)*bx(ip)*const*dtp(ip))
            ty = tan(gaminv(ip)*by(ip)*const*dtp(ip))
            tz = tan(gaminv(ip)*bz(ip)*const*dtp(ip))
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

      if (lw3dtimesubs) timebpusht3d = timebpusht3d + wtime() - substarttime
      return
      end
c=========================================================================== 
      subroutine xpush3d(np,xp,yp,zp,uxp,uyp,uzp,gaminv,ddt)
      use Subtimers
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

      if (lw3dtimesubs) timexpush3d = timexpush3d + wtime() - substarttime
      return
      end
c=========================================================================== 
      subroutine xpusht3d(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dtp)
      use Subtimers
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

      if (lw3dtimesubs) timexpusht3d = timexpusht3d + wtime() - substarttime
      return
      end
c============================================================================ 
      subroutine sete3d_relativity(np,ex,ey,vbeam)
      use Subtimers
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

      if (lw3dtimesubs) timesete3d_relativity = timesete3d_relativity + wtime() - substarttime
      return
      end
c============================================================================ 
      subroutine edamp(eomdamp,it,itdamp,center,np,
     &                 ex,ey,ez,exold,eyold,ezold,exlag,eylag,ezlag)
      use Subtimers
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
        if (lw3dtimesubs) timeedamp = timeedamp + wtime() - substarttime
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
         if (it .ge. itdamp) then
            if (eomdamp .gt. 0.) then
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
         if (eomdamp .gt. 0.) then
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
         if (it .lt. itdamp+1) then
            do ip = 1, np
               ex(ip) = exold(ip)
               ey(ip) = eyold(ip)
               ez(ip) = ezold(ip)
            enddo
         endif
      endif

      if (lw3dtimesubs) timeedamp = timeedamp + wtime() - substarttime
      return
      end
c============================================================================ 
      subroutine bendez3d(np,xp,zp,ez,bendres,bendradi,bends,bnezflag,linbend)
      use Subtimers
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
        if (lw3dtimesubs) timebendez3d = timebendez3d + wtime() - substarttime
        return
      endif

      do ip=1,np
         ez(ip) = ez(ip)*(1. - bendres(ip)*xp(ip)/(bendradi(ip) + xp(ip)))
      enddo

      if (lw3dtimesubs) timebendez3d = timebendez3d + wtime() - substarttime
      return
      end      
c=========================================================================== 
      subroutine zbendcor(np,xp,zp,uxp,uzp,gaminv,ddt,bendres,bendradi,bends,
     &                    linbend)
      use Subtimers
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
        if (lw3dtimesubs) timezbendcor = timezbendcor + wtime() - substarttime
        return
      endif

      do ip=1,np
        xprv = xp(ip) - ddt * uxp(ip) * gaminv(ip)
        xc = 0.5 * (xp(ip) + xprv)
        zp(ip) = zp(ip) + ddt*uzp(ip)*gaminv(ip)*bendres(ip)
     &                    *(bendradi(ip)/(bendradi(ip) + xc) - 1.)
      enddo

      if (lw3dtimesubs) timezbendcor = timezbendcor + wtime() - substarttime
      return
      end      
c=========================================================================== 
      subroutine sete3d(phi1d,selfe,np,xp,yp,zp,zgrid,xmmin,ymmin,zmmin,
     &                  dx,dy,dz,nx,ny,nz,efetch,ex,ey,ez,l2symtry,l4symtry)
      use Subtimers
      integer(ISZ):: np,nx,ny,nz
      real(kind=8):: zgrid,xmmin,ymmin,zmmin,dx,dy,dz
      real(kind=8):: phi1d(0:*),selfe(3,0:nx,0:ny,0:nz),xp(np),yp(np),zp(np)
      real(kind=8):: ex(np),ey(np),ez(np)
      integer(ISZ):: efetch
      logical(ISZ):: l2symtry,l4symtry

c  Sets self electric field for a group of up to NPARPGRP particles
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
      integer(ISZ):: indx(33,NPARPGRP)
      real(kind=8):: p(33,NPARPGRP), u1p(NPARPGRP), v1p(NPARPGRP), w1p(NPARPGRP)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      nnx  = nx + 1
      nnxy = (nx + 1) * (ny + 1)

c  Evaluation of offsets occurs on the first call only provided
c  local arrays are SAVE'd. nnxy is added to all offsets to account for
c  the fact that the phi1d passed begins at phi(0,0,-1), so the location
c  of phi(0,0,0) is equivalent to phi1d(nnxy).

      if (noff(1) .ne. -nnxy) then
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
      endif

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

           i =  abs(xp(ip) - xmmin)*dxi
           j =  abs(yp(ip) - ymmin)*dyi
           k =  (zp(ip) - zgrid - zmmin)*dzi
           ind0 = i + j*nnx + k*nnxy

           u1 = abs(xp(ip) - xmmin)*dxi - i
           v1 = abs(yp(ip) - ymmin)*dyi - j
           w1 = (zp(ip) - zgrid - zmmin)*dzi - k

           u0 = 1. - u1
           v0 = 1. - v1
           w0 = 1. - w1

c          --- Set offsets for points on symmetry axis of grid.  The offset
c          --- for points off the axis is zero.
           ox = 0
           oy = 0
           if (i == 0) ox = sox
           if (j == 0) oy = soy

c          --- Adjust sign of E field for appropriate quadrant.
           xsign = tdxi
           ysign = tdyi
           if (xp(ip) .lt. 0.) xsign = sx*tdxi
           if (yp(ip) .lt. 0.) ysign = sy*tdyi

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

            i = abs(xp(ip) - xmmin)*dxi
            j = abs(yp(ip) - ymmin)*dyi
            k = (zp(ip) - zgrid - zmmin)*dzi

            u1 = abs(xp(ip) - xmmin)*dxi - i
            v1 = abs(yp(ip) - ymmin)*dyi - j
            w1 = (zp(ip) - zgrid - zmmin)*dzi - k

            u0 = 1. - u1
            v0 = 1. - v1
            w0 = 1. - w1

c           --- Adjust sign of E field for appropriate quadrant.
            xsign = +1.
            ysign = +1.
            if (xp(ip) .lt. 0.) xsign = sx
            if (yp(ip) .lt. 0.) ysign = sy

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

             ez(ip) =       u0*v0*w0*selfe(3,i  ,j  ,k  )
     &                    + u1*v0*w0*selfe(3,i+1,j  ,k  )
     &                    + u0*v1*w0*selfe(3,i  ,j+1,k  )
     &                    + u1*v1*w0*selfe(3,i+1,j+1,k  )
     &                    + u0*v0*w1*selfe(3,i  ,j  ,k+1)
     &                    + u1*v0*w1*selfe(3,i+1,j  ,k+1)
     &                    + u0*v1*w1*selfe(3,i  ,j+1,k+1)
     &                    + u1*v1*w1*selfe(3,i+1,j+1,k+1)

          enddo

        endif

      endif

      if (lw3dtimesubs) timesete3d = timesete3d + wtime() - substarttime
      return
      end
c=========================================================================== 
      subroutine getselfe3d(phi,nx,ny,nz,selfe,nx_selfe,ny_selfe,nz_selfe,
     &                      dx,dy,dz)
      use Subtimers
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: nx_selfe,ny_selfe,nz_selfe
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),selfe(3,0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz

c Calculate the self-E via finite differences of phi.

      integer(ISZ):: ix,iy,iz
      real(kind=8):: tdxi,tdyi,tdzi
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Do some checking first.
      if (nx .ne. nx_selfe .or. ny .ne. ny_selfe .or. nz .ne. nz_selfe) then
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
      do iz=0,nz
        do iy=1,ny-1
          selfe(1, 0,iy,iz) = 0.
          selfe(2, 0,iy,iz) = tdyi*(phi( 0,iy-1,iz  ) - phi( 0,iy+1,iz  ))
          selfe(3, 0,iy,iz) = tdzi*(phi( 0,iy  ,iz-1) - phi( 0,iy  ,iz+1))
          selfe(1,nx,iy,iz) = 0.
          selfe(2,nx,iy,iz) = tdyi*(phi(nx,iy-1,iz  ) - phi(nx,iy+1,iz  ))
          selfe(3,nx,iy,iz) = tdzi*(phi(nx,iy  ,iz-1) - phi(nx,iy  ,iz+1))
        enddo
      enddo
      do iz=0,nz
        do ix=1,nx-1
          selfe(1,ix, 0,iz) = tdxi*(phi(ix-1, 0,iz  ) - phi(ix+1, 0,iz  ))
          selfe(2,ix, 0,iz) = 0.
          selfe(3,ix, 0,iz) = tdzi*(phi(ix  , 0,iz-1) - phi(ix  , 0,iz+1))
          selfe(1,ix,ny,iz) = tdxi*(phi(ix-1,ny,iz  ) - phi(ix+1,ny,iz  ))
          selfe(2,ix,ny,iz) = 0.
          selfe(3,ix,ny,iz) = tdzi*(phi(ix  ,ny,iz-1) - phi(ix  ,ny,iz+1))
        enddo
      enddo
      do iz=0,nz
          selfe(1, 0, 0,iz) = 0.
          selfe(2, 0, 0,iz) = 0.
          selfe(3, 0, 0,iz) = tdzi*(phi( 0, 0,iz-1) - phi( 0, 0,iz+1))
          selfe(1,nx, 0,iz) = 0.
          selfe(2,nx, 0,iz) = 0.
          selfe(3,nx, 0,iz) = tdzi*(phi(nx, 0,iz-1) - phi(nx, 0,iz+1))
          selfe(1, 0,ny,iz) = 0.
          selfe(2, 0,ny,iz) = 0.
          selfe(3, 0,ny,iz) = tdzi*(phi( 0,ny,iz-1) - phi( 0,ny,iz+1))
          selfe(1,nx,ny,iz) = 0.
          selfe(2,nx,ny,iz) = 0.
          selfe(3,nx,ny,iz) = tdzi*(phi(nx,ny,iz-1) - phi(nx,ny,iz+1))
      enddo

      if (lw3dtimesubs) timegetselfe3d = timegetselfe3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine stptcl3d
      use Subtimers
      use Beam_acc
      use Constant
      use InPart
      use InPart3d
      use InMesh3d
      use Particles
      use Setpwork3d
      use Picglb
#ifdef MPIPARALLEL
      use Parallel
#endif

c  Loads particles.
c
c  For the parallel version, each processor loads the fraction of the
c  longitudinally uniform beam of length (zimax-zimin) which is within its
c  region.  The cigar (or other longitudinal distribution) may move particles
c  out of that region so the routine reorg_particles is called which
c  reorganized particles according to their final longitudinal position.
c  The routine sw_globalsum recalculates the particle weights (sw)
c  appropriately.
c
c  Picglb is included only to use zbeam


      integer(ISZ):: npm,ip,is,i,j,k,i1,nxstri_s,nystri_s,nzstri_s,izstripe,icheck
      integer(ISZ):: nfib2,nfib3,nfib4,ii,ioff
      real(kind=8):: rnpmi,zmid,zlen,zleni,vtx,vty,vtz,rr,r0,rpp,phi1,phi2,rp,wz,wr
      real(kind=8):: vthr,vrbar
      real(kind=8):: zoff,vb,emx,emy,vtzz
      real(kind=8):: sumzdist,zdistmax
      real(kind=8):: sumrdist
      integer(ISZ):: iz,ipp,izz
      integer(ISZ):: ir,irr
      integer(ISZ):: envxport
      character(72):: errline
      real(kind=8):: wranf,rnorm,rnrev
#ifdef MPIPARALLEL
      real(kind=8):: zmin,zmax
#endif
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  If distrbtn is not set by the user, exit
      if (distrbtn == "none") then
        if (lw3dtimesubs) timestptcl3d = timestptcl3d + wtime() - substarttime
        return
      endif

c  Set npgrp so Setpwork3d arrays are correct size
      if (xrandom == "fibonacc" .or.
     &    xrandom == "digitrev" .or.
     &    xrandom == "pseudo") then
        npgrp = NPARPGRP
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
c  values for sq and sm.  For parallel version, also call reorg_particles to
c  get a value for sw.
      if (npm .gt. 0) then

c  Allocate scratch space for load
      call gallot("Setpwork3d",0)

c  Prepare for arbitrary particle distribution in z.  This coding will only
c  load a beam approximately matching the distribution.  A gathering of
c  the particles on a z mesh will not exactly reproduce the z distribution.
c  This problem is being studied.
      if (nzdist .gt. 0) then
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
      if (nrdist .gt. 0) then
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

c----------------------------------
c  Begin main loop over species
c----------------------------------

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

c       --- Inverse of the number of particles. When multiplied by the particle
c       --- index and added to the normalized value of zmin, the normalized
c       --- location of the particle in the beam, between 0 and 1, results.
        rnpmi = 1./dvnz(real(npm,kind=8))

c       --- Reset counter for arbitrary particle distribution in z.
c       --- It is reset for each species since each species has its own
c       --- initial uniform distribution.
c       --- Resetting izz here and starting each particle with its previous
c       --- value only works since the initial uniform z positions are
c       --- monotinically increasing.  Using this feature greatly increases
c       --- the speed of the code.
        if (nzdist .gt. 0) izz = 0

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
        call chckpart(is,0,npm,.false.)
#endif

c       --- Shift ioff by randoffset so that the seed for the quiet start
c       --- random number generators can be controlled.
        ioff = ioff + randoffset

c  MAIN INITIALIZATION LOOP

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
           if (nzdist .gt. 0) then
             do ipp=k,i1
c              --- Find what grid cell the particle is in.
c              --- Note that the previous value of izz is used (see note above)
               do while (zt(ipp-k+1) .gt. intzdist(izz+1))
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
         if (vthz_s(is) == 0. .and. k == 1) then
           do i=1,i1-k+1
              uzt(i) = 0.
           enddo
         endif

         if (cigarld) then
c          --- set vtz to vthz*2 since vtz is vz_max - vz_bar for cigar
           vtz = 2.*vthz_s(is)
c          --- load uzt with linear distribution (-.5 < uzt < .5)
           if (vtz .ne. 0.) then
             if (vzrandom == "pseudo") then
               do i=1,i1-k+1
                 uzt(i) = wranf() - .5
               enddo
             elseif (vzrandom == "digitrev") then
               do i=1,i1-k+1
                 uzt(i) = rnrev(i+k-1+ioff,dig7) - .5
               enddo
             else
               write(errline,'("ERROR: vzrandom has an improper value = ",a8)')
     &               vzrandom
               call remark (errline)
               call kaboom (0)
             endif
           endif
c          --- use cigar to adjust the z's
c          --- xt,yt,uxt,uyt passed as scratch arrays
           call cigar(i1-k+1,zt(1),uzt(1),zt(1),uzt(1),perpscal(1),
     &                straight_s(is),xt(1),yt(1),uxt(1),uyt(1))

c          --- Loading gaussian longitudinal velocity distribution.
           if (distr_l == "gaussian") then
             if (vzrandom == "pseudo") then
               do i=1,i1-k+1
                 uzt(i) = 0.5*perpscal(i)*rnorm()
               enddo
             elseif (vzrandom == "digitrev") then
               call rnormdig(k+ioff,i1-k+1,dig7,dig8,0.,uzt(1))
               do i=1,i1-k+1
                 uzt(i) = 0.5*perpscal(i)*uzt(i)
               enddo
             else
               write(errline,'("ERROR: vzrandom has an improper value = ",a8)')
     &               vzrandom
               call remark (errline)
               call kaboom (0)
             endif
           endif

         else
c          --- non-cigar load
c          --- set vtz to vthz
           vtz = vthz_s(is)
c          --- load uzt with gaussian distribution (mean 0, variance 1)
           if (vtz .ne. 0.) then
             if (vzrandom == "pseudo") then
               do i=1,i1-k+1
                 uzt(i) = rnorm()
               enddo
             elseif (vzrandom == "digitrev") then
               call rnormdig(k+ioff,i1-k+1,dig7,dig8,0.,uzt(1))
             else
               write(errline,'("ERROR: vzrandom has an improper value = ",a8)')
     &               vzrandom
               call remark (errline)
               call kaboom (0)
             endif
           endif
c          --- set perpscal to one first time through
           if (k == 1) then
             do i=1,i1-k+1
               perpscal(i) = 1.0
             enddo
           endif
         endif

c        --- unormalize z to fetch envelope
         do i=1,i1-k+1
           zt(i) = zbeam + zimin_s(is) + zlen*zt(i)
         enddo

c        --- fetch envelope
         if (nenvofz > 0) then
           do i=1,i1-k+1
             iz = (zt(i) - zimin_s(is) - zbeam)*zleni*nenvofz
             wz = (zt(i) - zimin_s(is) - zbeam)*zleni*nenvofz - iz
             at(i)  = aofz(iz)*(1.-wz)  + aofz(iz+1)*wz
             bt(i)  = bofz(iz)*(1.-wz)  + bofz(iz+1)*wz
             apt(i) = apofz(iz)*(1.-wz) + apofz(iz+1)*wz
             bpt(i) = bpofz(iz)*(1.-wz) + bpofz(iz+1)*wz
             xct(i)  = xofz(iz)*(1.-wz)  + xofz(iz+1)*wz
             yct(i)  = yofz(iz)*(1.-wz)  + yofz(iz+1)*wz
             xpct(i) = xpofz(iz)*(1.-wz) + xpofz(iz+1)*wz
             ypct(i) = ypofz(iz)*(1.-wz) + ypofz(iz+1)*wz
           enddo
         elseif (.not. cylinder) then
           icheck = envxport(i1-k+1,zt,at,apt,bt,bpt,xct,xpct,yct,ypct)
           if (icheck == 1) then
              call remark ("ERROR: out-of-range z sent to ENVXPORT")
              call kaboom (0)
           endif
         else
c          --- set only first time through
           if (k == 1) then
             do i=1,i1-k+1
               at(i)   = a0_s(is)
               bt(i)   = b0_s(is)
               apt(i)  = ap0_s(is)
               bpt(i)  = bp0_s(is)
               xct(i)  = x0_s(is)  + xcent_s(is)
               yct(i)  = y0_s(is)  + ycent_s(is)
               xpct(i) = xp0_s(is) + xpcent_s(is)
               ypct(i) = yp0_s(is) + ypcent_s(is)
             enddo
           endif
         endif


c  Load tranverse stuff
c  semi-gaussian distribution
         if (distrbtn == "semigaus") then

c          --- load normalized tranverse space variables
           if (xrandom == "pseudo") then
             do i=1,i1-k+1
               xt(i) = 2.*wranf() - 1.
               yt(i) = 2.*wranf() - 1.
             enddo
           elseif (xrandom == "fibonacc") then
             nfib2 = nfibgrps*fibg2
             nfib3 = nfibgrps*fibg3
             do i=1,i1-k+1
               xt(i) = 2.*mod((nfib2*(i+k-2.0+ioff) + 0.5)*rnpmi,1.0) - 1.
               yt(i) = 2.*mod((nfib3*(i+k-2.0+ioff) + 0.5)*rnpmi,1.0) - 1.
             enddo
           elseif (xrandom == "digitrev") then
             do i=1,i1-k+1
               xt(i) = 2.*rnrev(i+k-1+ioff,dig1) - 1.
               yt(i) = 2.*rnrev(i+k-1+ioff,dig2) - 1.
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
               if (nrdist .gt. 0) then
c                --- Find what grid cell the particle is in.
                 irr = 0
                 do while (rt(i) .gt. intrdist(irr+1))
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
               if (xt(i)**2 + yt(i)**2 .lt. 1.) then
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
                 if (rp .gt. 0.5) rpp = sqrt(1. - 0.5*sqrt(3. - 3.*rp**2))
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
     &                                .lt. (at(i)*bt(i)/r0)**2) then
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
             if (vtrandom == "pseudo") then
               do i=1,i1-k+1
                 uxt(i) = rnorm()
                 uyt(i) = rnorm()
               enddo
             elseif (vtrandom == "digitrev") then
               call rnormdig(k+ioff,i1-k+1,dig3,dig4,0.,uxt(1))
               call rnormdig(k+ioff,i1-k+1,dig5,dig6,0.,uyt(1))
             else
               write (errline,'("ERROR: vtrandom has an improper value = ",a8)')
     &             vtrandom
               call remark (errline)
               call kaboom (0)
             endif
           else if (distr_t == "uniform") then
             if (vtrandom == "pseudo") then
               do i=1,i1-k+1
                 uxt(i) = 2.*wranf() - 1.
                 uyt(i) = 2.*wranf() - 1.
               enddo
             elseif (vtrandom == "digitrev") then
               do i=1,i1-k+1
                 uxt(i) = 2.*rnrev(k+ioff-1+i,dig3) - 1.
                 uyt(i) = 2.*rnrev(k+ioff-1+i,dig5) - 1.
               enddo
             else
               write (errline,'("ERROR: vtrandom has an improper value = ",a8)')
     &             vtrandom
               call remark (errline)
               call kaboom (0)
             endif
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

c          --- unnormalize everything and load into particle arrays
           do i=1,j
             ii = indx(i)
c            --- put position into particle arrays
             xp(ip+i) = at(ii)*xt(ii)*perpscal(ii) + xct(ii)
             yp(ip+i) = bt(ii)*yt(ii)*perpscal(ii) + yct(ii)
             zp(ip+i) = zt(ii)
c            --- Use optional axially varying vbeam and emittance
             if (nvbeamofz == 0) then
               vb = vbeam_s(is)
             else
               iz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nvbeamofz
               wz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nvbeamofz - iz
               vb = vbeamofz(iz)*(1.-wz) + vbeamofz(iz+1)*wz
             endif
             if (nvthzofz == 0) then
               vtzz = vtz
             else
               iz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nvthzofz
               wz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nvthzofz - iz
               vtzz = vthzofz(iz)*(1.-wz) + vthzofz(iz+1)*wz
             endif
             if (nemitofz == 0) then
               emx = emitx_s(is)
               emy = emity_s(is) 
             else
               iz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nemitofz
               wz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nemitofz - iz
               emx = emitxofz(iz)*(1.-wz) + emitxofz(iz+1)*wz
               emy = emityofz(iz)*(1.-wz) + emityofz(iz+1)*wz
             endif
c            --- put velocity into particle arrays
             vtx = .5*vb*emx/at(ii)*perpscal(ii) + vthperp_s(is)
             vty = .5*vb*emy/bt(ii)*perpscal(ii) + vthperp_s(is)
             uxp(ip+i) = vb*apt(ii)*xt(ii)*perpscal(ii) + vtx*uxt(ii) +
     &                   vb*xpct(ii)
             uyp(ip+i) = vb*bpt(ii)*yt(ii)*perpscal(ii) + vty*uyt(ii) +
     &                   vb*ypct(ii)
             uzp(ip+i) = vb*(1. + vtilt_s(is)*(zmid - zt(ii))*zleni) +
     &                   vtzz*uzt(ii)
           enddo

c          --- add sinusoidal perturbation to uzp
           if (vzperamp .ne. 0.) then
             do i=ip+1,ip+j
               uzp(i) = uzp(i) + vzperamp*sin(2.*pi*zp(i)/vzperlam + vzperphs)
             enddo
           endif

c          --- increment number of particles by size of current group
           ip = ip + j

c  K-V distribution
         elseif (distrbtn == "K-V") then

c          --- fetch random numbers and put into xt,yt, & uxt  temporarily
           if (xrandom == "pseudo") then
             do i=1,i1-k+1
               xt(i)  = wranf()
               yt(i)  = wranf()
               uxt(i) = wranf()
             enddo
           elseif (xrandom == "fibonacc") then
             nfib2 = nfibgrps*fibg2
             nfib3 = nfibgrps*fibg3
             nfib4 = nfibgrps*fibg4
             do i=1,i1-k+1
               xt(i)  = mod((nfib2*(i+k-2.0+ioff) + 0.5)*rnpmi,1.0)
               yt(i)  = mod((nfib3*(i+k-2.0+ioff) + 0.5)*rnpmi,1.0)
               uxt(i) = mod((nfib4*(i+k-2.0+ioff) + 0.5)*rnpmi,1.0)
             enddo
           elseif (xrandom == "digitrev") then
             do i=1,i1-k+1
               xt(i)  = rnrev(i+k-1+ioff,dig1)
               yt(i)  = rnrev(i+k-1+ioff,dig2)
               uxt(i) = rnrev(i+k-1+ioff,dig3)
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
c            --- put position into particle arrays
             xp(ip+i) = rr*cos(phi1)*at(i)*perpscal(i) + xct(i)
             yp(ip+i) = rr*sin(phi1)*bt(i)*perpscal(i) + yct(i)
             zp(ip+i) = zt(i)
c            --- Use optional axially varying vbeam and emittance
             if (nvbeamofz == 0) then
               vb = vbeam_s(is)
             else
               iz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nvbeamofz
               wz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nvbeamofz - iz
               vb = vbeamofz(iz)*(1.-wz) + vbeamofz(iz+1)*wz
             endif
             if (nvthzofz == 0) then
               vtzz = vtz
             else
               iz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nvthzofz
               wz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nvthzofz - iz
               vtzz = vthzofz(iz)*(1.-wz) + vthzofz(iz+1)*wz
             endif
             if (nemitofz == 0) then
               emx = emitx_s(is)
               emy = emity_s(is) 
             else
               iz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nemitofz
               wz = (zp(ip+i) - zimin_s(is) - zbeam)*zleni*nemitofz - iz
               emx = emitxofz(iz)*(1.-wz) + emitxofz(iz+1)*wz
               emy = emityofz(iz)*(1.-wz) + emityofz(iz+1)*wz
             endif
c            --- put velocity into particle arrays
             rr = sqrt(1.-rr*rr)
             vtx = vb*emx/at(i)*perpscal(i) + vthperp_s(is)
             vty = vb*emy/bt(i)*perpscal(i) + vthperp_s(is)
             uxp(ip+i) = vb*xp(ip+i)*apt(i)/at(i) + vtx*rr*cos(phi2) +
     &                   vb*xpct(i)
             uyp(ip+i) = vb*yp(ip+i)*bpt(i)/bt(i) + vty*rr*sin(phi2) +
     &                   vb*ypct(i)
             uzp(ip+i) = vb*(1. + vtilt_s(is)*(zmid - zt(i))*zleni) +
     &                   vtzz*uzt(i)
           enddo

c          --- add sinusoidal perturbation to uzp
           if (vzperamp .ne. 0.) then
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

      enddo

c  Set particle number, indices, etc.

      nps(is) = ip - ins(is) + 1
      enddo

c--------------------------------
c  End main loop over species
c--------------------------------

        call gfree ("Setpwork3d")
c     --- end of if checking for positive npm
      endif

c     --- Find maximum of zdist, which is used to scale sw.  See below.
      if (nzdist .gt. 0) then
        zdistmax = nrmzdist(0)
        do iz=1,nzdist
          if (nrmzdist(iz) .gt. zdistmax) zdistmax = nrmzdist(iz)
        enddo
      endif

c     --- calculate values for sm, sq, and sw for each species.
      do is = 1, ns
        sq(is) = zion_s(is) * echarge
        sm(is) = aion_s(is) * amu
        zlen = zimax_s(is) - zimin_s(is)
        if (vbeam_s(is) .ne. 0. .and. nps(is) .ne. 0) then
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
        if (nzdist .gt. 0) then
          sw(is) = sw(is)/(zdistmax*nzdist)
        endif

      enddo

c     --- if slave, reorganize particles and recalculate sw
#ifdef MPIPARALLEL
      call reorg_particles
      call sw_globalsum(ns,sw(1))
#endif

      if (lw3dtimesubs) timestptcl3d = timestptcl3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine setrho3d(rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,depos)
      use Subtimers
      use InGen3d
      use InMesh3d
      use Fields3d
      use Picglb3d
      integer(ISZ):: np
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho1d(0:*)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      character(8):: depos

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
      integer(ISZ):: moff(0:7), indx(0:7,NPARPGRP)
      save moff
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ):: ii(0:NPARPGRP-1), jj(0:NPARPGRP-1), kk(0:NPARPGRP-1)
c     --- For both algorithms
      real(kind=8):: s(0:7,NPARPGRP)

      integer(ISZ):: ip,i,j,k,ind0,m,ir
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
      if (depos == "vector") then

c  Set up offset array for vectorized deposition:

      if (moff(4) .ne. (nx+1)*(ny+1)) then
        moff(0) = 0
        moff(1) = 1
        moff(2) = nx+1
        moff(3) = nx+2
        moff(4) = (nx+1)*(ny+1)
        moff(5) = (nx+1)*(ny+1)+1
        moff(6) = (nx+1)*(ny+1)+nx+1
        moff(7) = (nx+1)*(ny+1)+nx+2
      endif

c  Begin main loop over species, groups of NPARPGRP particles

      g = wght * q / (dx * dy * dz)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
c       --- special loop for 2-fold symmetry
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        g = g*0.5
        do ip = 1, np
           i  = int((xp(ip) - xmmin) * dxi)
           u1 = (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = int(abs(yp(ip))*dyi)
           v1 = abs(yp(ip))*dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 = (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ip) = ind0 + moff(0)
           indx(1,ip) = ind0 + moff(1)
           indx(2,ip) = ind0 + moff(2)
           indx(3,ip) = ind0 + moff(3)
           indx(4,ip) = ind0 + moff(4)
           indx(5,ip) = ind0 + moff(5)
           indx(6,ip) = ind0 + moff(6)
           indx(7,ip) = ind0 + moff(7)
           gyfact = 1.
           if (j == 0) gyfact = 2.
           s(0,ip) = u0 * v0 * w0 * g*gyfact
           s(1,ip) = u1 * v0 * w0 * g*gyfact
           s(2,ip) = u0 * v1 * w0 * g
           s(3,ip) = u1 * v1 * w0 * g
           s(4,ip) = u0 * v0 * w1 * g*gyfact
           s(5,ip) = u1 * v0 * w1 * g*gyfact
           s(6,ip) = u0 * v1 * w1 * g
           s(7,ip) = u1 * v1 * w1 * g
        enddo
      elseif (l4symtry) then
c       --- special loop for 4-fold symmetry
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        g = g*0.25
        do ip = 1, np
           i  = int(abs(xp(ip))*dxi)
           u1 = abs(xp(ip))*dxi - i
           u0 = 1. - u1
           j  = int(abs(yp(ip))*dyi)
           v1 = abs(yp(ip))*dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 = (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ip) = ind0 + moff(0)
           indx(1,ip) = ind0 + moff(1)
           indx(2,ip) = ind0 + moff(2)
           indx(3,ip) = ind0 + moff(3)
           indx(4,ip) = ind0 + moff(4)
           indx(5,ip) = ind0 + moff(5)
           indx(6,ip) = ind0 + moff(6)
           indx(7,ip) = ind0 + moff(7)
           gxfact = 1.
           gyfact = 1.
           if (i == 0) gxfact = 2.
           if (j == 0) gyfact = 2.
           s(0,ip) = u0 * v0 * w0 * g*gxfact*gyfact
           s(1,ip) = u1 * v0 * w0 * g*gyfact
           s(2,ip) = u0 * v1 * w0 * g*gxfact
           s(3,ip) = u1 * v1 * w0 * g
           s(4,ip) = u0 * v0 * w1 * g*gxfact*gyfact
           s(5,ip) = u1 * v0 * w1 * g*gyfact
           s(6,ip) = u0 * v1 * w1 * g*gxfact
           s(7,ip) = u1 * v1 * w1 * g
        enddo
      else
c       --- normal loop
        do ip = 1, np
           i  = int((xp(ip) - xmmin) * dxi)
           u1 = (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = int((yp(ip) - ymmin) * dyi)
           v1 = (yp(ip) - ymmin) * dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin) * dzi)
           w1 = (zp(ip) - zgrid - zmmin) * dzi - k
           w0 = 1. - w1
           ind0 = i + j*(nx+1) + k*(nx+1)*(ny+1)
           indx(0,ip) = ind0 + moff(0)
           indx(1,ip) = ind0 + moff(1)
           indx(2,ip) = ind0 + moff(2)
           indx(3,ip) = ind0 + moff(3)
           indx(4,ip) = ind0 + moff(4)
           indx(5,ip) = ind0 + moff(5)
           indx(6,ip) = ind0 + moff(6)
           indx(7,ip) = ind0 + moff(7)
           s(0,ip) = u0 * v0 * w0 * g
           s(1,ip) = u1 * v0 * w0 * g
           s(2,ip) = u0 * v1 * w0 * g
           s(3,ip) = u1 * v1 * w0 * g
           s(4,ip) = u0 * v0 * w1 * g
           s(5,ip) = u1 * v0 * w1 * g
           s(6,ip) = u0 * v1 * w1 * g
           s(7,ip) = u1 * v1 * w1 * g
        enddo
      endif
c     --- vectorized deposition over the 8 cells touched;
c     --- there'd be a hazard if we interchanged the loops.
      do ip = 1, np
         if (uzp(ip).gt.0) then
         do m = 0, 7
            rho1d(indx(m,ip)) = rho1d(indx(m,ip)) + s(m,ip)
         enddo
         endif
      enddo

c--------------------------------------
c  Begin scalar deposition loop
c--------------------------------------
      elseif (depos == "scalar") then

c  Begin main loop over species, groups of NPARPGRP particles

      g = wght * q / (dx * dy * dz)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        g = g*0.5
        do ip = 1, np
          ir = ip - 1
          ii(ir) = int((xp(ip) - xmmin) * dxi)
          u1     = (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int(abs(yp(ip))*dyi)
          v1     = abs(yp(ip))*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     = (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gyfact = 1.
          if (jj(ir) == 0) gyfact = 2.
          s(0,ip) = u0 * v0 * w0 * g*gyfact
          s(1,ip) = u1 * v0 * w0 * g*gyfact
          s(2,ip) = u0 * v1 * w0 * g
          s(3,ip) = u1 * v1 * w0 * g
          s(4,ip) = u0 * v0 * w1 * g*gyfact
          s(5,ip) = u1 * v0 * w1 * g*gyfact
          s(6,ip) = u0 * v1 * w1 * g
          s(7,ip) = u1 * v1 * w1 * g
        enddo
      elseif (l4symtry) then
        g = g*0.25
        do ip = 1, np
          ir = ip - 1
          ii(ir) = int(abs(xp(ip))*dxi)
          u1     = abs(xp(ip))*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int(abs(yp(ip))*dyi)
          v1     = abs(yp(ip))*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     = (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii(ir) == 0) gxfact = 2.
          if (jj(ir) == 0) gyfact = 2.
          s(0,ip) = u0 * v0 * w0 * g*gxfact*gyfact
          s(1,ip) = u1 * v0 * w0 * g*gyfact
          s(2,ip) = u0 * v1 * w0 * g*gxfact
          s(3,ip) = u1 * v1 * w0 * g
          s(4,ip) = u0 * v0 * w1 * g*gxfact*gyfact
          s(5,ip) = u1 * v0 * w1 * g*gyfact
          s(6,ip) = u0 * v1 * w1 * g*gxfact
          s(7,ip) = u1 * v1 * w1 * g
        enddo
      else
c       --- normal loop
        do ip = 1, np
          ir = ip - 1
          ii(ir) = int((xp(ip) - xmmin) * dxi)
          u1     = (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((yp(ip) - ymmin) * dyi)
          v1     = (yp(ip) - ymmin) * dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin) * dzi)
          w1     = (zp(ip) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          s(0,ip) = u0 * v0 * w0 * g
          s(1,ip) = u1 * v0 * w0 * g
          s(2,ip) = u0 * v1 * w0 * g
          s(3,ip) = u1 * v1 * w0 * g
          s(4,ip) = u0 * v0 * w1 * g
          s(5,ip) = u1 * v0 * w1 * g
          s(6,ip) = u0 * v1 * w1 * g
          s(7,ip) = u1 * v1 * w1 * g
        enddo
      endif
c     --- scalar loop does the actual deposition
      do ir = 0, np-1
         if (uzp(ir+1).gt.0) then
         rho(ii(ir)  ,jj(ir)  ,kk(ir)  )
     &  =rho(ii(ir)  ,jj(ir)  ,kk(ir)  ) + s(0,ir+1)
         rho(ii(ir)+1,jj(ir)  ,kk(ir)  )
     &  =rho(ii(ir)+1,jj(ir)  ,kk(ir)  ) + s(1,ir+1)
         rho(ii(ir)  ,jj(ir)+1,kk(ir)  )
     &  =rho(ii(ir)  ,jj(ir)+1,kk(ir)  ) + s(2,ir+1)
         rho(ii(ir)+1,jj(ir)+1,kk(ir)  )
     &  =rho(ii(ir)+1,jj(ir)+1,kk(ir)  ) + s(3,ir+1)
         rho(ii(ir)  ,jj(ir)  ,kk(ir)+1)
     &  =rho(ii(ir)  ,jj(ir)  ,kk(ir)+1) + s(4,ir+1)
         rho(ii(ir)+1,jj(ir)  ,kk(ir)+1)
     &  =rho(ii(ir)+1,jj(ir)  ,kk(ir)+1) + s(5,ir+1)
         rho(ii(ir)  ,jj(ir)+1,kk(ir)+1)
     &  =rho(ii(ir)  ,jj(ir)+1,kk(ir)+1) + s(6,ir+1)
         rho(ii(ir)+1,jj(ir)+1,kk(ir)+1)
     &  =rho(ii(ir)+1,jj(ir)+1,kk(ir)+1) + s(7,ir+1)
        endif
      enddo

      endif
      if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine loadrho3d(ins_i,nps_i,is_i,lzero)
      use Subtimers
      use InGen
      use InGen3d
      use InPart
      use InMesh3d
      use Picglb
      use Picglb3d
      use Particles
      use Fields3d
      use FRZmgrid
      integer(ISZ):: ins_i,nps_i,is_i
      logical(ISZ):: lzero

c --- This routine provides a simple call from the interpreter to load the
c --- rho array.  The value '-1' is used as a flag in the input to use
c --- all of the particles, otherwise the specified particles are loaded.


      integer(ISZ):: ins_u,nps_u
      integer(ISZ):: is1,is2
      integer(ISZ):: ip,ipmin,is
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- zero rho if requested
      if (lzero) call zeroarry(rho(0,0,0),(nx+1)*(ny+1)*(nz+1))
c     if (lzero) rho = 0.

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
            if(solvergeom==XYZgeom) then
              call setrho3d (rho(0,0,0), ip, xp(ipmin), yp(ipmin), zp(ipmin),
     &                       zgrid, uzp(ipmin), sq(is), sw(is), depos)
            elseif(solvergeom==RZgeom) then
              call rhoweightrz(xp(ipmin),yp(ipmin),zp(ipmin),ip,
     &                         sq(is)*sw(is),rho(0,0,0),nx,nz,dx,dz,zmmin)
            end if
         enddo
      enddo

c     --- enforce axial periodicity
      if (lzero) call perrho3d(rho(0,0,0),nx,ny,nz,periinz)

      if (lw3dtimesubs) timeloadrho3d = timeloadrho3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine stckxy3d(np,xp,xmmax,xmmin,dx,yp,ymmax,ymmin,dy,
     &                   zp,zmmin,dz,uxp,uyp,uzp,zgrid,zbeam,
     &                   nzzarr,prwallz,prwallxz,prwallyz,prwelips,
     &                   l2symtry,l4symtry,lostpars,zzmin,dzzi)
      use Subtimers
      use InGen3d, only: solvergeom, RZgeom
      integer(ISZ):: np,nzzarr
      real(kind=8):: xmmax,xmmin,dx,ymmax,ymmin,dy,zmmin,dz,zgrid,zbeam
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np)
      real(kind=8):: prwallz(0:nzzarr),prwallxz(0:nzzarr),prwallyz(0:nzzarr)
      real(kind=8):: prwelips(0:nzzarr)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: lostpars(0:nzzarr)
      real(kind=8):: zzmin,dzzi

c Enforces transverse sticky boundary conditions.
c Particles are considered stuck when outside a z-dependent aperture of
c radius prwallz, centered at (prwallxz, prwallyz). Particles are also
c lost if outside of the computational domain.
c Lost particles are moved to the edge transversely and to the left end of
c the grid longitudinally, as well as having the velocity set to zero. Zero
c longitudinal velocity is the flag for lost particles. Note that particles
c previously lost are continually moved to the left end of the grid (which
c is important with a moving grid).
c The array lostpars keeps track of where particles are lost (in the beam
c frame). A check is done so that already lost particles don't contribute to
c lostpars.

      integer(ISZ):: ip,iz
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

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
        ymin = -ymmax
      endif

      do ip=1,np
        iz = max(0,min(nzzarr,int((zp(ip) - zzmin - zbeam)*dzzi)))
        if (((xp(ip)-prwallxz(iz))**2 +
     &       (yp(ip)-prwallyz(iz)*prwelips(iz))**2) >= prwallz(iz)**2 .or.
     &      (xmin > xp(ip) .or. xp(ip) > xmax .or.
     &       ymin > yp(ip) .or. yp(ip) > ymax)) then
          if (uzp(ip) .ne. 0.) lostpars(iz) = lostpars(iz) + 1
          xp(ip) = xmmax - .1*dx
          yp(ip) = ymmax - .1*dy
          zp(ip) = zmmin - .9999*dz + zgrid
          uxp(ip) = 0.0
          uyp(ip) = 0.0
          uzp(ip) = 0.0
        endif
      enddo

      if (lw3dtimesubs) timestckxy3d = timestckxy3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine setrstar(rstar,nz,dz,zmmin,zgrid)
      use Subtimers
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

c     --- Get the bend radius of curvature and residence fraction at the
c     --- list of z points.
      call getbend(nz+3,nz+3,zz,vz,gi,bendres,rstar,-0.5*dz,0.5*dz,.false.)

c     --- Scale the radius by the reciprical of the residence fraction.
c     --- Also, make sure the rstar is non-zero.
      do iz = -1,nz+1
        if (bendres(iz) > 0.) rstar(iz) = rstar(iz)/bendres(iz)
        if (rstar(iz) == 0.) rstar(iz) = LARGEPOS
      enddo

      if (lw3dtimesubs) timesetrstar = timesetrstar + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine inject3d(itask)
      use Subtimers
      use Constant
      use Beam_acc
      use InGen
      use InGen3d
      use InPart
      use InPart3d
      use InMesh3d
      use Particles
      use Picglb
      use InjectVars
      use InjectVars3d
      use Picglb3d
      use Setpwork3d
      use Fields3d
#ifdef MPIPARALLEL
      use Parallel
#endif
      integer(ISZ):: itask

c This routine does the injection of particles.  It has three main parts, the
c first cleans up the particle arrays and makes sure that there is enough room
c for the new particles.  The second is the coding for constant current
c injection.  The third part is for space-charge limited injection.


      integer(ISZ):: is,ij,ii,ip,nn,ix,iy,iz,ipmin,ixx,iyy,ijz
      integer(ISZ):: inj_is,inj_ij
      real(kind=8):: rnpinjct,rnn,t
      real(kind=8):: rr,ehalf,vnorm,zz,wxx,wyy,aa,az,wz,xx,yy,xxi,xm,ym,zm
      real(kind=8):: rnorm,rnrev,qoverm,const,dxi,dyi,dzi,ez,vznorm,rp
      real(kind=8):: clightsqi
      real(kind=8):: xinj,yinj,ainj,ainji,ainjmini,binj,binji,binjmini
      real(kind=8):: apinj,bpinj,xpinj,ypinj,rinj,rinji
      real(kind=8):: sphericalcorrection
      real(kind=8):: vtx,vty
      real(kind=8):: zmid,zleni,ztilt,gamma
      real(kind=8):: zinjmin,zinjmax
      integer(ISZ):: load_fac,il,ilf
      real(kind=8):: signx(4),signy(4)
      real(kind=8):: p2x,p2y,p2z,w2x,w2y,w2z,area
      integer(ISZ):: i2x,i2y,i2z
      real(kind=8):: wranf
      integer(ISZ):: ith,nti
      real(kind=8):: cl_const(ns),jton(ns),vzconst(ns)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      if (lrelativ) clightsqi = 1./clight**2

c Print warning if vbeamfrm is not zero
      if (abs(vbeamfrm) .gt. 0.) then
       call remark("WARNING: vbeamfrm is not zero and injection is being used.")
       call remark("A core dump will result if the emitting surface moves off of the field grid.")
      endif

c Reset particle index so that dead particles at end of particle block are
c taken out of action.  Also zero out the array npinje_s.
      if (itask == 1) then
        do is=1,ns
          if (nps(is) == 0) cycle
          do while (uzp(ins(is) + nps(is)-1) == 0. .and. nps(is) > 0)
            nps(is) = nps(is) - 1
          enddo
        enddo
      endif

c Set inject z min and max. This is only really needed for
c the parallel version, but is useful for the serial version to
c prevent particles from being injected off of the grid.
#ifdef MPIPARALLEL
      zinjmin = zpslmin(my_index)
      zinjmax = zpslmax(my_index)
#else
      zinjmin = zmmin
      zinjmax = zmmax
#endif

c If itask = 1, zero npinject. This is done here since npinject is used both
c by normal and transverse space-charge limited injection. Both sections
c add onto npinject. Only do this if actually injecting particles.
      if (itask == 1 .and. (inject == 2 .or. inject == 3)) npinject = 0

c The check of whether or not there is enough space for the injected particles
c is now done in each injection method.

      if (inject == 1) then
c Injects particles with constant (or user varied) current.  It does a split
c leap frog advance over partial time steps, starting with all of the particles
c in the plane at zinject.  The particles are advanced uniformly distribuated
c fractions of a timestep to fill the injection region.  The E and B fields are
c first gathered (assuming that injection is not being done inside a bend) then
c a half advance is done on the velocity over the partial timestep.  The
c positions are then advanced for the full amount of the partial timestep and
c the particles are loaded onto the rho grid.  The routine then exits so the
c fields can be solved for.  The routine is called again to complete the
c velocity advance using the new fields including the injected particles.  The
c velocity ends up at time level it-1/2 to match the time level of the rest of
c the particles.
c
c One caveat is the gathering of the E fields (in the call to sete3d).  If the
c injection plane is at the left hand of the grid, the E gather will use the
c plane iz = -1 which is not necessarily correct.  For example, in injection
c off of a conducting plane, phi at iz=-1 is the same as phi at iz=0, this
c gives a E field that is too small (by about a factor of 2).

        if (itask == 1) then
c Start injected particles advance.  Particles end up with position
c at current time level and velocity at a time level of one half of the time
c each particle was advanced from the plane of injection (i.e. all different).

c         --- Zero out temp array.
          npinjtmp = 0

c         --- loop over injection sources
          do ij=1,ninject
            ainj = ainject(ij)
            binj = binject(ij)
            ainji = 1./ainject(ij)
            binji = 1./binject(ij)
            xinj = xinject(ij)
            yinj = yinject(ij)
            apinj = apinject(ij)
            bpinj = bpinject(ij)
            xpinj = xpinject(ij)
            ypinj = ypinject(ij)

c         --- Loop over species.
          do is=1,ns
            nn = npinje_s(is)*finject(ij,is)

c           --- Check if there is enough room for the new injected particles.
            call chckpart(is,nn,0,.true.)

c           --- Make sure there is enough temp space
            if (nn > npgrp) then
              npgrp = nn
              call gallot("Setpwork3d",0)
            endif

c           --- set random numbers
c           --- xt(ip) used to hold f(r**2)=constant
            if (xrandom == "pseudo") then
              do ip = 1,nn
                xt(ip) = wranf()*(1.-(ainjmin(ij)*ainji)**2)+
     &                               (ainjmin(ij)*ainji)**2
              enddo
            else
c             --- Default is to use digitrev
              do ip = 1,nn
                xt(ip) = rnrev(injctcnt+ip-1,dig1)*(1.-(ainjmin(ij)*ainji)**2)+
     &                                                 (ainjmin(ij)*ainji)**2
              enddo
            endif

c           --- or xt(ip) used to hold f(r)~(h+(1-h)r^2)
            if (hollow == 2) then
              do ip = 1,nn
                xt(ip) = (1 + hollow_h)*xt(ip)/
     &               (hollow_h + sqrt(hollow_h**2 + (1. - hollow_h**2)*xt(ip)))
              enddo
            endif

c           --- Convert r**2 in to r
            do ip = 1,nn
              xt(ip) = sqrt(xt(ip))
            enddo

c           --- yt(ip) used to hold f(theta)=constant
            if (xrandom == "pseudo") then
              do ip = 1,nn
                yt(ip) = 2.*pi*wranf()
              enddo
            else
c             --- Default is to use digitrev
              do ip = 1,nn
                yt(ip) = 2.*pi*rnrev(injctcnt+ip-1,dig2)
              enddo
            endif

c           --- Get velocity distribution.
            if (distrbtn == "K-V") then
c             --- For the K-V load, need another random number.
              if (vtrandom == "pseudo") then
                do ip = 1,nn
                  zt(ip) = wranf()
                enddo
              else
c               --- Default is to use digitrev
                do ip = 1,nn
                  zt(ip) = rnrev(injctcnt+ip-1,dig3)
                enddo
              endif
              do ip = 1,nn
                rr = sqrt(1. - xt(ip)*xt(ip))
                uxt(ip) = rr*cos(2.*pi*zt(ip))
                uyt(ip) = rr*sin(2.*pi*zt(ip))
              enddo
            else
c             --- Default is to use a Semi-Gaussian distribution
              if (vtrandom == "pseudo") then
                do ip = 1,nn
                  uxt(ip) = rnorm()
                  uyt(ip) = rnorm()
                enddo
              else
c               --- Default is to use digitrev
                call rnormdig(injctcnt,nn,dig3,dig4,0.,uxt(1))
                call rnormdig(injctcnt,nn,dig5,dig6,0.,uyt(1))
              endif
            endif

c           --- Axial thermal spread
            if (vzrandom == "pseudo") then
              do ip = 1,nn
                uzt(ip) = rnorm()
              enddo
            else
c             --- Default is to use digitrev
              call rnormdig(injctcnt,nn,dig7,dig8,0.,uzt(1))
            endif

c           --- Increment the particle counter for the digit reversed seed
            injctcnt = injctcnt + nn

c           --- Load particles into initial positions.
c           --- Shifts axial position by values in inj_grid.
c           --- uses beam envelope of ainject, binject, apinject and bpinject
c           --- (data is loaded into particle arrays)
            zmid = 0.5*(zimax_s(is) - zimin_s(is)) - vbeam_s(is)*dt*(it+1)
            zleni = 1./dvnz(zimax_s(is) - zimin_s(is))
            do ip=1,nn
              ii = ins(is) - 1 - npinjtmp(is,ij)
              xp(ii) = ainj*xt(ip)*cos(yt(ip)) + xinj
              yp(ii) = binj*xt(ip)*sin(yt(ip)) + yinj
c             --- find injection surface axial location via interpolation
c             --- The absolute values are needed for the cases of symmetry.
              ixx = abs(xp(ii) - xmmin)*dxi
              iyy = abs(yp(ii) - ymmin)*dyi
              wxx = abs(xp(ii) - xmmin)*dxi - ixx
              wyy = abs(yp(ii) - ymmin)*dyi - iyy
              zz = zmmin + dz*(inj_grid(ixx  ,iyy  )*(1. - wxx)*(1. - wyy) +
     &                         inj_grid(ixx+1,iyy  )*      wxx *(1. - wyy) +
     &                         inj_grid(ixx  ,iyy+1)*(1. - wxx)*      wyy  +
     &                         inj_grid(ixx+1,iyy+1)*      wxx *      wyy)
c             --- Only save particle if it is within the injection region.
              if (zinjmin <= zz .and. zz < zinjmax) then
                zp(ii) = zz
                vtx = vbeam_s(is)*0.5*emitx_s(is)*ainji + vthperp_s(is)
                vty = vbeam_s(is)*0.5*emity_s(is)*binji + vthperp_s(is)
                uxp(ii) = vbeam_s(is)*(apinj*(xp(ii)-xinj)*ainji + xpinj) +
     &                    vtx*uxt(ip)
                uyp(ii) = vbeam_s(is)*(bpinj*(yp(ii)-yinj)*binji + ypinj) +
     &                    vty*uyt(ip)
                uzp(ii) = vbeam_s(is) + 0.5*vthz_s(is)*uzt(ip)
                if (lrelativ) then
                  gaminv(ii) = sqrt(1.-
     &                           (uxp(ii)**2+uyp(ii)**2+uzp(ii)**2)*clightsqi)
                  gamma = 1./gaminv(ii)
                  uxp(ii) = uxp(ii)*gamma
                  uyp(ii) = uyp(ii)*gamma
                  uzp(ii) = uzp(ii)*gamma
                else
                  gaminv(ii) = 1.
                endif
c               --- Save the ID of the injection source.
                pid(ii) = ij
c               --- increment particle counter
                npinjtmp(is,ij) = npinjtmp(is,ij) + 1
              endif
            enddo

c           --- Get actual number of particles injected and ipmin
            nn = npinjtmp(is,ij)

c           --- make injected particles live particles
            ins(is) = ins(is) - nn
            nps(is) = nps(is) + nn

c           --- get E self-field at initial positions
            call inj_sete3d(nn,ins(is),xmmin,ymmin,zmmin,dx,dy,dz,
     &                      at(1),bt(1),apt(1))

c           --- get external fields at initial positions
c           --- (Assumes that injection will never be done inside a bend.)
            do ip=1,nn
              xt(ip) = 0.
              yt(ip) = 0.
              zt(ip) = 0.
              uxt(ip) = 0.
              uyt(ip) = 1.
            enddo
            do ipmin = ins(is),ins(is)+nn-1,NPARPGRP
              ip = min(NPARPGRP, ins(is)+nn-ipmin)
              ii = ipmin - ins(is) + 1
              call exteb3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                     gaminv(ipmin),0.,dt*0.5,
     &                     xt(ii),yt(ii),zt(ii),at(ii),bt(ii),apt(ii),
     &                     sm(is),sq(is),uxt(ii),uyt(ii),1.,dt)
            enddo

c           --- Set fractional time to advance particles
            rnpinjct = 1./dvnz(real(nn,kind=8))
            do ip=1,nn
              bpt(ip) = (ip - .5)*rnpinjct*dt
            enddo

c           --- Add in a velocity tilt (this is dependent on the fractional
c           --- time step and so must be done here).
            do ip=1,nn
              ii = ins(is) - 1 + ip
              ztilt = bpt(ip)*vbeam_s(is)
              uzp(ii) = uzp(ii)-vbeam_s(is)*vtilt_s(is)*(zmid+ztilt)*zleni
            enddo

c           --- do half velocity advance with E fields
            call epusht3d(nn,uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                    at(1),bt(1),apt(1),sq(is),sm(is),bpt(1),0.5)

c           --- do half velocity advance with B fields
            call bpusht3d(nn,uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                    gaminv(ins(is)),xt(1),yt(1),zt(1),sq(is),sm(is),
     &                    bpt(1),0.5,ibpush)

c           --- do full position advance
            call xpusht3d(nn,xp(ins(is)),yp(ins(is)),zp(ins(is)),
     &                    uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                    gaminv(ins(is)),bpt(1))

c           --- Calculate Gamma inverse
            if (lrelativ) then
              do ip=ins(is),ins(is)+nn-1
                gaminv(ip) = 1./sqrt(1. +
     &                            (uxp(ip)**2+uyp(ip)**2+uzp(ip)**2)*clightsqi)
              enddo
            else
              do ip=ins(is),ins(is)+nn-1
                gaminv(ip) = 1.
              enddo
            endif

c           --- setrho is now done by the setrho call in padvnc3d.

c         --- end of loop over species
          enddo

c         --- end of loop over injection sources
          enddo

        elseif (itask == 2) then
c Do second part of constant current injection: get new E fields (this is
c done after field solve including injected particles), synchronize velocity
c with position, and then move velocity one half timestep back to match time
c level of rest of particles.

c         loop over species and injection sources
          do is=1,ns
            ipmin = ins(is)
            nn = 0
          do ij=1,inject
            ipmin = ipmin + nn
c           --- Get number of particles injected
            nn = npinjtmp(is,ij)

c           --- Make sure there is enough temp space
            if (nn > npgrp) then
              npgrp = nn
              call gallot("Setpwork3d",0)
            endif

c           --- calculate new E self-fields
            call inj_sete3d(nn,ipmin,xmmin,ymmin,zmmin,dx,dy,dz,
     &                      at(1),bt(1),apt(1))

c           --- get external fields at current positions
            do ip=1,nn
              xt(ip) = 0.
              yt(ip) = 0.
              zt(ip) = 0.
              uxt(ip) = 0.
              uyt(ip) = 1.
            enddo
            do ip = 1,nn,NPARPGRP
              ii = min(NPARPGRP, nn+1-ip)
              call exteb3d(ii,xp(ip+ipmin-1),yp(ip+ipmin-1),zp(ip+ipmin-1),
     &                     uzp(ip+ipmin-1),gaminv(ip+ipmin-1),
     &                     -dt*0.5,0.,
     &                     xt(ip),yt(ip),zt(ip),at(ip),bt(ip),apt(ip),
     &                     sm(is),sq(is),uxt(ip),uyt(ip),1.,dt)
            enddo

c           --- Set fractional time to advance particles
            rnpinjct = 1./dvnz(real(nn,kind=8))
            do ip=1,nn
              bpt(ip) = (ip - .5)*rnpinjct*dt
            enddo

c           --- complete B advance
            call bpusht3d(nn,uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),xt(1),yt(1),zt(1),sq(is),sm(is),
     &                    bpt(1),0.5,ibpush)

c           --- complete the E advance
            call epusht3d(nn,uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    at(1),bt(1),apt(1),sq(is),sm(is),bpt(1),0.5)

c           --- Now, move velocites back one half a step
c           --- first half of a backward B advance
            call bpush3d(nn,uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                   gaminv(ipmin),xt(1),yt(1),zt(1),sq(is),sm(is),
     &                   -0.5*dt,ibpush)

c           --- then half of a backward E advance
            call epush3d(nn,uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                   at(1),bt(1),apt(1),sq(is),sm(is),-0.5*dt)

c        --- end of loop over species and injection sources
         enddo
         enddo

c Constant current injection is now complete for this time step.
        endif

c------------------------------------------------------------------------
c Space-Charge Limited injection

c Particles are initially placed on the emitting surface and given a
c time of emission uniformly distributed between 0 and dt. They are given
c an initial normal velocity equal to the Child-Langmuir velocity at dt
c minus the time of emission. This initial velocity gives the particles
c a small kick to get them moving which is needed since the normal field
c is zero at the surface, the particle starting location.  They are then
c advanced a partial isochronous leapfrog timestep, advancing dt minus the
c time of emission.  After the particles are injected, the position and
c velocity are at the same time level so, when this injection is used,
c every step is done with the split-leap frog algorithm.
c
c The voltage drop is calculated between the emitting surface and a
c secondary surface which is parallel too (is concentric too in the
c spherical case) the emitting surface.
c
c For inject == 2
c Space-charge limited injection, version 1. The number of particles injected is
c calculated from the field near the source using the Child-Langmuir relation.
c The number of particles loaded is J*dt*dx*dy/echarge/sw where J is the
c current density as calculated from the Child-Langmuir result.
c This algorithm works well for steady state flow, but does not give a
c good result for time-dependent flow.
c
c For inject == 3
c Space charge limited injection based off of the Gauss's Law.  This is
c the scheme to use with time-dependent injection.

      elseif (inject == 2 .or. inject == 3) then

        if (itask == 1) then
c         --- Zero npinje_s.
          npinje_s = 0

c         --- Copy inj_np to inj_prev, and zero inj_np
          inj_prev = inj_np
          inj_np = 0.

c         --- Calculate the charge density on the surface of the emitter.
          if (inject == 3) then
            if(solvergeom==XYZgeom) then
              call inj_setrho3d(rho(0,0,0),xmmin,ymmin,dx,dy,dz,nx,ny,nz,
     &                          l2symtry,l4symtry)
            elseif(solvergeom==RZgeom) then
              call inj_setrho3d_rz(rho(0,0,0),xmmin,dx,dz,nx,nz)
            endif
          endif

c         --- loop over injection sources
          do ij=1,ninject
            inj_ij = min(ij,inj_ninj)

c           --- Set some temporaries.
            ainj = ainject(ij)
            binj = binject(ij)
            ainji = 1./ainject(ij)
            binji = 1./binject(ij)
            ainjmini = 1./dvnz(ainjmin(ij))
            binjmini = 1./dvnz(binjmin(ij))
            xinj = xinject(ij)
            yinj = yinject(ij)
            apinj = apinject(ij)
            bpinj = bpinject(ij)
            xpinj = xpinject(ij)
            ypinj = ypinject(ij)
            rinji = 1./rinject(ij)

c           --- The potential drop in front of the surface is now
c           --- calculated by a call to getinj_phi from padvnc3d.

c           --- Set factor so correct number of particles are loaded when
c           --- symmetry is used.  With 2-fold, particles need to be loaded
c           --- on both sides of the symmetry axis.  With 4-fold, particles
c           --- need to be loaded in each of the four quadrants.
c           --- In the loop below, load_fac is adjust on the inner edge of the
c           --- grid so that particles arn't loaded into the same grid cell
c           --- multiple times.  The logic assumes that particles will not be
c           --- injected at edges of the grid that are not symmetry planes.
            load_fac = 1
            signx(1) = +1.
            signy(1) = +1.
            if (l2symtry) then
              load_fac = 2
              signx(2) = +1.
              signy(2) = -1.
            elseif (l4symtry .or. solvergeom == RZgeom) then
              load_fac = 4
              signx(2) = -1.
              signx(3) = -1.
              signx(4) = +1.
              signy(2) = -1.
              signy(3) = +1.
              signy(4) = -1.
            endif

            do is=1,ns

c             --- The factor to convert current density to particle number.
c             --- The absolute value is taken of sq so the factor is always
c             --- positive.
              jton(is) = dx*dy*dt/abs(sq(is))/sw(is)

c             --- This constant terms includes corrections for injection off
c             --- of a concentric spheres.  See I. Langmuir,
c             --- K. Blodgett, "Currents Limited by Space Charge Between
c             --- Concentric Spheres", PhysRev, 1924.
              if (inject == 2) then
                if (linj_spherical) then
                  zz = (dz*inj_d)*rinji
                  sphericalcorrection = (1. + 1.6*zz + 2.06*zz**2)
                else
                  sphericalcorrection = 1.
                endif
                cl_const(is) = 4./9.*eps0*(2.*abs(sq(is))/sm(is))**.5/
     &                        ((dz*inj_d)**2*sphericalcorrection)*jton(is)
              endif

c             --- The constant terms in the expression for the normal
c             --- velocity from the Child-Langmuir solution. Note that does
c             --- is only strictly correct for planar emission.
              vzconst(is)=2./9.*sqrt(2.)*(abs(sq(is))/sm(is))**1.5/(dz*inj_d)**2

            enddo

c           --- Load particles one grid cell at a time.
            do iy=0,inj_ny
              do ix=0,inj_nx

c               --- Mesh location
                xm = xmmin + ix*dx
                ym = ymmin + iy*dy

c               --- Only inject particles if within 2 grid cells of the source
c               --- and if the source is within the grid.
                if (((xm-xinj)*(binj+2.*dy))**2 +
     &              ((ym-yinj)*(ainj+2.*dx))**2 .gt.
     &              ((ainj+2.*dx)*(binj+2.*dy))**2 .or.
     &              (inj_grid(ix,iy) < -1 .or.
     &               inj_grid(ix,iy) >= nz)) cycle


c               --- Loop over the species
                do is=1,ns
                  inj_is = min(is,inj_ns)

c                 --- Skip this species if the fraction is zero.
                  if (finject(ij,is) == 0.) cycle

c                 --- number of particles injected in grid cell
                  if (inject == 2) then
c                   --- Child-Langmuir
                    rnn = cl_const(is)*abs(inj_phi(ix,iy))**1.5/
     &                    cos(inj_angl(ix,iy))
                    rnn = sign(rnn,sq(is)*inj_phi(ix,iy))
                  elseif (inject == 3) then
c                   --- Gauss's law
                    rnn = (eps0*inj_phi(ix,iy)/
     &                    (dz*inj_d) - inj_rho(ix,iy,inj_ij)*dz)*
     &                    dx*dy/sq(is)/sw(is)/cos(inj_angl(ix,iy))
                  endif

c                 --- Make sure that rnn is >= 0.
c                 --- Don't cycle at this point even if rnn is zero since
c                 --- averaging with the number from the previous step may
c                 --- give a positive value.
                  if (rnn < 0.) rnn = 0.

c                 --- Apply ad-hoc scaling factor
                  rnn = rnn*inj_f

c                 --- Average the number of particles for this step with
c                 --- that of the previous step. This helps the relaxation
c                 --- toward a steady state, especially in the Egun style
c                 --- iterative mode.
                  rnn = (inj_param*rnn +
     &                  (1. - inj_param)*inj_prev(ix,iy,inj_ij,inj_is))

c                 --- Force injected current to be between the range
c                 --- jmininj to jmaxinj.  The absolute values are used
c                 --- so that this works with either positive ions
c                 --- and electron.
                  rnn = min(rnn,jmaxinj*jton(is))
                  rnn = max(rnn,jmininj*jton(is))

c                 --- Scale the number of particles by the fraction of the
c                 --- total for this species.
                  rnn = rnn*finject(ij,is)

c                 --- Add a random number to the number of particles injected
c                 --- so that the average number of particles injected is
c                 --- correct.  For example, if rnn < 1., without the
c                 --- addition of the random number, no particles would ever
c                 --- be injected.  With the random number, particles will be
c                 --- injected and but the average number will be less than 1.
                  rnn = rnn + wranf()

c                 --- Save the number of particles injected
                  inj_np(ix,iy,inj_ij,inj_is) = inj_np(ix,iy,inj_ij,inj_is) + rnn

c                 --- Convert to integer.
                  nn = int(rnn)

c                 --- Only do following work if particles are to be injected
                  if (nn <= 0) cycle

c                 --- Make sure there is room for more particles.  Estimate
c                 --- total number of particles to be injected from the
c                 --- source. The second call makes it absolutely sure
c                 --- that there is enough room available in case the
c                 --- first estimate is wrong.
                  call chckpart(is,int(nn*ainj*binj*pi*dxi*dyi),0,.true.)
                  call chckpart(is,nn*load_fac+npinje_s(is),0,.true.)

c                 --- loop to load new particles
c                 --- In cases of symmetry, loop over the quadrants or
c                 --- half-planes.
                  ilf = load_fac
                  if (ix == 0) ilf = ilf/2
                  if (iy == 0) ilf = ilf/2
                  do il=1,ilf
                    do ip=1,nn
                      ii = ins(is) - npinje_s(is) - 1

c                     --- increment random number counter
                      injctcnt = injctcnt + 1

c                     --- calculate x, y and z of new particle
c                     --- if within injection source load it, else skip it
                      xp(ii)=signx(il)*xm+(rnrev(injctcnt,dig1)-.5)*dx*2.
                      yp(ii)=signy(il)*ym+(rnrev(injctcnt,dig2)-.5)*dy*2.

c                     --- Check if particle is within the elliptical annulus
c                     --- with outer major and minor radii of ainject and
c                     --- binject and inner radii of ainjmin and binjmin.
c                     --- The signs of x and y are adjusted to be in the
c                     --- first quadrant.
                      xx = xp(ii)*signx(il)-xinj
                      yy = yp(ii)*signy(il)-yinj
                      if ((xx*ainji)**2 + (yy*binji)**2 > 1. .or.
     &                    (xx*ainjmini)**2 + (yy*binjmini)**2 < 1.) cycle

c                     --- find injection surface axial location by interpolation
                      ixx = abs(xp(ii) - xmmin)*dxi
                      iyy = abs(yp(ii) - ymmin)*dyi
                      wxx = abs(xp(ii) - xmmin)*dxi - ixx
                      wyy = abs(yp(ii) - ymmin)*dyi - iyy
                      zp(ii)=zmmin+dz*(inj_grid(ixx  ,iyy  )*(1.-wxx)*(1.-wyy)+
     &                                 inj_grid(ixx+1,iyy  )*    wxx *(1.-wyy)+
     &                                 inj_grid(ixx  ,iyy+1)*(1.-wxx)*    wyy +
     &                                 inj_grid(ixx+1,iyy+1)*    wxx *    wyy)+
     &                       inj_zstart

c                     --- Only inject particles within the injection region.
                      if (zp(ii) < zinjmin .or. zp(ii) >= zinjmax) cycle

c                     --- find injection angle due to curvature by interpolation
                      az = inj_angl(ixx  ,iyy  )*(1. - wxx)*(1. - wyy) +
     &                     inj_angl(ixx+1,iyy  )*      wxx *(1. - wyy) +
     &                     inj_angl(ixx  ,iyy+1)*(1. - wxx)*      wyy  +
     &                     inj_angl(ixx+1,iyy+1)*      wxx *      wyy

c                     --- Calculate transverse angle
                      aa = atan2(yp(ii)-yinj,xp(ii)-xinj)

c                     --- Calculate and save the fraction of time step this
c                     --- particle is to be advanced. Also equal to dt minus
c                     --- the time of emission, all divided by dt.
c                     --- Note that the number is stored as a fraction between
c                     --- zero and one since pid does double duty, storing
c                     --- also the number of the injection source. Multiplying
c                     --- by dt would lose too many digits when added to an
c                     --- integer.
c                     --- The fractions are chosen using the digit reversed
c                     --- number generator and not particle number since they
c                     --- are correlated with transverse position.
                      pid(ii) = rnrev(injctcnt,dig5)

c                     --- Get the normal velocity for this particle. This is
c                     --- based on the Child-Langmuir solution, given a
c                     --- voltage drop, vnorm, and a time since emission,
c                     --- pid(ii)*dt.
                      vnorm = inj_phi(ixx  ,iyy  )*(1. - wxx)*(1. - wyy) +
     &                        inj_phi(ixx+1,iyy  )*      wxx *(1. - wyy) +
     &                        inj_phi(ixx  ,iyy+1)*(1. - wxx)*      wyy  +
     &                        inj_phi(ixx+1,iyy+1)*      wxx *      wyy
                      vznorm = vzconst(is)*vnorm**1.5*(pid(ii)*dt)**2

c                     --- Set transverse coordinates. The emittance term
c                     --- is probably not correct since it depends on vbeam.
c                     --- The correct way to add a transverse thermal spread
c                     --- is via vthperp_s. Also add in the transverse
c                     --- component of the normal velocity and the specified
c                     --- vbeam (which is also taken to be normal).
c                     --- Like the emittance term, the term with ap (and bp)
c                     --- are probably not correct.
                      vtx = vbeam_s(is)*0.5*emitx_s(is)*ainji + vthperp_s(is)
                      vty = vbeam_s(is)*0.5*emity_s(is)*binji + vthperp_s(is)
                      uxp(ii)=vbeam_s(is)*(signx(il)*apinj*(xp(ii)-xinj)*ainji
     &                        + xpinj) + vtx*rnorm()
     &                        - (vznorm + vbeam_s(is))*sin(az)*cos(aa)
                      uyp(ii)=vbeam_s(is)*(signy(il)*bpinj*(yp(ii)-yinj)*binji
     &                        + ypinj) + vty*rnorm()
     &                        - (vznorm + vbeam_s(is))*sin(az)*sin(aa)

c                     --- Include the user specified vbeam (assumed to be a
c                     --- normal), the thermal velocity, and the term from
c                     --- the Child-Langmuir solution.
c                     --- SMALLPOS is also added so that the field gathering
c                     --- routines don't treat the newly injected particles
c                     --- as dead particles when the other terms are zero.
                      uzp(ii) = vbeam_s(is)*cos(az) + vthz_s(is)*rnorm() +
     &                          vznorm*cos(az) + SMALLPOS

c                     --- Calculate gamma inverse and reset the particle
c                     --- massless momentums properly.
                      if (lrelativ) then
                        gaminv(ii) = sqrt(1.-
     &                          (uxp(ii)**2+uyp(ii)**2+uzp(ii)**2)*clightsqi)
                        gamma = 1./gaminv(ii)
                        uxp(ii) = uxp(ii)*gamma
                        uyp(ii) = uyp(ii)*gamma
                        uzp(ii) = uzp(ii)*gamma
                      else
                        gaminv(ii) = 1.
                      endif

c                     --- Save the ID of the injection source. Note that pid
c                     --- is doing double duty, saving the injection source
c                     --- number and the fraction of the timestep the
c                     --- particle is advanced.
                      pid(ii) = pid(ii) + ij

c                     --- increment particle counter
                      npinje_s(is) = npinje_s(is) + 1

                    enddo
c                   --- end of loop over particles
                  enddo
c                 --- end of loop for symmetries
                enddo
c               --- end loop over species
              enddo
c             --- end loop over ix
            enddo
c           --- end loop over iy
          enddo
c         --- end of loop over injection sources

          do is=1,ns
c           --- Sum npinje_s to get the total number of particles injected.
            npinject = npinject + npinje_s(is)

c           --- make injected particles live particles
            ins(is) = ins(is) - npinje_s(is)
            nps(is) = nps(is) + npinje_s(is)
          enddo

c         --- Now, advance the particles off of the emitting surface using
c         --- isochronous leapfrog.  Each particle has it's own time step
c         --- size, uniformly distributed between 0 and dt.
          do is=1,ns
c           --- Skip this if no particles of this species were injected.
            if(npinje_s(is)==0) cycle

c           --- Make sure that there is enough room in the temporary arrays.
            if (npinje_s(is) > npgrp) then
              npgrp = npinje_s(is)
              call gchange("Setpwork3d",0)
            endif

c           --- get E self-field at initial positions
            call inj_sete3d(npinje_s(is),ins(is),xmmin,ymmin,zmmin,
     &                      dx,dy,dz,at(1),bt(1),apt(1))

c           --- get external fields at initial positions
c           --- (Assumes that injection will never be done inside a bend.)
            do ip=1,npinje_s(is)
              xt(ip) = 0.
              yt(ip) = 0.
              zt(ip) = 0.
              uxt(ip) = 0.
              uyt(ip) = 1.
            enddo
            do ipmin = ins(is),ins(is)+npinje_s(is)-1,NPARPGRP
              ip = min(NPARPGRP, ins(is)+npinje_s(is)-ipmin)
              ii = ipmin - ins(is) + 1
              call exteb3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                     gaminv(ipmin),0.,dt*0.5,
     &                     xt(ii),yt(ii),zt(ii),at(ii),bt(ii),apt(ii),
     &                     sm(is),sq(is),uxt(ii),uyt(ii),1.,dt)
            enddo

c           --- Store the fractional timestep in bpt.
            do ip=ins(is),ins(is)+npinje_s(is)-1
              bpt(ip-ins(is)+1) = (pid(ip)-int(pid(ip)))*dt
            enddo

c           --- do half velocity advance with E fields
            call epusht3d(npinje_s(is),uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                    at(1),bt(1),apt(1),sq(is),sm(is),bpt(1),0.5)

c           --- do half velocity advance with B fields
            call bpusht3d(npinje_s(is),uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                    gaminv(ins(is)),xt(1),yt(1),zt(1),sq(is),sm(is),
     &                    bpt(1),0.5,ibpush)

c           --- do full position advance
            call xpusht3d(npinje_s(is),xp(ins(is)),yp(ins(is)),zp(ins(is)),
     &                    uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                    gaminv(ins(is)),bpt(1))

c           --- setrho is now done by the setrho call in padvnc3d.

          enddo

c         --- print warning if no particles were injected
          if (npinject == 0) then
            call remark("No particles were injected.")
          endif

        elseif (itask == 2) then
c Do second part of constant current injection: get new E fields (this is
c done after field solve including injected particles), synchronize velocity
c with position, and then move velocity one half timestep back to match time
c level of rest of particles.

c         --- loop over species
          do is=1,ns
            if(npinje_s(is)==0) cycle

c           --- Make sure that there is enough room in the temporary arrays.
            if (npinje_s(is) > npgrp) then
              npgrp = npinje_s(is)
              call gchange("Setpwork3d",0)
            endif

c           --- calculate new E self-fields
            call inj_sete3d(npinje_s(is),ins(is),xmmin,ymmin,zmmin,
     &                      dx,dy,dz,at(1),bt(1),apt(1))

c           --- Get external fields at current positions, assuming that there
c           --- is no bend.
            do ip=1,npinje_s(is)
              xt(ip) = 0.
              yt(ip) = 0.
              zt(ip) = 0.
              uxt(ip) = 0.
              uyt(ip) = 1.
            enddo
            do ip = 1,npinje_s(is),NPARPGRP
              ii = min(NPARPGRP, npinje_s(is)+1-ip)
              ipmin = ip+ins(is)-1
              call exteb3d(ii,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                     gaminv(ipmin),-dt*0.5,0.,
     &                     xt(ip),yt(ip),zt(ip),at(ip),bt(ip),apt(ip),
     &                     sm(is),sq(is),uxt(ip),uyt(ip),1.,dt)
            enddo

c           --- Store the fractional timestep in bpt.
            do ip=1,npinje_s(is)
              bpt(ip) = (pid(ip-1+ins(is)) - int(pid(ip-1+ins(is))))*dt
            enddo

c           --- complete B advance
            call bpusht3d(npinje_s(is),uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                    gaminv(ins(is)),xt(1),yt(1),zt(1),sq(is),sm(is),
     &                    bpt(1),0.5,ibpush)

c           --- complete the E advance
            call epusht3d(npinje_s(is),uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                    at(1),bt(1),apt(1),sq(is),sm(is),bpt(1),0.5)

c           --- Now, move velocites back one half a step
c           --- first half of a backward B advance
            call bpush3d(npinje_s(is),uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                   gaminv(ins(is)),xt(1),yt(1),zt(1),sq(is),sm(is),
     &                   -0.5*dt,ibpush)

c           --- then half of a backward E advance
            call epush3d(npinje_s(is),uxp(ins(is)),uyp(ins(is)),uzp(ins(is)),
     &                   at(1),bt(1),apt(1),sq(is),sm(is),-0.5*dt)

c        --- end of loop over species
         enddo

        endif

c Axially directed space-charge limited injection is now complete
c for this time step.
      endif

c------------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c Injection off of a tranverse facing surface.
c Uses any of the three methods implemented, based off the value of
c the variable 'inject'.

      if (ntinj .gt. 0 .and.
     &    1 .le. inject .and. inject .le. 3) then

        if (itask == 1) then

c         --- zero ntinject for inject = 2 or 3
          if (inject == 2 .or. inject == 3) then
            do is=1,ns
              ntinject(is) = 0
            enddo
          endif

c         --- Loop over injection sources
          do ij=1,ntinj

c           --- loop over species
            do is=1,ns
              nti = 0
              qoverm = sq(is)/sm(is)

c             --- Load particles one azimuthal section at a time,
c             --- only if there is a positive number of injected particles.
              do iz=nztinjmn(ij),nztinjmx(ij)-1
                zm = zmmin + iz*dz
                ijz = iz - nztinjmn(ij)
                rinj = abs(rtinject(ijz,ij))
                rinji = 1./abs(rtinject(ijz,ij))
                xx = sign(dx,rtinject(ijz,ij))
                xxi = 1./xx

c               --- Area of section of emitting surface.
                area = dz*rinj*2.*pi/nttinj(ij)

c               --- Const for inject = 2.  Denominator includes factor for
c               --- emitting on concentric cylinders.  See I. Langmuir,
c               --- K. Blodgett, "Currents Limited by Space Charge Between
c               --- Concentric Spheres", PhysRev, 1924.
                if (inject == 2) then
                  const = 4./9.*eps0*(2.*abs(qoverm))**.5*dxi**2*
     &             dt*area/echarge/sw(is)/
     &             (1. + 0.8*xx*rinji + 0.66*(xx*rinji)**2)
                endif

                do ith=0,nttinj(ij)-1

c                 --- angle of point in transverse plane
                  aa = ith*2.*pi/nttinj(ij)

c                 --- Find coordinates of the point a distance dx from
c                 --- the source along a line perpendicular to the
c                 --- emitting surface.
                  p2x = abs(xtinject(ij) + (rinj - xx)*cos(aa) - xmmin)*dxi
                  p2y = abs(ytinject(ij) + (rinj - xx)*sin(aa) - ymmin)*dyi
                  i2x = p2x
                  i2y = p2y
                  w2x = p2x - i2x
                  w2y = p2y - i2y

c                 --- Fetch difference between phi at that point and phi on
c                 --- the emitting surface.
                  ehalf = xxi*(vtinject(ij) -
     &                     phi(i2x  ,i2y  ,iz)*(1.-w2x)*(1.-w2y) -
     &                     phi(i2x+1,i2y  ,iz)*    w2x *(1.-w2y) -
     &                     phi(i2x  ,i2y+1,iz)*(1.-w2x)*    w2y  -
     &                     phi(i2x+1,i2y+1,iz)*    w2x *    w2y  )

c                 --- normal velocity
                  vnorm = -dt*qoverm*ehalf

c                 --- Calculate axial E field dx from surface.
                  ez = 0.5*dzi*(
     &     (phi(i2x  ,i2y  ,iz-1) - phi(i2x  ,i2y  ,iz+1))*(1.-w2x)*(1.-w2y) +
     &     (phi(i2x+1,i2y  ,iz-1) - phi(i2x+1,i2y  ,iz+1))*    w2x *(1.-w2y) +
     &     (phi(i2x  ,i2y+1,iz-1) - phi(i2x  ,i2y+1,iz+1))*(1.-w2x)*    w2y  +
     &     (phi(i2x+1,i2y+1,iz-1) - phi(i2x+1,i2y+1,iz+1))*    w2x *    w2y)

c                 --- Axial velocity
                  vznorm = dt*qoverm*ez

c                 --- Find coordinates of the point.
                  p2x = abs(xtinject(ij) + rinj*cos(aa) - xmmin)*dxi
                  p2y = abs(ytinject(ij) + rinj*sin(aa) - ymmin)*dyi
                  i2x = p2x
                  i2y = p2y
                  w2x = p2x - i2x
                  w2y = p2y - i2y

c                 --- number of particles injected in grid cell
                  if (inject == 1) then
                    rnn = ntinject(is)*ftinject(is,ij)/
     &                    real(nttinj(ij)*(nztinjmx(ij) - nztinjmn(ij)))
                  elseif (inject == 2) then
                    tinjprev(ith,ijz,ij) = ftinject(is,ij)*(
     &                                    inj_param*const*abs(ehalf*dx)**1.5 +
     &                                    (1. - inj_param)*tinjprev(ith,ijz,ij))
                    rnn = tinjprev(ith,ijz,ij)
                  elseif (inject == 3) then
c                   --- get charge density at injection point
                    rr = rho(i2x  ,i2y  ,iz)*(1.-w2x)*(1.-w2y) +
     &                   rho(i2x  ,i2y+1,iz)*(1.-w2x)*    w2y  +
     &                   rho(i2x+1,i2y  ,iz)*    w2x *(1.-w2y) +
     &                   rho(i2x+1,i2y+1,iz)*    w2x *    w2y
                    rnn=(eps0*ehalf - rr*dx)*area/sq(is)/sw(is)*ftinject(is,ij)
                  endif
                  rnn = min(rnn,
     &                 abs(jmaxtinj(ij)*area*dt/sq(is)/sw(is)*ftinject(is,ij)))
                  nn = int(rnn + wranf())

c                 --- Make sure there is room for more particles.  Estimate
c                 --- total number of particles to be injected from the
c                 --- source.
                  call chckpart(is,npinje_s(is) +
     &                             nn*nttinj(ij)*(nztinjmx(ij)-nztinjmn(ij)),
     &                          0,.true.)
                  call chckpart(is,nn,0,.true.)

                  do ip=1,nn
                    ii = ins(is) - npinje_s(is) - nti - 1
c                   --- calculate position of new particle
                    wz = rnrev(injctcnt,dig1)
                    zp(ii) = zm + wz*dz
                    rr = rnrev(injctcnt,dig2)
                    rp = rinj*(1.-wz) + abs(rtinject(ijz+1,ij))*wz - rr*xx
                    t = aa + (rnrev(injctcnt,dig3) - .5)*2.*pi/nttinj(ij)
                    xp(ii) = xtinject(ij) + rp*cos(t)
                    yp(ii) = ytinject(ij) + rp*sin(t)
                    vtx = 0.5*vbeam_s(is)*emitx_s(is)*rinji + vthperp_s(is)
                    vty = 0.5*vbeam_s(is)*emity_s(is)*rinji + vthperp_s(is)
                    uxp(ii) = vnorm*cos(t) + vtx*rnorm()
                    uyp(ii) = vnorm*sin(t) + vty*rnorm()
                    uzp(ii) = vthz_s(is)*rnorm() + vznorm*rr
                    if (lrelativ) then
                      gaminv(ii) = sqrt(1.-
     &                        (uxp(ii)**2+uyp(ii)**2+uzp(ii)**2)*clightsqi)
                      gamma = 1./gaminv(ii)
                      uxp(ii) = uxp(ii)*gamma
                      uyp(ii) = uyp(ii)*gamma
                      uzp(ii) = uzp(ii)*gamma
                    else
                      gaminv(ii) = 1.
                    endif
c                   --- increment particle counter
                    nti = nti + 1
c                   --- increment random number counter
                    injctcnt = injctcnt + 1
                  enddo
                enddo
              enddo

c             --- Set number of particles that were injected for inject = 2 or 3
              if (inject == 2 .or. inject == 3) then
                ntinject(is) = ntinject(is) + nti
              endif

c           --- end loop over species
            enddo

c         --- end loop over transverse injection sources
          enddo

c         --- Sum ntinject to get total number of particles injected.
          do is=1,ns
            npinject = npinject + ntinject(is)
          enddo


c         --- print warning if no particles were injected
          if (npinject == 0) then
            call remark("No particles were injected transversely.")
          endif

c         --- load new particles into rho
c         --- This must be done after all of the particles for this timestep
c         --- are injected.
          do is=1,ns
            do ipmin=ins(is)-npinje_s(is)-ntinject(is),
     &               ins(is)-npinje_s(is)-1,NPARPGRP
              ip = min(NPARPGRP, ins(is) - npinje_s(is) - ipmin)
              call setrho3d(rho(0,0,0),ip,xp(ipmin),yp(ipmin),zp(ipmin),zgrid,
     &                      uzp(ipmin),sq(is),sw(is),depos)
            enddo
          enddo

        elseif (itask == 3) then
c         --- Only include new particles now, after rest of particles have
c         --- gone through 'synchv' so position and velocity are at same time
c         --- level.

c         --- make injected particles live particles
          if (ntinj .gt. 0) then
            do is=1,ns
              ins(is) = ins(is) - ntinject(is)
              nps(is) = nps(is) + ntinject(is)
            enddo
          endif

        endif
      endif
c     --- End of transverse injection

      if (lw3dtimesubs) timeinject3d = timeinject3d + wtime() - substarttime

      return
      end
c=============================================================================
      subroutine injctint()
      use Subtimers
      use Constant
      use Beam_acc
      use InGen
      use InGen3d
      use InPart
      use InMesh3d
      use Particles
      use InjectVars
      use InjectVars3d
      use Setpwork3d
      use Picglb3d

      integer(ISZ):: is,ip,ij,nn,nnmax,ij1,ij2
      real(kind=8):: total_frac,dd,maxdxdy
      real(kind=8):: maxab1,maxab2,minab1,minab2
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c Make sure that space has been allocated for the arrays. This is needed
c in case the user changes ninject or ns but doesn't explicitly allocate
c the additional needed space.
      call gchange("InjectVars",0)

c Initialize beam size from initial envelope parameters if beam size not set
c One problem with this is that if the user ever wants apinject or bpinject
c to be zero while ap0 or bp0 not be zero.  The user would have to reset
c apinject or bpinject back to zero after the generate.
      do ij=1,ninject
        if (ainject(ij)  == 0.) ainject(ij)  = a0
        if (binject(ij)  == 0.) binject(ij)  = b0
        if (apinject(ij) == 0.) apinject(ij) = ap0
        if (bpinject(ij) == 0.) bpinject(ij) = bp0
      enddo

c Set the injection length and start, normally the length of one time step
c and at zimin, but can be set by user.
      if (leninjct == 0) leninjct = vbeam*dt

      if (inject == 1 .and. npinject == 0) then
c       --- Find the total number particles within the injection zone.
c       --- This gives the number of particles that will be injected
c       --- at each time step.
        do is=1,ns
          if (npinje_s(is) == 0) then
            do ip=ins(is),ins(is)+nps(is)-1
              if (zinject .lt. zp(ip) .and. zp(ip) .lt. zinject + leninjct)
     &          npinje_s(is) = npinje_s(is) + 1
            enddo
          endif
          npinject = npinject + npinje_s(is)
        enddo
      else
c       --- Calculate number of particles to inject for each species if
c       --- not set by the user.
        do is=1,ns
          if (npinje_s(is) == 0) then
            npinje_s(is) = npinject*sp_fract(is)
          endif
        enddo
      endif

c     --- Print warning if npinject=0
      if (npinject == 0) then
        call remark("ERROR: injctint, npinject is zero, no particles will
     &  be injected.")
      endif

c Each species is stored in seperate parts of the particle arrays, the size of
c which is based off the beam fraction of each species.  Check if there is
c enough space in each of the species blocks.
c Also set particle parameters: charge, mass and weight
      total_frac = 0.
      do is=1,ns
        total_frac = total_frac + sp_fract(is)
      enddo
      do is=1,ns
        call chckpart(is,int(injctspc*sp_fract(is)/total_frac+1),0,.true.)
        if (sq(is) == 0.) sq(is) = zion_s(is)*echarge
        if (sw(is) == 0.)
     &    sw(is)=(ibeam_s(is)*dt/(echarge*zion_s(is)))/npinje_s(is)*sp_fract(is)
        if (sm(is) == 0.) sm(is) = aion_s(is)*amu
      enddo

c Set up finject array if not set by the user (if sum of finject = 0).
c Give each injection source the same distribution of species.
      total_frac = 0.
      do is=1,ns
        do ij=1,ninject
          total_frac = total_frac + finject(ij,is)
        enddo
      enddo
      if (total_frac == 0.) then
        do is=1,ns
          do ij=1,ninject
            finject(ij,is) = sp_fract(is)
          enddo
        enddo
      endif

c Set the variable INJCTCNT
c This is used for loading the injection particles.  It ensures that the
c injected particles have new random numbers.
      injctcnt = nps(1) + 1 + randoffset

c Set up Setpwork3d arrays
c --- Estimate number of particles injected
      nnmax = 0
      do is=1,ns
        nn = 0
        do ij=1,ninject
          nn = nn + npinje_s(is)*finject(ij,is)
        enddo
        if (nn > nnmax) nnmax = nn
      enddo
      if (nnmax > npgrp) then
        npgrp = nnmax
        call gallot("Setpwork3d",0)
      endif

c Check for transverse space-charge limited injection.  Calculate number
c of azimuthal sections that the transverse emitting surface is broken into.
c That number is rounded to the nearest integer.  This must be done before
c the call to gchange("InjectVars3d",0) since ntting is used as an array
c dimension
c Also, find maximum length of surfaces, and set maximum current density
c to large number.
      do ij=1,ntinj
        nttinj(ij) = int(abs(rtinject(0,ij))*2.*pi/min(dx,dy) + .5)
        nttinjmax = max(nttinjmax,nttinj(ij))
        nztmax = max(nztmax, nztinjmx(ij) - nztinjmn(ij))
        if (jmaxtinj(ij) == 0.) jmaxtinj(ij) = LARGEPOS
      enddo

c Determine what the value of inj_ns should be. It can either be 1 or
c ns: if only one species will be injected from each source, the it is 1,
c otherwise, ns. This is needed in order that inj_prev and inj_np are set
c correctly. If multiple species are being injected from a single source,
c then inj_prev and inj_np must be saved for each species. In order to
c minimize the bookkeepping, if more than one species is injected from any
c source, than the two arrays are made big enough to save the data for all
c species for each source. If cases come up where this is too much wasted
c space, then the code can later be changed by only making inj_ns as big
c as is needed, i.e. the maximum number of species emitted from any one
c source.  
      inj_ns = 1
      do ij=1,ninject
        nn = 0
        do is=1,ns
          if (finject(ij,is) > 0.) nn = nn + 1
        enddo
        if (nn > 1) inj_ns = ns
      enddo

c Determine what the value if inj_ninj should be. It can either be 1 or
c ninject: if no sources are within 2 grid cells of another, then it is 1,
c otherwise ninject. This is needed since some parameters, inj_prev,
c inj_np, inj_area, and inj_rho, are stored on grids and must be known for
c each emitting source. If the sources are too close, then the data can
c overlap. In that case, the data is saved in seperate planes for each
c source. In order to minimize bookkeepping, if any two source (or more)
c sources are too close, then the arrays are made big enough so that all
c of the data is saved seperately for each source (whether or not they are
c too close to another). If cases comes up where this is too much wasted
c space, then the code can be changed so that inj_ninj is only made as big
c as needed, i.e. the maximum number of sources that are too close to each
c other.  
c The radius of each source is taken to be the maximum of a and b. This
c makes the checks much simpler but overly conservative.
      inj_ninj = 1
      maxdxdy = max(dx,dy)
      do ij1=1,ninject
        maxab1 = max(ainject(ij1),binject(ij1))
        minab1 = min(ainjmin(ij1),binjmin(ij1))
        do ij2=1,ninject
          if (ij1 == ij2) cycle
          dd = sqrt((xinject(ij1)-xinject(ij2))**2 +
     &              (yinject(ij1)-yinject(ij2))**2)
          maxab2 = max(ainject(ij2),binject(ij2))
          minab2 = min(ainjmin(ij2),binjmin(ij2))
c         --- Check if the outer edges are near each other. If not
c         --- go to next pair.
          if (maxab1 + maxab2 < dd - 2.*maxdxdy) cycle
c         --- Check if one is inside the other.
          if (dd + maxab1 < minab2 - 2.*maxdxdy) cycle
          if (dd + maxab2 < minab1 - 2.*maxdxdy) cycle
c         --- If it gets through the checks above, these two source are
c         --- too close
          inj_ninj = ninject
        enddo
      enddo

c Allocate inject arrays and fill with source location information.
c Note that gchange is called instead of gallot since some the of arrays
c in InjectVars3d may have been allocated by the user for specialized
c control of injection.
      if(solvergeom==XYZgeom) then
        inj_nx = nx
        inj_ny = ny
        call gchange("InjectVars3d",0)
        call fill_inj(dx,dy,dz,ix_axis,iy_axis)
      elseif(solvergeom==RZgeom) then
        inj_nx = nx
        inj_ny = nx
        call gchange("InjectVars3d",0)
        call fill_inj_rz(dx,dz,ix_axis,iy_axis)
      endif

c Allocate the particle index. The index keeps track of which source the
c particles are injected from.  It is also used to store the fractional
c timestep for each particle upon emission.
      npmaxi = npmax
      call gchange("Particles",0)

c Check if the voltages were specified.  If so, set so conductor will be
c set up for source.
      do ij=1,ninject
        if (vinject(ij) .ne. 0) lvinject = .true.
      enddo

c Set vbeamfrm to zero so that the grid does not move off of the source,
c causing a code crash.
      vbeamfrm = 0.

      if (lw3dtimesubs) timeinjctint = timeinjctint + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine fill_inj(dx,dy,dz,ix_axis,iy_axis)
      use Subtimers
      use InPart
      use InMesh3d
      use InjectVars
      use InjectVars3d
      real(kind=8):: dx,dy,dz
      integer(ISZ):: ix_axis,iy_axis

c     --- Fill the injection grid array.  A padding of 2 grid cells is also
c     --- set.
c     --- The array inj_grid holds the axial location of the injection source
c     --- as a function of the transverse coordinates.
c     --- The array inj_angl holds the angle due to the curvature of the source
c     --- as a function of the transverse coordinates.


      integer(ISZ):: ij,ix,iy
      real(kind=8):: xx,yy,zz,ainj,binj,ai_dx,bi_dy,dxi,dyi,dzi
      real(kind=8):: xl,yl,xu,yu,x1,x2,txl,tyl,txu,tyu
      real(kind=8):: ainji,binji
      real(kind=8):: area
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c Loop over the injection sources.
      do ij=1,ninject
c       --- Calculate the extent of the source plus guard cells.
        ainj    = ainject(ij)
        binj    = binject(ij)
        ainji   = 1./ainj
        binji   = 1./binj
        ai_dx = ainject(ij) + 2.*dx
        bi_dy = binject(ij) + 2.*dy
c       --- Calculate grid location of transverse center of source.
        zz = (zinject - zmmin)*dzi
c       --- Loop over transverse plane.
        do iy = 0,ny
          do ix = 0,nx
c           --- transverse coordinate relative to center of source
            xx = (ix - ix_axis)*dx - xinject(ij)
            yy = (iy - iy_axis)*dy - yinject(ij)
c           --- If the point is within source or guard cells, calculate both
c           --- axial distance and angle of curvature.
            if ((xx*bi_dy)**2 + (yy*ai_dx)**2 - (ai_dx*bi_dy)**2 <= 0.) then
              inj_id(ix,iy) = ij
c             --- Write the expression for inj_grid so as to avoid
c             --- subtraction of large, similar numbers when rinject is large.
c             inj_grid(ix,iy) = zz + (rinject(ij) -
c    &               sqrt(max(0.,rinject(ij)**2 - xx**2 - yy**2)))*dzi
              inj_grid(ix,iy) = zz + (xx**2 +  yy**2)/
     &         (rinject(ij) + sqrt(max(0.,rinject(ij)**2 - xx**2 - yy**2)))*dzi
c             inj_angl(ix,iy) = atan2(sqrt(xx**2 + yy**2),
c    &                             rinject(ij) - (inj_grid(ix,iy) - zz)*dz)
              inj_angl(ix,iy) = asin(sqrt(xx**2 + yy**2)/rinject(ij))

c             --- Calculate fraction of grid cell which is within the emitting
c             --- surface.  Used to scale the charge density on the emitting
c             --- surface.
c             --- The calculation is done with the emitting area scaled to
c             --- a unit circle.
c             --- The result is the fraction of the two by two block around
c             --- the grid point which is inside the emitting surface.

c             --- First, put the grid point in the first quadrant and scale it.
c             --- The max(0,) chop off any part that crosses the x or y axis.
              xl = max(0.,(abs(xx) - dx)*ainji)
              yl = max(0.,(abs(yy) - dy)*binji)
              xu =        (abs(xx) + dx)*ainji
              yu =        (abs(yy) + dy)*binji

c             --- If the two by two block is completely outside emitting
c             --- surface, set the area to zero.
              if (xl**2 + yl**2 > 1.) then
                area = 0.

c             --- If two by two block straddles edge of emitting surface, do
c             --- the calculation.
              elseif (xl**2 + yl**2 < 1. .and.
     &                                1. < xu**2 + yu**2) then
c               --- Calculate area of block within first quadrant.
                x1 = max(xl, sqrt(max(0.,1. - yu**2)))
                x2 = min(xu, sqrt(max(0.,1. - yl**2)))
                area = (x1 - xl)*(yu - yl) +
     &                 0.5*x2*sqrt(1. - x2**2) + 0.5*asin(x2) -
     &                 0.5*x1*sqrt(1. - x1**2) - 0.5*asin(x1) -
     &                 (x2 - x1)*yl
c               --- If block extended past x axis, calculate the area of
c               --- that piece.
                if (abs(xx) - dx < 0.) then
                  txu = - (abs(xx) - dx)*ainji
                  txl = 0.
                  x1 = max(txl, sqrt(max(0.,1. - yu**2)))
                  x2 = min(txu, sqrt(max(0.,1. - yl**2)))
                  area = area + (x1 - txl)*(yu - yl) +
     &                          0.5*x2*sqrt(1. - x2**2) + 0.5*asin(x2) -
     &                          0.5*x1*sqrt(1. - x1**2) - 0.5*asin(x1) -
     &                          (x2 - x1)*yl
                endif
c               --- If block extended past y axis, calculate the area of
c               --- that piece.
                if (abs(yy) - dy < 0.) then
                  tyu = - (abs(yy) - dy)*binji
                  tyl = 0.
                  x1 = max(xl, sqrt(max(0.,1. - tyu**2)))
                  x2 = min(xu, sqrt(max(0.,1. - tyl**2)))
                  area = area + (x1 - xl)*(tyu - tyl) +
     &                          0.5*x2*sqrt(1. - x2**2) + 0.5*asin(x2) -
     &                          0.5*x1*sqrt(1. - x1**2) - 0.5*asin(x1)
                endif

c               --- Normalize the area so that the total area of the block
c               --- would be 1 (since area was calculated above on a
c               --- unit circle).
                area = area*dxi*dyi*ainj*binj*0.25

c               --- If the block is completely within the emitting surface,
c               --- set the area to 1.
              else
                area = 1.
              endif

c             --- Take into account the tent function of the charge density,
c             --- i.e. convert the integral of the area to an integral
c             --- of the tent function over the part of the two by two block
c             --- within the emitting surface.  This is an approximate
c             --- calculation.  The conversion was obtained by completely
c             --- integrating over one dimension and partially intregrating
c             --- over the other.
              if (area <= 0.5) then
                area = 2.*area**2
              else
                area = 4.*area - 2.*area**2 -1.
              endif

c             --- Now assign to inj_area. This is done here to simplify the
c             --- coding above, primarily because of the third index here.
c             --- Note that the min works, since inj_ninj is either 1 or
c             --- ninject: if 1, then the min is always 1, if ninject, then
c             --- the min is always ij.
              inj_area(ix,iy,min(ij,inj_ninj)) = area

            endif
          enddo
        enddo
      enddo

      if (lw3dtimesubs) timefill_inj = timefill_inj + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine fill_inj_rz(dx,dz,ix_axis,iy_axis)
      use Subtimers
      use InPart
      use InMesh3d
      use InjectVars
      use InjectVars3d
      real(kind=8):: dx,dz
      integer(ISZ):: ix_axis,iy_axis

c     --- Fill the injection grid array.  A padding of 2 grid cells is also
c     --- set.
c     --- The array inj_grid holds the axial location of the injection source
c     --- as a function of the transverse coordinates.
c     --- The array inj_angl holds the angle due to the curvature of the source
c     --- as a function of the transverse coordinates.


      integer(ISZ):: ij,ix,iy
      real(kind=8):: xx,yy,zz,ainj,binj,ai_dx,bi_dy,dxi,dyi,dzi
      real(kind=8):: xl,yl,xu,yu,x1,x2
      real(kind=8):: ainji,binji
      real(kind=8):: area
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./dx
      dyi = 1./dx
      dzi = 1./dz

c Loop over the injection sources.
      do ij=1,ninject
c       --- Calculate the extent of the source plus guard cells.
        ainj    = ainject(ij)
        binj    = binject(ij)
        ainji   = 1./ainj
        binji   = 1./binj
        ai_dx = ainject(ij) + 2.*dx
        bi_dy = binject(ij) + 2.*dx
c       --- Calculate grid location of transverse center of source.
        zz = (zinject - zmmin)*dzi
c       --- Loop over transverse plane.
        do iy = 0,inj_ny
          do ix = 0,inj_nx
c           --- transverse coordinate relative to center of source
            xx = (ix - ix_axis)*dx - xinject(ij)
            yy = (iy - iy_axis)*dx - yinject(ij)
c           --- If the point is within source or guard cells, calculate both
c           --- axial distance and angle of curvature.
            if ((xx*bi_dy)**2 + (yy*ai_dx)**2 - (ai_dx*bi_dy)**2 .le. 0.) then
              inj_id(ix,iy) = ij
c             --- Write the expression for inj_grid so as to avoid
c             --- subtraction of large, similar numbers when rinject is large.
c             inj_grid(ix,iy) = zz + (rinject(ij) -
c    &               sqrt(max(0.,rinject(ij)**2 - xx**2 - yy**2)))*dzi
              inj_grid(ix,iy) = zz + (xx**2 +  yy**2)/
     &         (rinject(ij) + sqrt(max(0.,rinject(ij)**2 - xx**2 - yy**2)))*dzi
              inj_angl(ix,iy) = atan2(sqrt(xx**2 + yy**2),
     &                             rinject(ij) - (inj_grid(ix,iy) - zz)*dz)

c             --- Calculate fraction of grid cell which is within the emitting
c             --- surface.  Used to scale the charge density on the emitting
c             --- surface.
c             --- The calculation is done with the emitting area scaled to
c             --- a unit circle.
c             --- The result is the fraction of the two by two block around
c             --- the grid point which is inside the emitting surface.

c             --- First, put the grid point in the first quadrant and scale it.
c             --- The max(0,) chop off any part that crosses the x or y axis.
              xl = max(0.,(abs(xx) - dx)*ainji)
              yl = max(0.,(abs(yy) - dx)*binji)
              xu =        (abs(xx) + dx)*ainji
              yu =        (abs(yy) + dx)*binji

c             --- If the two by two block is completely outside emitting
c             --- surface, set the area to zero.
              if (xl**2 + yl**2 .gt. 1.) then
                area = 0.

c             --- If two by two block straddles edge of emitting surface, do
c             --- the calculation.
              elseif (xl**2 + yl**2 .lt. 1. .and.
     &                                   1. .lt. xu**2 + yu**2) then
c               --- Calculate area of block within first quadrant.
                x1 = max(xl, sqrt(max(0.,1. - yu**2)))
                x2 = min(xu, sqrt(max(0.,1. - yl**2)))
                area = (x1 - xl)*(yu - yl) +
     &                 0.5*x2*sqrt(1. - x2**2) + 0.5*asin(x2) -
     &                 0.5*x1*sqrt(1. - x1**2) - 0.5*asin(x1) -
     &                 (x2 - x1)*yl

c               --- Normalize the area so that the total area of the block
c               --- would be 1 (since area was calculated above on a
c               --- unit circle).
                area = area*dxi*dyi*ainj*binj*0.25

c               --- If the block is completely within the emitting surface,
c               --- set the area to 1.
              else
                area = 1.
              endif

c             --- Take into account the tent function of the charge density,
c             --- i.e. convert the integral of the area to an integral
c             --- of the tent function over the part of the two by two block
c             --- within the emitting surface.  This is an approximate
c             --- calculation.  The conversion was obtained by completely
c             --- integrating over one dimension and partially intregrating
c             --- over the other.
              if (area .le. 0.5) then
                area = 2.*area**2
              else
                area = 4.*area - 2.*area**2 -1.
              endif

c             --- Now assign to inj_area. This is done here to simplify the
c             --- coding above, primarily because of the third index here.
c             --- Note that the min works, since inj_ninj is either 1 or
c             --- ninject: if 1, then the min is always 1, if ninject, then
c             --- the min is always ij.
              inj_area(ix,iy,min(ij,inj_ninj)) = area

            endif
          enddo
        enddo
      enddo

      if (lw3dtimesubs) timefill_inj = timefill_inj + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine getinj_phi(nx,ny,nz,phi,dx,dy,dz,xmmin,ymmin)
      use InjectVars
      use InjectVars3d
#ifdef MPIPARALLEL
      use Parallel
#endif
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz,xmmin,ymmin

c Calculate potential drop from emitting surface at distance
c of dz from the surface.  This is only done for grid cells
c within two grid cells of the elliptical emitting surface, and
c within the axial extent of the grid. This is done over the
c full axial extent since points which are outside of the
c injection region maybe needed for the interlopation below
c to get zp.

      integer(ISZ):: ij,ix,iy,i2x,i2y,i2z
      real(kind=8):: dxi,dyi
      real(kind=8):: ainj,binj,ainjsqi,binjsqi,xinj,yinj
      real(kind=8):: rrsq,aa,p2x,p2y,p2z,w2x,w2y,w2z

#ifdef MPIPARALLEL
c     --- In the parallel version, only processor 0 calculates inj_phi.
      if (my_index == 0) then
#endif

      dxi = 1./dx
      dyi = 1./dy

c     --- loop over injection sources
      do ij=1,ninject

c       --- Set some temporaries.
        ainj = ainject(ij)
        binj = binject(ij)
        ainjsqi = 1./(ainject(ij) + 2.*dx)**2
        binjsqi = 1./(binject(ij) + 2.*dy)**2
        xinj = xinject(ij)
        yinj = yinject(ij)

        do iy = 0,ny
          do ix = 0,nx
            rrsq = (ix*dx + xmmin - xinj)**2*ainjsqi +
     &             (iy*dy + ymmin - yinj)**2*binjsqi
            if (rrsq < 1. .and.
     &          (-1 <= inj_grid(ix,iy) .and.
     &                 inj_grid(ix,iy) < nz)) then

c             --- angle of point in transverse plane
              aa = atan2(iy*dy+ymmin-yinj,dvnz(ix*dx+xmmin-xinj))

c             --- Find coordinates of the point a distance dz in front
c             --- of the source along a line perpendicular to the
c             --- emitting surface.
              p2x = ix-dz*cos(aa)*sin(inj_angl(ix,iy))*dxi*inj_d
              p2y = iy-dz*sin(aa)*sin(inj_angl(ix,iy))*dyi*inj_d
              p2z = inj_grid(ix,iy) + cos(inj_angl(ix,iy))*inj_d
              i2x = p2x
              i2y = p2y
              i2z = p2z
              w2x = p2x - i2x
              w2y = p2y - i2y
              w2z = p2z - i2z

c             --- Fetch difference between phi at that point and phi on
c             --- the emitting surface.
              inj_phi(ix,iy) = vinject(ij) -
     &                 phi(i2x  ,i2y,  i2z  )*(1.-w2x)*(1.-w2y)*(1.-w2z) -
     &                 phi(i2x+1,i2y,  i2z  )*    w2x *(1.-w2y)*(1.-w2z) -
     &                 phi(i2x  ,i2y+1,i2z  )*(1.-w2x)*    w2y *(1.-w2z) -
     &                 phi(i2x+1,i2y+1,i2z  )*    w2x *    w2y *(1.-w2z) -
     &                 phi(i2x  ,i2y,  i2z+1)*(1.-w2x)*(1.-w2y)*    w2z  -
     &                 phi(i2x+1,i2y,  i2z+1)*    w2x *(1.-w2y)*    w2z  -
     &                 phi(i2x  ,i2y+1,i2z+1)*(1.-w2x)*    w2y *    w2z  -
     &                 phi(i2x+1,i2y+1,i2z+1)*    w2x *    w2y *    w2z
            endif
          enddo
        enddo

c       --- Smooth out inj_phi for this emitting surface.
        if (inj_nsmooth > 0) then
          call inj_smoother(nx,ny,inj_phi,dx,dy,xmmin,ymmin,
     &                      xinj,yinj,ainj,binj,inj_nsmooth)
        endif

      enddo

#ifdef MPIPARALLEL
      endif
c     --- The result calculated by processor 0 is broadcast to the rest
c     --- of the processors.
      call parallelbroadcast(inj_phi,(nx+1)*(ny+1),0)
#endif

c     --- Calculate the transverse fields
      if (linj_eperp) then
        do iy = 1,ny-1
          do ix = 1,nx-1
            inj_ex(ix,iy) = (inj_phi(ix-1,iy  ) - inj_phi(ix+1,iy  ))*dxi*0.5
            inj_ey(ix,iy) = (inj_phi(ix  ,iy-1) - inj_phi(ix  ,iy+1))*dyi*0.5
          enddo
        enddo
        do ix = 1,nx-1
          inj_ex(ix,0) = (inj_phi(ix-1,0) - inj_phi(ix+1,0))*dxi*0.5
        enddo
        do iy = 1,ny-1
          inj_ey(0,iy) = (inj_phi(0,iy-1) - inj_phi(0,iy+1))*dyi*0.5
        enddo
      endif

      return
      end
c=============================================================================
      subroutine getinj_phi_rz(nx,nz,phi,dx,dz,xmmin,ymmin)
      use InjectVars
      use InjectVars3d
#ifdef MPIPARALLEL
      use Parallel
#endif
      integer(ISZ):: nx,nz
      real(kind=8):: phi(0:nx,0:nz)
      real(kind=8):: dx,dz,xmmin,ymmin

c Calculate potential drop from emitting surface at distance
c of dz from the surface.  This is only done for grid cells
c within two grid cells of the elliptical emitting surface, and
c within the axial extent of the grid. This is done over the
c full axial extent since points which are outside of the
c injection region maybe needed for the interlopation below
c to get zp.

      integer(ISZ):: ij,ix,iy,i2r,i2z
      real(kind=8):: dxi,dyi
      real(kind=8):: ainjsqi,binjsqi,xinj,yinj
      real(kind=8):: rrsq,aa,p2x,p2y,p2r,p2z,w2r,w2z

#ifdef MPIPARALLEL
c     --- In the parallel version, only processor 0 calculates inj_phi.
      if (my_index == 0) then
#endif

      dxi = 1./dx
      dyi = 1./dx

c     --- loop over injection sources
      do ij=1,ninject

c       --- Set some temporaries.
        ainjsqi = 1./(ainject(ij) + 2.*dx)**2
        binjsqi = 1./(binject(ij) + 2.*dx)**2
        xinj = xinject(ij)
        yinj = yinject(ij)

        do iy = 0,nx
          do ix = 0,nx
            rrsq = (ix*dx + xmmin - xinj)**2*ainjsqi +
     &             (iy*dx + ymmin - yinj)**2*binjsqi
            if (rrsq < 1. .and.
     &          (-1 <= inj_grid(ix,iy) .and.
     &                 inj_grid(ix,iy) < nz)) then

c             --- angle of point in transverse plane
              aa = atan2(iy*dx+ymmin-yinj,dvnz(ix*dx+xmmin-xinj))

c             --- Find coordinates of the point a distance dz in front
c             --- of the source along a line perpendicular to the
c             --- emitting surface.
              p2x = ix-dz*cos(aa)*sin(inj_angl(ix,iy))*dxi*inj_d
              p2y = iy-dz*sin(aa)*sin(inj_angl(ix,iy))*dyi*inj_d
              p2r = sqrt(p2x**2+p2y**2)
              p2z = inj_grid(ix,iy) + cos(inj_angl(ix,iy))*inj_d
              i2r = p2r
              i2z = p2z
              w2r = p2r - i2r
              w2z = p2z - i2z

c             --- Fetch difference between phi at that point and phi on
c             --- the emitting surface.
              inj_phi(ix,iy) = vinject(ij) -
     &                 phi(i2r  ,i2z  )*(1.-w2r)*(1.-w2z) -
     &                 phi(i2r+1,i2z  )*    w2r *(1.-w2z) -
     &                 phi(i2r  ,i2z+1)*(1.-w2r)*    w2z  -
     &                 phi(i2r+1,i2z+1)*    w2r *    w2z
            endif
          enddo
        enddo

c       --- Smooth out inj_phi for this emitting surface.
c        if (inj_nsmooth > 0) then
c          call inj_smoother(nx,ny,inj_phi,dx,dy,xmmin,ymmin,
c     &                      xinj,yinj,ainj,binj,inj_nsmooth)
c        endif

      enddo

#ifdef MPIPARALLEL
      endif
c     --- The result calculated by processor 0 is broadcast to the rest
c     --- of the processors.
      call parallelbroadcast(inj_phi,(nx+1)*(nx+1),0)
#endif

c     --- Calculate the transverse fields
      if (linj_eperp) then
        do iy = 1,nx-1
          do ix = 1,nx-1
            inj_ex(ix,iy) = (inj_phi(ix-1,iy  ) - inj_phi(ix+1,iy  ))*dxi*0.5
            inj_ey(ix,iy) = (inj_phi(ix  ,iy-1) - inj_phi(ix  ,iy+1))*dyi*0.5
          enddo
        enddo
        do ix = 1,nx-1
          inj_ex(ix,0) = (inj_phi(ix-1,0) - inj_phi(ix+1,0))*dxi*0.5
        enddo
        do iy = 1,nx-1
          inj_ey(0,iy) = (inj_phi(0,iy-1) - inj_phi(0,iy+1))*dyi*0.5
        enddo
      endif

      return
      end
c=============================================================================
      subroutine inj_smoother(nx,ny,inj_phi,dx,dy,xmmin,ymmin,
     &                        xinj,yinj,ainj,binj,inj_nsmooth)
      integer(ISZ):: nx,ny
      real(kind=8):: inj_phi(0:nx,0:ny)
      real(kind=8):: dx,dy,xmmin,ymmin,xinj,yinj,ainj,binj
      integer(ISZ):: inj_nsmooth

c This routine smoothes the normal electric field in front of the emitting
c surface, hopefully removing defects  at places where the emitting surface
c crosses grid lines.
c It holds the field at the center point fixed. For points outside of the
c emitter, a linear extrapolation is done using the results of a least
c square fit of the field as a function of radius for the points within
c a grid cell of the edge. The points outside are also held fixed.

      integer(ISZ):: ix,iy,ix0,iy0,is
      integer(ISZ):: ixm1,iym1,ixp1,iyp1
      real(kind=8):: rrsq1,rrsq2
      real(kind=8):: isum,fsum,xsum,xsqsum,fxsum,c1,c2
      real(kind=8):: s(0:nx,0:ny)
      real(kind=8):: ainjsqi1,binjsqi1
      real(kind=8):: ainjsqi2,binjsqi2

c     --- Calculate the parameters for the least square fit.
      isum = 0.
      fsum = 0.
      xsum = 0.
      xsqsum = 0.
      fxsum = 0.
      ainjsqi1 = 1./ainj**2
      binjsqi1 = 1./binj**2
      ainjsqi2 = 1./(ainj - dx)**2
      binjsqi2 = 1./(binj - dx)**2
      do iy = 0,ny
        do ix = 0,nx
          rrsq1 = (ix*dx + xmmin - xinj)**2*ainjsqi1 +
     &            (iy*dy + ymmin - yinj)**2*binjsqi1
          rrsq2 = (ix*dx + xmmin - xinj)**2*ainjsqi2 +
     &            (iy*dy + ymmin - yinj)**2*binjsqi2
          if (rrsq1 < 1. .and. rrsq2 > 1.) then
            isum = isum + 1.
            fsum = fsum + inj_phi(ix,iy)
            xsum = xsum + sqrt(rrsq1)
            xsqsum = xsqsum + rrsq1
            fxsum = fxsum + sqrt(rrsq1)*inj_phi(ix,iy)
          endif
        enddo
      enddo
      c1 = (fsum*xsqsum/xsum - fxsum)/(isum*xsqsum/xsum - xsum)
      c2 = (fsum - c1*isum)/xsum

c     --- Set values for points just outside of the emitter edge, using
c     --- the equation from the fit above.
      ainjsqi2 = 1./(ainj + 2.*dx)**2
      binjsqi2 = 1./(binj + 2.*dx)**2
      do iy = 0,ny
        do ix = 0,nx
          rrsq1 = (ix*dx + xmmin - xinj)**2*ainjsqi1 +
     &            (iy*dy + ymmin - yinj)**2*binjsqi1
          rrsq2 = (ix*dx + xmmin - xinj)**2*ainjsqi2 +
     &            (iy*dy + ymmin - yinj)**2*binjsqi2
          if (rrsq1 > 1. .and. rrsq2 < 1.) then
            inj_phi(ix,iy) = c1 + c2*sqrt(rrsq1)
          endif
        enddo
      enddo

c     --- Do the smoothing.
      ix0 = (xinj - xmmin)/dx
      iy0 = (yinj - ymmin)/dy
      do is = 1,inj_nsmooth
        s = inj_phi
        do iy = 0,ny
          do ix = 0,nx
            rrsq1 = (ix*dx + xmmin - xinj)**2*ainjsqi1 +
     &              (iy*dy + ymmin - yinj)**2*binjsqi1
            if (rrsq1 < 1. .and. (ix .ne. ix0 .or. iy .ne. iy0)) then
              ixm1 = ix - 1
              ixp1 = ix + 1
              iym1 = iy - 1
              iyp1 = iy + 1
              if (ix == 0) ixm1 = 1
              if (ix == nx) ixp1 = nx-1
              if (iy == 0) iym1 = 1
              if (iy == ny) iyp1 = ny-1
              inj_phi(ix,iy) = 0.0625*(s(ixm1,iym1) + s(ixm1,iyp1)  +
     &                                 s(ixp1,iym1) + s(ixp1,iyp1)) +
     &                         0.1250*(s(ixm1,iy  ) + s(ixp1,iy  )  +
     &                                 s(ix  ,iym1) + s(ix  ,iyp1)) +
     &                         0.2500*(s(ix  ,iy  ))
            endif
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine inj_sete3d(npart,ipmin,xmmin,ymmin,zmmin,dx,dy,dz,ex,ey,ez)
      use Subtimers
      use InPart
      use InjectVars
      use InjectVars3d
      use Particles
      integer(ISZ):: nx,ny,nz,npart,ipmin
      real(kind=8):: ex(npart),ey(npart),ez(npart)
      real(kind=8):: xmmin,ymmin,zmmin,dx,dy,dz

c Calculate the electric field from the grid for particles near the emitting
c surface.  The normal electric field is calculated from the potential drop
c across a length of 'dz' along a line normal to the emitting surface.  The
c field components are then obtained from the normal field. Optionally,
c the tangential fields near the emitting surface can also be included.
c Eventually this routine could be expanded to calculate the E field near
c any conductor surface.

      integer(ISZ):: ip,ix,iy,iz,ij
      real(kind=8):: xinj,yinj,zinj,rinj,rinji,zz,rr2,dd,ca,zii
      real(kind=8):: aa,az,atx,aty,px,py,pz,wx,wy,wz,en,etx,ety
      real(kind=8):: dxi,dyi,dzi,inj_di
      real(kind=8):: sx,sy
      real(kind=8):: fourthirds,onethird
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      zii = 1./(dz*inj_d)
      inj_di = 1./inj_d
      fourthirds = 4./3.
      onethird = 1./3.

c     --- Loop over particles
      do ip=ipmin,ipmin+npart-1

c       --- Get number of injection source of particle.  The particle id
c       --- has the number of the injection source the particle was
c       --- emitted from. The int is needed since it is doing double
c       --- duty and has other information stored in the fractional part.
        ij = int(pid(ip))

c       --- set temporaries
        xinj = xinject(ij)
        yinj = yinject(ij)
        zinj = zinject
        rinj = rinject(ij)
        rinji = 1./rinj

c       --- Calculate distance of particle from the emitting surface
        zz = zp(ip) - zinj
        rr2 = (xp(ip) - xinj)**2 + (yp(ip) - yinj)**2
        if (zz < rinj) then
c         --- The expression below for zz<rinj is mathematically identical
c         --- to the commented out expression below.  That form is
c         --- used so that for large values of rinj, i.e. a flat emitting
c         --- surface, the correct value of the distance is calculated,
c         --- namely zz.
c         dd = rinj - sqrt(rr2 + (rinj - zz)**2)
          dd = (2.*zz - (zz**2 + rr2)*rinji)/
     &         (1. + sqrt(rr2*rinji**2 + (1. - zz*rinji)**2))
        else
c         --- When zz>rinj, a seperate equation is needed.
          dd = rinj + sqrt(rr2 + (zz - rinj)**2)
        endif
        dd = dd*dzi

c       --- Only calculate E-field if particle close to emitting surface.
        if (0.0 <= dd .and. dd <= inj_d) then
c         --- The check if a particle is within the radial area of the
c         --- source is probably not needed. It has been removed since there
c         --- are some cases which this check is incorrect - for example
c         --- in a source with a diverging beam. The particles spread out
c         --- beyond the edge of the source and would not pass this check.
c    &    (((xp(ip)-xinj)*binject(ij))**2 + ((yp(ip)-yinj)*ainject(ij))**2 .le.
c    &    (ainject(ij)*binject(ij))**2) .and.
c    &    (((xp(ip)-xinj)*binjmin(ij))**2 + ((yp(ip)-yinj)*ainjmin(ij))**2 .ge.
c    &    (ainjmin(ij)*binjmin(ij))**2) ) then

c         --- angle of point in transverse plane
          aa = atan2(yp(ip)-yinj,dvnz(xp(ip)-xinj))

c         --- angle relative to z axis
          az = asin(sqrt(rr2)*rinji)

c         --- Find coordinates of the point on the secondary surface in front
c         --- of the emitting surface along a line perpendicular to the
c         --- emitting surface.
          px = abs(rinj*sin(az)*cos(aa) + xinj - xmmin)*dxi
          py = abs(rinj*sin(az)*sin(aa) + yinj - ymmin)*dyi
          ix = int(px)
          iy = int(py)
          wx = px - ix
          wy = py - iy

c         --- Calculate the normal electric field from the potential drop in
c         --- front of the emitting surface. The normal field given by the
c         --- Child-Langmuir solution (for planar surfaces) is used.
          en = zii*(inj_phi(ix  ,iy  )*(1.-wx)*(1.-wy) +
     &              inj_phi(ix+1,iy  )*    wx *(1.-wy) +
     &              inj_phi(ix  ,iy+1)*(1.-wx)*    wy  +
     &              inj_phi(ix+1,iy+1)*    wx *    wy   )
          if (linj_enormcl) then
c           --- Scale the normal E field to match the Child-Langmuir solution
            en = en*fourthirds*(dd*inj_di)**onethird
          endif

c         --- Set the particle's electric field based off of the normal
c         --- electric fields.
          ex(ip-ipmin+1) = -en*sin(az)*cos(aa)
          ey(ip-ipmin+1) = -en*sin(az)*sin(aa)
          ez(ip-ipmin+1) =  en*cos(az)

          if (linj_eperp) then
c           --- Add in the tangential fields if requested.
c           --- Fetch the transverse field components. The field is linearly
c           --- interpolated between the value of the emitting surface
c           --- (i.e. zero) and the value on the secondary surface.
            atx = asin((xp(ip)-xinj)/sqrt(rinj**2 - (yp(ip)-yinj)**2))
            aty = asin((yp(ip)-yinj)/sqrt(rinj**2 - (xp(ip)-xinj)**2))
            sx = 1.
            sy = 1.
            if (xp(ip) < xmmin) sx = -1.
            if (yp(ip) < ymmin) sy = -1.
            etx = (dd/inj_d)*sx*(inj_ex(ix  ,iy  )*(1.-wx)*(1.-wy) +
     &                           inj_ex(ix+1,iy  )*    wx *(1.-wy) +
     &                           inj_ex(ix  ,iy+1)*(1.-wx)*    wy  +
     &                           inj_ex(ix+1,iy+1)*    wx *    wy   )
            ety = (dd/inj_d)*sy*(inj_ey(ix  ,iy  )*(1.-wx)*(1.-wy) +
     &                           inj_ey(ix+1,iy  )*    wx *(1.-wy) +
     &                           inj_ey(ix  ,iy+1)*(1.-wx)*    wy  +
     &                           inj_ey(ix+1,iy+1)*    wx *    wy   )

            ex(ip-ipmin+1) = ex(ip-ipmin+1) + etx*cos(atx)
            ey(ip-ipmin+1) = ey(ip-ipmin+1) + ety*cos(aty)
            ez(ip-ipmin+1) = ez(ip-ipmin+1) + etx*sin(atx) + ety*sin(aty)
          endif

        endif
      enddo

      if (lw3dtimesubs) timeinj_sete3d = timeinj_sete3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine inj_setrho3d(rho,xmmin,ymmin,dx,dy,dz,nx,ny,nz,
     &                        l2symtry,l4symtry)
      use Subtimers
      use InPart
      use InjectVars
      use InjectVars3d
      use Particles
      integer(ISZ):: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: xmmin,ymmin,dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry

c Calculate the charge density on the emitting surface.  The particles are
c mapped onto the surface assuming a spherical emitter.  The charge density
c scale to account for symmetries and to account for the fraction of the
c grid cell within emitting surface.
c The scaling is done so that the charge density accurately represents
c the charge divided by the volume filled by that charge.  This only
c affects the edge of the emitting surface where only part of a grid cell
c will have charge in it, but the charge is divided by the area of all of
c the grid cell.

      integer(ISZ):: is,ip,ix,iy,ij,ij1
      real(kind=8):: xinj,yinj,zinj,rinj,rinji,zz,rr2,dd
      real(kind=8):: aa,az,px,py,wx1,wy1,wx0,wy0,ww
      real(kind=8):: dxi,dyi,dzi
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c     --- zero out the array
      call zeroarry(inj_rho(0,0,1),(1+nx)*(1+ny)*inj_ninj)
c     inj_rho = 0.

c     --- Loop over species and particles
      do is=1,ns
        ww = sq(is)*sw(is)*dxi*dyi*dzi
        do ip=ins(is),ins(is)+nps(is)-1
          ix = int(abs(xp(ip) - xmmin)*dxi)
          iy = int(abs(yp(ip) - ymmin)*dyi)

c         --- Get number of injection source of particle.  The particle id
c         --- has the number of the injection source the particle was
c         --- emitted from. The int is needed since it is doing double
c         --- duty and has other information stored in the fractional part.
          ij = int(pid(ip))

c         --- set temporaries
          rinji = 1./rinject(ij)
          rinj = rinject(ij)
          xinj = xinject(ij)
          yinj = yinject(ij)
          zinj = zinject

c         --- Calculate distance of particle from the emitting surface
          zz = zp(ip) - zinj
          rr2 = (xp(ip) - xinj)**2 + (yp(ip) - yinj)**2
          if (zz < rinj) then
c           --- The expression below for zz<rinj is mathematically identical
c           --- to the commented out expression below.  That form is
c           --- used so that for large values of rinj, i.e. a flat emitting
c           --- surface, the correct value of the distance is calculated,
c           --- namely zz.
c           dd = rinj - sqrt(rr2 + (rinj - zz)**2)
            dd = (2.*zz - (zz**2 + rr2)*rinji)/
     &           (1. + sqrt(rr2*rinji**2 + (1. - zz*rinji)**2))
          else
c           --- When zz>rinj, a seperate equation is needed.
            dd = rinj + sqrt(rr2 + (zz - rinj)**2)
          endif
          dd = dd*dzi

c         --- Only deposit charge of particles close to emitting surface.
          if (dd < 1. .and. uzp(ip) > 0.) then

c           --- angle of point in transverse plane
            aa = atan2(yp(ip)-yinj,dvnz(xp(ip)-xinj))

c           --- angle relative to z axis
            az = atan2(sqrt(rr2),rinj-zz)

c           --- Map the particle position onto the emitting surface.
            px = abs(rinj*sin(az)*cos(aa) + xinj - xmmin)*dxi
            py = abs(rinj*sin(az)*sin(aa) + yinj - ymmin)*dyi
            ix = int(px)
            iy = int(py)
            wx1 = px - ix
            wy1 = py - iy
            wx0 = 1. - wx1
            wy0 = 1. - wy1

c           --- Deposit the particle's charge onto the emitting surface.
            ij1 = min(ij,inj_ninj)
            inj_rho(ix  ,iy  ,ij1)=inj_rho(ix  ,iy  ,ij1) + ww*wx0*wy0*(1.-dd)
            inj_rho(ix+1,iy  ,ij1)=inj_rho(ix+1,iy  ,ij1) + ww*wx1*wy0*(1.-dd)
            inj_rho(ix  ,iy+1,ij1)=inj_rho(ix  ,iy+1,ij1) + ww*wx0*wy1*(1.-dd)
            inj_rho(ix+1,iy+1,ij1)=inj_rho(ix+1,iy+1,ij1) + ww*wx1*wy1*(1.-dd)

          endif
        enddo
      enddo

c     --- Loop over unique injection sources
      do ij1 = 1,min(ij,inj_ninj)

c     --- If two-fold symmetry, divide inj_rho by two except on x axis.
        if (l2symtry) then
          do iy=1,ny
            do ix=0,nx
              inj_rho(ix,iy,ij1) = inj_rho(ix,iy,ij1)*0.5
            enddo
          enddo
        endif

c     --- If four-fold symmetry, divide inj_rho by four, except on both axis
c     --- divide by two, except the origin.
        if (l4symtry) then
          do iy=1,ny
            do ix=1,nx
              inj_rho(ix,iy,ij1) = inj_rho(ix,iy,ij1)*0.25
            enddo
          enddo
          do ix=1,nx
            inj_rho(ix,0,ij1) = inj_rho(ix,0,ij1)*0.5
          enddo
          do iy=1,ny
            inj_rho(0,iy,ij1) = inj_rho(0,iy,ij1)*0.5
          enddo
        endif

c     --- Scale the charge density by the inverse of the fraction of the
c     --- contributing area which is within the emitting surface.
c     --- The factor of two is needed since particles are contributing to
c     --- rho only on the right side of the emitting surface.
        do iy=0,ny
          do ix=0,nx
            if (inj_area(ix,iy,ij1) > 0.)
     &        inj_rho(ix,iy,ij1) = 2.*inj_rho(ix,iy,ij1)/inj_area(ix,iy,ij1)
          enddo
        enddo

      enddo

      if (lw3dtimesubs) timeinj_setrho3d = timeinj_setrho3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine inj_setrho3d_rz(rho,xmmin,dx,dz,nx,nz)
      use Subtimers
      use InPart
      use InjectVars
      use InjectVars3d
      use Particles
      use Constant
      integer(ISZ):: nx,nz
      real(kind=8):: rho(0:nx,0:nz)
      real(kind=8):: xmmin,dx,dz

c Calculate the charge density on the emitting surface.  The particles are
c mapped onto the surface assuming a spherical emitter.  The charge density
c scale to account for symmetries and to account for the fraction of the
c grid cell within emitting surface.
c The scaling is done so that the charge density accurately represents
c the charge divided by the volume filled by that charge.  This only
c affects the edge of the emitting surface where only part of a grid cell
c will have charge in it, but the charge is divided by the area of all of
c the grid cell.

      integer(ISZ):: is,ip,ix,iy,ij,j,ij1
      real(kind=8):: xinj,yinj,zinj,rinj,rinji,zz,rr2,dd
      real(kind=8):: az,px,py,wx1,wx0,ww
      real(kind=8):: dxi,dyi,dzi
      real(kind=8):: substarttime,wtime
      real(kind=8):: invvol(0:nx)
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./dx
      dzi = 1./dz
      j = 0
      invvol(j) = 1._8 / (pi * (0.5_8*0.5_8*dx*dx)*dz)
      do j = 1, nx
        invvol(j) = 1._8 / (2._8 * pi * real(j,8) * dx * dx * dz)
      enddo

      iy = 0

c     --- zero out the array
      call zeroarry(inj_rho(0,0,1),(1+nx)*inj_ninj)
c      inj_rho = 0.

c     --- Loop over species and particles
      do is=1,ns
        ww = sq(is)*sw(is)*dxi*dzi
        do ip=ins(is),ins(is)+nps(is)-1
          ix = int(abs(sqrt(xp(ip)**2+yp(ip)**2) - xmmin)*dxi)

c         --- Get number of injection source of particle.  The particle id
c         --- has the number of the injection source the particle was
c         --- emitted from. The int is needed since it is doing double
c         --- duty and has other information stored in the fractional part.
          ij = int(pid(ip))

c         --- set temporaries
          rinji = 1./rinject(ij)
          rinj = rinject(ij)
          xinj = xinject(ij)
          yinj = yinject(ij)
          zinj = zinject

c         --- Calculate distance of particle from the emitting surface
          zz = zp(ip) - zinj
          rr2 = (xp(ip) - xinj)**2 + (yp(ip) - yinj)**2
          if (zz < rinj) then
c           --- The expression below for zz<rinj is mathematically identical
c           --- to the commented out expression below.  That form is
c           --- used so that for large values of rinj, i.e. a flat emitting
c           --- surface, the correct value of the distance is calculated,
c           --- namely zz.
c           dd = rinj - sqrt(rr2 + (rinj - zz)**2)
            dd = (2.*zz - (zz**2 + rr2)*rinji)/
     &           (1. + sqrt(rr2*rinji**2 + (1. - zz*rinji)**2))
          else
c           --- When zz>rinj, a seperate equation is needed.
            dd = rinj + sqrt(rr2 + (zz - rinj)**2)
          endif
          dd = dd*dzi

c         --- Only deposit charge of particles close to emitting surface.
          if (dd < 1. .and. uzp(ip) > 0.) then

c           --- angle relative to z axis
            az = atan2(sqrt(rr2),rinj-zz)

c           --- Map the particle position onto the emitting surface.
            px = abs(rinj*sin(az) + xinj - xmmin)*dxi
            ix = int(px)
            wx1 = px - ix
            wx0 = 1. - wx1

c           --- Deposit the particle's charge onto the emitting surface.
            ij1 = min(ij,inj_ninj)
            inj_rho(ix  ,iy,ij1)=inj_rho(ix  ,iy,ij1)+ww*wx0*(1.-dd)*invvol(ix)
            inj_rho(ix+1,iy,ij1)=inj_rho(ix+1,iy,ij1)+ww*wx1*(1.-dd)*invvol(ix+1)

          endif
        enddo
      enddo

c     --- Loop over unique injection sources
      do ij1 = 1,min(ij,inj_ninj)
c       --- Scale the charge density by the inverse of the fraction of the
c       --- contributing area which is within the emitting surface.
c       --- The factor of two is needed since particles are contributing to
c       --- rho only on the right side of the emitting surface.
        iy = 0
        do ix=0,nx
          if (inj_area(ix,iy,ij1) > 0.)
     &      inj_rho(ix,iy,ij1) = 2.*inj_rho(ix,iy,ij1)/inj_area(ix,iy,ij1)
        enddo

      enddo

      if (lw3dtimesubs) timeinj_setrho3d = timeinj_setrho3d + wtime() - substarttime
      return
      end
c===========================================================================
      subroutine sete3d_aperture(phi1d,np,xp,yp,zp,zgrid,xmmin,ymmin,zmmin,
     &                           zmmax,dx,dy,dz,nx,ny,nz,ex,ey,
     &                           l2symtry,l4symtry)
      use Subtimers
      use Apertures
      integer(ISZ):: np,nx,ny,nz
      real(kind=8):: zgrid,xmmin,ymmin,zmmin,zmmax,dx,dy,dz
      real(kind=8):: phi1d(0:*), xp(np), yp(np), zp(np)
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
                i =  abs(xp(ip) - xmmin)*dxi
                j =  abs(yp(ip) - ymmin)*dyi

                u1 = abs(xp(ip) - xmmin)*dxi - i
                v1 = abs(yp(ip) - ymmin)*dyi - j
                w1 = (zp(ip) - zgrid - zmmin)*dzi - k

                u0 = 1. - u1
                v0 = 1. - v1
                w0 = 1. - w1

c               --- Adjust sign of E field for approiate quadrant.
                xsign = 1.
                ysign = 1.
                if (xp(ip) .lt. 0.) xsign = sx
                if (yp(ip) .lt. 0.) ysign = sy

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

      if (lw3dtimesubs) timesete3d_aperture = timesete3d_aperture + wtime() - substarttime
      return
      end
c===========================================================================
      subroutine set_aperture_e()
      use Subtimers
      use InGen3d
      use Picglb
      use Picglb3d
      use InMesh3d
      use Fields3d
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
        izs = nz - int((zmmax - zs + zgrid) * dzi)
        ize = int((ze - zgrid - zmmin) * dzi)

c       --- Only do calculation if aperture is within the grid.
        if (0 < ize .and. izs < nz .and. izs <= ize) then

        midx = max(0,min(nx,int((aper_x(ia)                - xmmin)/dx + 0.5)))
        midy = max(0,min(ny,int((aper_y(ia)                - ymmin)/dy + 0.5)))
        minx = max(0,min(nx,int((aper_x(ia) - aper_rad(ia) - xmmin)/dx + 1  )))
        miny = max(0,min(ny,int((aper_y(ia) - aper_rad(ia) - ymmin)/dy + 1  )))
        maxx = max(0,min(nx,int((aper_x(ia) + aper_rad(ia) - xmmin)/dx      )))
        maxy = max(0,min(ny,int((aper_y(ia) + aper_rad(ia) - ymmin)/dy      )))
        minxm1 = max(0,minx-1)
        maxxp1 = min(nx,maxx+1)
        minym1 = max(0,miny-1)
        maxyp1 = min(ny,maxy+1)

c       -----------------------------------------------------------------------
c       --- Calculate Ex and Ey for all points using standard finite
c       --- difference.  For the planes before and after the aperture,
c       --- do the calculation over the full transverse extent of the grid.
c       --- Those points near the conductor will be recalculated below.
        do iz=izs-1,ize+1
          iza = iz - izs
          if (iz == izs-1 .or. iz == ize+1) then
            ixmin = 0
            ixmax = nx-1
            iymin = 0
            iymax = ny-1
          else
            ixmin = minxm1
            ixmax = maxxp1
            iymin = minym1
            iymax = maxyp1
          endif
          do iy=iymin,iymax
            do ix=ixmin,ixmax
              aper_ex(ix,iy,iza,ia)=tdxi*(phi(abs(ix-1),iy,iz)-phi(ix+1,iy,iz))
              aper_ey(ix,iy,iza,ia)=tdyi*(phi(ix,abs(iy-1),iz)-phi(ix,iy+1,iz))
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
                aper_ex(ix,iy,iza,ia) = tdxi*(phi(ix-1,iy,iz) - phi(ix,iy,iz) +
     &                                  (phi(ix  ,iy,iz) - aper_volt(ia))/wx)
              elseif (ix*dx == xx) then
                aper_ex(ix,iy,iza,ia) = tdxi*2.*(phi(ix-1,iy,iz)-aper_volt(ia))
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
                aper_ex(ix,iy,iza,ia) = tdxi*(phi(ix,iy,iz) - phi(ix+1,iy,iz) -
     &                                  (phi(ix,iy,iz) - aper_volt(ia))/wx)
              elseif (ix*dx == xx) then
                aper_ex(ix,iy,iza,ia) = tdxi*2.*(aper_volt(ia)-phi(ix+1,iy,iz))
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
                aper_ex(ix,iy,iza,ia) = tdxi*(phi(ix-1,iy,iz) - phi(ix,iy,iz) +
     &                                  (phi(ix  ,iy,iz) - aper_volt(ia))/wx)
              elseif (ix*dx == xx) then
                aper_ex(ix,iy,iza,ia) = tdxi*2.*(phi(ix-1,iy,iz)-aper_volt(ia))
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
                aper_ex(ix,iy,iza,ia) = tdxi*(phi(ix,iy,iz) - phi(ix+1,iy,iz) -
     &                                  (phi(ix,iy,iz) - aper_volt(ia))/wx)
              elseif (ix*dx == xx) then
                aper_ex(ix,iy,iza,ia) = tdxi*2.*(aper_volt(ia)-phi(ix+1,iy,iz))
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
                aper_ey(ix,iy,iza,ia) = tdyi*(phi(ix,iy-1,iz) - phi(ix,iy,iz) +
     &                                  (phi(ix,iy  ,iz) - aper_volt(ia))/wy)
              elseif (iy*dy == yy) then
                aper_ey(ix,iy,iza,ia) = tdyi*2.*(phi(ix,iy-1,iz)-aper_volt(ia))
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
                aper_ey(ix,iy,iza,ia) = tdyi*(phi(ix,iy-1,iz) - phi(ix,iy,iz) +
     &                                  (phi(ix,iy  ,iz) - aper_volt(ia))/wy)
              elseif (iy*dy == yy) then
                aper_ey(ix,iy,iza,ia) = tdyi*2.*(phi(ix,iy-1,iz)-aper_volt(ia))
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
                aper_ey(ix,iy,iza,ia) = tdyi*(phi(ix,iy,iz) - phi(ix,iy+1,iz) -
     &                                  (phi(ix,iy,iz) - aper_volt(ia))/wy)
              elseif (iy*dy == yy) then
                aper_ey(ix,iy,iza,ia) = tdyi*2.*(aper_volt(ia)-phi(ix,iy+1,iz))
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
                aper_ey(ix,iy,iza,ia) = tdyi*(phi(ix,iy,iz) - phi(ix,iy+1,iz) -
     &                                  (phi(ix,iy,iz) - aper_volt(ia))/wy)
              elseif (iy*dy == yy) then
                aper_ey(ix,iy,iza,ia) = tdyi*2.*(aper_volt(ia)-phi(ix,iy+1,iz))
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

      if (lw3dtimesubs) timeset_aperture_e = timeset_aperture_e + wtime() - substarttime
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
      call parallelmax(phiref,1)
#endif
      bndferr = bndferr / dvnz(phiref)
#ifdef MPIPARALLEL
      call parallelmax(bndferr,1)
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
      subroutine fieldsol3d(iwhich)
      use Subtimers
      use Timers
      use InGen
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Fields3d
      use Efields3d
      use LatticeInternal
      use InjectVars
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
      if (lw3dtimesubs) substarttime = wtime()
      timetemp = wtime()

c     --- If no field solve, return 
      if ( fstype == -1) return 

      if (fstype == 3 .or. fstype == 7) then
c       --- SOR and multigrid field solvers - deal with bends directly and set
c       --- axial boundary conditions (call to perphi3d is not needed)
        call vp3d (iwhich)

      else if (.not. linbend) then
c       --- If not in bend, call VP3D
        if (fstype == 0 .or. fstype == 1 .or. fstype == 2 .or.
     &      fstype == 4 .or. fstype == 5 .or. fstype == 6 .or.
     &      fstype == 8 .or. fstype == 9) then
          call copyarry (rho(0,0,0), phi(0,0,0), (nx+1)*(ny+1)*(nz+1))
        endif
        call vp3d (iwhich)
        if (periinz) call perphi3d(phi(0,0,-1),nx,ny,nz)

      else
c       --- Call perturbative bent beam field solver
        call bendfieldsol3d
        if (periinz) call perphi3d(phi(0,0,-1),nx,ny,nz)

      endif

#ifdef MPIPARALLEL
c     --- Distribute phi among the processes so each has phi in its
c     --- particle domain.
      call getphiforparticles(nx,ny,nz,phi(0,0,-1))
#endif

      if (lw3dtimesubs) timefieldsol3d = timefieldsol3d + wtime() - substarttime
      fstime = fstime + (wtime() - timetemp)
      return
      end



c=============================================================================
      subroutine vp3d (iwhich)
      use Subtimers
      use Constant
      use InGen
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Picglb
      use Fields3d
      use LatticeInternal
      use PSOR3d
      use FRZmgrid
      integer(ISZ):: iwhich
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

c Note that the select case statement was giving the f90 compiler on HPUX
c fits and so was replaced with if statements.

c     select case (fstype)
c       case (0)
      if (fstype == 0) then
          call vpois3d(iwhich,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq(0),kysq(0),kzsq(0),
     &                 attx(0),atty(0),attz(0),filt(1,1),xmmax-xmmin,ymmax-ymmin,
     &                 (izfsmax-izfsmin)*dz,nx,ny,izfsmax-izfsmin,
     &                 scrtch(0,0),xywork(0,0),zwork,0,l2symtry,l4symtry)
c       case (1)
      elseif (fstype == 1) then
          call pipe3df(iwhich,pipeshpe,rho(0,0,izfsmin),phi(0,0,izfsmin),
     &                kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xmmax-xmmin,ymmax-ymmin,(izfsmax-izfsmin)*dz,
     &                nx,ny,izfsmax-izfsmin,scrtch(0,0),xywork(0,0),
     &                zwork(1,0,0),l2symtry,l4symtry)
c       case (2)
      elseif (fstype == 2) then
          call vcap3d(iwhich,rho(0,0,izfsmin),phi(0,0,izfsmin),
     &                kxsq,kysq,kzsq,attx,atty,attz,filt,
     &                xmmax-xmmin,ymmax-ymmin,(izfsmax-izfsmin)*dz,
     &                nx,ny,izfsmax-izfsmin,scrtch(0,0),xywork(0,0),zwork,
     &                xmmax,zmmin+izfsmin*dz,zgrid,pipeshpe,periinz,
     &                l2symtry,l4symtry)
c       case (3)
      elseif (fstype == 3) then
c          --- Note that by default, izfsmax=nz, and izfsmin=0.
c          --- izfsmin and izfsmax allow field solves on a partial grid.
           call psor3df(iwhich,nx,ny,nz,phi(0,0,-1),rho(0,0,0),
     &                 phi(0,0,-1),rho(0,0,0),rstar, 
     &                 dx,dy,dz,xmmin,ymmin,zmmin,zbeam,zgrid,linbend,
     &                 l2symtry,l4symtry,scrtch(0,0),izfsmin,izfsmax)
c       case (4)
      elseif (fstype == 4) then
          if (iwhich == 1 .or. iwhich == 0) then
            call vpois3d(1,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq,kysq,kzsq,
     &                   attx,atty,attz,filt,xmmax-xmmin,ymmax-ymmin,
     &                   (izfsmax-izfsmin)*dz,nx,ny,izfsmax-izfsmin,
     &                   scrtch(0,0),xywork(0,0),zwork,0,l2symtry,l4symtry)
          endif
          if (iwhich == -1 .or. iwhich == 0) then
            call vpois3d(12,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq,kysq,kzsq,
     &                   attx,atty,attz,filt,xmmax-xmmin,ymmax-ymmin,
     &                   (izfsmax-izfsmin)*dz,nx,ny,izfsmax-izfsmin,
     &                   scrtch(0,0),xywork(0,0),zwork,izfsmin,l2symtry,l4symtry)
            call vpois3d(14,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq,kysq,kzsq,
     &                   attx,atty,attz,filt,xmmax-xmmin,ymmax-ymmin,
     &                   (izfsmax-izfsmin)*dz,nx,ny,izfsmax-izfsmin,
     &                   scrtch(0,0),xywork(0,0),zwork,0,l2symtry,l4symtry)
            call vpois3d(13,phi(0,0,izfsmin),phi(0,0,izfsmin),kxsq,kysq,kzsq,
     &                   attx,atty,attz,filt,xmmax-xmmin,ymmax-ymmin,
     &                   (izfsmax-izfsmin)*dz,nx,ny,izfsmax-izfsmin,
     &                   scrtch(0,0),xywork(0,0),zwork,0,l2symtry,l4symtry)
          endif

c       case (5)
      elseif (fstype == 5) then
c         --- General capacity matrix solver in kz space
          call capmatkz3d(iwhich,phi(0,0,0),rho,kxsq,kysq,kzsq,attx,atty,attz,
     &                    filt,xmmax-xmmin,ymmax-ymmin,(izfsmax-izfsmin)*dz,
     &                    nx,ny,izfsmax-izfsmin,dx,dy,dz,
     &                    xmmin,ymmin,zmmin+izfsmin*dz,scrtch(0,0),
     &                    xywork(0,0),zwork(1,0,0),l2symtry,l4symtry)
c       case (6)
      elseif (fstype == 6) then
c         --- General capacity matrix solver
          call capmat3df(iwhich,phi(0,0,0),rho,kxsq,kysq,kzsq,attx,atty,attz,
     &                  filt,xmmax-xmmin,ymmax-ymmin,(izfsmax-izfsmin)*dz,
     &                  nx,ny,izfsmax-izfsmin,dx,dy,dz,
     &                  xmmin,ymmin,zmmin+izfsmin*dz,scrtch(0,0),
     &                  xywork(0,0),zwork(1,0,0),l2symtry,l4symtry)

c       case (7)
      elseif (fstype == 7) then
          call multigrid3df(iwhich,nx,ny,nz,dx,dy,dz,phi(0,0,-1),rho(0,0,0),
     &                     rstar,linbend,l2symtry,l4symtry,scrtch(0,0),
     &                     xmmin,ymmin,zmmin,zbeam,zgrid)

#ifdef MPIPARALLEL
c       case (8)
      elseif (fstype == 8) then
c         --- Experimental parallel solver
c         --- Combines local tridiag solves with iteration to exchange
c         --- boundary data.
          call paralleltridiag(phi(0,0,izfsmin-1),kxsq,kysq,kzsq,
     &                         attx,atty,attz,filt,xmmax-xmmin,ymmax-ymmin,
     &                         (izfsmax-izfsmin)*dz,nx,ny,izfsmax-izfsmin,
     &                         scrtch(0,0),xywork(0,0),zwork,l2symtry,l4symtry)
c       case (9)
      elseif (fstype == 9) then
c         --- Experimental parallel solver based on the method
c         --- put forward by Lantz
          call lantzsolver(iwhich,phi(0,0,izfsmin-1),kxsq,kysq,kzsq,
     &                     attx,atty,attz,filt,xmmax-xmmin,ymmax-ymmin,
     $                     (izfsmax-izfsmin)*dz,nx,ny,izfsmax-izfsmin,
     &                     scrtch(0,0),xywork(0,0),zwork,l2symtry,l4symtry)

#endif
c       case (10)
      elseif (fstype == 10) then
c         --- RZ full-multigrid solver
          call multigridrzf(iwhich,phi(:,0,:),rho(:,0,:),nx,nz,dx,dz,
     &                      mgridrz_accuracy,
     &                      mgridrz_ncmax,
     &                      mgridrz_npre,
     &                      mgridrz_npost,
     &                      mgridrz_ncycles,
     &                      mgridrz_mgparam,
     &                      mgridrz_mgiters)
c       case (11)
      elseif (fstype == 11) then
c         --- Chombo AMR  full-multigrid solver
          call amrfieldsolve(iwhich)
c     endselect
      endif

      if (lw3dtimesubs) timevp3d = timevp3d + wtime() - substarttime
      return
      end
c=============================================================================
