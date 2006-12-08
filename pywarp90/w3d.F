#include "top.h"
c=============================================================================
c@(#) File W3D.F, version $Revision: 3.317 $, $Date: 2006/12/08 22:13:18 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package W3D of code WARP
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      module w3d_interfaces
      interface
      subroutine assignrhoandphi(rhopindts,phipindts)
      use Fields3dParticles
      real(kind=8),target:: rhopindts(0:nxp,0:nyp,0:nzp)
      real(kind=8),target:: phipindts(0:nxp,0:nyp,-1:nzp+1)
      end subroutine assignrhoandphi
      subroutine assignrhopandphip(rhopindts,phipindts)
      use Fields3dParticles
      real(kind=8),target:: rhopindts(0:nxp,0:nyp,0:nzp)
      real(kind=8),target:: phipindts(0:nxp,0:nyp,-1:nzp+1)
      end subroutine assignrhopandphip
      subroutine setrhoandphiforfieldsolve(rhopin,phipin)
      use Fields3dParticles
      real(kind=8),target:: rhopin(0:nxp,0:nyp,0:nzp)
      real(kind=8),target:: phipin(0:nxp,0:nyp,-1:nzp+1)
      end subroutine setrhoandphiforfieldsolve
      end interface
      end module w3d_interfaces
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
      use Subcycling,Only:nrhopndtscopies
      use Fields3dParticles
      use Efields3d
      use GridBoundary3d
      use Multipole
      use Io
      use Lattice
      use LatticeInternal
      use Particles, Only: pgroup,npmax,wpid,tpid,rpid,chdtspid,
     &                     xoldpid,yoldpid,zoldpid,lsaveoldpos
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
      use InDiag3d
      use ExtPart
      use w3d_interfaces
c     use ifcore

c  Invoked by the GENERATE command, it sets up the problem
c  This routine allots all of the neccesary dynamic arrays, calls the
c  particle loader and does the initial load onto the charge density
c  mesh, initializes arrays for the field solver and sets the mesh arrays,
c  does the initial field solve, and sets up other arrays that are needed.

      integer(ISZ):: i,j,k,ipmin,ip,iwin,nl
      integer(ISZ):: nextpid
c     integer(4):: fff
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
c     fff = for_set_fpe(14_4)

c  Announce that we're starting up

      call remark(" ***  particle simulation package W3D generating")

c  Put ny to zero for RZ geometry and setup a few flags

      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        ny = 0
        if (fstype .ne. -1 .and. fstype .ne. 12) fstype = 10
      elseif(solvergeom==XYgeom) then
        nz = 0
        if (fstype .ne. -1 .and. fstype .ne. 12) fstype = 10
      elseif(solvergeom==Zgeom) then
        nx = 0
        ny = 0
        if (fstype .ne. -1 .and. fstype .ne. 12) fstype = 10
      elseif(solvergeom==Rgeom) then
        ny = 0
        nz = 0
        if (fstype .ne. -1 .and. fstype .ne. 12) fstype = 10
      elseif(solvergeom==AMRgeom) then
        if (fstype .ne. -1 .and. fstype .ne. 12) fstype = 11
      endif

c  If a B field solver is being used, then turn on the fetching of the B field

      if (bfstype >= 0) lcallfetchb = .true.

c Set epflagpid if needed
      if (lepsaveonce) then
        epflagpid = nextpid()
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
      call setupSubcycling(pgroup)
      call setupSelfB(pgroup)
      call setupFields3dParticles()
      vbeamfrm = vbeam
      if(solvergeom==XYZgeom   .or. solvergeom==XYgeom) then
        if (l2symtry) then
          ymmin = 0.
        elseif (l4symtry) then
          xmmin = 0.
          ymmin = 0.
        endif
      else if(solvergeom==XZgeom) then
c        ymmin = 0.
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
      elseif(solvergeom==XYgeom) then
        dz = 1.
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
      do k = 0, nzfull
         zmesh(k) = k*dz + zmminglobal
      enddo
      do k = 0, nz
         zmeshlocal(k) = k*dz + zmmin
      enddo

c  Calculate location of axis in mesh, the term dx*1.e-5 acts as fuzp
      ix_axis = nint(-xmmin/dx)
      if(solvergeom==XYZgeom .or.solvergeom==AMRgeom .or. solvergeom==XYgeom) then
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
      if (fstype == 10) then
        if(solvergeom==RZgeom .or.
     &     solvergeom==XZgeom .or.
     &     solvergeom==Rgeom  .or.
     &     solvergeom==Zgeom) then
          call init_base(nx,nz,dx,dz,xmmin,zmmin)
        elseif(solvergeom==XYgeom) then
          call init_base(nx,ny,dx,dy,xmmin,ymmin)
        end if
      end if

c  Initialize the B field solver if it is being used

      call init_bfieldsolver(bfstype)

c  Initial call to fieldsolver in order to initialize attx, kxsq, etc.
      call fieldsol3d(1)

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
     &     solvergeom == XYgeom .or.
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

      call alotlostpart

c  Load the particles, calculate the charge density

      call remark(" ---  Loading particles")
      pgroup%ns = ns
      call setuppgroup(pgroup)
      call stptcl3d(pgroup)
      call setgamma(pgroup,lrelativ)
      call injctint(pgroup)
      call remark(" ---  Setting charge density")
      call loadrho3d(pgroup,-1,-1,-1,.true.)
      call loadj3d(pgroup,-1,-1,-1,.true.)
      if (chdtspid>0) then
        lsaveoldpos = .true.
        allspecl    = .true.
      endif
      if (lsaveoldpos) then
        xoldpid = nextpid()
        yoldpid = nextpid()
        zoldpid = nextpid()
      endif
      call remark(" ---  done")

c  Set up for injection

c  Create the dynamic arrays for the partcle qtys needed for the
c  damped mover.  They are always allocated, but with length 1 if not used.

      if (eomdamp /= 0.) then
        exeomoldpid = nextpid()
        eyeomoldpid = nextpid()
        ezeomoldpid = nextpid()
        exeomlagpid = nextpid()
        eyeomlagpid = nextpid()
        ezeomlagpid = nextpid()
      endif

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
      if (xplmin == 0.) xplmin = -xmmax
      if (xplmax == 0.) xplmax = xmmax
      if (yplmin == 0.) yplmin = -ymmax
      if (yplmax == 0.) yplmax = ymmax
      if (zplmin == 0.) zplmin = zmmin
      if (zplmax == 0.) zplmax = zmmax

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

      call prntpara(dx,dy,dz,lprntpara,pgroup)
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

c     if ( lfinishd ) then
c        call remark("w3dexe: problem completed.")
c$OMP MASTER
c        if (lw3dtimesubs) timew3dexe = timew3dexe + wtime() - substarttime
c$OMP END MASTER
c        return
c     elseif (nplive <= 0) then
c Check if there are any live particles left.
c If on first time step, set nplive to one and continue running since if using
c injection, there may not be any particles yet but still want to run.
c        if (it == 1) then
c          nplive = 1
c        else
c          print*," *** W3DEXE: stopping, nplive =",nplive
c$OMP MASTER
c          if (lw3dtimesubs) timew3dexe = timew3dexe + wtime() - substarttime
c$OMP END MASTER
c          return
c        endif
c     endif

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
      use Picglb
      use Picglb3d
      use Win_Moments
      use Z_Moments
      use Z_arrays
      use Hist
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
      use Particles, Only: pgroup,nplive
      use Picglb
      use Picglb3d
      use LatticeInternal
      use GridBoundary3d
      use Z_Moments,Only: tempmaxp,tempminp,tempzmmnts0,tempzmmnts
      use Timers
      character(*):: caller

c  When called by W3DEXE, step3d advances the system forward in time one 
c  timestep and gathers diagnostics.  When called by W3DGEN, step3d takes 
c  a step of zero size, to compute fields, and gather diagnostics at start 
c  of run.

      real(kind=8):: getbeamcom
      logical(ISZ):: thisstep,thiszbeam,dolabwn
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

c     --- This is a special routine needed when there is subcycling.
c     call setupevensubcyclingrho(it)

c  Main particle advance: x to t.l. it; v to t.l. it-1/2
c  Half-step in v from t.l. it-1   if last step was "special"
c  Full-step in v from t.l. it-3/2 if last step not "special"
c  No step at all if generating.

      if (caller == "w3dexe") then
        if (lspecial) then
          call padvnc3d ("halfv",pgroup)
        else
          call padvnc3d ("fullv",pgroup)
        endif

c       --- This is done here in case padvnc3d is called for multiple
c       --- pgroups
        if (lbeamcom) then
c         --- Set zbeam so that it follows the center of mass of the beam.
          zbeam = getbeamcom(pgroup) - zbeamcomoffset

        else
c         --- Advance beam frame location using the nominal beam frame velocity.
c          if (lgridqnt) then
c            zbeam = zbeam+vbeamfrm*dt
c          else
            zbeam = zgrid
c          endif
        endif

c       --- zgridprv needs to be updated for the "synchv" step
c       --- Note that zgridprv is also set at the beginning of w3dexe
        zgridprv = zgrid

        time = time + dt

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

c  Charge density contour plot diagnostics.  Note -- these diagnostics 
c  are done at this phase of the particle advance to allow for the eventual 
c  use of a single array for rho and phi.   

      if (lalways .or. lseldom) call pltfld3d("rho",ALWAYS)
      if (lseldom)            call pltfld3d("rho",SELDOM)

c  Set lattice; this is done just before field solve, and so is
c  relative to ZBEAM in the same way that self-fields are.

      call setlatt

c  Field-solve for potential.

      if (lbeforefs) call execuser("beforefs")
      call fieldsol3d(-1)
      call bfieldsol3d(-1)
      if (lafterfs) call execuser("afterfs")

c  Pre-calculate the self-E if it is needed for sete3d. This is done
c  after the call to afterfs in case some manipulation is done to phi.
      if (efetch == 3 .and.
     &    (solvergeom == XYZgeom .or.
     &     solvergeom == RZgeom .or.
     &     solvergeom == XZgeom .or.
     &     solvergeom == XYgeom .or.
     &     solvergeom == Rgeom  .or.
     &     solvergeom == Zgeom)) then
        if (maxval(pgroup%ndts) > 1) then
c         --- Note that with subcycling, this option will be incorrect since
c         --- the selfe for each ndts group is not (and will not) be
c         --- calculated.
          print*,"The efetch option 3 cannot be used with subcycling"
          stop
        endif
        call getselfe3d(phip(0,0,-1),nxp,nyp,nzp,selfe(1,0,0,0),
     &                  nx_selfe,ny_selfe,nz_selfe,dx,dy,dz,
     &                  boundxy,boundxy,boundxy,boundxy)
      endif

c  Set the potential near the emitting surface.
c     call getinj_phi()
c     call gettinj_phi()

c  Set the transverse E fields near any defined apertures.
      call set_aperture_e()

c  Complete constant current and axially directed space-charge limited
c  injection with new fields including injected particles.
      if (caller == "w3dexe") then
        call inject3d(2,pgroup)
      end if

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

      if (caller == "w3dexe" .and. lspecial) then
         call padvnc3d ("synchv",pgroup)
      elseif (caller == "w3dgen") then
         call padvnc3d ("gen",pgroup)
      endif

c  Finalize the moments calculation and do other diagnostics.
      call getzmmnt(1,0.,0.,0.,0.,0.,0.,0.,
     &              0.,0.,0.,0.,0.,3,nplive,0.,0.,0.,
     &              1,1,ns,tempmaxp,tempminp,tempzmmnts0,tempzmmnts)
      call getlabwn()
      call rhodia

c  Gather moments used in diagnostics at "special" timesteps only. 
c  Compute mean beam z velocity from current and line charge density 
c  on a 1-d mesh.  Also, calculate the electrostatic energy (getese), 
c  electrostatic potential on axis (sphiax), and the axial electric 
c  field (sezax).  

      if (lspecial) then
        if (lgetvzofz) call getvzofz
        call gtlchg
        call srhoax
        call getese
        call sphiax 
        call sezax
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
      call applyuniformfields(np,ex,ey,ez,bx,by,bz)

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

c     --- handle electric fields from 3-D grid
      call applyegrd(np,xp,yp,np,zp,.false.,ex,ey,ez)

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
      subroutine othere3d(np,xp,yp,zp,zbeam,zimax,zimin,straight,
     &                    ifeears,eears,eearsofz,dzzi,nzzarr,zzmin,
     &                    dedr,dexdx,deydy,dbdr,ex,ey,ez,bx,by,bz)
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
      subroutine seteears()
      use Subtimers3d
      use Beam_acc
      use InGen
      use InPart
      use InMesh3d
      use Picglb3d
      use Constant
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
     &              - (aion*amu)*vbeam**2*emitlong**2/((zion*echarge)*zm**4))

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
          eearsprs = - (aion*amu)*vbeam**2*emitlong**2/((zion*echarge)*zm**4)
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
      subroutine setuppadvncsubcyclingsampledaveraging(pgroup,it,is,
     &                                                 center,lcycle)
      use ParticleGroupmodule
      use Subcycling
      type(ParticleGroup):: pgroup
      integer(ISZ):: it,is
      character(*):: center
      logical(ISZ):: lcycle

c  Advance species according to subcycling rule governed by ndts.
c  The ndts groups are advanced as soon as all of the needed
c  information is available.
c
c  The halfv can be done on the first time step of each cycle.
c  The synchv is done on the last.
c
c  Note that (it-1) is used since it is advanced before this point.

      if (mod(it-1,pgroup%ndts(is-1)) == 0 .and.
     &    (center == "fullv" .or. center == "halfv")) then
        pgroup%ldts(is) = .true.
        ldts(ndtstorho(pgroup%ndts(is-1))) = pgroup%ldts(is)
        if (itndts(ndtstorho(pgroup%ndts(is-1))) == it-1) then
          itndts(ndtstorho(pgroup%ndts(is-1))) =
     &          itndts(ndtstorho(pgroup%ndts(is-1))) + pgroup%ndts(is-1)
        endif
      elseif ((mod(it-1,pgroup%ndts(is-1)) == pgroup%ndts(is-1) - 1 .and.
     &        center == "synchv") .or. center == "gen") then
        lcycle = .false.
      else
        pgroup%ldts(is) = .false.
        ldts(ndtstorho(pgroup%ndts(is-1))) = pgroup%ldts(is)
        lcycle = .true.
        if (center == "gen") lcycle = .false.
      endif
      
      return
      end
c=============================================================================
      subroutine setuppadvncsubcyclingfullvaveraging(pgroup,it,is,
     &                                               center,lcycle)
      use ParticleGroupmodule
      use Subcycling
      type(ParticleGroup):: pgroup
      integer(ISZ):: it,is
      character(*):: center
      logical(ISZ):: lcycle

c  Advance species according to subcycling rule governed by ndts.
c  If ndts is even, then the particles can be advanced when the faster
c  particles has gone half way through the larger time step. At that
c  point, there will be enough future charge densities from the faster
c  particles to form a time centered average.
c
c  If ndts is odd, then the particles are advanced on the time step that
c  is past the half way point, again so that all of the future charge
c  densities are known to form a time centered average.
c  The ndtstmp is just a clever way of finding that time step.
c
c  Note that (it-1) is used since it is advanced before this point.
c
c  The synchv is done one time step sooner since it is done after the
c  field solve which will include that time centered rho average.

      integer(ISZ):: ndtstmp

      ndtstmp = pgroup%ndts(is-1) - mod(pgroup%ndts(is-1),2)
      if (mod(it-1,pgroup%ndts(is-1)) == ndtstmp/2 .and.
     &        (center == "fullv" .or. center == "halfv")) then
        pgroup%ldts(is) = .true.
        ldts(ndtstorho(pgroup%ndts(is-1))) = pgroup%ldts(is)
        if (itndts(ndtstorho(pgroup%ndts(is-1))) < it) then
          itndts(ndtstorho(pgroup%ndts(is-1))) =
     &               itndts(ndtstorho(pgroup%ndts(is-1))) + pgroup%ndts(is-1)
        endif
      elseif ((pgroup%ndts(is-1) == 1 .or.
     &        (mod(it-1,pgroup%ndts(is-1)) == ndtstmp/2-1) .and.
     &         center == "synchv") .or. center == "gen") then
        lcycle = .false.
      else
        pgroup%ldts(is) = .false.
        ldts(ndtstorho(pgroup%ndts(is-1))) = pgroup%ldts(is)
        lcycle = .true.
        if (center == "gen") lcycle = .false.
      endif

      return
      end
c=============================================================================
      subroutine setuppadvncsubcyclinghalfvaveraging(pgroup,it,is,
     &                                               center,lcycle)
      use ParticleGroupmodule
      use Subcycling
      type(ParticleGroup):: pgroup
      integer(ISZ):: it,is
      character(*):: center
      logical(ISZ):: lcycle

      return
      end
c=============================================================================
      subroutine setuppadvncsubcyclingaveraging(pgroup,it,is,
     &                                          center,lcycle)
      use ParticleGroupmodule
      use Subcycling
      type(ParticleGroup):: pgroup
      integer(ISZ):: it,is
      character(*):: center
      logical(ISZ):: lcycle

      if (ndtsaveraging == 0) then
        call setuppadvncsubcyclingsampledaveraging(pgroup,it,is,center,lcycle)
      elseif (ndtsaveraging == 1) then
        call setuppadvncsubcyclingfullvaveraging(pgroup,it,is,center,lcycle)
      elseif (ndtsaveraging == 2) then
        call setuppadvncsubcyclinghalfvaveraging(pgroup,it,is,center,lcycle)
      endif

      return
      end
c=============================================================================
      subroutine positionadvance3d(pgroup,is,np,ipmin,interpdk,dt,
     &                             bendres,bendradi)
      use ParticleGroupmodule
      use Particles, Only: xoldpid,yoldpid,zoldpid
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,np,ipmin,interpdk
      real(kind=8):: dt
      real(kind=8):: bendres(np),bendradi(np)

      real(kind=8),pointer:: zpo(:)

      if (xoldpid > 0) pgroup%pid(ipmin:ipmin+np-1,xoldpid) = pgroup%xp(ipmin:ipmin+np-1)
      if (yoldpid > 0) pgroup%pid(ipmin:ipmin+np-1,yoldpid) = pgroup%yp(ipmin:ipmin+np-1)
      if (zoldpid > 0) pgroup%pid(ipmin:ipmin+np-1,zoldpid) = pgroup%zp(ipmin:ipmin+np-1)

      if (zoldpid > 0) then
        zpo => pgroup%pid(ipmin:ipmin+np-1,zoldpid)
      else
        allocate(zpo(np))
        zpo = pgroup%zp(ipmin:ipmin+np-1)
      endif

      if (interpdk .ne. 0) then
        call xpush3dintrp(pgroup,np,is,ipmin,dt)
      else
        call xpush3d(np,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &               pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),pgroup%gaminv(ipmin),dt)
      endif

c     --- Bend residence factor over POSITION advance step
      call getbend(np,np,pgroup%zp(ipmin),pgroup%uzp(ipmin),
     &             pgroup%gaminv(ipmin),
     &             bendres,bendradi,-dt,0.,.false.)

c     --- Correct position advance for warped mesh effect
      call zbendcor(pgroup,np,ipmin,dt,bendres,bendradi)

c     --- Correct position advance for slanted dipole entry/exit
      call sledgcor(pgroup,np,ipmin,zpo,0.,0.,0.,
     &              pgroup%sm(is),pgroup%sq(is),.false.)

      if (zoldpid == 0) then
        deallocate(zpo)
      endif

      return
      end
c=============================================================================
      subroutine padvnc3d(center,pgroup)
      use ParticleGroupmodule
      use GlobalVars
      use Subtimers3d
      use InMesh3d
      use InGen
      use InGen3d
      use InPart,Only: zimin,zimax
      use InGaps
      use InDiag
      use Lattice
      use LatticeInternal
      use Particles,Only: nplive,wpid,chdtspid,spid
      use Fields3d
      use Subcycling
      use Fields3dParticles
      use Efields3d
      use Picglb
      use Picglb3d
      use Beam_acc
      use Z_arrays
      use Z_Moments, only: nzmmnt,nszmmnt,tempmaxp,tempminp,tempzmmnts0,tempzmmnts
      use Damped_eom
      use Apertures
      use GridBoundary3d
      use DKInterp
      use DKInterptmp
      use FieldSolveAPI
      type(ParticleGroup):: pgroup
      character(*):: center

c  Advances the particles position and velocity according to CENTER,
c  and also loads RHO at the new time level.

c     --- Create local pointers to the arrays in pgroup.
      real(kind=8),pointer:: xp(:),yp(:),zp(:),uxp(:),uyp(:),uzp(:)
      real(kind=8),pointer:: ex(:),ey(:),ez(:),bx(:),by(:),bz(:)
      real(kind=8),pointer:: gaminv(:),pid(:,:)
      real(kind=8),pointer:: sm(:),sq(:),sw(:),dtscale(:)
      integer(ISZ),pointer:: ins(:),nps(:)

      integer(ISZ):: isid,is,ismax,ip,ipmin,i
      real(kind=8):: uxpadv,uypadv,uzpadv,halfdt
      real(kind=8):: fulldt_s,halfdt_s
      real(kind=8):: zgridprv_save
      real(kind=8),allocatable:: uxpo(:), uypo(:), uzpo(:)
      real(kind=8),allocatable:: exo(:),eyo(:),ezo(:),bxo(:),byo(:),bzo(:)
      real(kind=8),allocatable:: gaminvo(:)
      real(kind=8),allocatable:: bendres(:), bendradi(:)
#ifdef _OPENMP
      real(kind=8),allocatable:: threadmaxp(:,:,:),threadminp(:,:,:)
      real(kind=8),allocatable:: threadzmmnts0(:,:,:),threadzmmnts(:,:,:,:)
      integer(ISZ):: ithread,omp_get_thread_num
      integer(ISZ):: nthread,omp_get_num_threads
      integer(ISZ):: allocerror
#endif
      real(kind=8),pointer:: maxp(:,:),minp(:,:)
      real(kind=8),pointer:: zmmnts0(:,:),zmmnts(:,:,:)
      logical(ISZ):: lcycle
      real(kind=8):: substarttime,wtime

      if (lw3dtimesubs) substarttime = wtime()

      call setuppgroup(pgroup)

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
      if (pgroup%npid > 0) pid => pgroup%pid

      sm => pgroup%sm
      sq => pgroup%sq
      sw => pgroup%sw
      ins => pgroup%ins
      nps => pgroup%nps
      dtscale => pgroup%dtscale

      ismax = maxval(pgroup%sid)+1
      if (ismax == 0) return

      halfdt = 0.5*dt

!$OMP PARALLEL
!$OMP&PRIVATE(ip,xpo,ypo,zpo,uxpo,uypo,uzpo,bendres,bendradi,
!$OMP&        uxpadv,uypadv,uzpadv,maxp,minp,zmmnts0,zmmnts)

      allocate(uxpo(nparpgrp),uypo(nparpgrp),uzpo(nparpgrp))
      allocate(gaminvo(nparpgrp))
      allocate(bendres(nparpgrp), bendradi(nparpgrp))
      if (center == "synchv" .and. chdtspid>0) then
        allocate(exo(nparpgrp),eyo(nparpgrp),ezo(nparpgrp))
        allocate(bxo(nparpgrp),byo(nparpgrp),bzo(nparpgrp))
      endif

#ifdef _OPENMP
      allocate(maxp(6,0:nszmmnt),minp(6,0:nszmmnt),
     &         zmmnts0(NUMZMMNT,0:nszmmnt),
     &         zmmnts(0:nzmmnt,NUMZMMNT,0:nszmmnt),stat=allocerror)
      if (allocerror /= 0) then
        print*,"padvnc3d: allocation error ",allocerror,
     &         ": could not allocate temp arrays to shape ",nszmmnt
        stop
      endif
#endif

c  Zero the bend radius and residence arrays
c     call zeroarry (bendres,nparpgrp)
c     call zeroarry (bendradi,nparpgrp)
      bendres = 0.
      bendradi = 0.

c  Set lvdts array
c  This needs to be calculated ahead of time since the value of lvdts is needed
c  for speciess other that the one being looped over below.
      if (ndtsaveraging == 0) then
        do is=1,pgroup%ns
          if (pgroup%ndts(is-1) == 1) then
            pgroup%lvdts(is) = .true.
          elseif (mod(it-1,pgroup%ndts(is-1)) == pgroup%ndts(is-1) - 1 .and.
     &            center == "synchv") then
            pgroup%lvdts(is) = .true.
          else
            pgroup%lvdts(is) = .false.
          endif
        enddo
        pgroup%lvdts(pgroup%ns+1)=.false.
      endif
      
c  Loop over species
      do is=1,pgroup%ns
        isid = pgroup%sid(is-1) + 1
        lcycle = .false.
        call setuppadvncsubcyclingaveraging(pgroup,it,is,center,lcycle)
        if (lcycle) cycle

c Include the scaling factors in the time step size, both the one for
c subcycling and step size scaling for steady-state and slice modes.
        fulldt_s = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)
        halfdt_s = 0.5*fulldt_s

c If this species involves drift-kinetic interpolation, allocate
c  space for temporary arrays if sufficient space not already allocated
c Also set m_over_q, q_over_m**2
        if (interpdk(isid) .ne. 0) then
           if (npint < nparpgrp) then
              npint = nparpgrp
              call gchange("DKInterptmp",0)
           endif
           notusealphcalc(isid) = (1.-usealphacalc(isid))*alpha0(isid)
           if (sq(is) .ne. 0.) m_over_q(isid)=sm(is)/sq(is)
           if (sm(is) .ne. 0.) qovermsq(isid)= (sq(is)/sm(is))**2
        endif

c  Setup the moments if center is synchv or gen
c  If OPENMP, copies data to arrays which will be threadprivate
c  This needs to be fixed when there are multiple pgroups since the
c  species may be spread out among the pgroups.
        if (center == "synchv" .or. center == "gen") then
#ifdef _OPENMP
c         --- Same defaults as set in getzmmnt
          maxp(:,isid) = -LARGEPOS
          minp(:,isid) = +LARGEPOS
          zmmnts0(:,isid) = 0.
          zmmnts(:,:,isid) = 0.
#else
          maxp => tempmaxp
          minp => tempminp
          zmmnts0 => tempzmmnts0
          zmmnts => tempzmmnts
#endif
        endif
c  Loop over particle blocks; move each block separately
!$OMP DO
        do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
          ip = min(nparpgrp, ins(is)+nps(is)-ipmin)

c         --- Zero out B field arrays, but only if B is going
c         --- to be recomputed; B doesn't get recomputed if 
c         --- center is "synchv" and interpdk is nonzero and allspecl
c         --- is 0
          if (.not. (center .eq. "synchv" .and. interpdk(isid) .ne. 0 
     &        .and. .not. allspecl)) then
            bx(ipmin:ipmin+ip-1) = 0.
            by(ipmin:ipmin+ip-1) = 0.
            bz(ipmin:ipmin+ip-1) = 0.
          endif

c         --- Save old fields which are needed to check if particles
c         --- can have their ndts changed.
          if(center == "synchv" .and. chdtspid>0) then
            exo(1:ip) = ex(ipmin:ipmin+ip-1)
            eyo(1:ip) = ey(ipmin:ipmin+ip-1)
            ezo(1:ip) = ez(ipmin:ipmin+ip-1)
            bxo(1:ip) = bx(ipmin:ipmin+ip-1)
            byo(1:ip) = by(ipmin:ipmin+ip-1)
            bzo(1:ip) = bz(ipmin:ipmin+ip-1)
          endif

c         --- Obtain the self-fields
          call fetche3d(pgroup,ipmin,ip,is)
          call fetchb3d(pgroup,ipmin,ip,is)

c         --- Get E field for particles near the  injection surface.
          call inj_sete(pgroup,ipmin,ip,ex(ipmin),ey(ipmin),ez(ipmin))
          call inj_addtemp3d(pgroup,ip,ipmin,dz)

c         --- Get transverse E field for particles near any defined
c         --- apertures.
          call sete3d_aperture(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                         zgridprv,xmmin,ymmin,zmminp,zmmax,dx,dy,dz,
     &                         nxp,nyp,nzp,ex(ipmin),ey(ipmin),l2symtry,l4symtry)

c         --- Scale the self E-field to get the lowest order relativistic
c         --- correction.
          if (relativity == 1) then
            call sete3d_relativity(ip,ex(ipmin),ey(ipmin),vbeam)
          end if

c         --- Compute lag average for experimental damping algorithm
          if (eomdamp /= 0.) then
            call edamp(eomdamp,it,itdamp,center,ip,
     &                 ex(ipmin),ey(ipmin),ez(ipmin),
     &                 pid(ipmin,exeomoldpid),
     &                 pid(ipmin,eyeomoldpid),
     &                 pid(ipmin,ezeomoldpid),
     &                 pid(ipmin,exeomlagpid),
     &                 pid(ipmin,eyeomlagpid),
     &                 pid(ipmin,ezeomlagpid))
          end if
c         --- FULLV
c         Note if this is a species with interpolation between
c         drift kinetics and full PIC, and allspecl
c         is not 1, then ALWAYS do full-v steps.
c         That is, v is always interpreted as at a half step.
c         In this case fullv and halfv do the same thing and
c         synchv does nothing
          if (center == "fullv" .or. (center == "halfv" .and.
     &        interpdk(isid) == 1 .and. .not. allspecl)) then
c           --- Obtain bend radii and residence factors
            call getbend(ip,ip,pgroup%zp(ipmin),pgroup%uzp(ipmin),
     &                   pgroup%gaminv(ipmin),
     &                   bendres,bendradi,-halfdt_s,halfdt_s,.false.)
c           --- Correct Ez_self for warped mesh effect
            call bendez3d(ip,pgroup%xp(ipmin),pgroup%zp(ipmin),ez(ipmin),
     &                    bendres,bendradi,bends,bnezflag,linbend)
c           --- Add in Ez from axially-smoothed gaps 
            call gapfield(ip,pgroup%zp(ipmin),ez(ipmin),zbeam,zzmin,egap(0),dzz)

c           --- Add in ears and uniform focusing E field pieces
            call othere3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                    zbeam,zimax,zimin,straight,ifeears,eears,
     &                    eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,
     &                    ex(ipmin),ey(ipmin),ez(ipmin),
     &                    bx(ipmin),by(ipmin),bz(ipmin))
c           --- Set quad, dipole E and B; All: Bz
            call exteb3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                   gaminv(ipmin),-halfdt_s,halfdt_s,
     &                   bx(ipmin),by(ipmin),bz(ipmin),
     &                   ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),
     &                   bendres,bendradi,gammabar,fulldt_s)
c           --- Correction to z on entry/exit to accelerator gap
            call zgapcorr(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    -halfdt_s, halfdt_s, fulldt_s, sm(1), sq(1), time)
c           --- Apply cancellation of E+VxB correction
            if(pgroup%lebcancel) 
     &      call ebcancel3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin), 
     &                      ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin) )
c           --- First half-electric field increment to momenta
            call epush3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                   ex(ipmin),ey(ipmin),ez(ipmin),sq(is),sm(is),halfdt_s)
c           --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gamadv,lrelativ)
c           --- Magnetic field increment to momenta
c           --- If interpolating with drift kinetics, must do half
c           ---  a mu grad B correction before and half after Bpush
            if (interpdk(isid) .ne. 0) then
c               call setptrs(bx(ipmin),by(ipmin),bz(ipmin),
c    &                       ex(ipmin),ey(ipmin),ez(ipmin))
               call mugrdbpush(pgroup,ip,isid,ipmin,halfdt_s,fulldt_s,1)
            endif
            call bpush3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                   bx(ipmin),by(ipmin),bz(ipmin),sq(is),sm(is),
     &                   fulldt_s,ibpush)
            if (interpdk(isid) .ne. 0)
     &        call mugrdbpush(pgroup,ip,is,ipmin,halfdt_s,fulldt_s,0)
c           --- Final half-electric field increment to momenta
            call epush3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    ex(ipmin),ey(ipmin),ez(ipmin),sq(is),sm(is),halfdt_s)
c           --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gamadv,lrelativ)
c           --- Position advance
            call positionadvance3d(pgroup,is,ip,ipmin,interpdk(isid),
     &                             fulldt_s,bendres,bendradi)
c         --- HALFV
          elseif ((center == "halfv" .and. (interpdk(isid) .ne. 1 .or.
     &             allspecl)) .or. (center == "hbor2" .and.
     &             interpdk(isid) .ne. 0)) then
c           --- drift-kinetic-interpolated species are always treated
c           --- with full-v steps except when initalizing, which is
c           --- optionally taken care of by "hbor2" option
c           --- The hbor1 option does 1st half of a Boris step,
c           --- hbor2 does 2nd half.
c           --- Obtain bend radii and residence factors
            call getbend(ip,ip,pgroup%zp(ipmin),pgroup%uzp(ipmin),
     &                   pgroup%gaminv(ipmin),
     &                   bendres,bendradi,0.,halfdt_s,.false.)
c           --- Correct Ez_self for warped mesh effect
            call bendez3d(ip,pgroup%xp(ipmin),pgroup%zp(ipmin),ez(ipmin),
     &                    bendres,bendradi,bends,bnezflag,linbend)
c         --- Add in Ez from axially-smoothed gaps 
            call gapfield(ip,pgroup%zp(ipmin),ez(ipmin),zbeam,zzmin,egap(0),dzz)
c           --- Add in ears and uniform focusing E field pieces
            call othere3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                    zbeam,zimax,zimin,straight,ifeears,eears,
     &                    eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,
     &                    ex(ipmin),ey(ipmin),ez(ipmin),
     &                    bx(ipmin),by(ipmin),bz(ipmin))
c           --- Set quad, dipole E and B;  All: Bz
            call exteb3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                   gaminv(ipmin),0.,halfdt_s,
     &                   bx(ipmin),by(ipmin),bz(ipmin),
     &                   ex(ipmin),ey(ipmin),ez(ipmin),
     &                   sm(is),sq(is),bendres,bendradi,gammabar,fulldt_s)
c           --- Correction to z on entry/exit to accelerator gap
            call zgapcorr(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    0., halfdt_s, fulldt_s, sm(1), sq(1), time)
c           --- Apply cancellation of E+VxB correction
            if(pgroup%lebcancel) 
     &      call ebcancel3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin), 
     &                      ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin) )
c           --- Magnetic field increment to momenta
            if (interpdk(isid) .ne. 0) then
              call bpush3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                     bx(ipmin),by(ipmin),bz(ipmin),
     &                     sq(is),sm(is),fulldt_s,-1)
c              --- If interpolating with drift kinetics, must do 
c              ---  a mu grad B correction.  Since may be calling
c              ---  this on first timestep, must make sure grad B's are
c              ---  calculated, so last arg = 1.
               call mugrdbpush(pgroup,ip,is,ipmin,halfdt_s,fulldt_s,1)
            else
               call bpush3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                      bx(ipmin),by(ipmin),bz(ipmin),
     &                      sq(is),sm(is),halfdt_s,ibpush)
c           --- Final half-electric field increment to momenta
            endif
            call epush3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                   ex(ipmin),ey(ipmin),ez(ipmin),sq(is),sm(is),halfdt_s)
c           --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gamadv,lrelativ)
c           --- Position advance
            if (center .ne. "hbor2") then
c             --- hbor2 just does a v advance, no x push
              call positionadvance3d(pgroup,is,ip,ipmin,interpdk(isid),
     &                               fulldt_s,bendres,bendradi)
            endif

c         --- SYNCHV or GEN
c         Double "if" follows so that moments are calculated
c         for center = "synchv" or "gen" regardless of interpdk,
c         but the vpush is only done for interpdk=1 if allspecl=1.
          elseif (center == "synchv" .or. center == "gen"
     &            .or. center == "hbor1") then

c           --- Copy 'old' velocity into uxpo, uypo, and uzpo
c           --- This is needed for the call to getzmmnt below.
c           --- gaminvo is saved so it can be restored below
c           --- during a gen step.
            if (center == "gen" .or. center == "synchv") then
              do i=1,ip
                uxpo(i) = uxp(ipmin+i-1)
                uypo(i) = uyp(ipmin+i-1)
                uzpo(i) = uzp(ipmin+i-1)
                gaminvo(i) = gaminv(ipmin+i-1)
              enddo
            endif
            if ((center == "synchv" .and. (interpdk(isid) == 0
     &          .or. allspecl)) 
     &          .or. center == "gen" .or. (center == "hbor1"
     &          .and. interpdk(isid) .ne. 0)) then

c             --- drift-kinetic-interpolated species are always treated
c             --- with full-v steps except in initial setup, when we want
c             --- a half-v step packwards.  center = "hbor2" provides
c             --- the latter capability with dt set to be negative.
c             --- Exception, if allspecl = 1, then we synchronize every
c             --- step, using half Boris if interpdk = 1.
c             --- More generally hbor1 provides 1st half of a Boris
c             --- step, hbor2 provides 2nd half.

c             --- Obtain bend radii and residence factors
              call getbend(ip,ip,pgroup%zp(ipmin),pgroup%uzp(ipmin),
     &                     pgroup%gaminv(ipmin),bendres,
     &                     bendradi,-halfdt_s,0.,.false.)
c             --- Correct Ez_self for warped mesh effect
              call bendez3d(ip,pgroup%xp(ipmin),pgroup%zp(ipmin),ez(ipmin),
     &                     bendres,bendradi,bends,bnezflag,linbend)
c             --- Add in Ez from axially-smoothed gaps 
              call gapfield(ip,pgroup%zp(ipmin),ez(ipmin),zbeam,zzmin,egap(0),dzz)
c             --- Add in ears and uniform focusing E field pieces
              call othere3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                      zbeam,zimax,zimin,straight,ifeears,eears,
     &                      eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,
     &                      ex(ipmin),ey(ipmin),ez(ipmin),
     &                      bx(ipmin),by(ipmin),bz(ipmin))
c             --- Set quad, dipole E and B; All: Bz
              call exteb3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                     gaminv(ipmin),-halfdt_s,0.,
     &                     bx(ipmin),by(ipmin),bz(ipmin),
     &                     ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),
     &                     bendres,bendradi,gammabar,fulldt_s)
c             --- Apply cancellation of E+VxB correction
              if(pgroup%lebcancel) 
     &        call ebcancel3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin), 
     &                        ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin) )
c             --- Half electric field increment to momenta
              call epush3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                     ex(ipmin),ey(ipmin),ez(ipmin),sq(is),sm(is),halfdt_s)
c             --- Advance relativistic Gamma factor
              call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                      gamadv,lrelativ)
c             --- If interpolating with drift kinetics, must do half
c             ---  a mu grad B correction
c             --- Half magnetic field increment to momenta
              if (interpdk(isid) .ne. 0) then
                 call mugrdbpush(pgroup,ip,is,ipmin,halfdt_s,fulldt_s,1)
                 call bpush3d(ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                        bx(ipmin),by(ipmin),bz(ipmin),sq(is),sm(is),fulldt_s,-1)
c                Half a Boris push
              else
                 call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                      bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is), halfdt_s, ibpush)
              endif
              if (center == "gen") then
c               --- Reset uxp to uxpo, set uxpo to half step backward
c               --- for interpolation in moments calculation
c               --- Don't do this if we're doing drift-kinetic interp
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
                  gaminv(ipmin+i-1) = gaminvo(i)
                enddo
              endif
            endif
c           --- Calculate moments over particles, now that we're sync'd
c           --- THIS IS NOT DONE FOR center = hbor1.
            if (center == "gen" .or. center == "synchv") then
              if(wpid==0) then
               call getzmmnt(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                       uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),sq(is),sm(is),sw(is),
     &                       halfdt_s,dtscale(is),2,nplive,uxpo,uypo,uzpo,
     &                       is,isid,ismax,
     &                       maxp,minp,zmmnts0,zmmnts)
              else
               call getzmmnt_weights(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                       uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),pid(ipmin,wpid),
     &                       sq(is),sm(is),sw(is),halfdt_s,dtscale(is),
     &                       2,nplive,uxpo,uypo,uzpo,
     &                       is,isid,ismax,
     &                       maxp,minp,zmmnts0,zmmnts)
              endif
              call getextrapolatedparticles(ipmin,ip,pgroup,pgroup%npid,
     &                                      halfdt_s,uxpo,uypo,uzpo,
     &                                      isid,time,zbeam)
            endif
            if(center == "synchv" .and. chdtspid > 0) then
              call check_cc3d(pgroup,is,ipmin,ip,exo,eyo,ezo,bxo,byo,bzo)
            endif
          endif

c       --- End of loop over particle blocks for cache
        enddo
!$OMP END DO

c     --- End of loop over speices.
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

      deallocate(uxpo,uypo,uzpo)
      deallocate(gaminvo)
      deallocate(bendres, bendradi)

!$OMP END PARALLEL

c----------------------------------------------------------------------------
c     ---  Do final stuff for moments calculation
      if (center == "synchv" .or. center == "gen") then
#ifdef _OPENMP
        tempmaxp = max(tempmaxp,maxval(threadmaxp(:,1:nthread),3))
        tempminp = min(tempminp,minval(threadminp(:,1:nthread),3))
        tempzmmnts0 = tempzmmnts0 + sum(threadzmmnts0(:,1:nthread),3)
        tempzmmnts = tempzmmnts + sum(threadzmmnts(:,:,1:nthread),4)
        deallocate(threadmaxp,threadminp,threadzmmnts0,threadzmmnts)
#endif
      endif

      if (center == "synchv" .or. center == "gen" .or. ifgap) then

c       --- Calculate current
        call setcurr(pgroup,zbeam,ismax,wpid,lspeciesmoments,
     &               .not. laccumulate_rho,bound0)

      endif
c----------------------------------------------------------------------------
      if (center == "fullv" .or. center == "halfv") then

c       --- Inject more particles. This adds particles and alters ins and nps.
        call inject3d(1,pgroup)

c       --- Treat particles at boundaries
        call particleboundaries3d(pgroup)

c       --- Collect charge density and current
        call loadrho3d(pgroup,-1,-1,-1,.not. laccumulate_rho)
        call loadj3d(pgroup,-1,-1,-1,.not. laccumulate_rho)

      elseif (center == 'synchv') then

c       --- Rearrange particles group/species according to Courant condition.
c       --- Note that the particle's old ndts values are saved during
c       --- the rearrangement for the following synchv.
        if(chdtspid > 0) then
          call chgparticlesdts(pgroup)

c         --- Recollect charge density and current
c         --- Temporarily increment 'it' as if this was the start of the next
c         --- fullv or halfv step. Also set zgridprv to zgrid since the
c         --- particle positions are at zgrid and loadrho sets
c         --- zgrdndts = zgrid + (zgrid-zgridprv)*(pgroup%ndts(is-1)-1).
          it = it + 1
          zgridprv_save = zgridprv
          zgridprv = zgrid
          do is=1,pgroup%ns
            call setuppadvncsubcyclingaveraging(pgroup,it,is,'halfv',lcycle)
          enddo
          call loadrho3d(pgroup,-1,-1,-1,.not. laccumulate_rho)
c         call loadj3d(pgroup,-1,-1,-1,.not. laccumulate_rho)
          zgridprv = zgridprv_save
          it = it - 1
          do is=1,pgroup%ns
            call setuppadvncsubcyclingaveraging(pgroup,it,is,center,lcycle)
          enddo

        endif

      endif
c----------------------------------------------------------------------------

!$OMP MASTER
      if (lw3dtimesubs) timepadvnc3d = timepadvnc3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine getextpart(pgroup)
      use ParticleGroupmodule
      use GlobalVars
      use InGen
      use InPart
      use Picglb
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,ip,ipmin
      real(kind=8):: fulldt_s,halfdt_s

      do is=1,pgroup%ns

c       --- Get the extrpolated particles for the species.
c       --- This is only done after the positions have been advanced,
c       --- in which case ldts will be true.
c       --- Include time step size scaling for steady-state and slice modes.
        if (pgroup%ldts(is)) then
          fulldt_s = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)
          halfdt_s = 0.5*fulldt_s

          ipmin = pgroup%ins(is)
          call getextrapolatedparticles(pgroup%ins(is),pgroup%nps(is),pgroup,
     &               pgroup%npid,halfdt_s,
     &               pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     &               pgroup%sid(is-1)+1,time,zbeam)
         endif
       enddo

      return
      end
c=============================================================================
      subroutine perphi3d()
      use Subtimers3d
      use InGen3d
      use InMesh3d
      use Fields3d

c  Sets the slices on the exterior of phi for periodicity
c  sets slice at -1 equal to the slice at nz-1
c  sets slice at nz+1 equal to the slice at 1

      integer(ISZ):: ix,iy
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        call perphirz()
        return
      end if

#ifdef MPIPARALLEL
      call perpot3d_slave(phi,1,nx,ny,nz,0,0)
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
      subroutine perrho3d()
      use GlobalVars
      use Subtimers3d
      use InGen,Only: fstype
      use InGen3d
      use InMesh3d
      use Fields3d,Only: rho
      use GridBoundary3d

c  Sums the first and last slices of rho for periodicity
c  and puts the result into both slices.

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (fstype == 12) return
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
        call persource3d_slave(rho,1,nx,ny,nz)
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
      subroutine xpush3d(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dtp)
      use Subtimers3d
      integer(ISZ):: np
      real(kind=8):: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np),gaminv(np)
      real(kind=8):: dtp

c  Advance particle positions

      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      do ip=1,np
        xp(ip) = xp(ip) + uxp(ip)*gaminv(ip)*dtp
        yp(ip) = yp(ip) + uyp(ip)*gaminv(ip)*dtp
        zp(ip) = zp(ip) + uzp(ip)*gaminv(ip)*dtp
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timexpusht3d = timexpusht3d + wtime() - substarttime
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
c=============================================================================
      subroutine ebcancel3d(np,uxp,uyp,uzp,gaminv,exp,eyp,ezp,bxp,byp,bzp)
      use Subtimers3d
      integer(ISZ):: np
      real(kind=8):: uxp(np),uyp(np),uzp(np),gaminv(np)
      real(kind=8):: exp(np),eyp(np),ezp(np),bxp(np),byp(np),bzp(np)
      real(kind=8):: v,vx,vy,vz,e,ex,ey,ez,ve,vex,vey,vez
      real(kind=8):: eve,evex,evey,evez,b_e,b_ve,b_eve,v_eve,vcrossb_e,de,fact

c Apply E+VxB cancellation.
c Vectors are first decomposed on base B={E,VxE,Ex(VxE)}. 

      integer(ISZ):: ip
      real(kind=8):: const
      real(kind=8):: substarttime,wtime
!      if (lw3dtimesubs) substarttime = wtime()
      do ip=1,np
        ! get unit vectors e and v
        e  = sqrt(exp(ip)*exp(ip)+eyp(ip)*eyp(ip)+ezp(ip)*ezp(ip))
        if(e==0.) cycle
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)
        v = sqrt(vx*vx+vy*vy+vz*vz)
        if(v==0.) return
        ex = exp(ip)/e; ey = eyp(ip)/e; ez = ezp(ip)/e
        vx = vx/v; vy = vy/v; vz = vz/v
        ! get unit vectors ve
        vex = vy*ez-vz*ey
        vey = vz*ex-vx*ez
        vez = vx*ey-vy*ex
        ve = sqrt(vex*vex+vey*vey+vez*vez)
        vex = vex/ve; vey = vey/ve; vez = vez/ve
        if(ve==0.) cycle
        ! get unit vectors eve
        evex = ey*vez-ez*vey
        evey = ez*vex-ex*vez
        evez = ex*vey-ey*vex
        eve = sqrt(evex*evex+evey*evey+evez*evez)
        evex = evex/eve; evey = evey/eve; evez = evez/eve
        ! get B field in base {E,VxE,Ex(VxE)}
        b_e   = bxp(ip)*ex  +byp(ip)*ey  +bzp(ip)*ez
        b_ve  = bxp(ip)*vex +byp(ip)*vey +bzp(ip)*vez
        if(b_ve==0.) cycle
        b_eve = bxp(ip)*evex+byp(ip)*evey+bzp(ip)*evez
        ! get velocity V in base {E,VxE,Ex(VxE)}
!        v_e   = vx*ex  +vy*ey  +vz*ez  ! uninteresting
!        v_ve  = vx*vex +vy*vey +vz*vez ! should be 0
        v_eve = v*(vx*evex+vy*evey+vz*evez)
        if(v_eve==0.) cycle
        ! get VxB along E
        vcrossb_e = -v_eve*b_ve!+v_ve*b_eve
        if(vcrossb_e>=0.) cycle
        de = max(vcrossb_e,-e)
        b_ve = b_ve+de/v_eve
        fact = de/e
        exp(ip) = exp(ip)*(1.+fact)
        eyp(ip) = eyp(ip)*(1.+fact)
        ezp(ip) = ezp(ip)*(1.+fact)
        bxp(ip) = b_e*ex+b_ve*vex+b_eve*evex
        byp(ip) = b_e*ey+b_ve*vey+b_eve*evey
        bzp(ip) = b_e*ez+b_ve*vez+b_eve*evez
      enddo

!$OMP MASTER
!      if (lw3dtimesubs) timeepush3d = timeepush3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine check_cc3d(pgroup,is,ipmin,np,exo,eyo,ezo,bxo,byo,bzo)
      use Constant
      use InGen
      use Subtimers3d
      use Beam_acc
      use Picglb, only: time,zgrid
      use Picglb3d
      use InGen3d
      use InMesh3d
      use ParticleGroupmodule
      use Particles, toppgroup => pgroup
      use Fields3dParticles, only: nxc,nyc,nzc,isnearbycond
      implicit none
      type(ParticleGroup):: pgroup
      integer(ISZ):: np,is,ipmin
      real(kind=8):: exo(np),eyo(np),ezo(np),bxo(np),byo(np),bzo(np)

c Check if particle is under or over the Courant condition limits (ccmin and ccmax) 
c and register the ones that fall outside the limits for change to group (or species) 
c with appropriate time step.

      integer(ISZ):: ipmin0,ip0,ips,myndts
      real(kind=8):: tdx,tdy,tdz
      real(kind=8):: dd,denorm,dbnorm,denorm_up,dbnorm_up,dex,dey,dez,dbx,dby,dbz
      integer(ISZ),allocatable::pisnearbycond(:)
      logical(ISZ):: ldfield_below_threshold,ldfield_below_threshold_up,l_pisnearbycond
      integer(ISZ):: ip,ndtsmax
      real(kind=8):: substarttime,wtime,ddx,ddy,ddz,dtgaminv
      if (lw3dtimesubs) substarttime = wtime()

      tdx = dx
      tdy = dy
      tdz = dz
      ndtsmax = 2**(pgroup%lvdtsmax-1)
      myndts  = pgroup%ndts(is-1)
      ipmin0=ipmin-1
      if (nxc>0) then
        allocate(pisnearbycond(np))
        call getgridngp3di(np,pgroup%xp(ipmin:ipmin+np-1),pgroup%yp(ipmin:ipmin+np-1),pgroup%zp(ipmin:ipmin+np-1),
     &                     pisnearbycond,nxc,nyc,nzc,isnearbycond,xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,zgrid,
     &                     l2symtry,l4symtry)
      else
        l_pisnearbycond=.false.
      end if
      do ip=1,np
        ip0=ipmin0+ip
        if(tpid>0) then
          if ((time-pgroup%pid(ip0,tpid))<ndtsmax*dt) cycle
        endif
        if (dxpid > 0 .and. dypid > 0 .and. dzpid > 0) then
          tdx = pgroup%pid(ip0,dxpid)
          tdy = pgroup%pid(ip0,dypid)
          tdz = pgroup%pid(ip0,dzpid)
        endif
        ddx=abs(pgroup%xp(ip0)-pgroup%pid(ip0,xoldpid))/tdx
        ddy=abs(pgroup%yp(ip0)-pgroup%pid(ip0,yoldpid))/tdy
        ddz=abs(pgroup%zp(ip0)-pgroup%pid(ip0,zoldpid))/tdz
        dex=abs(pgroup%ex(ip0)-exo(ip))
        dey=abs(pgroup%ey(ip0)-eyo(ip))
        dez=abs(pgroup%ez(ip0)-ezo(ip))
        dbx=abs(pgroup%bx(ip0)-bxo(ip))
        dby=abs(pgroup%by(ip0)-byo(ip))
        dbz=abs(pgroup%bz(ip0)-bzo(ip))
        denorm = defieldmax*sqrt(exo(ip)*exo(ip)+eyo(ip)*eyo(ip)+ezo(ip)*ezo(ip))
        dbnorm = dbfieldmax*sqrt(bxo(ip)*bxo(ip)+byo(ip)*byo(ip)+bzo(ip)*bzo(ip))
        ldfield_below_threshold = dex<=denorm .and. dey<=denorm .and. dez<=denorm .and. 
     &                            dbx<=dbnorm .and. dby<=dbnorm .and. dbz<=dbnorm
        if (myndts<ndtsmax) then
c         --- the threshold on field change is lower for increasing time step in order to 
c         --- minimize particles changing time steps back and forth due to constant crossing of threshold
          denorm_up = denorm*real(pgroup%ndts(is-1))/pgroup%ndts(is+1-1)
          dbnorm_up = dbnorm*real(pgroup%ndts(is-1))/pgroup%ndts(is+1-1)
          ldfield_below_threshold_up = dex<=denorm_up .and. dey<=denorm_up .and. dez<=denorm_up .and. 
     &                                 dbx<=dbnorm_up .and. dby<=dbnorm_up .and. dbz<=dbnorm_up
        endif
        if (nxc>0) then
          if (pisnearbycond(ip)>0) then
            l_pisnearbycond=.true.
          else
            l_pisnearbycond=.false.
          endif
        end if
        dd = max(ddx,ddy,ddz)
c       --- if dd<ccmin or dd>ccmax, register particle for move to group/species with 
c       --- appropriate time step 
        if(myndts>1 .and. (dd>courantmax .or. (.not. ldfield_below_threshold) .or. l_pisnearbycond)) then
c         --- flag to move to group with smaller time step
          pgroup%pid(ip0,chdtspid) = -1.
        elseif(myndts<ndtsmax .and. pgroup%lvdts(is+1) .and. dd<courantmin .and. ldfield_below_threshold_up .and. 
     &         .not. l_pisnearbycond) then
c         --- flag to move to group with larger time step
          pgroup%pid(ip0,chdtspid) = +1.
        endif
      enddo

      if (nxc>0) deallocate(pisnearbycond)

!$OMP MASTER
      if (lw3dtimesubs) timecheck_cc3d = timecheck_cc3d + wtime() - substarttime
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
      subroutine sete3d(phi1d,selfe,np,xp,yp,zp,zgrid,xmmin,ymmin,zmmin,
     &                  dx,dy,dz,nx,ny,nz,efetch,ex,ey,ez,l2symtry,l4symtry,
     &                  lcylindrical)
      use Subtimers3d
      integer(ISZ):: np,nx,ny,nz
      real(kind=8):: zgrid,xmmin,ymmin,zmmin,dx,dy,dz
      real(kind=8):: phi1d(0:*),selfe(3,0:nx,0:ny,0:nz),xp(np),yp(np),zp(np)
      real(kind=8):: ex(np),ey(np),ez(np)
      integer(ISZ):: efetch
      logical(ISZ):: l2symtry,l4symtry,lcylindrical

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

      integer(ISZ):: nnx,nnxy,ip,i,j,k,ind0,m,inext,jnext,knext
      real(kind=8):: dxi,dyi,dzi,tdxi,tdyi,tdzi,u0,u1,v0,v1,w0,w1,ysign,xsign
      integer(ISZ):: ox,oy,sox,soy
      real(kind=8):: sx,sy
      real(kind=8):: xi,yj,zk,x,y
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

           if (lcylindrical) then
             x = sqrt(xp(ip)**2 + yp(ip)**2)
             y = 0.
           else
             x = xp(ip)
             y = yp(ip)
           endif

           i =  (x      - xmmin) * dxi
           j =  (y      - ymmin) * dyi
           k =  (zp(ip) - zgrid - zmmin) * dzi
           ind0 = i + j*nnx + k*nnxy

           u1 = (x      - xmmin) * dxi - i
           v1 = (y      - ymmin) * dyi - j
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

      else if (efetch == 5) then
c       --- This is the same as efetch==1, but has checks to skip particles that
c       --- are out of bounds.

        if (.not. (l2symtry .or. l4symtry)) then
          do ip = 1, np

            xi =  (xp(ip) - xmmin) * dxi
            yj =  (yp(ip) - ymmin) * dyi
            zk =  (zp(ip) - zgrid - zmmin) * dzi
            if (xi < 0. .or. xi > nx .or.
     &          yj < 0. .or. yj > ny .or.
     &          zk < 0. .or. zk > nz) then
              ex(ip) = 0.
              ey(ip) = 0.
              ez(ip) = 0.
              cycle
            endif
            i = int(xi)
            j = int(yj)
            k = int(zk)
            ind0 = i + j*nnx + k*nnxy

            u1 = xi - i
            v1 = yj - j
            w1 = zk - k

            u0 = 1. - u1
            v0 = 1. - v1
            w0 = 1. - w1

            ex(ip)=tdxi*(u0*v0*w0*(phi1d(noff( 7)+ind0) - phi1d(noff( 9)+ind0))
     &                 + u1*v0*w0*(phi1d(noff( 8)+ind0) - phi1d(noff(10)+ind0))
     &                 + u0*v1*w0*(phi1d(noff(11)+ind0) - phi1d(noff(13)+ind0))
     &                 + u1*v1*w0*(phi1d(noff(12)+ind0) - phi1d(noff(14)+ind0))
     &                 + u0*v0*w1*(phi1d(noff(19)+ind0) - phi1d(noff(21)+ind0))
     &                 + u1*v0*w1*(phi1d(noff(20)+ind0) - phi1d(noff(22)+ind0))
     &                 + u0*v1*w1*(phi1d(noff(23)+ind0) - phi1d(noff(25)+ind0))
     &                 + u1*v1*w1*(phi1d(noff(24)+ind0) - phi1d(noff(26)+ind0)))

            ey(ip)=tdyi*(u0*v0*w0*(phi1d(noff( 5)+ind0) - phi1d(noff(12)+ind0))
     &                 + u1*v0*w0*(phi1d(noff( 6)+ind0) - phi1d(noff(13)+ind0))
     &                 + u0*v1*w0*(phi1d(noff( 8)+ind0) - phi1d(noff(15)+ind0))
     &                 + u1*v1*w0*(phi1d(noff( 9)+ind0) - phi1d(noff(16)+ind0))
     &                 + u0*v0*w1*(phi1d(noff(17)+ind0) - phi1d(noff(24)+ind0))
     &                 + u1*v0*w1*(phi1d(noff(18)+ind0) - phi1d(noff(25)+ind0))
     &                 + u0*v1*w1*(phi1d(noff(20)+ind0) - phi1d(noff(27)+ind0))
     &                 + u1*v1*w1*(phi1d(noff(21)+ind0) - phi1d(noff(28)+ind0)))

            ez(ip)=tdzi*(u0*v0*w0*(phi1d(noff( 1)+ind0) - phi1d(noff(20)+ind0))
     &                 + u1*v0*w0*(phi1d(noff( 2)+ind0) - phi1d(noff(21)+ind0))
     &                 + u0*v1*w0*(phi1d(noff( 3)+ind0) - phi1d(noff(24)+ind0))
     &                 + u1*v1*w0*(phi1d(noff( 4)+ind0) - phi1d(noff(25)+ind0))
     &                 + u0*v0*w1*(phi1d(noff( 8)+ind0) - phi1d(noff(29)+ind0))
     &                 + u1*v0*w1*(phi1d(noff( 9)+ind0) - phi1d(noff(30)+ind0))
     &                 + u0*v1*w1*(phi1d(noff(12)+ind0) - phi1d(noff(31)+ind0))
     &                 + u1*v1*w1*(phi1d(noff(13)+ind0) - phi1d(noff(32)+ind0)))
          enddo

        else

c         --- Set offsets for indices on axis of symmetry.  The offsets change
c         --- the sign of the grid cells which are on the negative side
c         --- of the axis of symmetry.
          soy = 2*nnx
          sox = 0
          if (l4symtry) sox = 2

c         --- Set the signs of the E field for particles on negative side of
c         --- the axis of symmetry.
          sy = -1.
          sx = 1.
          if (l4symtry) sx = -1.

c         --- special loop symmetry is used
          do ip = 1, np

            xi =  (abs(xp(ip)) - xmmin)*dxi
            yj =  (abs(yp(ip)) - ymmin)*dyi
            zk =  (zp(ip) - zgrid - zmmin)*dzi
            if (xi < 0. .or. xi > nx .or.
     &          yj < 0. .or. yj > ny .or.
     &          zk < 0. .or. zk > nz) then
              ex(ip) = 0.
              ey(ip) = 0.
              ez(ip) = 0.
              cycle
            endif
            i = int(xi)
            j = int(yj)
            k = int(zk)
            ind0 = i + j*nnx + k*nnxy

            u1 = xi - i
            v1 = yj - j
            w1 = zk - k

            u0 = 1. - u1
            v0 = 1. - v1
            w0 = 1. - w1

c           --- Set offsets for points on symmetry axis of grid.  The offset
c           --- for points off the axis is zero.
            ox = 0
            oy = 0
            if (i == 0 .and. xmmin == 0.) ox = sox
            if (j == 0 .and. ymmin == 0.) oy = soy

c           --- Adjust sign of E field for appropriate quadrant.
            xsign = tdxi
            ysign = tdyi
            if (xp(ip) < 0.) xsign = sx*tdxi
            if (yp(ip) < 0.) ysign = sy*tdyi
           ex(ip)=xsign*(u0*v0*w0*(phi1d(noff( 7)+ind0+ox)-phi1d(noff( 9)+ind0))
     &                 + u1*v0*w0*(phi1d(noff( 8)+ind0   )-phi1d(noff(10)+ind0))
     &                 + u0*v1*w0*(phi1d(noff(11)+ind0+ox)-phi1d(noff(13)+ind0))
     &                 + u1*v1*w0*(phi1d(noff(12)+ind0   )-phi1d(noff(14)+ind0))
     &                 + u0*v0*w1*(phi1d(noff(19)+ind0+ox)-phi1d(noff(21)+ind0))
     &                 + u1*v0*w1*(phi1d(noff(20)+ind0   )-phi1d(noff(22)+ind0))
     &                 + u0*v1*w1*(phi1d(noff(23)+ind0+ox)-phi1d(noff(25)+ind0))
     &                 + u1*v1*w1*(phi1d(noff(24)+ind0   )-phi1d(noff(26)+ind0)))

           ey(ip)=ysign*(u0*v0*w0*(phi1d(noff( 5)+ind0+oy)-phi1d(noff(12)+ind0))
     &                 + u1*v0*w0*(phi1d(noff( 6)+ind0+oy)-phi1d(noff(13)+ind0))
     &                 + u0*v1*w0*(phi1d(noff( 8)+ind0   )-phi1d(noff(15)+ind0))
     &                 + u1*v1*w0*(phi1d(noff( 9)+ind0   )-phi1d(noff(16)+ind0))
     &                 + u0*v0*w1*(phi1d(noff(17)+ind0+oy)-phi1d(noff(24)+ind0))
     &                 + u1*v0*w1*(phi1d(noff(18)+ind0+oy)-phi1d(noff(25)+ind0))
     &                 + u0*v1*w1*(phi1d(noff(20)+ind0   )-phi1d(noff(27)+ind0))
     &                 + u1*v1*w1*(phi1d(noff(21)+ind0   )-phi1d(noff(28)+ind0)))

            ez(ip)=tdzi*(u0*v0*w0*(phi1d(noff( 1)+ind0) - phi1d(noff(20)+ind0))
     &                 + u1*v0*w0*(phi1d(noff( 2)+ind0) - phi1d(noff(21)+ind0))
     &                 + u0*v1*w0*(phi1d(noff( 3)+ind0) - phi1d(noff(24)+ind0))
     &                 + u1*v1*w0*(phi1d(noff( 4)+ind0) - phi1d(noff(25)+ind0))
     &                 + u0*v0*w1*(phi1d(noff( 8)+ind0) - phi1d(noff(29)+ind0))
     &                 + u1*v0*w1*(phi1d(noff( 9)+ind0) - phi1d(noff(30)+ind0))
     &                 + u0*v1*w1*(phi1d(noff(12)+ind0) - phi1d(noff(31)+ind0))
     &                 + u1*v1*w1*(phi1d(noff(13)+ind0) - phi1d(noff(32)+ind0)))
          enddo

        endif

      else if (efetch == 6) then

c     --- Same as efetch=1 but the i,j,k and are precalculated for the next
c     --- particle since the float to integer conversion can be slow

      if (.not. (l2symtry .or. l4symtry)) then
        inext = (xp(1) - xmmin) * dxi
        jnext = (yp(1) - ymmin) * dyi
        knext = (zp(1) - zgrid - zmmin) * dzi
        do ip = 1, np

           i = inext
           j = jnext
           k = knext
           if (ip < np) then
             inext = (xp(ip+1) - xmmin) * dxi
             jnext = (yp(ip+1) - ymmin) * dyi
             knext = (zp(ip+1) - zgrid - zmmin) * dzi
           endif
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
        inext = (abs(xp(1)) - xmmin)*dxi
        jnext = (abs(yp(1)) - ymmin)*dyi
        knext = (zp(1) - zgrid - zmmin)*dzi
        do ip = 1, np

           i = inext
           j = jnext
           k = knext
           if (ip < np) then
             inext = (abs(xp(ip+1)) - xmmin)*dxi
             jnext = (abs(yp(ip+1)) - ymmin)*dyi
             knext = (zp(ip+1) - zgrid - zmmin)*dzi
           endif
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
      subroutine stptcl3d(pgroup)
      use ParticleGroupmodule
      use GlobalVars
      use Subtimers3d
      use Beam_acc
      use Constant
      use InGen
      use InGen3d
      use InPart
      use InPart3d
      use InMesh3d
      use Particles,Only: npmax
      use Picglb
      use Picglb3d
#ifdef MPIPARALLEL
      use Parallel
#endif
      type(ParticleGroup):: pgroup

c  Loads particles.
c  At end, particle boundary conditions are enforced.
c
c  For the parallel version, each processor loads the fraction of the
c  longitudinally uniform beam of length (zimax-zimin) which is within its
c  region. The routine sw_globalsum recalculates the particle weights (sw)
c  appropriately.
c
c  Picglb is included only to use zbeam

      integer(ISZ):: npm,ip,js,isid,i,j,k,i1
      integer(ISZ):: nxstri_s,nystri_s,nzstri_s,izstripe,icheck
      integer(ISZ):: nfib2,nfib3,nfib4,ii,ioff
      real(kind=8):: rnpmi,zmid,zlen,zleni,vtx,vty,vtz,rr,r0
      real(kind=8):: rpp,phi1,phi2,rp,wz,wr
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
      integer(ISZ),allocatable,dimension(:):: indx,nptot
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

c  Set npm to be greater than zero if any particles are to be loaded.
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
      if (npm > 0 .and. distrbtn .ne. "preload") then

      call setuppgroup(pgroup)

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

      do js=0,pgroup%ns-1
        if (pgroup%sid(js) == -1) pgroup%sid(js) = js
        isid = pgroup%sid(js)+1

c       --- Get number of particles to be loaded for this species.
        if (np_s(isid) == 0) then
c         --- Only set np_s is it had not been set be the user
          if (xrandom == "fibonacc") then
            np_s(isid) = nfibgrps*fibg1*sp_fract(isid)
          elseif (xrandom == "digitrev") then
            np_s(isid) = npmax*sp_fract(isid)
          elseif (xrandom == "pseudo") then
            np_s(isid) = npmax*sp_fract(isid)
          elseif (xrandom == "grid") then
            np_s(isid) = nxstripe*nystripe*nzstripe*sp_fract(isid)
          endif
        endif
        npm = np_s(isid)

        zmid = .5 * (zimax_s(isid) + zimin_s(isid)) + zbeam
        zlen = zimax_s(isid) - zimin_s(isid)
        zleni = 1./zlen
        if (zlen == 0) zleni = 0.

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
        zmin = max(zimin_s(isid),zpslmin(my_index))
        zmax = min(zimax_s(isid),zpslmax(my_index))
c       --- calculate number of particles that are loaded by processors
c       --- to the left
        ioff = int(npm*(zmin - zimin_s(isid))*zleni + 1.e-5)
c       --- calculate axial offset so the first particle for this processor
c       --- is the correct distance from the last particle of the processor
c       --- to the left.
        zoff = ioff*rnpmi
c       --- number of particles to be loaded in this processor: the number
c       --- of particles between zmax and zimin minus the number between
c       --- zmin and zimin.
        npm = int(npm*(zmax - zimin_s(isid))*zleni + 1.e-5) - ioff
c       --- Make sure that there is enough room to load the particles.
c       --- This is only needed here in the parallel version since it is
c       --- difficult to gaurantee that there is room when an arbitrary
c       --- domain decomposition is allowed.
#endif

!$OMP MASTER
c       --- Make sure there is space allocated for the new particles.
        call chckpart(pgroup,js+1,0,npm,.false.)
!$OMP END MASTER

        ip = pgroup%ins(js+1) - 1
!$OMP SINGLE
#ifdef _OPENMP
        iptotal = pgroup%ins(js+1) - 1
#endif
!$OMP END SINGLE

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
           nzstri_s = nzstripe*sp_fract(isid)**(1./3.)
           izstripe = (k-1)/npgrp + 1
           do i=k,i1
              zt(i-k+1) = (izstripe-.5)/nzstri_s + zjig*(wranf()-.5)/nzstri_s +
     &                    zoff
           enddo
         endif

c        --- zero vz if vthz is zero only the first time through
         if (vthz_s(isid) == 0. .and. firstk) then
           do i=1,npgrp
              uzt(i) = 0.
           enddo
         endif

         if (cigarld) then
c          --- set vtz to vthz*2 since vtz is vz_max - vz_bar for cigar
           vtz = 2.*vthz_s(isid)
c          --- load uzt with linear distribution (-.5 < uzt < .5)
           if (vtz /= 0.) then
             do i=1,i1-k+1
               uzt(i) = wrandom(vzrandom,i+k-1+ioff,dig7,fibg1,npm) - .5
             enddo
           endif
c          --- use cigar to adjust the z's
c          --- xt,yt,uxt,uyt passed as scratch arrays
           call cigar(i1-k+1,zt(1),uzt(1),zt(1),uzt(1),perpscal(1),
     &                straight_s(isid),xt(1),yt(1),uxt(1),uyt(1))

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
           vtz = vthz_s(isid)
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
         zt = zbeam + zimin_s(isid) + zlen*zt

c        --- fetch envelope, centroid, and emittance.  Reset values 
c        --- after initial set if axial varying emittance, axial velocity, 
c        --- and axial thermal velocity are input
         if (nenvofz > 0) then
           do i=1,i1-k+1
c            --- envelope interpolations 
             iz = (zt(i) - zimin_s(isid) - zbeam)*zleni*nenvofz
             wz = (zt(i) - zimin_s(isid) - zbeam)*zleni*nenvofz - iz
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
             vtzt(i) = vthz_s(isid)
           enddo 
         else
c          --- set only first time through
           if (firstk) then
             do i=1,npgrp
               at(i)   = a0_s(isid)
               bt(i)   = b0_s(isid)
               apt(i)  = ap0_s(isid)
               bpt(i)  = bp0_s(isid)
c
               xct(i)  = x0_s(isid)  + xcent_s(isid)
               yct(i)  = y0_s(isid)  + ycent_s(isid)
               xpct(i) = xp0_s(isid) + xpcent_s(isid)
               ypct(i) = yp0_s(isid) + ypcent_s(isid)
c
               epsxt(i) = emitx_s(isid) 
               epsyt(i) = emity_s(isid) 
c
               vzt(i)  = vbeam_s(isid) 
               vtzt(i) = vthz_s(isid) 
c
               ibeamt(i) = ibeam_s(isid)
             enddo
           endif
         endif
c        --- reset emittance interpolations               
         if (nemitofz > 0) then 
           do i=1,i1-k+1
             iz = (zt(i) - zimin_s(isid) - zbeam)*zleni*nemitofz
             wz = (zt(i) - zimin_s(isid) - zbeam)*zleni*nemitofz - iz
c
             epsxt(i) = emitxofz(iz)*(1.-wz) + emitxofz(iz+1)*wz
             epsyt(i) = emityofz(iz)*(1.-wz) + emityofz(iz+1)*wz
           enddo
         endif  
c        --- reset axial velocity interpolations 
         if (nvbeamofz > 0) then 
           do i=1,i1-k+1
c
             iz = (zt(i) - zimin_s(isid) - zbeam)*zleni*nvbeamofz
             wz = (zt(i) - zimin_s(isid) - zbeam)*zleni*nvbeamofz - iz
c
             vzt(i) = vbeamofz(iz)*(1.-wz) + vbeamofz(iz+1)*wz
           enddo 
         endif 
c        --- reset axial thermal velocity interpolations 
         if (nvthzofz > 0) then 
           do i=1,i1-k+1
c
             iz = (zt(i) - zimin_s(isid) - zbeam)*zleni*nvthzofz
             wz = (zt(i) - zimin_s(isid) - zbeam)*zleni*nvthzofz - iz
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
         if (distrbtn(:8) == "semigaus" .or. distrbtn == "SG" .or. 
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
             nxstri_s = nxstripe*sp_fract(isid)**(1./3.)
             nystri_s = nystripe*sp_fract(isid)**(1./3.)
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
             r0 = max(a0_s(isid),b0_s(isid))
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
c            --- Force the first value to be done in such a way so that the
c            --- calculation of the gaussian random numbers is reset (since
c            --- they are calculated two at a time. This is only really
c            --- needed for the parallel version to get the same numbers
c            --- as the serial case.
             uyt(1) = wrandomgauss(vtrandom,1+k-1+ioff,dig5,dig6,1,1,.false.)
             do i=2,i1-k+1
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
             pgroup%xp(ip+i) = at(ii)*xt(ii)*perpscal(ii) + xct(ii)
             pgroup%yp(ip+i) = bt(ii)*yt(ii)*perpscal(ii) + yct(ii)
             pgroup%zp(ip+i) = zt(ii)
c            --- put velocity into particle arrays
             vtx = .5*vzt(ii)*epsxt(ii)/at(ii)*perpscal(ii) + vthperp_s(isid)
             vty = .5*vzt(ii)*epsyt(ii)/bt(ii)*perpscal(ii) + vthperp_s(isid)
             pgroup%uxp(ip+i) = vzt(ii)*apt(ii)*xt(ii)*perpscal(ii) +
     &                   vtx*uxt(ii) + vzt(ii)*xpct(ii)
             pgroup%uyp(ip+i) = vzt(ii)*bpt(ii)*yt(ii)*perpscal(ii) +
     &                   vty*uyt(ii) + vzt(ii)*ypct(ii)
             pgroup%uzp(ip+i) = vzt(ii)*(1.+vtilt_s(isid)*(zmid - zt(ii))*zleni) +
     &                   vtzt(ii)*uzt(ii)
             if (pgroup%npid > 0) pgroup%pid(ip+i,:) = 0.
           enddo

c          --- add sinusoidal perturbation to uzp
           if (vzperamp /= 0.) then
             do i=ip+1,ip+j
               pgroup%uzp(i) = pgroup%uzp(i) +
     &                   vzperamp*sin(2.*pi*pgroup%zp(i)/vzperlam + vzperphs)
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
             pgroup%xp(ip+i) = rr*cos(phi1)*at(i)*perpscal(i)
             pgroup%yp(ip+i) = rr*sin(phi1)*bt(i)*perpscal(i)
             pgroup%zp(ip+i) = zt(i)
c            --- put velocity into particle arrays
             rr = sqrt(1.-rr*rr)
             vtx = vzt(i)*epsxt(i)/at(i)*perpscal(i) + vthperp_s(isid)
             vty = vzt(i)*epsyt(i)/bt(i)*perpscal(i) + vthperp_s(isid)
             pgroup%uxp(ip+i) = vzt(i)*pgroup%xp(ip+i)*apt(i)/at(i) +
     &                   vtx*rr*cos(phi2) + vzt(i)*xpct(i)
             pgroup%uyp(ip+i) = vzt(i)*pgroup%yp(ip+i)*bpt(i)/bt(i) +
     &                   vty*rr*sin(phi2) + vzt(i)*ypct(i)
             pgroup%uzp(ip+i) = vzt(i)*(1.+vtilt_s(isid)*(zmid - zt(i))*zleni) +
     &                   vtzt(i)*uzt(i)
c            --- The centroids are added to the transverse positions
c            --- after the variables are used above to calculate the
c            --- coherent velocity component.
             pgroup%xp(ip+i) = pgroup%xp(ip+i) + xct(i)
             pgroup%yp(ip+i) = pgroup%yp(ip+i) + yct(i)
           enddo

c          --- add sinusoidal perturbation to uzp
           if (vzperamp /= 0.) then
             do i=ip+1,ip+i1-k+1
               pgroup%uzp(i) = pgroup%uzp(i) +
     &               vzperamp*sin(2.*pi*pgroup%zp(i)/vzperlam + vzperphs)
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
             pgroup%xp(ip+i) = xt(i) + xct(i)
             pgroup%yp(ip+i) = yt(i) + yct(i)
c            --- set the longitudinal particle coordinate 
             pgroup%zp(ip+i) = zt(i)
c            --- set transverse particle gamma*velocity   
c                  1st term: rescale uxp returned in x' units to uxp  
c                  2nd term: centroid component 
             pgroup%uxp(ip+i) = vzt(i)*uxt(i) + vzt(i)*xpct(i)
             pgroup%uyp(ip+i) = vzt(i)*uyt(i) + vzt(i)*ypct(i)
c            --- set longitudinal particle gamma*velocity 
             pgroup%uzp(ip+i) = vzt(i)*(1. + vtilt_s(isid)*(zmid - zt(i))*zleni) +
     &                   vtzt(i)*uzt(i)
             if (pgroup%npid > 0) pgroup%pid(ip+i,:) = 0.
           enddo

c          --- add sinusoidal perturbation to uzp
           if (vzperamp /= 0.) then
             do i=ip+1,ip+i1-k+1
               pgroup%uzp(i) = pgroup%uzp(i) +
     &               vzperamp*sin(2.*pi*pgroup%zp(i)/vzperlam + vzperphs)
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
      pgroup%nps(js+1) = ip - pgroup%ins(js+1) + 1
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

c     --- end of if checking for positive npm and distrbtn is not preload
      endif

c     --- Make sure that sid is set if it is not set by the user.
      do js=0,pgroup%ns-1
        if (pgroup%sid(js) == -1) pgroup%sid(js) = js
      enddo

c     --- Find maximum of zdist, which is used to scale sw.  See below.
      if (nzdist > 0) then
        zdistmax = nrmzdist(0)
        do iz=1,nzdist
          if (nrmzdist(iz) > zdistmax) zdistmax = nrmzdist(iz)
        enddo
      endif

c     --- Get the global number of particles
      allocate(nptot(0:pgroup%ns-1))
      nptot = pgroup%nps
#ifdef MPIPARALLEL
      call parallelsumintegerarray(nptot,pgroup%ns)
#endif

c     --- calculate values for sm, sq, and sw for each species.
      do js = 0,pgroup%ns-1
        isid = pgroup%sid(js) + 1
        if (pgroup%sq(js+1) == 0) pgroup%sq(js+1) = zion_s(isid) * echarge
        if (pgroup%sm(js+1) == 0) pgroup%sm(js+1) = aion_s(isid) * amu
        if (vbeam_s(isid) /= 0. .and. nptot(js) /= 0 .and.
     &      zion_s(isid) /= 0.  .and. pgroup%sw(js+1) == 0) then
          zlen = zimax_s(isid) - zimin_s(isid)
          pgroup%sw(js+1)=ibeam_s(isid)*zlen*sp_fract(isid)/
     &                   (vbeam_s(isid)*echarge*zion_s(isid)*nptot(js))
        endif
c       --- Don't zero out sw in case the user has set it.

c       --- Adjust weighting of particles since cigar() makes beam more dense
        if (cigarld) then
          pgroup%sw(js+1) = pgroup%sw(js+1)*
     &                        (straight_s(isid) + (1. - straight_s(isid))*2./3.)
        endif

c       --- Scale sw when using zdist so that ibeam/vbeam gives the maximum
c       --- linecharge of the distribution.  sw is scaled by the ratio of the
c       --- sum of zdist and the the product of the max of zdist and the
c       --- number of grid points in zdist.  That is the ratio of the loaded
c       --- linecharge and what the linecharge would have been with a uniform
c       --- distribution.
        if (nzdist > 0) then
          pgroup%sw(js+1) = pgroup%sw(js+1)/(zdistmax*nzdist)
        endif

      enddo

      deallocate(nptot)

c     --- This is not needed since the global number of particles is used
c     --- above to calculate sw.
c     --- if slave, recalculate sw
c#ifdef MPIPARALLEL
c      call sw_globalsum(pgroup%ns,pgroup%sw)
c#endif

c     --- Apply particle boundary conditions
      call zpartbnd(pgroup,zmmax,zmmin,dz,zgrid)
      do js=0,pgroup%ns-1
        if (pgroup%nps(js+1) > 0) then
           ip = pgroup%ins(js+1)
           call stckxy3d(pgroup%nps(js+1),
     &                   pgroup%xp(ip),xmmax,xmmin,dx,
     &                   pgroup%yp(ip),ymmax,ymmin,dy,
     &                   pgroup%zp(ip),zmmin,dz,pgroup%uxp(ip),
     &                   pgroup%uyp(ip),pgroup%uzp(ip),pgroup%gaminv(ip),
     &                   zgrid,zbeam,l2symtry,l4symtry,pboundxy,.true.)
        endif
        call processlostpart(pgroup,js+1,clearlostpart,time,zbeam)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timestptcl3d = timestptcl3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setrho3d(rho,rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,depos,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry,
     &                    lcylindrical)
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
      logical(ISZ):: l2symtry,l4symtry,lcylindrical

c  Sets charge density using various algorithms

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c--------------------------------------
c  Vectorized deposition loop
c--------------------------------------
      if (depos == "vector") then

        call setrho3dvector(rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                      nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                      l2symtry,l4symtry)

c--------------------------------------
c  Scalar deposition loop
c--------------------------------------
      elseif (depos == "scalar") then

        call setrho3dscalar(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                      nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                      l2symtry,l4symtry)

c--------------------------------------
c  Direct deposition loop
c--------------------------------------
      elseif (depos == "direct") then

        call setrho3ddirect(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                      nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                      l2symtry,l4symtry)

c--------------------------------------
c  Direct deposition loop with precalculated integer conversions
c--------------------------------------
      elseif (depos == "direct1") then

        if (.not. lcylindrical) then
        call setrho3ddirect1(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                       nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                       l2symtry,l4symtry)

        else
        call setrho3ddirect2(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                       nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                       l2symtry,l4symtry,lcylindrical)
        endif

c--------------------------------------
c  Vectorized deposition loop with precalculated integer conversions
c--------------------------------------
      else if (depos == "vector1") then

        call setrho3dvector1(rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                       nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                       l2symtry,l4symtry)

c--------------------------------------
c  Direct deposition using a second order spline
c--------------------------------------
      else if (depos == "dspline2") then

        call setrho3ddirectspline2(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                             nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,
     &                             l2symtry,l4symtry)

      endif

!$OMP END PARALLEL

      if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine setrho3dvector(rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry)
      use GlobalVars
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho1d(0:(1+nx)*(1+ny)*(1+nz)-1)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c  Sets charge density

c  Algorithm notes: rho array is dimensioned (0:nx,0:ny,0:nz) outside,
c  but is made one dimensional in this routine
c  so cell index into 1d rho array for vectorized deposition is:
c     i + j*(nx+1) + k*(nx+1)*(ny+1)
c  In each case,
c     rho(i  ,j  ,k  ) = rho(i  ,j  ,k  ) + u0*v0*w0*g
c     rho(i+1,j  ,k  ) = rho(i+1,j  ,k  ) + u1*v0*w0*g
c  Note that many changes are possible; for example, we might define
c  ind0(ir) and not use indx; this saves some store operations but
c  leads to a more complicated indirect address for the vectorized
c  gather-add-scatter loop.  It seems about 3% slower than the present way.
c  RHO must be zeroed in ZERORHO since it is not zeroed here (to allow
c  handling of blocks of particles at a time)

      integer(ISZ):: moff(0:7)
      integer(ISZ),allocatable:: indx(:,:)
      real(kind=8),allocatable:: s(:,:)

      integer(ISZ):: ipmin,nptmp,ip,i,j,k,ind0,m,ir
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact

c     --- Set up offset array for vectorized deposition:
      moff(0) = 0
      moff(1) = 1
      moff(2) = nx+1
      moff(3) = nx+2
      moff(4) = (nx+1)*(ny+1)
      moff(5) = (nx+1)*(ny+1)+1
      moff(6) = (nx+1)*(ny+1)+nx+1
      moff(7) = (nx+1)*(ny+1)+nx+2

      g = wght*q/(dx*dy*dz)
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
!$OMP&                 gyfact,gxfact,m)
!$OMP&FIRSTPRIVATE(np)

      allocate(indx(0:7,0:nparpgrp-1),s(0:7,0:nparpgrp-1))

!$OMP DO
      do ipmin = 1,np,nparpgrp
        nptmp = min(nparpgrp, np+1-ipmin)

c       --- vectorized loop to compute indices, weights
        if (l2symtry) then
c         --- special loop for 2-fold symmetry
c         --- The particle weight is reduced by a factor of 2 except near the
c         --- transverse boundaries.
          do ip = ipmin,ipmin+nptmp-1
            i  = int((xp(ip) - xmmin)*dxi)
            u1 =     (xp(ip) - xmmin)*dxi - i
            u0 = 1. - u1
            j  = int((abs(yp(ip)) - ymmin)*dyi)
            v1 =     (abs(yp(ip)) - ymmin)*dyi - j
            v0 = 1. - v1
            k  = int((zp(ip) - zgrid - zmmin)*dzi)
            w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
            s(0,ir) = u0*v0*w0*g*gyfact
            s(1,ir) = u1*v0*w0*g*gyfact
            s(2,ir) = u0*v1*w0*g
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g*gyfact
            s(5,ir) = u1*v0*w1*g*gyfact
            s(6,ir) = u0*v1*w1*g
            s(7,ir) = u1*v1*w1*g
          enddo
        elseif (l4symtry) then
c         --- special loop for 4-fold symmetry
c         --- The particle weight is reduced by a factor of 4 except near the
c         --- transverse boundaries.
          do ip = ipmin,ipmin+nptmp-1
            i  = int((abs(xp(ip)) - xmmin)*dxi)
            u1 =     (abs(xp(ip)) - xmmin)*dxi - i
            u0 = 1. - u1
            j  = int((abs(yp(ip)) - ymmin)*dyi)
            v1 =     (abs(yp(ip)) - ymmin)*dyi - j
            v0 = 1. - v1
            k  = int((zp(ip) - zgrid - zmmin)*dzi)
            w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
            s(0,ir) = u0*v0*w0*g*gxfact*gyfact
            s(1,ir) = u1*v0*w0*g*gyfact
            s(2,ir) = u0*v1*w0*g*gxfact
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g*gxfact*gyfact
            s(5,ir) = u1*v0*w1*g*gyfact
            s(6,ir) = u0*v1*w1*g*gxfact
            s(7,ir) = u1*v1*w1*g
          enddo
        else
c         --- normal loop
          do ip = ipmin,ipmin+nptmp-1
            i  = int((xp(ip) - xmmin)*dxi)
            u1 =     (xp(ip) - xmmin)*dxi - i
            u0 = 1. - u1
            j  = int((yp(ip) - ymmin)*dyi)
            v1 =     (yp(ip) - ymmin)*dyi - j
            v0 = 1. - v1
            k  = int((zp(ip) - zgrid - zmmin)*dzi)
            w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
            s(0,ir) = u0*v0*w0*g
            s(1,ir) = u1*v0*w0*g
            s(2,ir) = u0*v1*w0*g
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g
            s(5,ir) = u1*v0*w1*g
            s(6,ir) = u0*v1*w1*g
            s(7,ir) = u1*v1*w1*g
          enddo
        endif

c       --- vectorized deposition over the 8 cells touched;
c       --- there'd be a hazard if we interchanged the loops.
!$OMP CRITICAL (CRITICAL_SETRHO3DVECTOR)
        do ir = 0,nptmp-1
          if (uzp(ir+ipmin) /= 0.) then
            do m = 0, 7
              rho1d(indx(m,ir)) = rho1d(indx(m,ir)) + s(m,ir)
            enddo
          endif
        enddo
!$OMP END CRITICAL (CRITICAL_SETRHO3DVECTOR)

      enddo
!$OMP END DO

      deallocate(indx,s)
!$OMP END PARALLEL

      return
      end
c=============================================================================
      subroutine setrho3dscalar(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry)
      use GlobalVars
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c  Sets charge density
c  Similar to vector, but rho is treated as a 3d array rather than a 1d array

      integer(ISZ),allocatable:: ii(:), jj(:), kk(:)
      real(kind=8),allocatable:: s(:,:)

      integer(ISZ):: ipmin,nptmp,ip,i,j,k,ind0,ir
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact

      g = wght*q/(dx*dy*dz)
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
!$OMP PARALLEL PRIVATE(ipmin,nptmp,u1,u0,v1,v0,w1,w0,ir,ip,
!$OMP&                 gyfact,gxfact,s,ii,jj,kk)
!$OMP&FIRSTPRIVATE(np)

      allocate(ii(0:nparpgrp-1), jj(0:nparpgrp-1), kk(0:nparpgrp-1))
      allocate(s(0:7,0:nparpgrp-1))

!$OMP DO
      do ipmin = 1,np,nparpgrp
        nptmp = min(nparpgrp, np+1-ipmin)

c       --- vectorized loop to compute indices, weights
        if (l2symtry) then
          do ip = ipmin,ipmin+nptmp-1
            ir = ip - ipmin
            ii(ir) = int((xp(ip) - xmmin)*dxi)
            u1     =     (xp(ip) - xmmin)*dxi - ii(ir)
            u0     = 1. - u1
            jj(ir) = int((abs(yp(ip)) - ymmin)*dyi)
            v1     =     (abs(yp(ip)) - ymmin)*dyi - jj(ir)
            v0     = 1. - v1
            kk(ir) = int((zp(ip) - zgrid - zmmin)*dzi)
            w1     =     (zp(ip) - zgrid - zmmin)*dzi - kk(ir)
            w0     = 1. - w1
            gyfact = 1.
            if (jj(ir) == 0 .and. ymmin == 0.) gyfact = 2.
            s(0,ir) = u0*v0*w0*g*gyfact
            s(1,ir) = u1*v0*w0*g*gyfact
            s(2,ir) = u0*v1*w0*g
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g*gyfact
            s(5,ir) = u1*v0*w1*g*gyfact
            s(6,ir) = u0*v1*w1*g
            s(7,ir) = u1*v1*w1*g
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
            kk(ir) = int((zp(ip) - zgrid - zmmin)*dzi)
            w1     =     (zp(ip) - zgrid - zmmin)*dzi - kk(ir)
            w0     = 1. - w1
            gxfact = 1.
            gyfact = 1.
            if (ii(ir) == 0 .and. xmmin == 0.) gxfact = 2.
            if (jj(ir) == 0 .and. ymmin == 0.) gyfact = 2.
            s(0,ir) = u0*v0*w0*g*gxfact*gyfact
            s(1,ir) = u1*v0*w0*g*gyfact
            s(2,ir) = u0*v1*w0*g*gxfact
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g*gxfact*gyfact
            s(5,ir) = u1*v0*w1*g*gyfact
            s(6,ir) = u0*v1*w1*g*gxfact
            s(7,ir) = u1*v1*w1*g
          enddo
        else
c         --- normal loop
          do ip = ipmin,ipmin+nptmp-1
            ir = ip - ipmin
            ii(ir) = int((xp(ip) - xmmin)*dxi)
            u1     = (xp(ip) - xmmin)*dxi - ii(ir)
            u0     = 1. - u1
            jj(ir) = int((yp(ip) - ymmin)*dyi)
            v1     = (yp(ip) - ymmin)*dyi - jj(ir)
            v0     = 1. - v1
            kk(ir) = int((zp(ip) - zgrid - zmmin)*dzi)
            w1     = (zp(ip) - zgrid - zmmin)*dzi - kk(ir)
            w0     = 1. - w1
            s(0,ir) = u0*v0*w0*g
            s(1,ir) = u1*v0*w0*g
            s(2,ir) = u0*v1*w0*g
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g
            s(5,ir) = u1*v0*w1*g
            s(6,ir) = u0*v1*w1*g
            s(7,ir) = u1*v1*w1*g
          enddo
        endif
c       --- scalar loop does the actual deposition
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

      enddo
!$OMP END DO

      deallocate(ii,jj,kk,s)

!$OMP END PARALLEL

      return
      end
c=============================================================================
      subroutine setrho3ddirect(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry)
      use GlobalVars
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c  Sets charge density
c  No particle blocks are used (since there are no temporary arrays).

      integer(ISZ):: ip,ii,jj,kk
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact

      g = wght*q/(dx*dy*dz)
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
!$OMP PARALLEL PRIVATE(ii,jj,kk,u1,u0,v1,v0,w1,w0,ip,gyfact,gxfact)
!$OMP&FIRSTPRIVATE(np)

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = int((xp(ip) - xmmin)*dxi)
          u1 =     (xp(ip) - xmmin)*dxi - ii
          u0 = 1. - u1
          jj = int((abs(yp(ip)) - ymmin)*dyi)
          v1 =     (abs(yp(ip)) - ymmin)*dyi - jj
          v0 = 1. - v1
          kk = int((zp(ip) - zgrid - zmmin)*dzi)
          w1 =     (zp(ip) - zgrid - zmmin)*dzi - kk
          w0 = 1. - w1
          gyfact = 1.
          if (jj == 0 .and. ymmin == 0.) gyfact = 2.
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0*g*gyfact
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0*g*gyfact
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0*g
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0*g
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1*g*gyfact
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1*g*gyfact
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1*g
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1*g
        enddo
!$OMP END DO
      elseif (l4symtry) then
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = int((abs(xp(ip)) - xmmin)*dxi)
          u1 =     (abs(xp(ip)) - xmmin)*dxi - ii
          u0 = 1. - u1
          jj = int((abs(yp(ip)) - ymmin)*dyi)
          v1 =     (abs(yp(ip)) - ymmin)*dyi - jj
          v0 = 1. - v1
          kk = int((zp(ip) - zgrid - zmmin)*dzi)
          w1 =     (zp(ip) - zgrid - zmmin)*dzi - kk
          w0 = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii == 0 .and. xmmin == 0.) gxfact = 2.
          if (jj == 0 .and. ymmin == 0.) gyfact = 2.
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0*g*gxfact*gyfact
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0*g*gyfact
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0*g*gxfact
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0*g
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1*g*gxfact*gyfact
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1*g*gyfact
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1*g*gxfact
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1*g
        enddo
!$OMP END DO
      else
c       --- normal loop
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = int((xp(ip) - xmmin)*dxi)
          u1 = (xp(ip) - xmmin)*dxi - ii
          u0 = 1. - u1
          jj = int((yp(ip) - ymmin)*dyi)
          v1 = (yp(ip) - ymmin)*dyi - jj
          v0 = 1. - v1
          kk = int((zp(ip) - zgrid - zmmin)*dzi)
          w1 = (zp(ip) - zgrid - zmmin)*dzi - kk
          w0 = 1. - w1
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0*g
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0*g
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0*g
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0*g
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1*g
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1*g
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1*g
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1*g
        enddo
!$OMP END DO
      endif

!$OMP END PARALLEL

      return
      end
c=============================================================================
      subroutine setrho3ddirect1(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry)
      use GlobalVars
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c  Sets charge density
c  No particle blocks are used (since there are no temporary arrays).
c  Also, the float to integer conversions are precalculated.
c  This seems to be the fastest version.

      integer(ISZ):: ip,ii,jj,kk
      integer(ISZ):: iinext,jjnext,kknext
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact
      real(kind=8):: u1next,v1next,w1next

      g = wght*q/(dx*dy*dz)
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
!$OMP PARALLEL PRIVATE(ii,jj,kk,u1,u0,v1,v0,w1,w0,ip,gyfact,gxfact,
!$OMP&                 iinext,jjnext,kknext,u1next,v1next,w1next)
!$OMP&FIRSTPRIVATE(np)

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        iinext = int((xp(1) - xmmin)*dxi)
        jjnext = int((abs(yp(1)) - ymmin)*dyi)
        kknext = int((zp(1) - zgrid - zmmin)*dzi)
        u1next = (xp(1) - xmmin)*dxi - iinext
        v1next = (abs(yp(1)) - ymmin)*dyi - jjnext
        w1next = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = iinext
          jj = jjnext
          kk = kknext
          u1 = u1next
          v1 = v1next
          w1 = w1next
          if (ip < np) then
            iinext = int((xp(ip+1) - xmmin)*dxi)
            jjnext = int((abs(yp(ip+1)) - ymmin)*dyi)
            kknext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            u1next = (xp(ip+1) - xmmin)*dxi - iinext
            v1next = (abs(yp(ip+1)) - ymmin)*dyi - jjnext
            w1next = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1
          gyfact = 1.
          if (jj == 0 .and. ymmin == 0.) gyfact = 2.
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0*g*gyfact
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0*g*gyfact
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0*g
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0*g
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1*g*gyfact
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1*g*gyfact
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1*g
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1*g
        enddo
!$OMP END DO
      elseif (l4symtry) then
        iinext = int((abs(xp(1)) - xmmin)*dxi)
        jjnext = int((abs(yp(1)) - ymmin)*dyi)
        kknext = int((zp(1) - zgrid - zmmin)*dzi)
        u1next = (abs(xp(1)) - xmmin)*dxi - iinext
        v1next = (abs(yp(1)) - ymmin)*dyi - jjnext
        w1next = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = iinext
          jj = jjnext
          kk = kknext
          u1 = u1next
          v1 = v1next
          w1 = w1next
          if (ip < np) then
            iinext = int((abs(xp(ip+1)) - xmmin)*dxi)
            jjnext = int((abs(yp(ip+1)) - ymmin)*dyi)
            kknext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            u1next = (abs(xp(ip+1)) - xmmin)*dxi - iinext
            v1next = (abs(yp(ip+1)) - ymmin)*dyi - jjnext
            w1next = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii == 0 .and. xmmin == 0.) gxfact = 2.
          if (jj == 0 .and. ymmin == 0.) gyfact = 2.
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0*g*gxfact*gyfact
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0*g*gyfact
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0*g*gxfact
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0*g
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1*g*gxfact*gyfact
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1*g*gyfact
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1*g*gxfact
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1*g
        enddo
!$OMP END DO
      else
c       --- normal loop
        iinext = int((xp(1) - xmmin)*dxi)
        jjnext = int((yp(1) - ymmin)*dyi)
        kknext = int((zp(1) - zgrid - zmmin)*dzi)
        u1next = (xp(1) - xmmin)*dxi - iinext
        v1next = (yp(1) - ymmin)*dyi - jjnext
        w1next = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = iinext
          jj = jjnext
          kk = kknext
          u1 = u1next
          v1 = v1next
          w1 = w1next
          if (ip < np) then
            iinext = int((xp(ip+1) - xmmin)*dxi)
            jjnext = int((yp(ip+1) - ymmin)*dyi)
            kknext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            u1next = (xp(ip+1) - xmmin)*dxi - iinext
            v1next = (yp(ip+1) - ymmin)*dyi - jjnext
            w1next = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0*g
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0*g
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0*g
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0*g
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1*g
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1*g
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1*g
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1*g
        enddo
!$OMP END DO
      endif

!$OMP END PARALLEL

      return
      end
c=============================================================================
      subroutine setrho3ddirect2(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry,
     &                    lcylindrical)
      use GlobalVars
      use Constant,Only: pi
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry,lcylindrical

c  Sets charge density
c  No particle blocks are used (since there are no temporary arrays).
c  Also, the float to integer conversions are precalculated.
c  This seems to be the fastest version.

      integer(ISZ):: ip,ii,jj,kk,ix
      integer(ISZ):: iinext,jjnext,kknext
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact
      real(kind=8):: u1next,v1next,w1next

c     --- Work array holding q/cell volume, the charge density per
c     --- real particle. This is primarily needed for the RZ version
c     --- since the cell volume there has radial dependence.
      real(kind=8):: cdens(0:nx)
c     --- Temp arrays to hold particle data
c     --- These are needed when lcylindrical is true, in which case x=r, y=0.
      real(kind=8):: x,y

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      if (lcylindrical) then
        if (xmmin == 0.) then
c         --- The factor 0.75 corrects for overdeposition due to linear
c         --- weighting (for uniform distribution)
c         --- see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
c         --- and Verboncoeur, J. of Comp. Phys.,
          cdens(0) = 0.75/(pi*(0.5*0.5*dx*dx*dz))
        else
          cdens(0) = 1./(2.*pi*(xmmin)*dx*dz)
        endif
        do ix = 1,nx
          cdens(ix) = 1./(2.*pi*(ix*dx+xmmin)*dx*dz)
        enddo
        cdens = cdens*wght*q
      else
        g = wght*q/(dx*dy*dz)
        if (l2symtry) then
c         --- The particle weight is reduced by a factor of 2 except near the
c         --- transverse boundaries.
          g = g*0.5
        elseif (l4symtry) then
c         --- The particle weight is reduced by a factor of 4 except near the
c         --- transverse boundaries.
          g = g*0.25
        endif
        cdens = g
      endif

! np was made FIRSTPRIVATE to get around a bug when the expression
! np+1-ipmin was evaluating to 1-ipmin (as if np was zero).
! I don't know why it works, but it does.
!$OMP PARALLEL PRIVATE(ii,jj,kk,u1,u0,v1,v0,w1,w0,ip,gyfact,gxfact,
!$OMP&                 iinext,jjnext,kknext,u1next,v1next,w1next)
!$OMP&FIRSTPRIVATE(np)

c     --- vectorized loop to compute indices, weights
      if (l2symtry .and. .not. lcylindrical) then
        iinext = int((xp(1) - xmmin)*dxi)
        jjnext = int((abs(yp(1)) - ymmin)*dyi)
        kknext = int((zp(1) - zgrid - zmmin)*dzi)
        u1next = (xp(1) - xmmin)*dxi - iinext
        v1next = (abs(yp(1)) - ymmin)*dyi - jjnext
        w1next = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = iinext
          jj = jjnext
          kk = kknext
          u1 = u1next
          v1 = v1next
          w1 = w1next
          if (ip < np) then
            iinext = int((xp(ip+1) - xmmin)*dxi)
            jjnext = int((abs(yp(ip+1)) - ymmin)*dyi)
            kknext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            u1next = (xp(ip+1) - xmmin)*dxi - iinext
            v1next = (abs(yp(ip+1)) - ymmin)*dyi - jjnext
            w1next = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1
          gyfact = 1.
          if (jj == 0 .and. ymmin == 0.) gyfact = 2.
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0*g*gyfact
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0*g*gyfact
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0*g
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0*g
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1*g*gyfact
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1*g*gyfact
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1*g
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1*g
        enddo
!$OMP END DO
      elseif (l4symtry .and. .not. lcylindrical) then
        iinext = int((abs(xp(1)) - xmmin)*dxi)
        jjnext = int((abs(yp(1)) - ymmin)*dyi)
        kknext = int((zp(1) - zgrid - zmmin)*dzi)
        u1next = (abs(xp(1)) - xmmin)*dxi - iinext
        v1next = (abs(yp(1)) - ymmin)*dyi - jjnext
        w1next = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = iinext
          jj = jjnext
          kk = kknext
          u1 = u1next
          v1 = v1next
          w1 = w1next
          if (ip < np) then
            iinext = int((abs(xp(ip+1)) - xmmin)*dxi)
            jjnext = int((abs(yp(ip+1)) - ymmin)*dyi)
            kknext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            u1next = (abs(xp(ip+1)) - xmmin)*dxi - iinext
            v1next = (abs(yp(ip+1)) - ymmin)*dyi - jjnext
            w1next = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii == 0 .and. xmmin == 0.) gxfact = 2.
          if (jj == 0 .and. ymmin == 0.) gyfact = 2.
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0*g*gxfact*gyfact
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0*g*gyfact
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0*g*gxfact
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0*g
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1*g*gxfact*gyfact
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1*g*gyfact
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1*g*gxfact
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1*g
        enddo
!$OMP END DO
      else
c       --- normal loop
        if (lcylindrical) then
          x = sqrt(xp(1)**2 + yp(1)**2)
          y = 0.
        else
          x = xp(1)
          y = yp(1)
        endif
        iinext = int((x - xmmin)*dxi)
        jjnext = int((y - ymmin)*dyi)
        kknext = int((zp(1) - zgrid - zmmin)*dzi)
        u1next = (x - xmmin)*dxi - iinext
        v1next = (y - ymmin)*dyi - jjnext
        w1next = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          if (uzp(ip) == 0) cycle
          ii = iinext
          jj = jjnext
          kk = kknext
          u1 = u1next
          v1 = v1next
          w1 = w1next
          if (ip < np) then
            if (lcylindrical) then
              x = sqrt(xp(ip+1)**2 + yp(ip+1)**2)
              y = 0.
            else
              x = xp(ip+1)
              y = yp(ip+1)
            endif
            iinext = int((x - xmmin)*dxi)
            jjnext = int((y - ymmin)*dyi)
            kknext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            u1next = (x - xmmin)*dxi - iinext
            v1next = (y - ymmin)*dyi - jjnext
            w1next = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1
          u0 = u0*cdens(ii  )
          u1 = u1*cdens(ii+1)
          rho(ii  ,jj  ,kk  ) = rho(ii  ,jj  ,kk  ) + u0*v0*w0
          rho(ii+1,jj  ,kk  ) = rho(ii+1,jj  ,kk  ) + u1*v0*w0
          rho(ii  ,jj+1,kk  ) = rho(ii  ,jj+1,kk  ) + u0*v1*w0
          rho(ii+1,jj+1,kk  ) = rho(ii+1,jj+1,kk  ) + u1*v1*w0
          rho(ii  ,jj  ,kk+1) = rho(ii  ,jj  ,kk+1) + u0*v0*w1
          rho(ii+1,jj  ,kk+1) = rho(ii+1,jj  ,kk+1) + u1*v0*w1
          rho(ii  ,jj+1,kk+1) = rho(ii  ,jj+1,kk+1) + u0*v1*w1
          rho(ii+1,jj+1,kk+1) = rho(ii+1,jj+1,kk+1) + u1*v1*w1
        enddo
!$OMP END DO
      endif

!$OMP END PARALLEL

      return
      end
c=============================================================================
      subroutine setrho3dvector1(rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                    nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry)
      use GlobalVars
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho1d(0:(1+nx)*(1+ny)*(1+nz)-1)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c  Sets charge density
c  Same as vector except that the float to integer conversions are
c  precalculated for the next particle since the conversion can be a bottleneck.

      integer(ISZ):: moff(0:7)
      integer(ISZ),allocatable:: indx(:,:)
      real(kind=8),allocatable:: s(:,:)

      integer(ISZ):: ipmin,nptmp,ip,i,j,k,ind0,m,ir
      integer(ISZ):: inext,jnext,knext
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact

c  Set up offset array for vectorized deposition:

      moff(0) = 0
      moff(1) = 1
      moff(2) = nx+1
      moff(3) = nx+2
      moff(4) = (nx+1)*(ny+1)
      moff(5) = (nx+1)*(ny+1)+1
      moff(6) = (nx+1)*(ny+1)+nx+1
      moff(7) = (nx+1)*(ny+1)+nx+2

      g = wght*q/(dx*dy*dz)
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
!$OMP&                 gyfact,gxfact,m,inext,jnext,knext)
!$OMP&FIRSTPRIVATE(np)

      allocate(indx(0:7,0:nparpgrp-1),s(0:7,0:nparpgrp-1))

!$OMP DO
      do ipmin = 1,np,nparpgrp
        nptmp = min(nparpgrp, np+1-ipmin)

c       --- Here, the i,j,k are precalculated for the next particle
c       --- since the float to integer conversion can be expensive.
c       --- Tests however only show about a 5% speed up.

c       --- vectorized loop to compute indices, weights
        if (l2symtry) then
c         --- special loop for 2-fold symmetry
c         --- The particle weight is reduced by a factor of 2 except near the
c         --- transverse boundaries.
          inext = int((xp(ipmin) - xmmin)*dxi)
          jnext = int((abs(yp(ipmin)) - ymmin)*dyi)
          knext = int((zp(ipmin) - zgrid - zmmin)*dzi)
          do ip = ipmin,ipmin+nptmp-1
            i = inext
            j = jnext
            k = knext
            if (ip < ipmin+nptmp-1) then
              inext = int((xp(ip+1) - xmmin)*dxi)
              jnext = int((abs(yp(ip+1)) - ymmin)*dyi)
              knext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            endif
            u1 =     (xp(ip) - xmmin)*dxi - i
            u0 = 1. - u1
            v1 =     (abs(yp(ip)) - ymmin)*dyi - j
            v0 = 1. - v1
            w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
            s(0,ir) = u0*v0*w0*g*gyfact
            s(1,ir) = u1*v0*w0*g*gyfact
            s(2,ir) = u0*v1*w0*g
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g*gyfact
            s(5,ir) = u1*v0*w1*g*gyfact
            s(6,ir) = u0*v1*w1*g
            s(7,ir) = u1*v1*w1*g
          enddo
        elseif (l4symtry) then
c         --- special loop for 4-fold symmetry
c         --- The particle weight is reduced by a factor of 4 except near the
c         --- transverse boundaries.
          inext  = int((abs(xp(ipmin)) - xmmin)*dxi)
          jnext  = int((abs(yp(ipmin)) - ymmin)*dyi)
          knext  = int((zp(ipmin) - zgrid - zmmin)*dzi)
          do ip = ipmin,ipmin+nptmp-1
            i = inext
            j = jnext
            k = knext
            if (ip < ipmin+nptmp-1) then
               inext = int((abs(xp(ip+1)) - xmmin)*dxi)
               jnext = int((abs(yp(ip+1)) - ymmin)*dyi)
               knext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            endif
            u1 =     (abs(xp(ip)) - xmmin)*dxi - i
            u0 = 1. - u1
            v1 =     (abs(yp(ip)) - ymmin)*dyi - j
            v0 = 1. - v1
            w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
            s(0,ir) = u0*v0*w0*g*gxfact*gyfact
            s(1,ir) = u1*v0*w0*g*gyfact
            s(2,ir) = u0*v1*w0*g*gxfact
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g*gxfact*gyfact
            s(5,ir) = u1*v0*w1*g*gyfact
            s(6,ir) = u0*v1*w1*g*gxfact
            s(7,ir) = u1*v1*w1*g
          enddo
        else
c         --- normal loop
          inext = int((xp(ipmin) - xmmin)*dxi)
          jnext = int((yp(ipmin) - ymmin)*dyi)
          knext = int((zp(ipmin) - zgrid - zmmin)*dzi)
          do ip = ipmin,ipmin+nptmp-1
            i = inext
            j = jnext
            k = knext
            if (ip < ipmin+nptmp-1) then
              inext = int((xp(ip+1) - xmmin)*dxi)
              jnext = int((yp(ip+1) - ymmin)*dyi)
              knext = int((zp(ip+1) - zgrid - zmmin)*dzi)
            endif
            u1 =     (xp(ip) - xmmin)*dxi - i
            u0 = 1. - u1
            v1 =     (yp(ip) - ymmin)*dyi - j
            v0 = 1. - v1
            w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
            s(0,ir) = u0*v0*w0*g
            s(1,ir) = u1*v0*w0*g
            s(2,ir) = u0*v1*w0*g
            s(3,ir) = u1*v1*w0*g
            s(4,ir) = u0*v0*w1*g
            s(5,ir) = u1*v0*w1*g
            s(6,ir) = u0*v1*w1*g
            s(7,ir) = u1*v1*w1*g
          enddo
        endif

c       --- vectorized deposition over the 8 cells touched;
c       --- there'd be a hazard if we interchanged the loops.
!$OMP CRITICAL (CRITICAL_SETRHO3D1)
        do ir = 0,nptmp-1
          if (uzp(ir+ipmin) /= 0.) then
            do m = 0, 7
              rho1d(indx(m,ir)) = rho1d(indx(m,ir)) + s(m,ir)
            enddo
          endif
        enddo
!$OMP END CRITICAL (CRITICAL_SETRHO3D1)

      enddo
!$OMP END DO

      deallocate(indx,s)

!$OMP END PARALLEL

      return
      end
c=============================================================================
      subroutine setrho3ddirectspline2(rho,np,xp,yp,zp,zgrid,uzp,q,wght,
     &                   nx,ny,nz,dx,dy,dz,xmmin,ymmin,zmmin,l2symtry,l4symtry)
      use GlobalVars
      use InGen,Only: pbound0,pboundnz,pboundxy
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: zgrid,q,wght
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np)
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry

c  Sets charge density
c  Uses a second order spline with direct deposition (no particle blocks)

      integer(ISZ):: ip,ii,jj,kk
      integer(ISZ):: iinext,jjnext,kknext
      real(kind=8):: g,dxi,dyi,dzi
      real(kind=8):: wx,wy,wz
      real(kind=8):: wxnext,wynext,wznext
      real(kind=8):: u0,u1,u2,v0,v1,v2,w0,w1,w2
      real(kind=8):: gxf,gyf,gxfm1,gyfm1
      real(kind=8),allocatable:: rhotemp(:,:,:)

      if (np == 0) return

c     --- This is very memory inefficient, but make dealing with boundary
c     --- conditions much easier.
      allocate(rhotemp(-1:nx+1,-1:ny+1,-1:nz+1))

      g = wght*q/(dx*dy*dz)
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
!$OMP PARALLEL PRIVATE(ii,jj,kk,wx,wy,wz,u2,u1,u0,v2,v1,v0,w2,w1,w0,ip,
!$OMP&                 gyfact,gxfact)
!$OMP&FIRSTPRIVATE(np)

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        iinext = nint((xp(1) - xmmin)*dxi)
        jjnext = nint((abs(yp(1)) - ymmin)*dyi)
        kknext = nint((zp(1) - zgrid - zmmin)*dzi)
        wxnext = (xp(1) - xmmin)*dxi - iinext
        wynext = (abs(yp(1)) - ymmin)*dyi - jjnext
        wznext = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          ii = iinext
          jj = jjnext
          kk = kknext
          wx = wxnext
          wy = wynext
          wz = wznext
          if (ip < np) then
            iinext = nint((xp(ip+1) - xmmin)*dxi)
            jjnext = nint((abs(yp(ip+1)) - ymmin)*dyi)
            kknext = nint((zp(ip+1) - zgrid - zmmin)*dzi)
            wxnext = (xp(ip+1) - xmmin)*dxi - iinext
            wynext = (abs(yp(ip+1)) - ymmin)*dyi - jjnext
            wznext = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          if (uzp(ip) == 0.) cycle

          u0 = 0.5*(0.5 - wx)**2
          u1 = (0.75 - wx**2)
          u2 = 0.5*(0.5 + wx)**2
          v0 = 0.5*(0.5 - wy)**2
          v1 = (0.75 - wy**2)
          v2 = 0.5*(0.5 + wy)**2
          w0 = 0.5*(0.5 - wz)**2
          w1 = (0.75 - wz**2)
          w2 = 0.5*(0.5 + wz)**2

          rhotemp(ii-1,jj-1,kk-1) = rhotemp(ii-1,jj-1,kk-1) + u0*v0*w0*g
          rhotemp(ii  ,jj-1,kk-1) = rhotemp(ii  ,jj-1,kk-1) + u1*v0*w0*g
          rhotemp(ii+1,jj-1,kk-1) = rhotemp(ii+1,jj-1,kk-1) + u2*v0*w0*g
          rhotemp(ii-1,jj  ,kk-1) = rhotemp(ii-1,jj  ,kk-1) + u0*v1*w0*g
          rhotemp(ii  ,jj  ,kk-1) = rhotemp(ii  ,jj  ,kk-1) + u1*v1*w0*g
          rhotemp(ii+1,jj  ,kk-1) = rhotemp(ii+1,jj  ,kk-1) + u2*v1*w0*g
          rhotemp(ii-1,jj+1,kk-1) = rhotemp(ii-1,jj+1,kk-1) + u0*v2*w0*g
          rhotemp(ii  ,jj+1,kk-1) = rhotemp(ii  ,jj+1,kk-1) + u1*v2*w0*g
          rhotemp(ii+1,jj+1,kk-1) = rhotemp(ii+1,jj+1,kk-1) + u2*v2*w0*g

          rhotemp(ii-1,jj-1,kk  ) = rhotemp(ii-1,jj-1,kk  ) + u0*v0*w1*g
          rhotemp(ii  ,jj-1,kk  ) = rhotemp(ii  ,jj-1,kk  ) + u1*v0*w1*g
          rhotemp(ii+1,jj-1,kk  ) = rhotemp(ii+1,jj-1,kk  ) + u2*v0*w1*g
          rhotemp(ii-1,jj  ,kk  ) = rhotemp(ii-1,jj  ,kk  ) + u0*v1*w1*g
          rhotemp(ii  ,jj  ,kk  ) = rhotemp(ii  ,jj  ,kk  ) + u1*v1*w1*g
          rhotemp(ii+1,jj  ,kk  ) = rhotemp(ii+1,jj  ,kk  ) + u2*v1*w1*g
          rhotemp(ii-1,jj+1,kk  ) = rhotemp(ii-1,jj+1,kk  ) + u0*v2*w1*g
          rhotemp(ii  ,jj+1,kk  ) = rhotemp(ii  ,jj+1,kk  ) + u1*v2*w1*g
          rhotemp(ii+1,jj+1,kk  ) = rhotemp(ii+1,jj+1,kk  ) + u2*v2*w1*g

          rhotemp(ii-1,jj-1,kk+1) = rhotemp(ii-1,jj-1,kk+1) + u0*v0*w2*g
          rhotemp(ii  ,jj-1,kk+1) = rhotemp(ii  ,jj-1,kk+1) + u1*v0*w2*g
          rhotemp(ii+1,jj-1,kk+1) = rhotemp(ii+1,jj-1,kk+1) + u2*v0*w2*g
          rhotemp(ii-1,jj  ,kk+1) = rhotemp(ii-1,jj  ,kk+1) + u0*v1*w2*g
          rhotemp(ii  ,jj  ,kk+1) = rhotemp(ii  ,jj  ,kk+1) + u1*v1*w2*g
          rhotemp(ii+1,jj  ,kk+1) = rhotemp(ii+1,jj  ,kk+1) + u2*v1*w2*g
          rhotemp(ii-1,jj+1,kk+1) = rhotemp(ii-1,jj+1,kk+1) + u0*v2*w2*g
          rhotemp(ii  ,jj+1,kk+1) = rhotemp(ii  ,jj+1,kk+1) + u1*v2*w2*g
          rhotemp(ii+1,jj+1,kk+1) = rhotemp(ii+1,jj+1,kk+1) + u2*v2*w2*g

        enddo
!$OMP END DO
      elseif (l4symtry) then
        iinext = nint((abs(xp(1)) - xmmin)*dxi)
        jjnext = nint((abs(yp(1)) - ymmin)*dyi)
        kknext = nint((zp(1) - zgrid - zmmin)*dzi)
        wxnext = (abs(xp(1)) - xmmin)*dxi - iinext
        wynext = (abs(yp(1)) - ymmin)*dyi - jjnext
        wznext = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          ii = iinext
          jj = jjnext
          kk = kknext
          wx = wxnext
          wy = wynext
          wz = wznext
          if (ip < np) then
            iinext = nint((abs(xp(ip+1)) - xmmin)*dxi)
            jjnext = nint((abs(yp(ip+1)) - ymmin)*dyi)
            kknext = nint((zp(ip+1) - zgrid - zmmin)*dzi)
            wxnext = (abs(xp(ip+1)) - xmmin)*dxi - iinext
            wynext = (abs(yp(ip+1)) - ymmin)*dyi - jjnext
            wznext = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          if (uzp(ip) == 0.) cycle

          u0 = 0.5*(0.5 - wx)**2
          u1 = (0.75 - wx**2)
          u2 = 0.5*(0.5 + wx)**2
          v0 = 0.5*(0.5 - wy)**2
          v1 = (0.75 - wy**2)
          v2 = 0.5*(0.5 + wy)**2
          w0 = 0.5*(0.5 - wz)**2
          w1 = (0.75 - wz**2)
          w2 = 0.5*(0.5 + wz)**2

          rhotemp(ii-1,jj-1,kk-1) = rhotemp(ii-1,jj-1,kk-1) + u0*v0*w0*g
          rhotemp(ii  ,jj-1,kk-1) = rhotemp(ii  ,jj-1,kk-1) + u1*v0*w0*g
          rhotemp(ii+1,jj-1,kk-1) = rhotemp(ii+1,jj-1,kk-1) + u2*v0*w0*g
          rhotemp(ii-1,jj  ,kk-1) = rhotemp(ii-1,jj  ,kk-1) + u0*v1*w0*g
          rhotemp(ii  ,jj  ,kk-1) = rhotemp(ii  ,jj  ,kk-1) + u1*v1*w0*g
          rhotemp(ii+1,jj  ,kk-1) = rhotemp(ii+1,jj  ,kk-1) + u2*v1*w0*g
          rhotemp(ii-1,jj+1,kk-1) = rhotemp(ii-1,jj+1,kk-1) + u0*v2*w0*g
          rhotemp(ii  ,jj+1,kk-1) = rhotemp(ii  ,jj+1,kk-1) + u1*v2*w0*g
          rhotemp(ii+1,jj+1,kk-1) = rhotemp(ii+1,jj+1,kk-1) + u2*v2*w0*g

          rhotemp(ii-1,jj-1,kk  ) = rhotemp(ii-1,jj-1,kk  ) + u0*v0*w1*g
          rhotemp(ii  ,jj-1,kk  ) = rhotemp(ii  ,jj-1,kk  ) + u1*v0*w1*g
          rhotemp(ii+1,jj-1,kk  ) = rhotemp(ii+1,jj-1,kk  ) + u2*v0*w1*g
          rhotemp(ii-1,jj  ,kk  ) = rhotemp(ii-1,jj  ,kk  ) + u0*v1*w1*g
          rhotemp(ii  ,jj  ,kk  ) = rhotemp(ii  ,jj  ,kk  ) + u1*v1*w1*g
          rhotemp(ii+1,jj  ,kk  ) = rhotemp(ii+1,jj  ,kk  ) + u2*v1*w1*g
          rhotemp(ii-1,jj+1,kk  ) = rhotemp(ii-1,jj+1,kk  ) + u0*v2*w1*g
          rhotemp(ii  ,jj+1,kk  ) = rhotemp(ii  ,jj+1,kk  ) + u1*v2*w1*g
          rhotemp(ii+1,jj+1,kk  ) = rhotemp(ii+1,jj+1,kk  ) + u2*v2*w1*g

          rhotemp(ii-1,jj-1,kk+1) = rhotemp(ii-1,jj-1,kk+1) + u0*v0*w2*g
          rhotemp(ii  ,jj-1,kk+1) = rhotemp(ii  ,jj-1,kk+1) + u1*v0*w2*g
          rhotemp(ii+1,jj-1,kk+1) = rhotemp(ii+1,jj-1,kk+1) + u2*v0*w2*g
          rhotemp(ii-1,jj  ,kk+1) = rhotemp(ii-1,jj  ,kk+1) + u0*v1*w2*g
          rhotemp(ii  ,jj  ,kk+1) = rhotemp(ii  ,jj  ,kk+1) + u1*v1*w2*g
          rhotemp(ii+1,jj  ,kk+1) = rhotemp(ii+1,jj  ,kk+1) + u2*v1*w2*g
          rhotemp(ii-1,jj+1,kk+1) = rhotemp(ii-1,jj+1,kk+1) + u0*v2*w2*g
          rhotemp(ii  ,jj+1,kk+1) = rhotemp(ii  ,jj+1,kk+1) + u1*v2*w2*g
          rhotemp(ii+1,jj+1,kk+1) = rhotemp(ii+1,jj+1,kk+1) + u2*v2*w2*g

        enddo
!$OMP END DO
      else
c       --- normal loop
        iinext = nint((xp(1) - xmmin)*dxi)
        jjnext = nint((yp(1) - ymmin)*dyi)
        kknext = nint((zp(1) - zgrid - zmmin)*dzi)
        wxnext = (xp(1) - xmmin)*dxi - iinext
        wynext = (yp(1) - ymmin)*dyi - jjnext
        wznext = (zp(1) - zgrid - zmmin)*dzi - kknext
!$OMP DO
        do ip = 1,np
          ii = iinext
          jj = jjnext
          kk = kknext
          wx = wxnext
          wy = wynext
          wz = wznext
          if (ip < np) then
            iinext = nint((xp(ip+1) - xmmin)*dxi)
            jjnext = nint((yp(ip+1) - ymmin)*dyi)
            kknext = nint((zp(ip+1) - zgrid - zmmin)*dzi)
            wxnext = (xp(ip+1) - xmmin)*dxi - iinext
            wynext = (yp(ip+1) - ymmin)*dyi - jjnext
            wznext = (zp(ip+1) - zgrid - zmmin)*dzi - kknext
          endif
          if (uzp(ip) == 0.) cycle

          u0 = 0.5*(0.5 - wx)**2
          u1 = (0.75 - wx**2)
          u2 = 0.5*(0.5 + wx)**2
          v0 = 0.5*(0.5 - wy)**2
          v1 = (0.75 - wy**2)
          v2 = 0.5*(0.5 + wy)**2
          w0 = 0.5*(0.5 - wz)**2
          w1 = (0.75 - wz**2)
          w2 = 0.5*(0.5 + wz)**2

          rhotemp(ii-1,jj-1,kk-1) = rhotemp(ii-1,jj-1,kk-1) + u0*v0*w0*g
          rhotemp(ii  ,jj-1,kk-1) = rhotemp(ii  ,jj-1,kk-1) + u1*v0*w0*g
          rhotemp(ii+1,jj-1,kk-1) = rhotemp(ii+1,jj-1,kk-1) + u2*v0*w0*g
          rhotemp(ii-1,jj  ,kk-1) = rhotemp(ii-1,jj  ,kk-1) + u0*v1*w0*g
          rhotemp(ii  ,jj  ,kk-1) = rhotemp(ii  ,jj  ,kk-1) + u1*v1*w0*g
          rhotemp(ii+1,jj  ,kk-1) = rhotemp(ii+1,jj  ,kk-1) + u2*v1*w0*g
          rhotemp(ii-1,jj+1,kk-1) = rhotemp(ii-1,jj+1,kk-1) + u0*v2*w0*g
          rhotemp(ii  ,jj+1,kk-1) = rhotemp(ii  ,jj+1,kk-1) + u1*v2*w0*g
          rhotemp(ii+1,jj+1,kk-1) = rhotemp(ii+1,jj+1,kk-1) + u2*v2*w0*g

          rhotemp(ii-1,jj-1,kk  ) = rhotemp(ii-1,jj-1,kk  ) + u0*v0*w1*g
          rhotemp(ii  ,jj-1,kk  ) = rhotemp(ii  ,jj-1,kk  ) + u1*v0*w1*g
          rhotemp(ii+1,jj-1,kk  ) = rhotemp(ii+1,jj-1,kk  ) + u2*v0*w1*g
          rhotemp(ii-1,jj  ,kk  ) = rhotemp(ii-1,jj  ,kk  ) + u0*v1*w1*g
          rhotemp(ii  ,jj  ,kk  ) = rhotemp(ii  ,jj  ,kk  ) + u1*v1*w1*g
          rhotemp(ii+1,jj  ,kk  ) = rhotemp(ii+1,jj  ,kk  ) + u2*v1*w1*g
          rhotemp(ii-1,jj+1,kk  ) = rhotemp(ii-1,jj+1,kk  ) + u0*v2*w1*g
          rhotemp(ii  ,jj+1,kk  ) = rhotemp(ii  ,jj+1,kk  ) + u1*v2*w1*g
          rhotemp(ii+1,jj+1,kk  ) = rhotemp(ii+1,jj+1,kk  ) + u2*v2*w1*g

          rhotemp(ii-1,jj-1,kk+1) = rhotemp(ii-1,jj-1,kk+1) + u0*v0*w2*g
          rhotemp(ii  ,jj-1,kk+1) = rhotemp(ii  ,jj-1,kk+1) + u1*v0*w2*g
          rhotemp(ii+1,jj-1,kk+1) = rhotemp(ii+1,jj-1,kk+1) + u2*v0*w2*g
          rhotemp(ii-1,jj  ,kk+1) = rhotemp(ii-1,jj  ,kk+1) + u0*v1*w2*g
          rhotemp(ii  ,jj  ,kk+1) = rhotemp(ii  ,jj  ,kk+1) + u1*v1*w2*g
          rhotemp(ii+1,jj  ,kk+1) = rhotemp(ii+1,jj  ,kk+1) + u2*v1*w2*g
          rhotemp(ii-1,jj+1,kk+1) = rhotemp(ii-1,jj+1,kk+1) + u0*v2*w2*g
          rhotemp(ii  ,jj+1,kk+1) = rhotemp(ii  ,jj+1,kk+1) + u1*v2*w2*g
          rhotemp(ii+1,jj+1,kk+1) = rhotemp(ii+1,jj+1,kk+1) + u2*v2*w2*g

        enddo
!$OMP END DO
      endif

!$OMP END PARALLEL

      rho(0:nx,0:ny,0:nz) = rhotemp(0:nx,0:ny,0:nz)

c     --- Note that the commented out lines below are taken care of during
c     --- the call to perrho3d, but are listed here for completeness.
      if (l4symtry .or. pboundxy == neumann) then
        rho(0,:,:) = 2.*rho(0,:,:)
        rho(1,:,:) = rho(1,:,:) + rhotemp(-1,0:ny,0:nz)
        if (pboundxy == neumann) then
          rho(nx,:,:) = 2.*rho(nx,:,:)
          rho(nx-1,:,:) = rho(nx-1,:,:) + rhotemp(nx+1,0:ny,0:nz)
        endif
      else if (pboundxy == periodic) then
c       rho(0,:,:) = rho(0,:,:) + rho(nx,:,:)
c       rho(nx,:,:) = rho(0,:,:)
        rho(1,:,:) = rho(1,:,:) + rhotemp(nx+1,0:ny,0:nz)
        rho(nx-1,:,:) = rho(nx-1,:,:) + rhotemp(-1,0:ny,0:nz)
      endif

      if (l4symtry .or. l2symtry .or. pboundxy == neumann) then
        rho(:,0,:) = 2.*rho(:,0,:)
        rho(:,1,:) = rho(:,1,:) + rhotemp(0:nx,-1,0:nz)
        if (pboundxy == neumann) then
          rho(:,ny,:) = 2.*rho(:,ny,:)
          rho(:,ny-1,:) = rho(:,ny-1,:) + rhotemp(0:nx,ny+1,0:nz)
        endif
      else if (pboundxy == periodic) then
c       rho(:,0,:) = rho(:,0,:) + rho(:,ny,:)
c       rho(:,ny,:) = rho(:,0,:)
        rho(:,1,:) = rho(:,1,:) + rhotemp(0:nx,ny+1,0:nz)
        rho(:,ny-1,:) = rho(:,ny-1,:) + rhotemp(0:nx,-1,0:nz)
      endif

      if (pbound0 == periodic .or. pboundnz == periodic) then
c       rho(:,:,0) = rho(:,:,0) + rho(:,:,nz)
c       rho(:,:,nz) = rho(:,:,0)
        rho(:,:,1) = rho(:,:,1) + rhotemp(0:nx,0:ny,nz+1)
        rho(:,:,nz-1) = rho(:,:,nz-1) + rhotemp(0:nx,0:ny,-1)
      else
        if (pbound0 == neumann) then
          rho(:,:,0) = 2.*rho(:,:,0)
          rho(:,:,1) = rho(:,:,1) + rhotemp(0:nx,0:ny,-1)
        endif
        if (pboundnz == neumann) then
          rho(:,:,nz) = 2.*rho(:,:,nz)
          rho(:,:,nz-1) = rho(:,:,nz-1) + rhotemp(0:nx,0:ny,nz+1)
        endif
      endif

      deallocate(rhotemp)

      return
      end
c=============================================================================
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
c     rho(i  ,j  ,k  ) = rho(i  ,j  ,k  ) + u0*v0*w0*g
c     rho(i+1,j  ,k  ) = rho(i+1,j  ,k  ) + u1*v0*w0*g
c  Note that many changes are possible; for example, we might define
c  ind0(ir) and not use indx; this saves some store operations but
c  leads to a more complicated indirect address for the vectorized
c  gather-add-scatter loop.  It seems about 3% slower than the present way.
c  RHO must be zeroed in ZERORHO since it is not zeroed here (to allow
c  handling of blocks of particles at a time)

c     --- For vectorized algorithm
      integer(ISZ):: moff(0:7)
      integer(ISZ),allocatable:: indx(:,:)
      save moff
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ),allocatable:: ii(:), jj(:), kk(:)
c     --- For both algorithms
      real(kind=8),allocatable:: s(:,:)

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

      g = wght*q / (dx*dy*dz)
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

      if (depos == "vector") then 
        allocate(indx(0:7,0:nparpgrp-1),s(0:7,0:nparpgrp-1))
      elseif (depos == "scalar") then
        allocate(ii(0:nparpgrp-1), jj(0:nparpgrp-1), kk(0:nparpgrp-1),
     &           s(0:7,0:nparpgrp-1))
      endif

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
           i  = int((xp(ip) - xmmin)*dxi)
           u1 =     (xp(ip) - xmmin)*dxi - i
           u0 = 1. - u1
           j  = int((abs(yp(ip)) - ymmin)*dyi)
           v1 =     (abs(yp(ip)) - ymmin)*dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin)*dzi)
           w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
           s(0,ir) = u0*v0*w0*g*gyfact
           s(1,ir) = u1*v0*w0*g*gyfact
           s(2,ir) = u0*v1*w0*g
           s(3,ir) = u1*v1*w0*g
           s(4,ir) = u0*v0*w1*g*gyfact
           s(5,ir) = u1*v0*w1*g*gyfact
           s(6,ir) = u0*v1*w1*g
           s(7,ir) = u1*v1*w1*g
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
           k  = int((zp(ip) - zgrid - zmmin)*dzi)
           w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
           s(0,ir) = u0*v0*w0*g*gxfact*gyfact
           s(1,ir) = u1*v0*w0*g*gyfact
           s(2,ir) = u0*v1*w0*g*gxfact
           s(3,ir) = u1*v1*w0*g
           s(4,ir) = u0*v0*w1*g*gxfact*gyfact
           s(5,ir) = u1*v0*w1*g*gyfact
           s(6,ir) = u0*v1*w1*g*gxfact
           s(7,ir) = u1*v1*w1*g
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
           if (uzp(ip) == 0.) cycle
           i  = int((xp(ip) - xmmin)*dxi)
           u1 =     (xp(ip) - xmmin)*dxi - i
           u0 = 1. - u1
           j  = int((yp(ip) - ymmin)*dyi)
           v1 =     (yp(ip) - ymmin)*dyi - j
           v0 = 1. - v1
           k  = int((zp(ip) - zgrid - zmmin)*dzi)
           w1 =     (zp(ip) - zgrid - zmmin)*dzi - k
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
           s(0,ir) = u0*v0*w0*g
           s(1,ir) = u1*v0*w0*g
           s(2,ir) = u0*v1*w0*g
           s(3,ir) = u1*v1*w0*g
           s(4,ir) = u0*v0*w1*g
           s(5,ir) = u1*v0*w1*g
           s(6,ir) = u0*v1*w1*g
           s(7,ir) = u1*v1*w1*g
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
          ii(ir) = int((xp(ip) - xmmin)*dxi)
          u1     =     (xp(ip) - xmmin)*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(yp(ip)) - ymmin)*dyi)
          v1     =     (abs(yp(ip)) - ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin)*dzi)
          w1     =     (zp(ip) - zgrid - zmmin)*dzi - kk(ir)
          w0     = 1. - w1
          gyfact = 1.
          if (jj(ir) == 0 .and. ymmin == 0.) gyfact = 2.
          s(0,ir) = u0*v0*w0*g*gyfact
          s(1,ir) = u1*v0*w0*g*gyfact
          s(2,ir) = u0*v1*w0*g
          s(3,ir) = u1*v1*w0*g
          s(4,ir) = u0*v0*w1*g*gyfact
          s(5,ir) = u1*v0*w1*g*gyfact
          s(6,ir) = u0*v1*w1*g
          s(7,ir) = u1*v1*w1*g
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
          kk(ir) = int((zp(ip) - zgrid - zmmin)*dzi)
          w1     =     (zp(ip) - zgrid - zmmin)*dzi - kk(ir)
          w0     = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii(ir) == 0 .and. xmmin == 0.) gxfact = 2.
          if (jj(ir) == 0 .and. ymmin == 0.) gyfact = 2.
          s(0,ir) = u0*v0*w0*g*gxfact*gyfact
          s(1,ir) = u1*v0*w0*g*gyfact
          s(2,ir) = u0*v1*w0*g*gxfact
          s(3,ir) = u1*v1*w0*g
          s(4,ir) = u0*v0*w1*g*gxfact*gyfact
          s(5,ir) = u1*v0*w1*g*gyfact
          s(6,ir) = u0*v1*w1*g*gxfact
          s(7,ir) = u1*v1*w1*g
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
          if (uzp(ip) == 0.) cycle
          ir = ip - ipmin
          ii(ir) = int((xp(ip) - xmmin)*dxi)
          u1     = (xp(ip) - xmmin)*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((yp(ip) - ymmin)*dyi)
          v1     = (yp(ip) - ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - zmmin)*dzi)
          w1     = (zp(ip) - zgrid - zmmin)*dzi - kk(ir)
          w0     = 1. - w1
          s(0,ir) = u0*v0*w0*g
          s(1,ir) = u1*v0*w0*g
          s(2,ir) = u0*v1*w0*g
          s(3,ir) = u1*v1*w0*g
          s(4,ir) = u0*v0*w1*g
          s(5,ir) = u1*v0*w1*g
          s(6,ir) = u0*v1*w1*g
          s(7,ir) = u1*v1*w1*g
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

      if (depos == "vector") then 
        deallocate(indx,s)
      elseif (depos == "scalar") then
        deallocate(ii,jj,kk,s)
      endif

!$OMP END PARALLEL

      if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine zerorhowithsampledsubcycling()
      use Subcycling
      use Fields3dParticles

      integer(ISZ):: in1

c     --- Zero the rho copy for species when the positions
c     --- are advanced.
c     --- rhopndts(:,:,:,1,in1) holds the old rho which is still needed
c     --- for the faster advanced groups.
      do in1=0,nsndts-1
        if (ldts(in1)) then
          if (nrhopndtscopies == 2) then
            rhopndts(:,:,:,1,in1) = rhopndts(:,:,:,0,in1)
          endif
          call zeroarry(rhopndts(:,:,:,0,in1),(nxp+1)*(nyp+1)*(nzp+1))
        endif
      enddo

      return
      end
c=============================================================================
      subroutine averagerhowithsampledsubcycling()
      use Picglb,Only: it
      use Subcycling
      use Fields3dParticles

      integer(ISZ):: in1,in2,in2tmp
      real(kind=8):: wranf,oldornew

      if (ndtsmax == 1) return

c     --- During the generate, do the copy of the new rhop to the old rhop
c     --- for group 0, which is normally done during the zerorhop call.
      if (it == 0) rhopndts(:,:,:,1,0) = rhopndts(:,:,:,0,0)

c     --- Save the rho where the fastest particle's rho is. For now,
c     --- assume that this is in1=0
      do in1=1,nsndts-1
        if (it == 0) then
c         --- At it==0, before the first step, always add the new rho.
          rhopndts(:,:,:,1,0) = rhopndts(:,:,:,1,0) + rhopndts(:,:,:,0,in1)
        elseif (mod(rhotondts(in1),2) == 1) then
c         --- Use the rho that is closest in time to the current time.
          if (mod(it-1,rhotondts(in1)) > rhotondts(in1)/2. .or.
     &        mod(it-1,rhotondts(in1)) == 0) then
            rhopndts(:,:,:,1,0) = rhopndts(:,:,:,1,0) + rhopndts(:,:,:,0,in1)
          else
            rhopndts(:,:,:,1,0) = rhopndts(:,:,:,1,0) + rhopndts(:,:,:,1,in1)
          endif
        else
c         --- When ndts is even, at the mid point of the step, take the
c         --- average of the old and the new
c         --- Otherwise, use the rho that is closest in time to the currentc
c         --- time.
          if (mod(it-1,rhotondts(in1)) == rhotondts(in1)/2) then
            rhopndts(:,:,:,1,0) = rhopndts(:,:,:,1,0) +
     &       0.5*(rhopndts(:,:,:,0,in1) + rhopndts(:,:,:,1,in1))
          else if (mod(it-1,rhotondts(in1)) > rhotondts(in1)/2. .or.
     &             mod(it-1,rhotondts(in1)) == 0) then
            rhopndts(:,:,:,1,0) = rhopndts(:,:,:,1,0) + rhopndts(:,:,:,0,in1)
          else
            rhopndts(:,:,:,1,0) = rhopndts(:,:,:,1,0) + rhopndts(:,:,:,1,in1)
          endif
        endif
      enddo

      return
      end
c=============================================================================
      subroutine zerorhowithfullvsubcycling()
      use Subcycling
      use Fields3dParticles

      integer(ISZ):: in1

c     --- Zero the rho copy for species when the positions
c     --- are advanced.
      do in1=0,nsndts-1
        if (mod(rhotondts(in1),2) == 1) then
c         --- rho is zeroed for all of the odd cases
          if (ldts(in1)) then
            call zeroarry(rhopndts(:,:,:,0,in1),(nxp+1)*(nyp+1)*(nzp+1))
          endif
        else
          print*,"This version of subcyling is not supported with even step sizes."
          stop
c         --- rho is zeroed for even cases only when there is no group
c         --- with half its time step size (since those cases already
c         --- have a contribution added.
c         if (ldts(in1)) then
c           lzeroeven = .true.
c           do in2 = 0,nsndts-1
c             if (rhotondts(in1) > rhotondts(in2)) then
c               if (mod(rhotondts(in1),rhotondts(in2)) == 0) then
c                 lzeroeven = .false.
c                 exit
c               endif
c             endif
c           enddo
c           if (lzeroeven) rhopndts(:,:,:,0,in1) = 0.
c         endif
        endif
      enddo

      return
      end
c=============================================================================
      subroutine zerorhowithselfb()
      use SelfB
      use Fields3dParticles

      if (nsselfb3d > 0) rhopselfb = 0.

      return
      end
c=============================================================================
      subroutine averagerhowithfullvsubcycling()
      use Subcycling
      use Fields3dParticles
      use Picglb,Only: it

      integer(ISZ):: in1,in2,ndtstmp,in2tmp
      real(kind=8):: rndts

      if (ndtsmax == 1) return

c     --- The copy is only needed if there are multiple groups of ndts.
      if (nrhopndtscopies == 2) then
        do in1 = 0,nsndts-1
          if (ldts(in1)) then
c           --- Copy the individual speices rho to its array holding the
c           --- total rho.
            if (mod(rhotondts(in1),2) == 1) then
              rhopndts(:,:,:,1,in1) = rhopndts(:,:,:,0,in1)
            else
c             --- For ndts even, some contribution may already have been made
c             --- in setupevensubcyclingrho.
              rhopndts(:,:,:,1,in1) = rhopndts(:,:,:,1,in1) +
     &                                rhopndts(:,:,:,0,in1)
            endif
          endif
        enddo
      endif
      do in1 = 0,nsndts-1
        do in2 = 0,nsndts-1
          if (in1 == in2) cycle
          ndtstmp = rhotondts(in1) - mod(rhotondts(in1),2)
          if (rhotondts(in1) < rhotondts(in2) .and. ldts(in1)) then
c           --- For groups with larger time steps, add in their rho
c           --- now since it only needs to be done once.
c           --- Also, if only sampling instead of averaging.
            rhopndts(:,:,:,1,in1) = rhopndts(:,:,:,1,in1) +
     &                              rhopndts(:,:,:,0,in2)
          else if (rhotondts(in1) > rhotondts(in2)) then
c           --- For groups with smaller times steps, add in the
c           --- fractional contribution.
            rndts = 1./rhotondts(in1)
c           --- Check the special case when the rho for one group is
c           --- exactly half way between two of the larger time steps.
c           --- This can happen when the larger is an integer multiple
c           --- of the smaller.
c           --- Only half is added in ala the trapezoidal rule.
            in2tmp = rhotondts(in2) - mod(rhotondts(in2),2)
            if (mod(rhotondts(in1),2) == 0 .and. ldts(in2) .and.
     &          rhotondts(in1) == 2*mod(itndts(in2),rhotondts(in1))) then
              rndts = rndts*0.5
            endif
c           --- The first time step is special since the past histories
c           --- of the faster moving groups is not known.
c           --- It is assumed that the starting time step is it=0,
c           --- which means that it = 0 now. As an estimate, the initial
c           --- charge density is multiplied by a scale factor to replace
c           --- the unknown values from the past.
            if (it == 0) then
              rndts = ((rhotondts(in1) + 1)/2.)/rhotondts(in1)
            endif
c           --- Note that time could be saved since rhotondts(in1) > 1
c           --- and its rho won't change every step and its contribution
c           --- doesn't needed to be added in every step. Right now
c           --- though I'm too lazy to figure out the logic and the
c           --- weighting factors.
            rhopndts(:,:,:,1,in1) = rhopndts(:,:,:,1,in1) +
     &                              rhopndts(:,:,:,0,in2)*rndts
          endif
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine zerorhowithhalfvsubcycling()
      use Subcycling
      use Fields3dParticles

      if (ndtsmax == 1) return

      print*,"Version 2 of averaging for subcyling is not yet supported"
      stop

      return
      end
c=============================================================================
      subroutine averagerhowithhalfvsubcycling()
      use Subcycling
      use Fields3dParticles

      return
      end
c=============================================================================
      subroutine loadrho3d(pgroup,ins_i,nps_i,is_i,lzero)
      use ParticleGroupmodule
      use GlobalVars
      use Subtimers3d
      use InGen
      use InGen3d
      use InPart
      use InMesh3d
      use Picglb
      use Picglb3d
      use Particles,Only: wpid
      use Fields3d
      use Subcycling
      use Fields3dParticles
      use GridBoundary3d
      use FieldSolveAPI, Only: lzerorhofsapi
      type(ParticleGroup):: pgroup
      integer(ISZ):: ins_i,nps_i,is_i
      logical(ISZ):: lzero

c --- This routine provides a simple call from the interpreter to load the
c --- rhop array.  The value '-1' is used as a flag in the input to use
c --- all of the particles, otherwise the specified particles are loaded.

      integer(ISZ):: ins_u,nps_u
      integer(ISZ):: is1,is2,isid
      integer(ISZ):: ip,ipmin,is,indts
      real(kind=8):: swtmp,zgrdndts
      real(kind=8):: substarttime,wtime
      integer(ISZ):: allocerror
      logical(ISZ):: lzeroeven

      if (lw3dtimesubs) substarttime = wtime()

#ifdef MPIPARALLEL
c     --- Set the flag signalling that rho has been changed and that
c     --- the parallel sums on boundaries will need to be done again.
      lrhosummedonboundaries = .false.
#endif

      if (lbeforelr) call execuser("beforelr")

c     --- Ensure that subcycling is setup correctly (in case the user had
c     --- made changes.
      call setupSubcycling(pgroup)
      call setupSelfB(pgroup)

      if(depos == 'none') return

      if (fstype == 12) then
c       --- Use user registered solver to load rho
        lzerorhofsapi = lzero
        call execuser("loadrhoregistered")
      else

c       --- zero rhop if requested
        if (lzero) then
          if (solvergeom == XYZgeom) then
            call setupFields3dParticles()
c           --- Zero the rho copy for species when the positions
c           --- are advanced.
            if (ndtsaveraging == 0) then
              call zerorhowithsampledsubcycling()
            elseif (ndtsaveraging == 1) then
              call zerorhowithfullvsubcycling()
            elseif (ndtsaveraging == 2) then
              call zerorhowithhalfvsubcycling()
            endif
            call zerorhowithselfb()
          endif
          if(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==XYgeom .or.
     &       solvergeom==Rgeom  .or. solvergeom==Zgeom) call reset_rzmgrid_rho()
        end if

c       --- change AMR grid if necessary
        if(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==XYgeom .or.
     &     solvergeom==Rgeom  .or. solvergeom==Zgeom) call change_loc_part()

c       --- set limits on loop over species
        if (is_i == -1) then
          is1 = 1
          is2 = pgroup%ns
        else
          is1 = is_i
          is2 = is_i
        endif

c       --- set initial limits from input
c       --- (will be changed if necessary in the loop)
        ins_u = ins_i
        nps_u = nps_i

c       --- loop over species
        do is=is1,is2
           isid = pgroup%sid(is-1) + 1
           if (isid == 0) cycle

c          --- For the AMR/Chombo version, all of the lost particles must be
c          --- cleared out.
           if (solvergeom == AMRgeom) call clearpart(pgroup,is,1)

c          --- get loop limits for particles if needed
           if (ins_i == -1) ins_u = pgroup%ins(is)
           if (nps_i == -1) nps_u = pgroup%nps(is)
           if (nps_u == 0) cycle

c          --- Scale the weight, sw, by the time step scale size. This only
c          --- makes sense for steady-state and slice modes. In time-dependent
c          --- mode, it is assumed that dtscale has not been changed from 1.
           swtmp = pgroup%sw(is)*pgroup%dtscale(is)
  
           if(solvergeom==XYZgeom) then
             ipmin = ins_u
             if (pgroup%ldts(is)) then
c              --- Only load the rho if the positions have been advanced.
               indts = ndtstorho(pgroup%ndts(is-1))
               zgrdndts = zgrid + (zgrid-zgridprv)*(pgroup%ndts(is-1)-1)
               call setrho3d(rhopndts(:,:,:,0,indts),
     &                       rhopndts(:,:,:,0,indts),nps_u,
     &                    pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                    zgrdndts,pgroup%uzp(ipmin),pgroup%sq(is),swtmp,depos,
     &                    nxp,nyp,nzp,dx,dy,dz,xmmin,ymmin,zmminp,
     &                    l2symtry,l4symtry,solvergeom==RZgeom)
             endif
           elseif(solvergeom==RZgeom .or. solvergeom==XZgeom) then
c            --- loop over particle blocks
             do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
               ip = min(nparpgrp, ins_u+nps_u-ipmin)
               if(wpid==0) then
                 call rhoweightrz(pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                            pgroup%zp(ipmin),ip,
     &                            pgroup%sq(is)*swtmp,nx,nz,dx,dz,xmmin,zgrid)
               else
                 call rhoweightrz_weights(pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                            pgroup%zp(ipmin),pgroup%pid(ipmin,wpid),ip,
     &                            pgroup%sq(is)*swtmp,nx,nz,dx,dz,xmmin,zgrid)
               end if
             enddo
           elseif(solvergeom==XYgeom) then
c            --- loop over particle blocks
             do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
               ip = min(nparpgrp, ins_u+nps_u-ipmin)
               if(wpid==0) then
                 call rhoweightrz(pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                            pgroup%yp(ipmin),ip,
     &                            pgroup%sq(is)*swtmp,nx,ny,dx,dy,xmmin,ymmin)
               else
                 call rhoweightrz_weights(pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                            pgroup%yp(ipmin),pgroup%pid(ipmin,wpid),ip,
     &                            pgroup%sq(is)*swtmp,nx,ny,dx,dy,xmmin,ymmin)
               end if
             enddo
           elseif(solvergeom==Zgeom) then
c            --- loop over particle blocks
             do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
                ip = min(nparpgrp, ins_u+nps_u-ipmin)
                call rhoweightz(pgroup%zp(ipmin),ip,pgroup%sq(is)*swtmp,
     &                          nz,dz,zgrid)
             enddo
           elseif(solvergeom==Rgeom) then
c            --- loop over particle blocks
             do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
                ip = min(nparpgrp, ins_u+nps_u-ipmin)
                if(wpid==0) then
                  call rhoweightr(pgroup%xp(ipmin),pgroup%yp(ipmin),ip,
     &                            pgroup%sq(is)*swtmp,nx,dx,xmmin)
                else
                  call rhoweightr_weights(pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                                    pgroup%pid(ipmin,wpid),ip,
     &                                    pgroup%sq(is)*swtmp,nx,dx,xmmin)
                end if
             enddo
           elseif(solvergeom==AMRgeom) then
c            --- loop over particle blocks
             do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
                ip = min(nparpgrp, ins_u+nps_u-ipmin)
                call cho_setrho3d(ip,pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                            pgroup%zp(ipmin),zgrid,
     &                            pgroup%sq(is),swtmp,isid,
     &                            (ipmin-pgroup%ins(is))/nparpgrp)
             enddo
           endif
        enddo

        if (lzero) then
c         --- Add together the rhop from the species with different time steps
c         --- This must be done after all of the species has been deposited.
          if(solvergeom==XYZgeom) then
            if (ndtsaveraging == 0) then
              call averagerhowithsampledsubcycling()
            elseif (ndtsaveraging == 1) then
              call averagerhowithfullvsubcycling()
            elseif (ndtsaveraging == 2) then
              call averagerhowithhalfvsubcycling()
            endif
          endif
        endif

      endif

c     --- Do any loading that is needed for injection
      call loadrho_inject(pgroup)

!$OMP MASTER
      if (lw3dtimesubs) timeloadrho3d = timeloadrho3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setupevensubcyclingrho(it)
c     use Subcycling
c     use Fields3dParticles
c     use InGen3d,Only:solvergeom,XYZgeom
      integer(ISZ):: it

c     --- When, sometimes rho is know at a time that is exactly half way
c     --- between two larger times steps. When this happens, the trapezoidal
c     --- rule is used and one half of that rho is added to the larger time
c     --- steps rho.

c     integer(ISZ):: iz,in1,in2,ndtstmp
c     real(kind=8):: rndts

c     --- This code is not supported

c     if (solvergeom/=XYZgeom) return

c     --- First, zero out the rho or those cases
c     do in1 = 0,nsndts-1
c       do in2 = 0,nsndts-1
c         if (in1 == in2) cycle
c         if (rhotondts(in1) > rhotondts(in2)) then
c           if (mod(rhotondts(in1),rhotondts(in2)) == 0) then
c             ndtstmp = rhotondts(in1) - mod(rhotondts(in1),2)
c             if (mod(it-1,rhotondts(in1)) == ndtstmp/2) then
c               rhopndts(:,:,:,1,in1) = 0.
c               exit
c             endif
c           endif
c         endif
c       enddo
c     enddo

c     --- Then add in the half contribution
c     do in1 = 0,nsndts-1
c       do in2 = 0,nsndts-1
c         if (in1 == in2) cycle
c         if (rhotondts(in1) > rhotondts(in2)) then
c           if (mod(rhotondts(in1),rhotondts(in2)) == 0) then
c             ndtstmp = rhotondts(in1) - mod(rhotondts(in1),2)
c             if (mod(it-1,rhotondts(in1)) == ndtstmp/2) then
c               rndts = 0.5/rhotondts(in1)
c               rhopndts(:,:,:,1,in1) = rhopndts(:,:,:,1,in1) +
c    &                                  rhopndts(:,:,:,0,in2)*rndts
c             endif
c           endif
c         endif
c       enddo
c     enddo

      return
      end
c=============================================================================
      subroutine assignrhoandphi(rhopindts,phipindts)
      use InGen3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles
      real(kind=8),target:: rhopindts(0:nxp,0:nyp,0:nzp)
      real(kind=8),target:: phipindts(0:nxp,0:nyp,-1:nzp+1)

c     --- The assignments must be done in a subroutine like this so that
c     --- the lower and upper bounds of rho and phi are correct.

#ifndef MPIPARALLEL

      if (solvergeom==XYZgeom) then
c        --- Point rho to the block for the currently active ndts group.
         rho => rhopindts
         phi => phipindts
      endif

#else

#endif

      return
      end
c=============================================================================
      subroutine assignrhopandphip(rhopindts,phipindts)
      use InGen3d
      use Fields3d
      use Fields3dParticles
      real(kind=8),target:: rhopindts(0:nxp,0:nyp,0:nzp)
      real(kind=8),target:: phipindts(0:nxp,0:nyp,-1:nzp+1)

c     --- The assignments must be done in a subroutine like this so that
c     --- the lower and upper bounds of rho and phi are correct.

      if (solvergeom==XYZgeom) then
c        --- Point rho to the block for the currently active ndts group.
c        --- Note that nrhopndtscopies is used in case ndts is only 1
c        --- (and nrhopndtscopies is 1) since the extra rho copy is not saved.
         rhop => rhopindts
         phip => phipindts
      endif

      return
      end
c=============================================================================
      subroutine setrhoandphiforfieldsolve(rhopin,phipin)
      use InGen
      use InGen3d
      use InMesh3d
      use Picglb3d, only: dz
      use Fields3d
      use Subcycling
      use Fields3dParticles
#ifdef MPIPARALLEL
      use Parallel
#endif
      real(kind=8),target:: rhopin(0:nxp,0:nyp,0:nzp)
      real(kind=8),target:: phipin(0:nxp,0:nyp,-1:nzp+1)

#ifndef MPIPARALLEL

c      --- Point rho to the block for the currently active ndts group.
       rho => rhopin
       phi => phipin

#else

c     --- The rho and phi arrays need to be allocated if they havn't already
c     --- been.
      if (.not. associated(rho)) allocate(rho(0:nx,0:ny,0:nz))
      if (.not. associated(phi)) allocate(phi(0:nx,0:ny,-1:nz+izextra))

c  For parallel version, each processor sends rho to neighboring processors
c  whose field solve region overlap its particle region.
      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        call getrhoforfieldsolverz(nx,nz,rho)
      else if(solvergeom==XYgeom) then
        call remark("XYgeom not tested in parallel")
        stop
      else if(solvergeom==Zgeom) then
        call getrhoforfieldsolvez(nz,rho)
      elseif (solvergeom==XYZgeom) then
        call setrhoforfieldsolve3d(nx,ny,nz,rho,nxp,nyp,nzp,rhopin,nzpguard,
     &                             my_index,nslaves,
     &                             izpslave,nzpslave,izfsslave,nzfsslave)
      end if

#endif

      return
      end
c=============================================================================
      subroutine setrhoforfieldsolve3d(nx,ny,nz,rho,nxp,nyp,nzp,rhop,nzpguard,
     &                                 my_index,nslaves,izpslave,nzpslave,
     &                                 izfsslave,nzfsslave)
      use Fields3dParticles,Only: lrhosummedonboundaries
      integer(ISZ):: nx,ny,nz,nxp,nyp,nzp,nzpguard
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

#ifdef MPIPARALLEL

      call parallelbarrier()

c     --- If rho has already been summed on the boundaries since the last
c     --- call to loadrho, then don't do it again.
      if (.not. lrhosummedonboundaries) then
        call sumsourcepondomainboundaries(rhop,1,nxp,nyp,nzp,nzpguard,
     &                                    my_index,nslaves,izpslave,nzpslave)
      endif
      call parallelbarrier()
      call setsourceforfieldsolve3d_parallel(1,nx,ny,nz,rho,nxp,nyp,nzp,rhop,
     &                                       my_index,nslaves,
     &                                       izpslave,nzpslave,izfsslave,nzfsslave)

#endif

      return
      end
c=============================================================================
      subroutine setupFields3dParticles()
      use InGen,Only: fstype
      use Picglb3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles
      use Subcycling,Only:nrhopndtscopies,nsndts,nsndtsphi
      use SelfB,Only: nsselfb
      use w3d_interfaces
#ifdef MPIPARALLEL
      use Parallel
#endif

c Ensures that the phip and rhop arrays are setup properly.

      if (fstype == 12) return

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

      nrhopndtscopies3d = nrhopndtscopies
      nsndts3d = nsndts
      nsndtsphi3d = nsndtsphi
      nsselfb3d = nsselfb
      call gchange("Fields3dParticles",0)

c     --- If either rhop of phip are not associated or either does not have
c     --- the correct dimensions as given by nxp, etc, something 
c     --- needs to be done for fix them.
c     if (.not. associated(rhop) .or. .not. associated(phip) .or.
c    &    ANY(ubound(phip) .ne. (/ nxp,nyp,nzp+1 /) ) .or.
c    &    ANY(ubound(rhop) .ne. (/ nxp,nyp,nzp   /) )) then

c       --- Free whatever the rhop and phip arrays are now refering to.
c       call gfree('Fields3dParticles')

c       if (ALL( (/ nxp,nyp,nzp /) == (/ nx,ny,nzfull /) )) then
c         --- If the sizes of rhop and rho, and phip and phi, are the same,
c         --- then associate rhop with rho and phip with phi
c         rhop => rho
c         phip => phi
c       else
c         --- If the sizes are not the same, then allocate the group
c         call gallot('Fields3dParticles',0)
c       endif
c     endif

c     --- rhop and phip have the special values for the species which are
c     --- incremented every time step.
      call assignrhopandphip(rhopndts(:,:,:,nrhopndtscopies-1,0),
     &                       phipndts(:,:,:,0))

c     --- Note that special coding is needed to clean up the rhopndts array
c     --- for cases when there are no longer certain values of ndts. This
c     --- code will need the full list of particle groups so it can check
c     --- all values of ndts to find any that are no longer used. This could
c     --- be fairly important to conserve memory since having extra copies of
c     --- rhop is expensive.

      return
      end
c=============================================================================
      subroutine fetche3d(pgroup,ipmin,ip,is)
      use ParticleGroupmodule
      use Particles, toppgroup => pgroup
      use GlobalVars
      use Subtimers3d
      use Picglb
      use Picglb3d
      use InGen
      use InGen3d
      use InMesh3d
      use Fields3d
      use Subcycling
      use Fields3dParticles
      use Efields3d
      use FieldSolveAPI
      type(ParticleGroup):: pgroup
      integer(ISZ):: ipmin,ip,is

c     --- Obtain the self-field from the electrostatic potential

      real(kind=8),pointer:: ex(:),ey(:),ez(:)
      integer(ISZ):: indts,tmpnsndts,getnsndtsforsubcycling
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Don't do anything if there are no particles input.
      if (ip == 0) return

      ex => pgroup%ex(ipmin:ipmin+ip-1)
      ey => pgroup%ey(ipmin:ipmin+ip-1)
      ez => pgroup%ez(ipmin:ipmin+ip-1)

      if(fstype == 12) then
        pgroupfsapi = pgroup
        jsfsapi = is - 1
        ipminfsapi = ipmin
        npfsapi = ip
        call execuser("fetcheregistered")
        jsfsapi = -1
        npfsapi = 0
      elseif(solvergeom==XYZgeom) then
        tmpnsndts = getnsndtsforsubcycling()
        indts = min(tmpnsndts-1,ndtstorho(pgroup%ndts(is-1)))
        call sete3d(phipndts(:,:,:,indts),selfe,ip,
     &              pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &              zgridprv,xmmin,ymmin,zmminp,dx,dy,dz,nxp,nyp,nzp,efetch,
     &              ex,ey,ez,l2symtry,l4symtry,solvergeom==RZgeom)
      elseif(solvergeom==RZgeom) then
        call setemgridrz(ipmin,ip,is,ex,ey,ez,pgroup)
      elseif(solvergeom==XZgeom) then
         call fieldweightxz(pgroup%xp(ipmin),pgroup%zp(ipmin),ex,ez,ip,
     &                      zgridprv)
      elseif(solvergeom==XYgeom) then
         call fieldweightxz(pgroup%xp(ipmin),pgroup%yp(ipmin),ex,ey,ip,0.)
      elseif(solvergeom==Zgeom) then
         call fieldweightz(pgroup%zp(ipmin),ez,ip,zgridprv)
      elseif(solvergeom==Rgeom) then
         call fieldweightr(pgroup%xp(ipmin),pgroup%yp(ipmin),ex,ey,ip)
      elseif(solvergeom==AMRgeom) then
        call cho_gete3d(ip,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                  zgridprv,ex,ey,ez,is,(ipmin-pgroup%ins(is))/nparpgrp)
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
      subroutine fetche3dfrompositions(is,ndts,n,x,y,z,ex,ey,ez)
      use GlobalVars
      use Subtimers3d
      use Picglb
      use Picglb3d
      use InGen
      use InGen3d
      use InMesh3d
      use Fields3d
      use Subcycling
      use Fields3dParticles
      use Efields3d
      use FieldSolveAPI
      integer(ISZ):: is,ndts,n
      real(kind=8),target:: x(n),y(n),z(n),ex(n),ey(n),ez(n)

c     --- Obtain the self-field from the electrostatic potential
c     --- This is nearly the same as fetche except that the positions
c     --- are passed in via the argument list rather than through the
c     --- particles group.

      integer(ISZ):: indts,tmpnsndts,getnsndtsforsubcycling
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Don't do anything if there are no particles input.
      if (n == 0) return

      if(fstype == 12) then
#if defined(XLF) || defined(XLF2)
c       --- The pointer assignment breaks with the xlf compilers, so instead
c       --- new arrays are allocated and those are used as the work space.
c       --- This is less efficient since there is overhead in the allocation
c       --- and the data needs to be copied around.
        allocate(xfsapi(n),yfsapi(n),zfsapi(n))
        allocate(exfsapi(n),eyfsapi(n),ezfsapi(n))
        xfsapi = x
        yfsapi = y
        zfsapi = z
#else
        xfsapi => x
        yfsapi => y
        zfsapi => z
        exfsapi => ex
        eyfsapi => ey
        ezfsapi => ez
#endif
        npfsapi = n
        jsfsapi = is - 1
        call execuser("fetcheregistered")
        npfsapi = 0
        jsfsapi = -1
        npfsapi = 0
#if defined(XLF) || defined(XLF2)
        ex(1:n) = exfsapi
        ey(1:n) = eyfsapi
        ez(1:n) = ezfsapi
        deallocate(xfsapi,yfsapi,zfsapi)
        deallocate(exfsapi,eyfsapi,ezfsapi)
#else
        nullify(xfsapi)
        nullify(yfsapi)
        nullify(zfsapi)
#endif
        nullify(exfsapi)
        nullify(eyfsapi)
        nullify(ezfsapi)
      elseif(solvergeom==XYZgeom) then
        tmpnsndts = getnsndtsforsubcycling()
        indts = min(tmpnsndts-1,ndtstorho(ndts))
        if (indts < 0 .or. indts > nsndts-1) then
          print*,"Error: fetche3dfrompositions: An improper value on indts has been passed in"
          return
        endif
        call sete3d(phipndts(:,:,:,indts),selfe,n,
     &              x,y,z,
     &              zgridprv,xmmin,ymmin,zmminp,dx,dy,dz,nxp,nyp,nzp,
     &              efetch,ex,ey,ez,l2symtry,l4symtry,solvergeom==RZgeom)
      elseif(solvergeom==RZgeom) then
        call fieldweightrz(x,y,z,ex,ey,ez,n,zgridprv)
      elseif(solvergeom==XZgeom) then
        call fieldweightxz(x,z,ex,ez,n,zgridprv)
      elseif(solvergeom==XYgeom) then
        call fieldweightxz(x,y,ex,ey,n,0.)
      elseif(solvergeom==Zgeom) then
        call fieldweightz(z,ez,n,zgridprv)
      elseif(solvergeom==Rgeom) then
        call fieldweightr(x,y,ex,ey,n)
      elseif(solvergeom==AMRgeom) then
        call cho_gete3d(n,x,y,z,zgridprv,ex,ey,ez,is,-1)
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
      subroutine fetchb3d(pgroup,ipmin,ip,is)
      use ParticleGroupmodule
      type(ParticleGroup):: pgroup
      integer(ISZ):: ipmin,ip,is

c     --- Don't do anything if there are no particles input.
      if (ip == 0) return

      call fetchb3dfrompositions(is,ip,
     &                           pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                           pgroup%zp(ipmin),
     &                           pgroup%bx(ipmin),pgroup%by(ipmin),
     &                           pgroup%bz(ipmin))
      return
      end
c=============================================================================
      subroutine particleboundaries3d(pgroup)
      use ParticleGroupmodule
      use Subtimers3d
      use InGen
      use InGen3d
      use InPart
      use InMesh3d
      use Picglb
      use Picglb3d
      type(ParticleGroup):: pgroup

c Does everything needed to apply particle boundary conditions.

      integer(ISZ):: is,ipmin
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Apply particle boundary conditions
      call zpartbnd(pgroup,zmmax,zmmin,dz,zgrid)

c     --- absorb particles passing through semitransparent disc
      call semitransparent_disc(pgroup,dz)

      do is=1,pgroup%ns
        if (pgroup%ldts(is) .and. pgroup%nps(is) > 0) then
           ipmin = pgroup%ins(is)
           call stckxy3d(pgroup%nps(is),
     &                   pgroup%xp(ipmin),xmmax,xmmin,dx,
     &                   pgroup%yp(ipmin),ymmax,ymmin,dy,
     &                   pgroup%zp(ipmin),zmmin,dz,pgroup%uxp(ipmin),
     &                   pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     &                   pgroup%gaminv(ipmin),
     &                   zgrid,zbeam,l2symtry,l4symtry,pboundxy,.true.)
        endif
      enddo

      if (lcallscraper) call execuser("callscraper")

      do is=1,ns
c       --- do it only if species was advanced on this time step
        if (pgroup%ldts(is)) then
          call processlostpart(pgroup,is,clearlostpart,
     &                         time+dt*pgroup%ndts(is-1),zbeam)
        endif
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
      elseif(solvergeom==XYgeom) then
        nz = 0
      elseif(solvergeom==Zgeom) then
        nx = 0
        ny = 0
      elseif(solvergeom==Rgeom) then
        ny = 0
        nz = 0
      endif

c     --- Put grid min's to zero as appropriate for the symmetry
      if(solvergeom==XYZgeom) then
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
      use InGen
      use InGen3d
      use Picglb
      use FieldSolveAPI
      integer(ISZ):: n
      real(kind=8),target:: x(n),y(n),z(n),p(n)

c     --- Don't do anything if there are no particles input.
      if (n == 0) return

c     --- Call the appropriate routine to get the phi
      if (fstype == 12) then
        jsfsapi = -1
        npfsapi = n
        xfsapi => x
        yfsapi => y
        zfsapi => z
        phifsapi => p
        call execuser("fetchphiregistered")
        npfsapi = 0
        nullify(xfsapi)
        nullify(yfsapi)
        nullify(zfsapi)
        nullify(phifsapi)
      elseif (solvergeom==XYZgeom) then
        call fetchphi3d(n,x,y,z,p)
      elseif (solvergeom==RZgeom) then
        call setphirz(n,x,y,z,p,zgrid)
      elseif (solvergeom==XZgeom) then
        call setphixz(n,x,y,z,p,zgrid)
      elseif (solvergeom==XYgeom) then
        call remark("fetchphi not implemented in XY.")
        stop
      elseif (solvergeom==Zgeom) then
        call setphiz(n,z,p,zgrid)
      elseif (solvergeom==AMRgeom) then
        call cho_getphi3d(n,x,y,z,0.,p,-1,-1)
      else
        p = 0.
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
      integer(ISZ):: ixp1,iyp1,izp1
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
          ixp1 = ix + 1
          iyp1 = iy + 1
          izp1 = iz + 1
          if (abs(x(ii)) == xmmax) ixp1 = ix
          if (abs(y(ii)) == ymmax) iyp1 = iy
          if (abs(z(ii)) == zmmaxp+dz) izp1 = iz
          if (xmmin  <= abs(x(ii)) .and. abs(x(ii)) <= xmmax .and.
     &        ymmin  <= abs(y(ii)) .and. abs(y(ii)) <= ymmax .and.
     &        zmminp-dz <=     z(ii)  .and.     z(ii)  <= zmmaxp+dz) then
            p(ii) = phip(ix  ,iy,  iz  )*(1.-wx)*(1.-wy)*(1.-wz) +
     &              phip(ixp1,iy,  iz  )*    wx *(1.-wy)*(1.-wz) +
     &              phip(ix  ,iyp1,iz  )*(1.-wx)*    wy *(1.-wz) +
     &              phip(ixp1,iyp1,iz  )*    wx *    wy *(1.-wz) +
     &              phip(ix  ,iy,  izp1)*(1.-wx)*(1.-wy)*    wz  +
     &              phip(ixp1,iy,  izp1)*    wx *(1.-wy)*    wz  +
     &              phip(ix  ,iyp1,izp1)*(1.-wx)*    wy *    wz  +
     &              phip(ixp1,iyp1,izp1)*    wx *    wy *    wz
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
          ixp1 = ix + 1
          iyp1 = iy + 1
          izp1 = iz + 1
          if (abs(x(ii)) == xmmax) ixp1 = ix
          if (abs(y(ii)) == ymmax) iyp1 = iy
          if (abs(z(ii)) == zmmaxp+dz) izp1 = iz
          if (xmmin  <=     x(ii)  .and.     x(ii)  <= xmmax .and.
     &        ymmin  <= abs(y(ii)) .and. abs(y(ii)) <= ymmax .and.
     &        zmminp-dz <=     z(ii)  .and.     z(ii)  <= zmmaxp+dz) then
            p(ii) = phip(ix  ,iy,  iz  )*(1.-wx)*(1.-wy)*(1.-wz) +
     &              phip(ixp1,iy,  iz  )*    wx *(1.-wy)*(1.-wz) +
     &              phip(ix  ,iyp1,iz  )*(1.-wx)*    wy *(1.-wz) +
     &              phip(ixp1,iyp1,iz  )*    wx *    wy *(1.-wz) +
     &              phip(ix  ,iy,  izp1)*(1.-wx)*(1.-wy)*    wz  +
     &              phip(ixp1,iy,  izp1)*    wx *(1.-wy)*    wz  +
     &              phip(ix  ,iyp1,izp1)*(1.-wx)*    wy *    wz  +
     &              phip(ixp1,iyp1,izp1)*    wx *    wy *    wz
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
          ixp1 = ix + 1
          iyp1 = iy + 1
          izp1 = iz + 1
          if (abs(x(ii)) == xmmax) ixp1 = ix
          if (abs(y(ii)) == ymmax) iyp1 = iy
          if (abs(z(ii)) == zmmaxp+dz) izp1 = iz
          if (xmmin  <= x(ii) .and. x(ii) <= xmmax .and.
     &        ymmin  <= y(ii) .and. y(ii) <= ymmax .and.
     &        zmminp-dz <= z(ii) .and. z(ii) <= zmmaxp+dz) then
            p(ii) = phip(ix  ,iy,  iz  )*(1.-wx)*(1.-wy)*(1.-wz) +
     &              phip(ixp1,iy,  iz  )*    wx *(1.-wy)*(1.-wz) +
     &              phip(ix  ,iyp1,iz  )*(1.-wx)*    wy *(1.-wz) +
     &              phip(ixp1,iyp1,iz  )*    wx *    wy *(1.-wz) +
     &              phip(ix  ,iy,  izp1)*(1.-wx)*(1.-wy)*    wz  +
     &              phip(ixp1,iy,  izp1)*    wx *(1.-wy)*    wz  +
     &              phip(ix  ,iyp1,izp1)*(1.-wx)*    wy *    wz  +
     &              phip(ixp1,iyp1,izp1)*    wx *    wy *    wz
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

      if (napertures == 0) return

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

      if (napertures == 0) return

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
      subroutine getphipforparticles(indts)
      use InGen3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles
      integer(ISZ):: indts

#ifdef MPIPARALLEL

c     --- Distribute phi among the processes so each has phi in its
c     --- particle domain.
      if(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==Zgeom) then
        call getphiforparticlesrz()
      elseif (solvergeom==XYgeom) then
        call remark("getphipforparticles not implemented in parallel for XYgeom.")
        stop
c      else if(solvergeom==Zgeom) then
c        call getphipforparticlesz()
      elseif (solvergeom==XYZgeom) then
        call getphipforparticles3d(nx,ny,nz,phi,
     &                             nxp,nyp,nzp,phipndts(:,:,:,indts))
      end if
#endif

      return
      end
c===========================================================================
      subroutine getphipforparticles3d(nx,ny,nz,phi,
     &                                nxp,nyp,nzp,phip)
      integer(ISZ):: nx,ny,nz,nxp,nyp,nzp
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: phip(0:nxp,0:nyp,-1:nzp+1)

#ifdef MPIPARALLEL
        call getphipforparticles3d_parallel(nx,ny,nz,phi,nxp,nyp,nzp,phip)
#endif

      return
      end
c===========================================================================
      subroutine getphiforfields()
      use InGen3d
      use InMesh3d
      use Fields3d
      use Fields3dParticles

c This gets the phi at iz=nz and in the guard cells. It is only needed for
c the 3d FFT and tridiagonal based solvers for the parallel version.

#ifdef MPIPARALLEL

c     --- Distribute phi among the processes so each has phi in its
c     --- particle domain.
      if(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==Zgeom) then
c       call getphiforparticlesrz()
c      else if(solvergeom==Zgeom) then
c        call getphiforparticlesz()
      elseif (solvergeom==XYZgeom) then
        call getphiforfields3d(nx,ny,nz,phi)
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
      use Picglb, only: zgrid
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Fields3d
      use Subcycling
      use Fields3dParticles
      use Efields3d
      use LatticeInternal
      use InjectVars
      use GridBoundary3d
      use Parallel,Only:my_index,nslaves
      use w3d_interfaces
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
      integer(ISZ):: izmin,izmax,tnz,indts,tmpnsndts,getnsndtsforsubcycling
      integer(ISZ):: isndts
      if (lw3dtimesubs) substarttime = wtime()
      timetemp = wtime()

c     --- If no field solve, return 
      if ( fstype == -1) return 

c     --- If a field solver is registered, call it and return.
      if (fstype == 12) then
        call execuser('fieldsolregistered')

#ifdef MPIPARALLEL
c       --- Set the flag signalling that rho has been parallel summed
c       --- on the boundaries. This is done here assuming that
c       --- after the field solve, the sums have been done.
c       --- This also avoids having the registered solvers having
c       --- to set the flag.
        lrhosummedonboundaries = .true.
#endif

!$OMP MASTER
        if (lw3dtimesubs) timefieldsol3d = timefieldsol3d + wtime() - substarttime
        fstime = fstime + (wtime() - timetemp)
!$OMP END MASTER
        return
      endif

c     --- Calculate rstar here since it is only used by the field solver and
c     --- to ensure that it is set when the field solver is called.
      if (linbend) call setrstar(rstar(-1),nz,dz,zmmin,zgrid)

c     --- Loop over the subcyling groups and do any field solves that
c     --- are necessary.
c     --- Do loop in reverse order so that rho and phi end up with the arrays
c     --- for the speices with the smallest timestep.
      tmpnsndts = getnsndtsforsubcycling()
      do indts=tmpnsndts-1,0,-1
        isndts = min(indts,nsndtsphi3d-1)

        if (.not. ldts(indts) .and.
     &      (ndtsaveraging == 0 .and. .not. any(ldts))) cycle

c       --- For serial version, the arrays rho and phi are pointed to the
c       --- appropriate arrays for the current ndts group.
c       --- For parallel version, each processor sends rho to neighboring
c       --- processors whose field solve region overlap its particle region.
        call setrhoandphiforfieldsolve(rhopndts(:,:,:,nrhopndtscopies-1,indts),
     &                                 phipndts(:,:,:,isndts))

c       --- enforce periodicity
c       --- This is done on the rho array, rather than rhop, since that will
c       --- be used for the field solve. Also, in some cases, the rhop array
c       --- will not cover the full axial extent of the system, so 
c       --- periodicity can never be enforced on rho which does cover the full
c       --- system.
        call perrho3d()

        if (fstype == 3 .or. fstype == 7 .or. fstype == 11 .or. fstype == 13) then
c         --- SOR and multigrid field solvers - deal with bends directly and
c         --- set axial boundary conditions (call to perphi3d is not needed)
c         --- Chombo knows nothing of bends and does its own b.c.'s
          call vp3d(iwhich)

        else if (.not. linbend) then
c         --- If not in bend, call VP3D
          if (fstype == 0 .or. fstype == 1 .or. fstype == 2 .or.
     &        fstype == 4 .or. fstype == 5 .or. fstype == 6 .or.
     &        fstype == 8 .or. fstype == 9) then
            izmin = 0
            izmax = nz
            if (bound0  == dirichlet .and. my_index == 0) izmin = 1
            if (boundnz == dirichlet .and. my_index >= nslaves-1) izmax = nz-1
            tnz = izmax - izmin + 1
#ifndef MPIPARALLEL
            call copyarry(rhopndts(:,:,izmin,nrhopndtscopies-1,indts),
     &                    phipndts(:,:,izmin,isndts),
     &                    (nx+1)*(ny+1)*(tnz))
#else
            call copyarry(rho(:,:,izmin),phi(:,:,izmin),(nx+1)*(ny+1)*(tnz))
#endif
          endif
          call vp3d(iwhich)
          call getphiforfields()
          if (bound0==periodic) call perphi3d()

        else
c         --- Call perturbative bent beam field solver
          call bendfieldsol3d()
          call getphiforfields()
          if (bound0==periodic) call perphi3d()

        endif

c       --- Distribute phi among the processes so each has phi in its
c       --- particle domain.
        call getphipforparticles(indts)

      enddo
c     --- End loop over ndts groups

c     --- Points to the values that are used by the most frequently time
c     --- advanced group.
      call assignrhoandphi(rhopndts(:,:,:,nrhopndtscopies-1,0),
     &                     phipndts(:,:,:,0))

#ifdef MPIPARALLEL
c     --- Set the flag signalling that rho has been parallel summed
c     --- on the boundaries. This is done here only after rho has been
c     --- setup for all subcycling and selfb groups.
      lrhosummedonboundaries = .true.
#endif

!$OMP MASTER
      if (lw3dtimesubs) timefieldsol3d = timefieldsol3d + wtime() - substarttime
      fstime = fstime + (wtime() - timetemp)
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
      use Picglb3d
      use Fields3d

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
      call vp3d(-1)

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
c=============================================================================
      subroutine vp3d(iwhich)
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
      use BoltzmannElectrons
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

      if (fstype .ne. -1 .and. fstype .ne. 12) then
        if(solvergeom==RZgeom .or.
     &     solvergeom==XZgeom .or.
     &     solvergeom==XYgeom .or.
     &     solvergeom==Rgeom  .or.
     &     solvergeom==Zgeom) fstype=10
        if(solvergeom==AMRgeom) fstype = 11
      endif

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
           call psor3df(iwhich,nx,ny,nz,phi,rho,
     &                 phi,rho,rstar, 
     &                 dx,dy,dz,xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid,linbend,
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
        if (minval(electrontemperature) == 0) then
          call multigrid3df(iwhich,nx,ny,izfsmax-izfsmin,nzfull,dx,dy,dz,
     &            phi(0,0,izfsmin-1),rho(0,0,izfsmin),rstar(izfsmin-1),
     &            linbend,bound0,boundnz,boundxy,l2symtry,l4symtry,
     &            xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid)
        else
          call multigridbe3df(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                        rstar,linbend,
     &                        bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                        xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid)
        endif

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
c       --- RZ full-multigrid solver
       if (solvergeom==RZgeom .or. solvergeom==XZgeom) then
         call multigridrzf(iwhich,phi,rho,nx,nz)
       elseif (solvergeom==XYgeom) then
         call multigridxyf2(iwhich,phi(:,:,0),rho,nx,ny)
       endif
      elseif (fstype == 11) then
c         --- Chombo AMR  full-multigrid solver
        call cho_solve3d(iwhich,nx,ny,nz,nzfull,dx,dy,dz,
     &                   l2symtry,l4symtry,xmmin,ymmin,zmmin,zmminglobal)
      elseif (fstype == 12) then
        if (iwhich <= 0) then
          call execuser('fieldsolregistered')
        endif

      elseif (fstype == 13) then

        call multigridbe3df(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                      rstar,linbend,
     &                      bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                      xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid)

      endif

!$OMP MASTER
      if (lw3dtimesubs) timevp3d = timevp3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
