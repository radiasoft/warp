#include "top.h"
c=======================================================================
c@(#) File ENV.M, version $Revision: 3.21 $, $Date: 2003/04/10 16:21:55 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  It is the source file for the package ENV of the PIC code WARP,
c  but it may be useful by itself.  It consists of an simple 
c  envelope code and a BASIS interface.  
c
c=======================================================================
      subroutine envinit
      use ENVversion

c  Called at first reference to package (not nec. a "run" etc.).


      call envvers (STDOUT)

      return
      end
c=============================================================================
      subroutine envvers (iout)
      use ENVversion
      integer(ISZ):: iout
c  Echoes code version,etc. to output files when they're created
      call printpkgversion(iout,"Envelope solver",versenv)
      return
      end
c=============================================================================
      subroutine envgen()
      use Constant
      use ENVvars
      use Beam_acc
      use InGen
      use Lattice
      use LatticeInternal

c  Invoked by the GENERATE command, it allocates storage
c  for the envelope code.

c  Send a startup message to the standard output

      if (lenvout) then
        call remark(" ***  envelope solver package ENV generating")
      endif

c  Calculate derived quantities

      call derivqty

      nenv = (zu - zl)/dzenv + 0.5

c     --- allocate dynamic arrays ("group allot")
      call gallot("ENVvars",0)

c     --- reset Lattice arrays to final size
c     --- This is not really needed since it is also done in envexe below.
c     --- The absolute value is taken since dzenv may be negative (when
c     --- the envelope is running backward).
      zlatbuffer = abs(dzenv)
      call resetlat
      nzl = 0
      zlmin = 0.
      zlmin = 0.

c     --- Set default values for region over which the tune is calculated.
      if (tunezs == tuneze) then
        tunezs = zl
        if (tunelen == 0.) then
          tuneze = zu
        else
          tuneze = tunezs + tunelen
        endif
      endif

      return
      end
c=============================================================================
      subroutine envexe()
      use ENVvars
      use OutParams
      use Parallel
      use Lattice
      use LatticeInternal

c  Invoked by the STEP command, it runs the envelope code.
c  passes sigma and sigma0 back to TOP through OutParams


      real(kind=8):: wtimeoff

c  Send a startup message to the standard output

      if (lenvout) then
        call remark(" ***  envelope solver package ENV running")
      endif

c     --- reset Lattice arrays to final size
c     --- The absolute value is taken since dzenv may be negative (when
c     --- the envelope is running backward).
      zlatbuffer = abs(dzenv)
      if (nzl .ne. 0) then
        call resetlat
        nzl = 0
        zlmin = 0.
        zlmin = 0.
      endif

c  Execute evelope solve, and time it
      call wtimeon
      call envx
      envtime = wtimeoff()
 
c     --- Only print the output if requested.
      if (lenvout .and. my_index == 0) then
        print*," ***  envelope solver: time =",envtime," milliseconds"

c       --- Print the phase advance if it was calculated.
        if (zl <= tunezs .and. tunezs <  zu .and.
     &      zl <  tuneze .and. tuneze <= zu) then
          write(STDOUT,99) sig0x, sigx, sigx/dvnz(sig0x), 
     &                     sig0y, sigy, sigy/dvnz(sig0y), 
     &                     deltaa, deltab, deltaap, deltabp
   99     format("       sigma0x = ",f13.5,"    sigmax = ",f13.5,/,
     &           "       sigmax/sigma0x = ",f10.7,/, 
     &           "       sigma0y = ",f13.5,"    sigmay = ",f13.5,/,
     &           "       sigmay/sigma0y = ",f10.7,/, 
     &           "       delta a = ",e13.5,"   delta b = ",e13.5,/,
     &           "      delta ap = ",e13.5,"  delta bp = ",e13.5)
        endif
      endif

      sigmax  = sigx
      sigma0x = sig0x

      sigmay  = sigy
      sigma0y = sig0y

      return
      end
c=======================================================================
      subroutine envfin()
      use ENVvars

c  Deallocates dynamic storage used by test driver


      call gfree ("ENVvars")

      return
      end
c=======================================================================
      integer function envxport(np,z,a,ap,b,bp,x,xp,y,yp)
      use ENVvars
      integer(ISZ):: np
      real(kind=8):: z(np)
      real(kind=8):: a(np),ap(np),b(np),bp(np)
      real(kind=8):: x(np),xp(np),y(np),yp(np)

c  Export routine from envelope package:
c      envelope parameters at a set of z values passed in
c  Called during the particle loading process of the particle code.
c  Returns 1 if any input z is outside the range over which the envelope
c      is defined.
c  Inputs:  np, z(1:np)
c  Outputs: a(1:np), ap(1:np), b(1:np), bp(1:np)
c           x(1:np), xp(1:np), y(1:np), yp(1:np)


      integer(ISZ):: i,ip
      real(kind=8):: zpmin,zpmax,w0,w1

c  Validity check on input

      zpmin = 1.e6
      zpmax = -1.e6
      do 500 ip=1,np
         zpmin = min(zpmin, z(ip))
         zpmax = max(zpmax, z(ip))
  500 continue
      if (zpmin < zl .or. zpmax > zu) then
         envxport = 1
         return
      endif

c  Interpolate to get envelope values

      do ip=1,np
        i = int((z(ip)-zl)/dzenv)
        w1 = (z(ip)-zl)/dzenv - i
        w0 = 1. - w1
        a(ip)  = w0*aenv(i)  + w1*aenv(i+1)
        ap(ip) = w0*apenv(i) + w1*apenv(i+1)
        b(ip)  = w0*benv(i)  + w1*benv(i+1)
        bp(ip) = w0*bpenv(i) + w1*bpenv(i+1)
        x(ip)  = w0*xenv(i)  + w1*xenv(i+1)
        xp(ip) = w0*xpenv(i) + w1*xpenv(i+1)
        y(ip)  = w0*yenv(i)  + w1*yenv(i+1)
        yp(ip) = w0*ypenv(i) + w1*ypenv(i+1)
      enddo

      envxport = 0
      return
      end
c=======================================================================
      subroutine envx
      use ENVvars
      use Beam_acc
      use InGen

c  Interface to ENVELOPE, using variables from the database for
c  package ENV, and global variables from package TOP.


c  Work around compiler bug on Sun that leads to missing VAR_SEG7
      integer(ISZ):: junk
      junk = 1

      call envelatt (
     & nenv, lrelativ, gammabar, genprv, dzenv,
     & aion, zion, ekin, vbeam, ibeam, a0, b0, ap0, bp0, x0, y0, xp0, yp0,
     & emitx, emity, 
     & dedr, dexdx, deydy, dbdr, zl, zu,
     & tunezs, tuneze, 
     & sig0x, sigx, sig0y, sigy, deltaa, deltaap, deltab, deltabp,
     & aenv(0), apenv(0), benv(0), bpenv(0), vzenv(0), xorb(0), xporb(0), yorb(0), yporb(0),
     & fqxenv(0), fqyenv(0), fuxenv(0), fuyenv(0), xenv(0), yenv(0), xpenv(0), ypenv(0), zenv(0))
 
      return
      end
c=======================================================================
      subroutine envelatt (
     & nenv, lrelativ, gammabar, genprv, dzenv,
     & aion, zion, ekin, vbeam, ibeam, a0, b0, ap0, bp0, x0, y0, xp0, yp0,
     & emitx, emity, 
     & dedr, dexdx, deydy, dbdr, zl, zu,
     & tunezs, tuneze, 
     & sig0x, sigx, sig0y, sigy, deltaa, deltaap, deltab, deltabp,
     & aenv, apenv, benv, bpenv, vzenv, xorb, xporb, yorb, yporb,
     & fqxenv, fqyenv, fuxenv, fuyenv, xenv, yenv, xpenv, ypenv, zenv)
      use ENVtune
      use Constant
      use InGen
      integer(ISZ):: nenv
      logical(ISZ):: lrelativ
      real(kind=8):: ibeam
      real(kind=8):: mx11,mx22,mx110,mx220
      real(kind=8):: my11,my22,my110,my220
      real(kind=8):: gammabar,genprv,dzenv
      real(kind=8):: aion,zion,ekin,vbeam,a0,b0,ap0,bp0,x0,y0,xp0,yp0
      real(kind=8):: emitx,emity 
      real(kind=8):: dedr,dexdx,deydy,dbdr,zl,zu
      real(kind=8):: tunezs, tuneze 
      real(kind=8):: sig0x,sigx,sig0y,sigy,deltaa,deltaap,deltab,deltabp
      real(kind=8):: aenv(0:nenv),apenv(0:nenv),benv(0:nenv),bpenv(0:nenv)
      real(kind=8):: xorb(0:nenv),xporb(0:nenv),yorb(0:nenv),yporb(0:nenv)
      real(kind=8):: vzenv(0:nenv)
      real(kind=8):: fqxenv(0:nenv),fqyenv(0:nenv),fuxenv(0:nenv),fuyenv(0:nenv)
      real(kind=8):: xenv(0:nenv), yenv(0:nenv), xpenv(0:nenv), ypenv(0:nenv)
      real(kind=8):: zenv(0:nenv)

c  This routine does the actual envelope calculation for the general lattice. 
c  The input arguments are described in the variable description files
c  of packages ENV and TOP.  The package takes linear focusing components 
c  from various lattice and multipole arrays.  Any linear field componets 
c  imported as gridded data will be missed.  

      integer(ISZ):: istep
      real(kind=8):: p,gamma,emitnx_l, emitny_l
      real(kind=8):: dz,z,vz,ek,time
      real(kind=8):: a,ap,b,bp,fa,fb
      real(kind=8):: x,xp,fx,y,yp,fy
      real(kind=8):: x1,xp1,fx1,x2,xp2,fx2
      real(kind=8):: y1,yp1,fy1,y2,yp2,fy2
      real(kind=8):: x01,x0p1,fx01,x02,x0p2,fx02
      real(kind=8):: y01,y0p1,fy01,y02,y0p2,fy02
      real(kind=8):: faccl,fqxplot,fqyplot,fuxplot,fuyplot
      real(kind=8):: cossig0x,cossigx,tansig0x,tansigx
      real(kind=8):: cossig0y,cossigy,tansig0y,tansigy
      real(kind=8):: wz

c  For convenience, this routine accessesd the lattice database and is not
c  linked to the universe solely by its calling sequence.


c  Set fundamental constants, derived quantities, initial conditions

c     --- beam axial velocity and energy
      vz = vbeam
      ek = ekin
c     --- gammabar is a database qty passed in thru calling sequence
      gamma = gammabar
c     --- generalized perveance; load database variable, too
      p = ibeam * zion*echarge / (2.*pi*eps0 * aion*amu * (vz*gamma)**3)
      genprv = p
c     --- get value of normalized emittance and store in local variable
c     --- which can be changed
      emitnx_l = emitx*vz*gamma/clight
      emitny_l = emity*vz*gamma/clight
c     --- step size from database
      dz = dzenv
c     Initial condtions:
      time = 0.
c     --- for envelope
      z    = zl
      a    = a0
      ap   = ap0
      b    = b0
      bp   = bp0
c     --- envelope centroid (single particle)
      x    = x0
      xp   = xp0
      y    = y0
      yp   = yp0
c     --- x-plane orbit: for depressed tune, principal  
      x1   = a
      xp1  = ap
c     --- x-plane orbit: for depressed tune, starts at x = 0
      x2   = 0.
      xp2  = 1.
c     --- x-plane orbit: for undepressed tune, doesn't feel s.c. 
      x01  = a
      x0p1 = ap
c     --- x-plane orbit: for undepressed tune, doesn't feel s.c., starts at x=0
      x02  = 0.
      x0p2 = 1.
c     --- y-plane orbit: for depressed tune, principal  
      y1   = b
      yp1  = bp
c     --- y-plane orbit: for depressed tune, starts at x = 0
      y2   = 0.
      yp2  = 1.
c     --- y-plane orbit: for undepressed tune, doesn't feel s.c. 
      y01  = b
      y0p1 = bp
c     --- y-plane orbit: for undepressed tune, doesn't feel s.c., starts at x=0
      y02  = 0.
      y0p2 = 1.

c  Initialize facility and update any variables which are functions of z.

      call envfofzinit()
      call envsetfofz(z,p,vz,emitnx_l,emitny_l,gamma)

c  Obtain initial force

       call envflatt (z, time, zion, aion, lrelativ, gamma, vz, ek, dz,
     & dedr,dexdx,deydy,dbdr,p, 
     & a, b, x, y, x1, x2, x01, x02, y1, y2, y01, y02, emitnx_l, emitny_l,  
     & fa, fb, fx, fy, fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     & fqxplot, fqyplot, fuxplot, fuyplot)

c  Main loop of envelope code

      do istep = 0,nenv
c        --- store data for this step
         aenv(istep)   = a
         apenv(istep)  = ap
         benv(istep)   = b
         bpenv(istep)  = bp
         xenv(istep)   = x
         xpenv(istep)  = xp
         yenv(istep)   = y
         ypenv(istep)  = yp
         vzenv(istep)  = vz
         zenv(istep)   = z
         xorb(istep)   = x1
         xporb(istep)  = xp1
         yorb(istep)   = y1
         yporb(istep)  = yp1
         fqxenv(istep) = fqxplot
         fqyenv(istep) = fqyplot
         fuxenv(istep) = fuxplot
         fuyenv(istep) = fuyplot 
c        --- Check if just before starting point of region over which the tune
c        --- is to be calculated.
c        --- NOTE: this calculation is not correct when dzenv < 0.
         if (tunezs-dz < z .and. z <= tunezs) then
           astart    = a
           apstart   = ap
           bstart    = b
           bpstart   = bp
c
           x1start   = x1
           xp1start  = xp1
           x2start   = x2
           xp2start  = xp2
           x01start  = x01
           x0p1start = x0p1
           x02start  = x02
           x0p2start = x0p2
c
           y1start   = y1
           yp1start  = yp1
           y2start   = y2
           yp2start  = yp2
           y01start  = y01
           y0p1start = y0p1
           y02start  = y02
           y0p2start = y0p2
         endif
c        --- Check if just after starting point of region over which the tune
c        --- is to be calculated.  Linear interpolate to find values at the
c        --- starting point.
c        --- NOTE: this calculation is not correct when dzenv < 0.
         if (tunezs < z .and. z < tunezs+dz) then
           wz = 1. - (z - tunezs)/dz
           astart    = astart   *(1. - wz) + a   *wz
           apstart   = apstart  *(1. - wz) + ap  *wz
           bstart    = bstart   *(1. - wz) + b   *wz
           bpstart   = bpstart  *(1. - wz) + bp  *wz
c
           x1start   = x1start  *(1. - wz) + x1  *wz
           xp1start  = xp1start *(1. - wz) + xp1 *wz
           x2start   = x2start  *(1. - wz) + x2  *wz
           xp2start  = xp2start *(1. - wz) + xp2 *wz
           x01start  = x01start *(1. - wz) + x01 *wz
           x0p1start = x0p1start*(1. - wz) + x0p1*wz
           x02start  = x02start *(1. - wz) + x02 *wz
           x0p2start = x0p2start*(1. - wz) + x0p2*wz
c
           y1start   = y1start  *(1. - wz) + y1  *wz
           yp1start  = yp1start *(1. - wz) + yp1 *wz
           y2start   = y2start  *(1. - wz) + y2  *wz
           yp2start  = yp2start *(1. - wz) + yp2 *wz
           y01start  = y01start *(1. - wz) + y01 *wz
           y0p1start = y0p1start*(1. - wz) + y0p1*wz
           y02start  = y02start *(1. - wz) + y02 *wz
           y0p2start = y0p2start*(1. - wz) + y0p2*wz
         endif
c        --- Check if just before ending point of region over which the tune
c        --- is to be calculated.
c        --- NOTE: this calculation is not correct when dzenv < 0.
         if (tuneze-dz < z .and. z <= tuneze) then
           aend    = a
           apend   = ap
           bend    = b
           bpend   = bp
c
           x1end   = x1
           xp1end  = xp1
           x2end   = x2
           xp2end  = xp2
           x01end  = x01
           x0p1end = x0p1
           x02end  = x02
           x0p2end = x0p2
c
           y1end   = y1
           yp1end  = yp1
           y2end   = y2
           yp2end  = yp2
           y01end  = y01
           y0p1end = y0p1
           y02end  = y02
           y0p2end = y0p2
         endif
c        --- Check if just after ending point of region over which the tune
c        --- is to be calculated.  Linear interpolate to find values at the
c        --- ending point.
c        --- NOTE: this calculation is not correct when dzenv < 0.
         if (tuneze < z .and. z < tuneze+dz) then
           wz = 1. - (z - tuneze)/dz
           aend    = aend   *(1. - wz) + a   *wz
           apend   = apend  *(1. - wz) + ap  *wz
           bend    = bend   *(1. - wz) + b   *wz
           bpend   = bpend  *(1. - wz) + bp  *wz
c
           x1end   = x1end  *(1. - wz) + x1  *wz
           xp1end  = xp1end *(1. - wz) + xp1 *wz
           x2end   = x2end  *(1. - wz) + x2  *wz
           xp2end  = xp2end *(1. - wz) + xp2 *wz
           x01end  = x01end *(1. - wz) + x01 *wz
           x0p1end = x0p1end*(1. - wz) + x0p1*wz
           x02end  = x02end *(1. - wz) + x02 *wz
           x0p2end = x0p2end*(1. - wz) + x0p2*wz
c
           y1end   = y1end  *(1. - wz) + y1  *wz
           yp1end  = yp1end *(1. - wz) + yp1 *wz
           y2end   = y2end  *(1. - wz) + y2  *wz
           yp2end  = yp2end *(1. - wz) + yp2 *wz
           y01end  = y01end *(1. - wz) + y01 *wz
           y0p1end = y0p1end*(1. - wz) + y0p1*wz
           y02end  = y02end *(1. - wz) + y02 *wz
           y0p2end = y0p2end*(1. - wz) + y0p2*wz
         endif
c        Advance the envelope and orbits; 
c        --- envelope a in x-z plane
         ap   = ap   + .5*dz*fa   - 0.5*dz*faccl*ap
         a    = a    + ap*dz
c        --- envelope b in y-z plane
         bp   = bp   + .5*dz*fb   - 0.5*dz*faccl*bp
         b    = b    + bp*dz
c        --- envelope x centroid
         xp   = xp   + .5*dz*fx   - 0.5*dz*faccl*xp
         x    = x    + xp*dz
c        --- envelope y centroid
         yp   = yp   + .5*dz*fy   - 0.5*dz*faccl*yp
         y    = y    + yp*dz
c        --- x-orbit, 1 
         xp1  = xp1  + .5*dz*fx1  - 0.5*dz*faccl*xp1
         x1   = x1   + xp1*dz
c        --- x-orbit, 2 
         xp2  = xp2  + .5*dz*fx2  - 0.5*dz*faccl*xp2
         x2   = x2   + xp2*dz
c        --- x-orbit, 3 
         x0p1 = x0p1 + .5*dz*fx01 - 0.5*dz*faccl*x0p1
         x01  = x01  + x0p1*dz
c        --- x-orbit, 4 
         x0p2 = x0p2 + .5*dz*fx02 - 0.5*dz*faccl*x0p2
         x02  = x02  + x0p2*dz
c        --- y-orbit, 1 
         yp1  = yp1  + .5*dz*fy1  - 0.5*dz*faccl*yp1
         y1   = y1   + yp1*dz
c        --- y-orbit, 2 
         yp2  = yp2  + .5*dz*fy2  - 0.5*dz*faccl*yp2
         y2   = y2   + yp2*dz
c        --- y-orbit, 3 
         y0p1 = y0p1 + .5*dz*fy01 - 0.5*dz*faccl*y0p1
         y01  = y01  + y0p1*dz
c        --- y-orbit, 4 
         y0p2 = y0p2 + .5*dz*fy02 - 0.5*dz*faccl*y0p2
         y02  = y02  + y0p2*dz
c        --- z coordinate
         z = z + dz
         time = time + dz/vz
c        --- Update any variables which are functions of z.
         call envsetfofz(z,p,vz,emitnx_l,emitny_l,gamma)
c        --- get force at new position
         call envflatt (z, time, zion, aion, lrelativ, gamma, vz, ek, dz,
     &    dedr,dexdx,deydy,dbdr,p,
     &    a, b, x, y, x1, x2, x01, x02, y1, y2, y01, y02, emitnx_l, emitny_l,
     &    fa, fb, fx, fy, fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     &    fqxplot, fqyplot, fuxplot, fuyplot)
c        --- complete the "velocity" advance steps
         ap   = ap   + .5*dz*fa   - 0.5*dz*faccl*ap
         bp   = bp   + .5*dz*fb   - 0.5*dz*faccl*bp
         xp   = xp   + .5*dz*fx   - 0.5*dz*faccl*xp
         yp   = yp   + .5*dz*fy   - 0.5*dz*faccl*yp

         xp1  = xp1  + .5*dz*fx1  - 0.5*dz*faccl*xp1
         xp2  = xp2  + .5*dz*fx2  - 0.5*dz*faccl*xp2
         x0p1 = x0p1 + .5*dz*fx01 - 0.5*dz*faccl*x0p1
         x0p2 = x0p2 + .5*dz*fx02 - 0.5*dz*faccl*x0p2

         yp1  = yp1  + .5*dz*fy1  - 0.5*dz*faccl*yp1
         yp2  = yp2  + .5*dz*fy2  - 0.5*dz*faccl*yp2
         y0p1 = y0p1 + .5*dz*fy01 - 0.5*dz*faccl*y0p1
         y0p2 = y0p2 + .5*dz*fy02 - 0.5*dz*faccl*y0p2
      enddo 

c     --- Only do the following if tunezs and tuneze and within the
c     --- envelope calculation range zl and zu.
      if (zl <= tunezs .and. tunezs <  zu .and.
     &    zl <  tuneze .and. tuneze <= zu) then

c       --- compute the undepressed x-tune (phase adv / cell)
        mx110 = (x0p1start*x02end -   x01end*  x0p2start)/
     &      dvnz(x0p1start*x02start - x01start*x0p2start)
        mx220 = (x01start*x0p2end -   x0p1end*  x02start)/
     &      dvnz(x01start*x0p2start - x0p1start*x02start)
        cossig0x = .5 * (mx110 + mx220)
        cossig0x = max(-1., min(1.,cossig0x))
        tansig0x = sqrt(1. - cossig0x**2)/dvnz(cossig0x)
        sig0x = (180./pi) * atan(tansig0x)
        if (sig0x < 0.) sig0x = sig0x + 180.

c       --- compute the depressed x-tune (phase adv / cell)
        mx11 = (xp1start*x2end -   x1end*  xp2start)/
     &     dvnz(xp1start*x2start - x1start*xp2start)
        mx22 = (x1start*xp2end -   xp1end*  x2start)/
     &     dvnz(x1start*xp2start - xp1start*x2start)
        cossigx = .5 * (mx11 + mx22)
        cossigx = max(-1., min(1.,cossigx))
        tansigx = sqrt(1. - cossigx**2)/dvnz(cossigx)
        sigx = (180./pi) * atan(tansigx)
        if (sigx < 0.) sigx = sigx + 180.

c       --- compute the undepressed y-tune (phase adv / cell)
        my110 = (y0p1start*y02end -   y01end*  y0p2start)/
     &      dvnz(y0p1start*y02start - y01start*y0p2start)
        my220 = (y01start*y0p2end -   y0p1end*  y02start)/
     &      dvnz(y01start*y0p2start - y0p1start*y02start)
        cossig0y = .5 * (my110 + my220)
        cossig0y = max(-1., min(1.,cossig0y))
        tansig0y = sqrt(1. - cossig0y**2)/dvnz(cossig0y)
        sig0y = (180./pi) * atan(tansig0y)
        if (sig0y < 0.) sig0y = sig0y + 180.

c       --- compute the depressed y-tune (phase adv / cell)
        my11 = (yp1start*y2end -   y1end*  yp2start)/
     &     dvnz(yp1start*y2start - y1start*yp2start)
        my22 = (y1start*yp2end -   yp1end*  y2start)/
     &     dvnz(y1start*yp2start - yp1start*y2start)
        cossigy = .5 * (my11 + my22)
        cossigy = max(-1., min(1.,cossigy))
        tansigy = sqrt(1. - cossigy**2)/dvnz(cossigy)
        sigy = (180./pi) * atan(tansigy)
        if (sigy < 0.) sigy = sigy + 180.

c       --- compute quantities that serve as measures of "matchedness"
        deltaa  = aend  - astart
        deltaap = apend - apstart
        deltab  = bend  - bstart
        deltabp = bpend - bpstart

      endif

      return
      end
c=======================================================================
      subroutine envflatt(z,time,zion,aion,lrelativ,gamma,vz,ek,dz,
     & dedr,dexdx,deydy,dbdr,p,
     & a, b, x, y, x1, x2, x01, x02, y1, y2, y01, y02, emitnx_l, emitny_l,
     & fa, fb, fx, fy, fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     & fqxplot, fqyplot, fuxplot, fuyplot)
      use Constant
      use InGen
      use Lattice
      use LatticeInternal
      use Mult_data
      real(kind=8):: z,time,zion,aion,gamma,vz,ek,dz
      real(kind=8):: dedr,dexdx,deydy,dbdr,p,a,b,x,y,x1,x2,x01,x02,y1,y2,y01,y02
      real(kind=8):: emitnx_l,emitny_l 
      real(kind=8):: fa,fb,fx,fy,fx1,fx2,fx01,fx02,fy1,fy2,fy01,fy02,faccl
      real(kind=8):: fqxplot,fqyplot,fuxplot,fuyplot
      logical(ISZ):: lrelativ

c Computes the forces needed for the general-lattice envelope calculation.
c Inlcudes both residence corrected hard edged quadrupole elements 
c and hard-edge multipole elements, and elements defined by their 
c multipoles.  Only linear forces consistent with a truncated second-order 
c envelope formulation are calculated.  These include a linear applied 
c forces (electric and magnetic quadrupole fields, and a radial 
c E-field proportional to the radial coodrinate r), and 
c linear space-charge forces.  

c Note that this routine makes calls directly to getquadid etc instead of
c relying on a call to setlattzt. This is done since it is significantly
c faster, though this only really matters in cases where the envelope code is
c called many times. The increase is speed is essentially non-noticable
c if the envelope code is called once. It still uses the celemid arrays from
c LatticeInternal since the id's must be saved and they must be known for
c each level of overlap. It is easier to use those existing arrays rather
c than to create new local ones.

      real(kind=8):: z1,z2,zl,fracl,zr,fracr
      real(kind=8):: fuzz
      real(kind=8):: frac,ez,ke,de,db
      real(kind=8):: dedrr,dbdxq,dedxq,xdbdx,xdedx,ydbdx,ydedx,xex,yey,alpha
      real(kind=8):: fux,fuy,fq,feself,fexself,feyself 
      real(kind=8):: wie,wim
      real(kind=8):: offset
      integer(ISZ):: io
      integer(ISZ):: ii,iie,iim,iem,imm
      integer(ISZ):: iq,ih,ie,im,ia

      z1 = z - 0.5 * dz
      z2 = z + 0.5 * dz
      zl = min(z1,z2)
      zr = max(z1,z2)

c     --- Initialize variables to the values from the uniform fields.
      dedrr = dedr
      dbdxq = 0.
      dedxq = 0.
      xdbdx = dbdr*vz*x
      xdedx = (dedr + dexdx)*x
      ydbdx = dbdr*vz*y
      ydedx = (dedr + deydy)*y
      xex   = 0.
      yey   = 0.
      ez = 0.
      faccl = 0.

c     --- Hard edge quadrupoles

      if (quads) then

        do io=1,nquadol
          call getquadid(z,offset,cquadid(0,io),io)
          iq = cquadid(0,io)

c         --- Find quad at left edge of velocity advance step.
          fracl = 0.
          if (quadzs(iq)+offset < zl .and. zl < quadze(iq)+offset) fracl = 1.

c         --- Find quad at right edge of velocity advance step.
          fracr = 0.
          if (quadzs(iq)+offset < zr .and. zr < quadze(iq)+offset) fracr = 1.

c         --- Calculate residence fraction.
          frac = fracl
          if (fracl > fracr) frac = (quadze(iq)+offset - zl) / abs(dz)
          if (fracr > fracl) frac = (zr - quadzs(iq)-offset) / abs(dz)

c         --- quadrupole fields, including residence fraction
          if (frac > 0.) then
            dbdxq = dbdxq + quaddb(iq)*frac 
            dedxq = dedxq + quadde(iq)*frac 
            xdbdx = xdbdx + quaddb(iq)*frac*(x - qoffx(iq))
            xdedx = xdedx + quadde(iq)*frac*(x - qoffx(iq))
            ydbdx = ydbdx + quaddb(iq)*frac*(y - qoffy(iq))
            ydedx = ydedx + quadde(iq)*frac*(y - qoffy(iq))
          endif

        enddo

      endif 

c     --- Hard edge multipoles, can include linearly focusing 
c     --- quadrupole or azimuthally symmetric field 
 
      if (heles) then

        fuzz = 1.e-3
        do io=1,nheleol
          call getheleid(z,offset,cheleid(0,io),io)
          ih = cheleid(0,io)

c         --- Find multipole at left edge of velocity advance step.
          fracl = 0.
          if (zl < heleze(ih)+offset .and. helezs(ih)+offset < zl) fracl = 1.
          fracr = 0.
          if (zr < heleze(ih)+offset .and. helezs(ih)+offset < zr) fracr = 1.

c         --- Calculate residence fraction.
          frac = fracl
          if (fracl > fracr) frac = (heleze(ih)+offset - zl) / abs(dz)
          if (fracr > fracl) frac = (zr - helezs(ih)-offset) / abs(dz)
          if (frac == 0.) cycle

c         --- Get only the quadrupole and uniform components
          do ii = 1, max(helenm(ih),helene(ih)) 
            if ((nint(hele_n(ii,ih)) == 2) .and. 
     &          (nint(hele_v(ii,ih)) == 0)) then
c             --- quadrupole field  
              if ( abs(cos(helepm(ii,ih))*sin(helepm(ii,ih))) > fuzz ) 
     &         call remark("Warning: magnetic multipole structure incompatible")
              if ( abs(cos(helepe(ii,ih))*sin(helepe(ii,ih))) > fuzz ) 
     &         call remark("Warning: electric multipole structure incompatible")
              db = heleam(ii,ih)*cos(helepm(ii,ih))*frac
              de = heleae(ii,ih)*cos(helepe(ii,ih))*frac
              dbdxq = dbdxq + db
              dedxq = dedxq + de
              xdbdx = xdbdx + db*(x - heleox(ih))
              xdedx = xdedx + de*(x - heleox(ih))
              ydbdx = ydbdx + db*(y - heleoy(ih))
              ydedx = ydedx + de*(y - heleoy(ih))
            endif 
            if ((nint(hele_n(ii,ih)) == 0) .and.
     &          (nint(hele_v(ii,ih)) == 0)) then
c             --- uniform field
              dedrr = dedrr + 0.5*heleep(ii,ih)*frac
              xdedx = xdedx + 0.5*heleep(ii,ih)*frac*(x - heleox(ih))
              ydedx = ydedx + 0.5*heleep(ii,ih)*frac*(y - heleoy(ih))
            endif
          enddo

        enddo

      endif 

c     --- Extract any linear focusing field component from the multipole 
c     --- elements.

      if (emlts) then

        do io=1,nemltol
          call getemltid(z,offset,cemltid(0,io),io)
          ie = cemltid(0,io)

c         --- If z is within an element, fetch the field.
          if (emltzs(ie)+offset <= z .and. z < emltze(ie)+offset) then

c           --- Make sure the element was defined.
            iem = emltid(ie)
            if (iem == 0) cycle

c           --- Find location of z in the electric multipole grid.
            iie = int((z - emltzs(ie)-offset)/dzemlt(iem))
            wie =     (z - emltzs(ie)-offset)/dzemlt(iem) - iie

c           --- Lookup field strengths with linear interpolation.
            do ii=1,nesmult
              if ((nint(emlt_n(ii)) == 2) .and. 
     &            (nint(emlt_v(ii)) == 0)) then
                de = (esemlt(iie  ,ii,iem)*(1.-wie)+ esemlt(iie+1,ii,iem)*wie)*
     &               (emltsc(ie) + emltsf(ie))
                dedxq = dedxq + de
                xdedx = xdedx + de*(x - emltox(ie))
                ydedx = ydedx + de*(y - emltoy(ie))
              endif
              if ((nint(emlt_n(ii)) == 0) .and. 
     &            (nint(emlt_v(ii)) == 0)) then
                de = 0.5*(esemltp(iie  ,ii,iem)*(1.-wie) +
     &                    esemltp(iie+1,ii,iem)*wie)*
     &                   (emltsc(ie) + emltsf(ie))
                dedrr = dedrr - de
                xdedx = xdedx - de*(x - emltox(ie))
                ydedx = ydedx - de*(y - emltoy(ie))
                ez = ez + (esemlt(iie  ,ii,iem)*(1.-wie) +
     &                     esemlt(iie+1,ii,iem)*wie)*
     &                    (emltsc(ie) + emltsf(ie))
              endif
              if ((nint(emlt_n(ii)) == 1) .and.
     &            (nint(emlt_v(ii)) == 0)) then
                de = (esemlt(iie  ,ii,iem)*(1.-wie)+ esemlt(iie+1,ii,iem)*wie)*
     &               (emltsc(ie) + emltsf(ie))
                alpha = emltph(ie) +
     &                  esemltph(iie  ,ii,iem)*(1. - wie) +
     &                  esemltph(iie+1,ii,iem)*      wie
                xex = xex + de*cos(alpha)
                yey = yey + de*sin(alpha)
              endif
            enddo

          endif

        enddo

      endif

      if (mmlts) then

        do io=1,nmmltol
          call getmmltid(z,offset,cmmltid(0,io),io)
          im = cmmltid(0,io)

c         --- If z is within an element, fetch the field.
          if (mmltzs(im)+offset <= z .and. z < mmltze(im)+offset) then

c           --- Make sure the element was defined.
            imm = mmltid(im)
            if (imm == 0) cycle

c           --- Find location of z in the electric multipole grid.
            iim = int((z - mmltzs(im)-offset)/dzmmlt(imm))
            wim =     (z - mmltzs(im)-offset)/dzmmlt(imm) - iim

c           --- Lookup field strengths with linear interpolation.
            do ii=1,nmsmult
              if ((nint(mmlt_n(ii)) == 2) .and. 
     &            (nint(mmlt_v(ii)) == 0)) then
                db = (msmmlt(iim  ,ii,imm)*(1.-wim)+ msmmlt(iim+1,ii,imm)*wim)*
     &               (mmltsc(im) + mmltsf(im))
                dbdxq = dbdxq + db
                xdbdx = xdbdx + db*(x - mmltox(im))
                ydbdx = ydbdx + db*(y - mmltoy(im))
              endif
            enddo

          endif

        enddo

      endif

c     --- Get accelerating gradient from accl gaps.
      if (accls) then

        do io=1,nacclol
          call getacclid(z,offset,cacclid(0,io),io)
          ia = cacclid(0,io)

c         --- Find accl at left edge of velocity advance step.
          fracl = 0.
          if (zl < acclze(ia)+offset .and. acclzs(ia)+offset < zl) fracl = 1.

c         --- Find accl at right edge of velocity advance step.
          fracr = 0.
          if (zr < acclze(ia)+offset .and. acclzs(ia)+offset < zr) fracr = 1.

c         --- Calculate residence fraction.
          frac = fracl
          if (fracl > fracr) frac = (acclze(ia)+offset - zl) / abs(dz)
          if (fracr > fracl) frac = (zr - acclzs(ia)-offset) / abs(dz)

          ez = ez + acclez(ia)*frac

        enddo
      endif

c     --- Add in the acceleration (or deceleration)
      if (ez .ne. 0) then
        ek = ek + ez*dz

c       --- Remove old vz and gamma from the perveance
c       --- (Probably should completely recalculate the perveance, but
c       --- that requires adding several more variables to the argument
c       --- list.)
        p = p*(vz*gamma)**3

c       --- Calculate the new vz and gamma
        if (lrelativ) then
          ke = jperev*ek/dvnz(aion*amu*clight**2)
          gamma = 1. + ke
          vz = clight*sqrt((2*ke+ke**2)/gamma**2)
        else
          vz = sqrt(2.*ek*jperev/dvnz(aion)/amu)
        endif

c       --- Recalculate the perveance using update quantities
        p = p/(vz*gamma)**3

c       --- "Force" to account for the dvz/dz term which appears in the
c       --- conversion from d/dt to d/dz
        faccl = 0.5*ez/ek

      endif

c     --- Compute external uniform focusing (proportional to r and 
c     --- azimuthally symmetric) and quadrupole focusing strengths
      fux = zion*echarge*(dedrr + dexdx - dbdr*vz)/(aion*amu*gamma*vz**2)
      fuy = zion*echarge*(dedrr + deydy - dbdr*vz)/(aion*amu*gamma*vz**2)
      fq  = zion*echarge*(dedxq - dbdxq*vz)/(aion*amu*gamma*vz**2)

c     --- Calculate forces for advancing envelope and test particles
      feself = 2.*p/(a + b)
      fa = ( fq + fux)*a + feself + (emitnx_l*clight/(vz*gamma))**2/a**3
      fb = (-fq + fuy)*b + feself + (emitny_l*clight/(vz*gamma))**2/b**3
      fx = zion*echarge/(aion*amu*gamma*vz**2)*(xdedx - xdbdx*vz + xex)
      fy = zion*echarge/(aion*amu*gamma*vz**2)*(ydedx - ydbdx*vz + yey)

      fexself = feself/a
      fx1  = (fq + fux + fexself)*x1
      fx2  = (fq + fux + fexself)*x2
      fx01 = (fq + fux)*x01
      fx02 = (fq + fux)*x02

      feyself = feself/b
      fy1  = (-fq + fuy + feyself)*y1
      fy2  = (-fq + fuy + feyself)*y2
      fy01 = (fq + fuy)*y01
      fy02 = (fq + fuy)*y02

      fqxplot =  fq
      fqyplot = -fq
      fuxplot =  fux
      fuyplot =  fuy 

      return
      end
c=======================================================================
      subroutine envfofzinit()
      use ENVvars
      use ENVfofz
c Initialize variables, checking if any quantities are set to vary.
c Also, do some consistency checks.

c     --- Check if any quantities are set to vary.
c     --- Note that the exact same line should appear at the end of this
c     --- subroutine since there may be problems and some of the variations
c     --- may have been turned off.
      lefofz = libeame_z .or. lemitne_z

c     --- Do some checks on consistency.
      if (lefofz) then
c       --- First set mesh range to the default if it hasn't been set yet.
        if (zlefofz == zuefofz) then
          zlefofz = zl
          zuefofz = zu
        endif

c       --- Check value of dzefofz.  If not set, and if nzefofz is not
c       --- set, get value from dzenv, otherwize calculate new value.
        if (dzefofz <= 0.) then
          if (nzefofz <= 0) then
            dzefofz = dzenv
          else
            dzefofz = (zuefofz - zlefofz)/nzefofz
          endif
        endif

c       --- Make sure that nzefofz is set correctly.
        if (nzefofz <= 0) then
          nzefofz = nenv
        endif

c       --- Make sure that the arrays have been allocated.  If not,
c       --- print warning and turn variation off.

        if (libeame_z) then
          if (.not. ASSOCIATED(ibeame_z)) then
            call remark("ENV warning: ibeame_z has not been allocated")
            call remark("             current will not be varied")
            libeame_z = .false.
          endif
        endif

        if (lemitne_z) then
          if (.not. (ASSOCIATED(emitnxe_z) .or. ASSOCIATED(emitnye_z))) then
            call remark("ENV warning: emitnxe_z or emitnye_z not allocated")
            call remark("             normalized emittance will not be varied")
            lemitne_z = .false.
          endif
        endif

      endif

c     --- Repeat of first line of this subroutine.  See comment there.
      lefofz = libeame_z .or. lemitne_z

      return
      end
c=======================================================================
      subroutine envsetfofz(z,p,vz,emitnx_l,emitny_l,gamma)
      use Constant
      use Beam_acc
      use InGen
      use ENVvars
      use ENVfofz
      real(kind=8):: z,p,vz,emitnx_l,emitny_l,gamma

c Updates quantities which are specified to vary with z.


      integer(ISZ):: iz
      real(kind=8):: wz
      real(kind=8):: ibeam_z

c     --- If nothing is to vary, return.
      if (.not. lefofz) return

c     --- Get location in mesh
      iz = (z - zlefofz)/dzefofz
      wz = (z - zlefofz)/dzefofz - iz

c     --- If current location outside of mesh, return
      if (iz < 0 .or. nzefofz < iz) return

c     --- Update the current.
      ibeam_z = ibeam
      if (libeame_z) then
        if (iz < nzefofz) then
          ibeam_z = ibeame_z(iz)*(1.-wz)+ibeame_z(iz+1)*wz
        elseif (iz == nzefofz) then
          ibeam_z = ibeame_z(iz)
        endif
      endif

c     --- Update the normalized emittance
      if (lemitne_z) then
        if (iz < nzefofz) then
          emitnx_l = emitnxe_z(iz)*(1.-wz)+emitnxe_z(iz+1)*wz
          emitny_l = emitnye_z(iz)*(1.-wz)+emitnye_z(iz+1)*wz
        elseif (iz == nzefofz) then
          emitnx_l = emitnxe_z(iz)
          emitny_l = emitnye_z(iz)
        endif
      endif

c     --- Recalculate the perveance using update quantities
      p = ibeam_z*zion*echarge/(2.*pi*eps0*aion*amu*(vz*gamma)**3)

      return
      end
c=======================================================================
c Envelope match search routines
c=========================================================================
      subroutine mtchinit
      use Beam_acc
      use ENVvars
      use Match
c Initialize the envelope match search
      integer(ISZ):: i

c set up match arrays with initial deltas
      do i=1,7
        mtch(i,1) = a0
        mtch(i,2) = ap0
        mtch(i,3) = b0
        mtch(i,4) = bp0
        mtch(i,5) = ibeam
        mtch(i,6) = emit
      enddo
      mtch(2,1) = mtch(2,1) + mtch(2,1)*mtchdlta(1)
      mtch(3,2) = mtch(3,2) + mtch(3,2)*mtchdlta(2)
      mtch(4,3) = mtch(4,3) + mtch(4,3)*mtchdlta(3)
      mtch(5,4) = mtch(5,4) + mtch(5,4)*mtchdlta(4)
      mtch(6,5) = mtch(6,5) + mtch(6,5)*mtchdlta(5)
      mtch(7,6) = mtch(7,6) + mtch(7,6)*mtchdlta(6)

c calculate initial errors
c (uses same error expression as in the function mtchfunc)
      do i=1,7
        a0    = mtch(i,1)
        ap0   = mtch(i,2)
        b0    = mtch(i,3)
        bp0   = mtch(i,4)
        ibeam = mtch(i,5)
        emit  = mtch(i,6)
        call envexe()
        mtcherrs(i) = abs(deltaa)*wdeltaa + abs(deltaap)*wdeltaap +
     &                abs(deltab)*wdeltab + abs(deltabp)*wdeltabp +
     &                abs(sigx - sig_desr)*wsig
      enddo

      return
      end
c=========================================================================
      subroutine envmatch()
      use Match
c user interface to the envelope match search routine
      call amoeba(mtch,mtcherrs,7,6,6,mtch_tol,mtchiter)
      return
      end
c=========================================================================
      real(kind=8) function mtchfunc(xx)
      use Beam_acc
      use ENVvars
      use Match
c function that is called by the minimization routine, ameoba, to find
c a matched envelope solution
      real(kind=8):: xx(6)

c set envelope values and print them
      a0    = xx(1)
      ap0   = xx(2)
      b0    = xx(3)
      bp0   = xx(4)
      ibeam = xx(5)
      emit  = xx(6)
      print*, ' a0 = ', a0, '  b0 = ', b0
      print*, 'ap0 = ', ap0,' bp0 = ', bp0
      print*, 'ibeam = ', ibeam, ' emit = ', emit

c do envelope calculation
      call envexe()

c calculate error
      mtchfunc = abs(deltaa)*wdeltaa + abs(deltaap)*wdeltaap +
     &           abs(deltab)*wdeltab + abs(deltabp)*wdeltabp +
     &           abs(sigx - sig_desr)*wsig

      return
      end
c=========================================================================
c routine that does the envelope match search
c This is a multidimensional function minimization routine that was
c taken from Numerical Recipes, William H. Press et. al., Cambridge University
c Press, 1986, pg 292, Ch 10
      subroutine amoeba(p,y,mp,np,ndim,ftol,itmax)
      integer(ISZ):: mp,np,ndim,itmax
      real(kind=8):: ftol
      real(kind=8):: p(mp,np),y(mp)
      integer(ISZ):: nmax
      real(kind=8):: alpha,beta,gamma
      parameter (nmax=20,alpha=1.0,beta=0.5,gamma=2.0)
      real(kind=8):: pr(nmax),prr(nmax),pbar(nmax)
      integer(ISZ):: mpts,iter,ilo,ihi,inhi,i,j
      real(kind=8):: rtol,ypr,yprr
      real(kind=8):: mtchfunc
      mpts=ndim+1
      iter=0
1     ilo=1
      if(y(1) > y(2))then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 11 i=1,mpts
        if(y(i) < y(ilo)) ilo=i
        if(y(i) > y(ihi))then
          inhi=ihi
          ihi=i
        else if(y(i) > y(inhi))then
          if(i.ne.ihi) inhi=i
        endif
11    continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if(rtol < ftol) then
        print*,'Envelope match search reached tolerence'
        return
      endif
      if(iter == itmax) then
        print*,'Envelope match reached maximum iterations.'
        return
      endif
      iter=iter+1
      do 12 j=1,ndim
        pbar(j)=0.
12    continue
      do 14 i=1,mpts
        if(i.ne.ihi)then
          do 13 j=1,ndim
            pbar(j)=pbar(j)+p(i,j)
13        continue
        endif
14    continue
      do 15 j=1,ndim
        pbar(j)=pbar(j)/ndim
        pr(j)=(1.+alpha)*pbar(j)-alpha*p(ihi,j)
15    continue
      ypr=mtchfunc(pr)
      if(ypr <= y(ilo))then
        do 16 j=1,ndim
          prr(j)=gamma*pr(j)+(1.-gamma)*pbar(j)
16      continue
        yprr=mtchfunc(prr)
        if(yprr < y(ilo))then
          do 17 j=1,ndim
            p(ihi,j)=prr(j)
17        continue
          y(ihi)=yprr
        else
          do 18 j=1,ndim
            p(ihi,j)=pr(j)
18        continue
          y(ihi)=ypr
        endif
      else if(ypr >= y(inhi))then
        if(ypr < y(ihi))then
          do 19 j=1,ndim
            p(ihi,j)=pr(j)
19        continue
          y(ihi)=ypr
        endif
        do 21 j=1,ndim
          prr(j)=beta*p(ihi,j)+(1.-beta)*pbar(j)
21      continue
        yprr=mtchfunc(prr)
        if(yprr < y(ihi))then
          do 22 j=1,ndim
            p(ihi,j)=prr(j)
22        continue
          y(ihi)=yprr
        else
          do 24 i=1,mpts
            if(i.ne.ilo)then
              do 23 j=1,ndim
                pr(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=pr(j)
23            continue
              y(i)=mtchfunc(pr)
            endif
24        continue
        endif
      else
        do 25 j=1,ndim
          p(ihi,j)=pr(j)
25      continue
        y(ihi)=ypr
      endif
      go to 1
      end
c=========================================================================
      subroutine envtypetest(n1,n2)
      use TTester
      integer(ISZ):: n1,n2
      real(kind=8),pointer:: bp(:)
      TYPE(ENVtype),POINTER:: newe
      e1%oo = 1
      e2%oo = 1
      e1%nn = n1
      e2%nn = n2
      call gallot("TTester",0)
      e1%aa(1:3) = 2.
      e2%ee(2:4) = 1.
      allocate(bp(0:e1%nn))
      bp = 3.1416
      e1%bp => bp
      newe => NewENVType()
      newe%oo = 1
      newe%nn = 21
      call ENVtypeallot(newe)
      e1%ack => newe
      xxx => NewENVType()
      return
      end
      subroutine printxxx()
      use TTester
      print*,"xxx%oo = ",xxx%oo
      print*,"xxx%nn = ",xxx%nn
      print*,"xxx%ee = ",xxx%ee
      return
      end
      subroutine setzzz()
      use TTester
      real(kind=8),allocatable,target:: newzzz(:)
      allocate(newzzz(10))
      newzzz(3:5) = 18.
      zzz => newzzz
      return
      end
      subroutine printzzz()
      use TTester
      print*,"zzz = ",zzz
      NULLIFY(zzz)
      return
      end
      subroutine iqtest(iq)
      use ATester
      integer(ISZ):: iq
      real(kind=8),pointer:: nn(:)
      iqiq = iq
      call gsetdims("ATester")
      allocate(nn(iq))
      nn(iq) = 12345.
      xqxq => nn
      return
      end
