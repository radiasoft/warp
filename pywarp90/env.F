#include "top.h"
c=======================================================================
c@(#) File ENV.M, version $Revision: 3.39 $, $Date: 2008/07/07 19:59:51 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  It is the source file for the package ENV of the PIC code WARP,
c  but it may be useful by itself.  It consists of an simple 
c  envelope code and a Python interface.  
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
      integer function envxport(np,z,a,ap,b,bp,x,xp,y,yp,vz,
     &                          epsx,epsy,ib)
      use ENVvars
      use Beam_acc
      use ENVfofz
      use Constant
      integer(ISZ):: np
      real(kind=8):: z(np)
      real(kind=8):: a(np),ap(np),b(np),bp(np)
      real(kind=8):: x(np),xp(np),y(np),yp(np)
      real(kind=8):: vz(np),epsx(np),epsy(np),ib(np) 

c  Export routine from envelope package:
c      envelope parameters at a set of z values passed in
c  Called during the particle loading process of the particle code.
c  Returns 1 if any input z is outside the range over which the envelope
c      is defined.
c  Inputs:  np, z(1:np)
c  Outputs: a(1:np), ap(1:np), b(1:np), bp(1:np)
c           x(1:np), xp(1:np), y(1:np), yp(1:np)
c           vz(1:np), epsx(1:np), epsy(1:np), ib(1:np) 


      integer(ISZ):: i,ip
      real(kind=8):: zzmin,zzmax,w0,w1,gammabeta

c  Check if there is any envelope data
      if (.not. ASSOCIATED(aenv)) then
        call kaboom("Error: envelope code has not been intialized")
        return
      endif

c  Validity check on input

      zzmin = 1.e6
      zzmax = -1.e6
      do ip = 1,np
         zzmin = min(zzmin, z(ip))
         zzmax = max(zzmax, z(ip))
      enddo
      if (zzmin < zl .or. zzmax > zu) then
         envxport = 1
         return
      endif

c  Interpolate to get envelope values and as a first pass load constant 
c  values for the beam emittance and current 

      do ip=1,np
        i = int((z(ip)-zl)/dzenv)
        w1 = (z(ip)-zl)/dzenv - i
        w0 = 1. - w1
        a(ip)     = w0*aenv(i)  + w1*aenv(i+1)
        ap(ip)    = w0*apenv(i) + w1*apenv(i+1)
        b(ip)     = w0*benv(i)  + w1*benv(i+1)
        bp(ip)    = w0*bpenv(i) + w1*bpenv(i+1)
        x(ip)     = w0*xenv(i)  + w1*xenv(i+1)
        xp(ip)    = w0*xpenv(i) + w1*xpenv(i+1)
        y(ip)     = w0*yenv(i)  + w1*yenv(i+1)
        yp(ip)    = w0*ypenv(i) + w1*ypenv(i+1)
        vz(ip)    = w0*vzenv(i) + w1*vzenv(i+1)
        epsx(ip) = emitx 
        epsy(ip) = emity 
        ib(ip)   = ibeam 
      enddo 

c If the beam emittance and current is varied in s, reset the 
c exported values consistently

      if (lefofz) then
        do ip=1,np 
          i = int((z(ip)-zlefofz)/dzefofz)
          w1 = (z(ip)-zlefofz)/dzefofz - i
          w0 = 1. - w1
          if (lemitne_z) then 
            if (lrelativ) then
              gammabeta = (vz(ip)/clight)/sqrt( 1. - (vz(ip)/clight)**2 )
            else
              gammabeta = vz(ip)/clight 
            endif 
            epsx(ip) = ( w0*emitnxe_z(i) + w1*emitnxe_z(i+1) )/gammabeta
            epsy(ip) = ( w0*emitnye_z(i) + w1*emitnye_z(i+1) )/gammabeta
          endif 
          if (libeame_z) ib(ip) = w0*ibeame_z(i) + w1*ibeame_z(i+1) 
        enddo
      endif 

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
     & lrelativ, gammabar,
     & aion, zion, ekin, vbeam, ibeam, a0, b0, ap0, bp0, x0, y0, xp0, yp0,
     & xx0, xxp0, xpxp0, yy0, yyp0, ypyp0, xy0, xpy0, xyp0, xpyp0, bph0, 
     & emitx, emity, dedr, dexdx, deydy, dbdr)
 
      return
      end
c=======================================================================
      subroutine envelatt (
     & lrelativ, gammabar,
     & aion, zion, ekin, vbeam, ibeam, a0, b0, ap0, bp0, x0, y0, xp0, yp0,
     & xx0, xxp0, xpxp0, yy0, yyp0, ypyp0, xy0, xpy0, xyp0, xpyp0, bph0, 
     & emitx, emity, dedr, dexdx, deydy, dbdr)
      use ENVvars
      use ENVtune
      use Constant
      use InGen
      logical(ISZ):: lrelativ
      real(kind=8):: ibeam
      real(kind=8):: gammabar
      real(kind=8):: aion,zion,ekin,vbeam,a0,b0,ap0,bp0,x0,y0,xp0,yp0
      real(kind=8):: xx0,xxp0,xpxp0,yy0,yyp0,ypyp0,xy0,xpy0,xyp0,xpyp0,bph0
      real(kind=8):: emitx,emity 
      real(kind=8):: dedr,dexdx,deydy,dbdr

c  This routine does the actual envelope calculation for the general lattice. 
c  The input arguments are described in the variable description files
c  of packages ENV and TOP.  The package takes linear focusing components 
c  from various lattice and multipole arrays.  

      integer(ISZ):: istep
      real(kind=8):: mx11,mx22,mx110,mx220
      real(kind=8):: my11,my22,my110,my220
      real(kind=8):: p,gamma,emitnx_l, emitny_l
      real(kind=8):: dz,z,vz,ek,time
      real(kind=8):: a,ap,b,bp,fa,fb
      real(kind=8):: x,xp,fx,y,yp,fy,fxtrans,fytrans,flong
      real(kind=8):: x1,xp1,fx1,x2,xp2,fx2
      real(kind=8):: y1,yp1,fy1,y2,yp2,fy2
      real(kind=8):: x01,x0p1,fx01,x02,x0p2,fx02
      real(kind=8):: y01,y0p1,fy01,y02,y0p2,fy02
      real(kind=8):: xx,xxp,xpxp,yy,yyp,ypyp,xy,xpy,xyp,xpyp
      real(kind=8):: bph,bphtemp
      real(kind=8):: xpxpint,ypypint,xpypint
      real(kind=8):: faccl,fqxplot,fqyplot,fuxplot,fuyplot
      real(kind=8):: kfxx,kfxy,kfyx,kfyy,kxx,kxy,kyx,kyy
      real(kind=8):: cossig0x,cossigx,tansig0x,tansigx
      real(kind=8):: cossig0y,cossigy,tansig0y,tansigy
      real(kind=8):: wz
      real(kind=8):: t,s,xpm,ypm,xps,xpp,ypp
      real(kind=8):: somvmat(7,7),somvvec(7),somovec(7),somuvec(7)
      integer(ISZ):: indx(7)

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
c     --- Second-order beam moments w.r.t. beam centroid
      xx   = xx0
      xxp  = xxp0
      xpxp = xpxp0
      yy   = yy0
      yyp  = yyp0
      ypyp = ypyp0
      xy   = xy0
      xpy  = xpy0
      xyp  = xyp0
      xpyp = xpyp0
c     --- These integrals of second-order moments don't really have any 
c         physical significance.  They're just here to serve as analogues 
c         to position in the synchronized leapfrog advance.
      xpxpint = 0.
      ypypint = 0.
      xpypint = 0.
c     --- Transverse beam tilt angle
      bph  = bph0
c     --- Make sure that the beam tilt angle bph is consistent with 
c         xx, yy, and xy, and that the consistent angle bph is as close as 
c         possible to the initially specified bph0.
c         When running the program, try to ensure that the difference between 
c         bph0 and bphtemp is not an odd integer multiple of pi/4. radians, 
c         as that makes for ambiguous logic
      if (2.*xy*cos(2.*bph) .ne. (xx-yy)*sin(2.*bph)) then
        if (xx == yy) then
          if (xy == 0.) then
            bphtemp =  0.
          else if (xy > 0.) then
            bphtemp =  pi/4.
          else if (xy < 0.) then
            bphtemp = -pi/4.
          endif
        else 
          bphtemp = 0.5*atan(2.*xy/(xx-yy))
        endif
c       Make sure that the angle bph consistent with xx, yy, and xy is as 
c       close to the initially specified bph as possible.  The pi/2 terms 
c       appear because tan(2*x) has a periodicity of pi/2.
        bph = bphtemp + (pi/2.)*nint((bph-bphtemp)/(pi/2.))
      endif

c  Initialize facility and update any variables which are functions of z.

      call envfofzinit()
      call envsetfofz(z,p,vz,emitnx_l,emitny_l,gamma)

c  Obtain initial force

       call envflatt (z, time, zion, aion, lrelativ, gamma, vz, ek, dz,
     & dedr,dexdx,deydy,dbdr,p, 
     & a, b, x, y, xp, yp, x1, x2, x01, x02, y1, y2, y01, y02, xx, yy, xy, 
     & bph, emitnx_l, emitny_l,  
     & fa, fb, fx, fy, fxtrans, fytrans, flong, 
     & fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     & fqxplot, fqyplot, fuxplot, fuyplot,   
     & kfxx, kfxy, kfyx, kfyy, kxx, kxy, kyx, kyy )

c  Main loop of envelope code

      do istep = 0,nenv
c        --- store data for this step
         aenv(istep)    = a
         apenv(istep)   = ap
         benv(istep)    = b
         bpenv(istep)   = bp
         xenv(istep)    = x
         xpenv(istep)   = xp
         yenv(istep)    = y
         ypenv(istep)   = yp
         vzenv(istep)   = vz
         zenv(istep)    = z
         xorb(istep)    = x1
         xporb(istep)   = xp1
         yorb(istep)    = y1
         yporb(istep)   = yp1
         xxenv(istep)   = xx
         xxpenv(istep)  = xxp
         xpxpenv(istep) = xpxp
         yyenv(istep)   = yy
         yypenv(istep)  = yyp
         ypypenv(istep) = ypyp
         xyenv(istep)   = xy
         xpyenv(istep)  = xpy
         xypenv(istep)  = xyp
         xpypenv(istep) = xpyp
         fqxenv(istep)  = fqxplot
         fqyenv(istep)  = fqyplot
         fuxenv(istep)  = fuxplot
         fuyenv(istep)  = fuyplot
         bphenv(istep)  = bph
c        --- Calculate more complicated stored quantities
         rxenv(istep)     = 2.*sqrt(   xx*cos(bph)**2 + yy*sin(bph)**2 + 
     &                              2.*xy*sin(bph)*cos(bph)  )
         ryenv(istep)     = 2.*sqrt(   xx*sin(bph)**2 + yy*cos(bph)**2 - 
     &                              2.*xy*sin(bph)*cos(bph)  )
         emitxenv(istep)  = 4.*sqrt(xx*xpxp - xxp**2)
         emityenv(istep)  = 4.*sqrt(yy*ypyp - yyp**2)
         emitnxenv(istep) = 4.*sqrt(xx*xpxp - xxp**2)*gamma*vz/clight
         emitnyenv(istep) = 4.*sqrt(yy*ypyp - yyp**2)*gamma*vz/clight
         emitngenv(istep) = 4.*sqrt(0.5*(xx*xpxp - xxp**2 + yy*ypyp - yyp**2) +
     &                              (xy*xpyp - xyp*xpy))*gamma*vz/clight
         emitnhenv(istep) = 4.*( (xx*xpxp - xxp**2)*(yy*ypyp - yyp**2) + 
     &                           (xy*xpyp)**2 + (xyp*xpy)**2 - xx*yy*xpyp**2 - 
     &                           xx*ypyp*xpy**2 - xpxp*yy*xyp**2 - 
     &                           xpxp*ypyp*xy**2 - 2.*xy*xyp*xpy*xpyp + 
     &                           2.*xxp*ypyp*xy*xpy - 2.*xxp*yyp*xy*xpyp - 
     &                           2.*xxp*yyp*xyp*xpy + 2.*xpxp*yyp*xy*xyp + 
     &                           2.*xx*yyp*xpy*xpyp + 2.*xxp*yy*xpyp*xyp 
     &                         )**(1./4.)*gamma*vz/clight

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
c
c        For the case of solenoidal magnetic fields, all orbits and 
c        centroid advances are carried out implicitly in the Larmor 
c        frame.  This gives correct results for a uniform density beam 
c        with zero canonical angular momentum.  If WARP inputs are 
c        generalized to work with nonzero canonical angular momentum, then 
c        both the envelope and centroid advance equations will need to be 
c        generalized.  In this case it would probably be easier to just 
c        advance first and second order moments for a uniform density 
c        elliptical beam and then project out all quantities needed (envelope, 
c        centroid, emittance, etc.) from the statistical moments.   
c
c        --- envelope a in x-z plane
         ap   = ap   + .5*dz*fa   - 0.5*dz*faccl*ap
         a    = a    + ap*dz
c        --- envelope b in y-z plane
         bp   = bp   + .5*dz*fb   - 0.5*dz*faccl*bp
         b    = b    + bp*dz
c        --- envelope x and y centroids
c        Incorporate 1/4 of the transverse and acceleration terms
         xpm  = xp   + .25*dz*(fxtrans - faccl*xp)
         ypm  = yp   + .25*dz*(fytrans - faccl*yp)
c        Rotate angles via Boris method
         t    =        .25*dz*flong
         s    = 2.*t/(1. + t**2)
         xps  = xpm  + ypm*t
         ypp  = ypm  - xps*s
         xpp  = xps  + ypp*t
c        Incorporate 1/4 of the transverse and acceleration terms 
c        by inverting implicit equations
         xp   = (xpp + .25*dz*fxtrans)/(1. + 0.25*dz*faccl)
         yp   = (ypp + .25*dz*fytrans)/(1. + 0.25*dz*faccl)
c        Advance positions a full step
         x    = x    + xp*dz
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
c        --- Second-order moments
c        See Sven Chilton's MS thesis, "Implementation of an iterative 
c        matching scheme for the Kapchinskij-Vladimirskij envelope equations 
c        in the WARP code," for the single particle equations of motion from 
c        which the second order moment equations are derived.
c        Set up a matrix for advancing second order moment "velocities", somvmat
         somvmat(1,:) = (/ 1.-0.5*dz*faccl, 0.5*dz, 0., 0., 0., 0.5*dz*flong,
     &                     0. /)
         somvmat(2,:) = (/ dz*kxx, 1.-dz*faccl, 0., 0., dz*kxy, 0., dz*flong /)
         somvmat(3,:) = (/ 0., 0., 1.-0.5*dz*faccl, 0.5*dz, -0.5*dz*flong, 
     &                     0., 0. /)
         somvmat(4,:) = (/ 0., 0., dz*kyy, 1.-dz*faccl, 0., dz*kyx, -dz*flong /)
         somvmat(5,:) = (/ 0., 0., 0.5*dz*flong, 0., 1.-0.5*dz*faccl, 0., 
     &                     0.5*dz /)
         somvmat(6,:) = (/-0.5*dz*flong, 0., 0., 0., 0., 1.-0.5*dz*faccl, 
     &                     0.5*dz /)
         somvmat(7,:) = (/ 0.5*dz*kyx, -0.5*dz*flong, 0.5*dz*kxy, 0.5*dz*flong, 
     &                     0.5*dz*kyy, 0.5*dz*kxx, 1.-dz*faccl/)
c        Set up a second order moment "velocity" vector at step i, somvvec
         somvvec = (/ xxp, xpxp, yyp, ypyp, xpy, xyp, xpyp /)
c        Set up a second order moment "velocity" offset vector at step i, 
c        somovec
         somovec = 0.5*dz*(/ kxx*xx + kxy*xy, 0., kyx*xy + kyy*yy, 0., 
     &                       kxx*xy + kxy*yy, kyy*xy + kyx*xx, 0. /)
c        Multiply the velocity vector by the acceleration matrix
         somuvec = 0.
         call matvecprod(somvmat,somvvec,somuvec,7,7)
c        Update the velocity vector, i.e. advance it to step i+1/2
         somvvec = somuvec + somovec
c        Extract the second order moment velocity analogues from somvvec
         xxp  = somvvec(1)
         xpxp = somvvec(2)
         yyp  = somvvec(3)
         ypyp = somvvec(4)
         xpy  = somvvec(5)
         xyp  = somvvec(6)
         xpyp = somvvec(7)
c        Advance the second order position analogues a whole step 
c        (from step i to i+1)
         xx      = xx + 2.*xxp*dz
         xpxpint = xpxpint + xpxp*dz
         yy      = yy + 2.*yyp*dz
         ypypint = ypypint + ypyp*dz
         xy      = xy + (xpy+xyp)*dz
         xpypint = xpypint + xpyp*dz
c        --- Recalculate beam tilt angle bph
         if (xx == yy) then
           if (xy == 0.) then
             bphtemp =  0.
           else if (xy > 0.) then
             bphtemp =  pi/4.
           else if (xy < 0.) then
             bphtemp = -pi/4.
           endif
         else 
           bphtemp = 0.5*atan(2.*xy/(xx-yy))
         endif
         bph = bphtemp + (pi/2.)*nint((bph-bphtemp)/(pi/2.))
c        --- z coordinate
         z = z + dz
         time = time + dz/vz
c        --- Update any variables which are functions of z.
         call envsetfofz(z,p,vz,emitnx_l,emitny_l,gamma)
c        --- get force at new position
         call envflatt (z, time, zion, aion, lrelativ, gamma, vz, ek, dz,
     &    dedr,dexdx,deydy,dbdr,p,
     &    a, b, x, y, xp, yp, x1, x2, x01, x02, y1, y2, y01, y02, xx, yy, xy, 
     &    bph, emitnx_l, emitny_l, 
     &    fa, fb, fx, fy, fxtrans, fytrans, flong, 
     &    fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     &    fqxplot, fqyplot, fuxplot, fuyplot,  
     &    kfxx, kfxy, kfyx, kfyy, kxx, kxy, kyx, kyy )
c        --- complete the "velocity" advance steps
c        --- envelope angles
         ap   = (ap   + .5*dz*fa )/( 1. + 0.5*dz*faccl*ap)
         bp   = (bp   + .5*dz*fb )/( 1. + 0.5*dz*faccl*bp)
c        --- envelope x and y centroid angles
c        Incorporate 1/4 of the transverse and acceleration terms
         xpm  = xp   + .25*dz*(fxtrans - faccl*xp)
         ypm  = yp   + .25*dz*(fytrans - faccl*yp)
c        Rotate angles via Boris method
         t    =        .25*dz*flong
         s    = 2.*t/(1. + t**2)
         xps  = xpm  + ypm*t
         ypp  = ypm  - xps*s
         xpp  = xps  + ypp*t
c        Incorporate 1/4 of the transverse and acceleration terms 
c        by inverting implicit equations
         xp   = (xpp + .25*dz*fxtrans)/(1. + 0.25*dz*faccl)
         yp   = (ypp + .25*dz*fytrans)/(1. + 0.25*dz*faccl)
c        --- x-plane orbit angles with and without space-charge
         xp1  = (xp1  + .5*dz*fx1 )/(1. + 0.5*dz*faccl*xp1 )
         xp2  = (xp2  + .5*dz*fx2 )/(1. + 0.5*dz*faccl*xp2 )
         x0p1 = (x0p1 + .5*dz*fx01)/(1. + 0.5*dz*faccl*x0p1)
         x0p2 = (x0p2 + .5*dz*fx02)/(1. + 0.5*dz*faccl*x0p2)         
c        --- y-plane orbit angles with and without space-charge
         yp1  = (yp1  + .5*dz*fy1 )/(1. + 0.5*dz*faccl*yp1 )
         yp2  = (yp2  + .5*dz*fy2 )/(1. + 0.5*dz*faccl*yp2 )
         y0p1 = (y0p1 + .5*dz*fy01)/(1. + 0.5*dz*faccl*y0p1)
         y0p2 = (y0p2 + .5*dz*fy02)/(1. + 0.5*dz*faccl*y0p2)
c        --- Second-order moments
c        Set up a matrix, somvmat, to advance second order moment "velocities"
         somvmat(1,:) = (/ 1.+0.5*dz*faccl, -0.5*dz, 0., 0., 0., -0.5*dz*flong,
     &                     0. /)
         somvmat(2,:) = (/-dz*kxx, 1.+dz*faccl, 0., 0., -dz*kxy, 0., 
     &                    -dz*flong /)
         somvmat(3,:) = (/ 0., 0., 1.+0.5*dz*faccl, -0.5*dz, 0.5*dz*flong, 
     &                     0., 0. /)
         somvmat(4,:) = (/ 0., 0., -dz*kyy, 1.+dz*faccl, 0., -dz*kyx, 
     &                     dz*flong /)
         somvmat(5,:) = (/ 0., 0., -0.5*dz*flong, 0., 1.+0.5*dz*faccl, 0., 
     &                    -0.5*dz/)
         somvmat(6,:) = (/ 0.5*dz*flong, 0., 0., 0., 0., 1.+0.5*dz*faccl, 
     &                    -0.5*dz /)
         somvmat(7,:) = (/-0.5*dz*kyx, 0.5*dz*flong, -0.5*dz*kxy, 
     &                    -0.5*dz*flong, -0.5*dz*kyy, -0.5*dz*kxx, 1.+dz*faccl/)
c        Set up a second order moment "velocity" offset vector at step i+1, 
c        somovec
         somovec = 0.5*dz*(/ kxx*xx + kxy*xy, 0., kyx*xy + kyy*yy, 0., 
     &                       kxx*xy + kxy*yy, kyy*xy + kyx*xx, 0. /)
c        Offset the velocity vector somvvec.  Note that somvvec on the right 
c        side of the following equation is at step i+1/2, while somovec is at 
c        step i+1
         somvvec = somvvec + somovec
c        Solve for somvvec at step i+1 by performing an LU decomposition on 
c        somvmat and an LU backsubstitution on somvvec (offset)
         indx = 0
         call ludcmp(somvmat,7,7,indx,0.)
         call lubksb(somvmat,7,7,indx,somvvec)
c        The vector somvvec has now been advanced to step i+1
c        Extract the second order moment velocity analogues from somvvec
         xxp  = somvvec(1)
         xpxp = somvvec(2)
         yyp  = somvvec(3)
         ypyp = somvvec(4)
         xpy  = somvvec(5)
         xyp  = somvvec(6)
         xpyp = somvvec(7)
         
      enddo 

c     Compute the undepressed and depressed phase advances from 
c     the relevant transfer matrices, i.e. 
c     cos(sig) = 0.5*Trace(m) = 0.5*(C+S), 
c     where sig is the phase advance, m is the transfer matrix 
c     evaluated at the end of the lattice period, and C and S 
c     are, respectively, the cosine-like and sine-like principal 
c     orbit functions evaluated at the end of the lattice period.  

c     Notation: the transfer matrix element names obey the form 
c     mj##?, with 
c     m  = transfer matrix 
c     j  = x or y: denotes the plane
c     ## = 11 or 22: denotes the element of the transfer matrix
c     ?  = 0 or (blank): 0 denotes undepressed, (blank) denotes depressed
c     Phase advances are calculated over the interval [tunezs,tuneze].  
c     If [tunezs,tuneze] in contained within [zl,zu] AND tunezs == zl, 
c     mj11? = Cj? and mj22? = Sj?.  Otherwise, the transfer matrix 
c     elements are not as easily related to the standard cosine-like and 
c     sine-like principal orbits, but either way, 
c     cos(sigj?) = 0.5*(mj11? + mj22?) = 0.5*(Cj? + Sj?).

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
     & a, b, x, y, xp, yp, x1, x2, x01, x02, y1, y2, y01, y02, xx, yy, xy, 
     & bph, emitnx_l, emitny_l, 
     & fa, fb, fx, fy, fxtrans, fytrans, flong, 
     & fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     & fqxplot, fqyplot, fuxplot, fuyplot, 
     & kfxx, kfxy, kfyx, kfyy, kxx, kxy, kyx, kyy )
      use Constant
      use InGen
      use Lattice
      use LatticeInternal
      use Mult_data
      use SemiTransparentDisc
      use EGRDdata
      use BGRDdata
      use PGRDdata
      use ENVvars,Only:llarmorframe
      real(kind=8):: z,time,zion,aion,gamma,vz,ek,dz
      real(kind=8):: dedr,dexdx,deydy,dbdr,p
      real(kind=8):: a,b,x,y,xp,yp,x1,x2,x01,x02,y1,y2,y01,y02
      real(kind=8):: xx,xxp,xpxp,yy,yyp,ypyp,xy,xpy,xyp,xpyp,bph
      real(kind=8):: emitnx_l,emitny_l
      real(kind=8):: fa,fb,fx,fy,fxtrans,fytrans,flong
      real(kind=8):: fx1,fx2,fx01,fx02,fy1,fy2,fy01,fy02,faccl
      real(kind=8):: fqxplot,fqyplot,fuxplot,fuyplot
      real(kind=8):: kfxx,kfxy,kfyx,kfyy,kxx,kxy,kyx,kyy
      logical(ISZ):: lrelativ

c Computes the forces needed for the general-lattice envelope calculation.
c Inlcudes both residence corrected hard edged quadrupole elements 
c and hard-edge multipole elements, elements defined by their 
c multipoles, hard-edged dipoles, and gridded field elements.  Only linear 
c forces consistent with a truncated second-order 
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

c Dipole forces are not consistently implemented for all elements in the 
c present version.  This would require some work to make consistent with 
c the presence of any bends. 

c Sign conventions:
c    quadde positive => defocusing in x,   focusing in y 
c    quaddb positive =>   focusing in x, defocusing in y 
c 
c In centroid force accumulations, note that offsets to each element 
c may be different and thus the x and y coordinates with shifts must be 
c accumulated in each element.  Also for the CENTROID:
c   electric:
c    xdedx =  quadde*x = d E_x/dx * x     d E_x/dx = - d E_y/dy 
c    ydedy = -quadde*y = d E_y/dy * y 
c   magnetic:
c    xdbdx =  quaddb*x = d B_y/dx * x     d B_y/dx = d B_x/dy 
c    ydbdy =  quaddb*y = d B_x/dy * y     
c
c   F_x = const*( xdedx - ydbdy*v_z )    const = q/(m*gamma*v_z^2)
c   F_y = const*( ydedy + xdbdx*v_z ) 


      real(kind=8):: z1,z2,zl,fracl,zr,fracr
      real(kind=8):: fuzz
      real(kind=8):: frac,ez,bz,ke,de,db
      real(kind=8):: dedrr,dbdrr,dbdxq,dedxq
      real(kind=8):: xdbdx,xdedx,ydbdy,ydedy,bx,by,ex,ey,alpha
      real(kind=8):: exd,eyd,bxd,byd
      real(kind=8):: fux,fuy,fq,feself,fexself,feyself,fsol
      real(kind=8):: fexxself,feyyself
      real(kind=8):: wie,wim
      real(kind=8):: offset
      integer(ISZ):: io
      integer(ISZ):: ii,iie,iim,iem,imm
      integer(ISZ):: id,iq,ih,ie,im,ia
      real(kind=8):: pscale,maxpscale,sscale,maxsscale,maxzscale,rscale
      real(kind=8):: xpg(5),ypg(5),exg(5),eyg(5),ezg(5),bxg(5),byg(5),bzg(5)
      real(kind=8):: dxegrd,dyegrd,dxbgrd,dybgrd,dxpgrd,dypgrd
      real(kind=8):: dbdxqc2p,dbdxqs2p,dedxqc2p,dedxqs2p
      real(kind=8):: rx,ry,ksxx,ksxy,ksyx,ksyy
      real(kind=8):: lph

      z1 = z - 0.5 * dz
      z2 = z + 0.5 * dz
      zl = min(z1,z2)
      zr = max(z1,z2)

c     --- Initialize variables to the values from the uniform fields.
c     --- xdbdx and ydbdy had erroneous factors of vz
      dedrr = dedr
      dbdrr = dbdr
      dbdxq = 0.
      dedxq = 0.
      xdbdx = dbdr*x
      xdedx = (dedr + dexdx)*x
      ydbdy = dbdr*y
      ydedy = (dedr + deydy)*y
      bx = 0.
      by = 0.
      ex = 0.
      ey = 0.
      ez = 0.
      bz = bz0
      faccl = 0.
c     --- These quantities are related to the standard quadrupole field 
c         gradients.  The suffix c2p means the gradient is multiplied by the 
c         cosine of twice the skew coupling phase angle, while the suffix s2p 
c         means the gradient is multiplied by the sine of twice the skew 
c         coupling phase angle.
      dbdxqc2p = 0.
      dbdxqs2p = 0.
      dedxqc2p = 0.
      dedxqs2p = 0.

c     --- Hard edge dipoles

      if (dipos) then

        do io=1,ndipool
          call getdipoid(z,offset,cdipoid(0,io),io)
          id = cdipoid(0,io)

c         --- Find dipole at left edge of velocity advance step.
          fracl = 0.
          if (dipozs(id)+offset < zl .and. zl < dipoze(id)+offset) fracl = 1.

c         --- Find dipole at right edge of velocity advance step.
          fracr = 0.
          if (dipozs(id)+offset < zr .and. zr < dipoze(id)+offset) fracr = 1.

c         --- Calculate residence fraction.
          frac = fracl
          if (fracl > fracr) frac = (dipoze(id)+offset - zl) / abs(dz)
          if (fracr > fracl) frac = (zr - dipozs(id)-offset) / abs(dz)

c         --- dipole fields, including residence fraction
          if (frac > 0.) then
            bx = bx + dipobx(id)*frac
            ex = ex + dipoex(id)*frac
            by = by + dipoby(id)*frac
            ey = ey + dipoey(id)*frac
          endif

        enddo

      endif 

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
            dbdxqc2p = dbdxqc2p + quaddb(iq)*frac*cos(2.*quadph(iq))
            dbdxqs2p = dbdxqs2p + quaddb(iq)*frac*sin(2.*quadph(iq))
            dedxqc2p = dedxqc2p + quadde(iq)*frac*cos(2.*quadph(iq))
            dedxqs2p = dedxqs2p + quadde(iq)*frac*sin(2.*quadph(iq))
            
            dbdxq = dbdxq + quaddb(iq)*frac 
            dedxq = dedxq + quadde(iq)*frac
c            xdbdx = xdbdx + quaddb(iq)*frac*(x - qoffx(iq))
c            xdedx = xdedx + quadde(iq)*frac*(x - qoffx(iq))
c            xdbdx = xdbdx + quaddb(iq)*frac*(y - qoffx(iq))
c            xdedx = xdedx - quadde(iq)*frac*(y - qoffx(iq))
c           Electric and magnetic fields in lab frame: skew coupling is 
c           taken into account
            xdbdx = xdbdx + 
     &                quaddb(iq)*frac*( (x - qoffx(iq))*cos(2.*quadph(iq)) + 
     &                                  (y - qoffy(iq))*sin(2.*quadph(iq)) )
            xdedx = xdedx + 
     &                quadde(iq)*frac*( (x - qoffx(iq))*cos(2.*quadph(iq)) + 
     &                                  (y - qoffy(iq))*sin(2.*quadph(iq)) )
            ydbdy = ydbdy + 
     &                quaddb(iq)*frac*( (y - qoffy(iq))*cos(2.*quadph(iq)) - 
     &                                  (x - qoffx(iq))*sin(2.*quadph(iq)) )
            ydedy = ydedy - 
     &                quadde(iq)*frac*( (y - qoffy(iq))*cos(2.*quadph(iq)) - 
     &                                  (x - qoffx(iq))*sin(2.*quadph(iq)) )
          endif

        enddo

      endif 

c     --- Hard edge multipoles, can include linearly focusing 
c     --- quadrupole, solenoid, or azimuthally symmetric field 
 
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

c         --- Get only the quadrupole, solenoid, and uniform components
          do ii = 1, max(helenm(ih),helene(ih)) 
            if ((nint(hele_n(ii,ih)) == 2) .and. 
     &          (nint(hele_v(ii,ih)) == 0)) then
c             --- quadrupole field  
              if ( abs(cos(helepm(ii,ih))*sin(helepm(ii,ih))) > fuzz ) 
     &         call remark("Warning: magnetic multipole structure incompatible")
              if ( abs(cos(helepe(ii,ih))*sin(helepe(ii,ih))) > fuzz ) 
     &         call remark("Warning: electric multipole structure incompatible")
c              db = heleam(ii,ih)*cos(helepm(ii,ih))*frac
c              de = heleae(ii,ih)*cos(helepe(ii,ih))*frac
              db = heleam(ii,ih)*frac
              de = heleae(ii,ih)*frac
              dbdxqc2p = dbdxqc2p + db*cos(2.*helepm(ii,ih))
              dbdxqs2p = dbdxqs2p + db*sin(2.*helepm(ii,ih))
              dedxqc2p = dedxqc2p + de*cos(2.*helepe(ii,ih))
              dedxqs2p = dedxqs2p + de*sin(2.*helepe(ii,ih))
              dbdxq = dbdxq + db
              dedxq = dedxq + de
c              xdbdx = xdbdx + db*(x - heleox(ih))
c              xdedx = xdedx + de*(x - heleox(ih))
c              ydbdy = ydbdy + db*(y - heleoy(ih))
c              ydedy = ydedy - de*(y - heleoy(ih))
c             Electric and magnetic fields in lab frame: skew coupling is 
c             taken into account
              xdbdx = xdbdx + db*( (x - heleox(ih))*cos(2.*helepm(ii,ih)) + 
     &                             (y - heleoy(ih))*sin(2.*helepm(ii,ih)) )
              xdedx = xdedx + de*( (x - heleox(ih))*cos(2.*helepe(ii,ih)) + 
     &                             (y - heleoy(ih))*sin(2.*helepe(ii,ih)) )
              ydbdy = ydbdy + db*( (y - heleoy(ih))*cos(2.*helepm(ii,ih)) - 
     &                             (x - heleox(ih))*sin(2.*helepm(ii,ih)) )
              ydedy = ydedy - de*( (y - heleoy(ih))*cos(2.*helepe(ii,ih)) - 
     &                             (x - heleox(ih))*sin(2.*helepe(ii,ih)) )
            endif 
            if ((nint(hele_n(ii,ih)) == 1) .and. 
     &          (nint(hele_v(ii,ih)) == 0)) then
c             --- dipole field
              db  = heleam(ii,ih)*cos(helepm(ii,ih))*frac
              de  = heleae(ii,ih)*cos(helepe(ii,ih))*frac
              bx = bx + db
              ex = ex + de
              by = by + db
              ey = ey - de
            endif
            if ((nint(hele_n(ii,ih)) == 0) .and.
     &          (nint(hele_v(ii,ih)) == 0)) then
c             --- solenoidal field 
              bz = bz + heleam(ii,ih)*frac
            endif
            if ((nint(hele_n(ii,ih)) == 0) .and.
     &          (nint(hele_v(ii,ih)) == 0)) then
c             --- uniform field
              dedrr = dedrr - 0.5*heleep(ii,ih)*frac
              xdedx = xdedx - 0.5*heleep(ii,ih)*frac*(x - heleox(ih))
              ydedy = ydedy - 0.5*heleep(ii,ih)*frac*(y - heleoy(ih))
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
c               --- quadrupole field
                de = (esemlt(iie  ,ii,iem)*(1.-wie)+ esemlt(iie+1,ii,iem)*wie)*
     &               (emltsc(ie) + emltsf(ie))
                dedxq = dedxq + de
c                xdedx = xdedx + de*(x - emltox(ie))
c                ydedy = ydedy - de*(y - emltoy(ie))
                xdedx = xdedx + de*( (x - emltox(ih))*cos(2.*emltph(ie)) + 
     &                               (y - emltoy(ih))*sin(2.*emltph(ie)) )
                ydedy = ydedy - de*( (y - emltoy(ih))*cos(2.*emltph(ie)) - 
     &                               (x - emltox(ih))*sin(2.*emltph(ie)) )
                dedxqc2p = dedxqc2p + de*cos(2.*emltph(ie))
                dedxqs2p = dedxqs2p + de*sin(2.*emltph(ie))
              endif
              if ((nint(emlt_n(ii)) == 0) .and. 
     &            (nint(emlt_v(ii)) == 0)) then
c               --- uniform transverse field
                de = 0.5*(esemltp(iie  ,ii,iem)*(1.-wie) +
     &                    esemltp(iie+1,ii,iem)*wie)*
     &                   (emltsc(ie) + emltsf(ie))
                dedrr = dedrr - de
                xdedx = xdedx - de*(x - emltox(ie))
                ydedy = ydedy - de*(y - emltoy(ie))
c               --- longitudinal field
                ez = ez + (esemlt(iie  ,ii,iem)*(1.-wie) +
     &                     esemlt(iie+1,ii,iem)*wie)*
     &                    (emltsc(ie) + emltsf(ie))
              endif
              if ((nint(emlt_n(ii)) == 1) .and.
     &            (nint(emlt_v(ii)) == 0)) then
c               --- dipole field
                de = (esemlt(iie  ,ii,iem)*(1.-wie)+ esemlt(iie+1,ii,iem)*wie)*
     &               (emltsc(ie) + emltsf(ie))
                alpha = emltph(ie) +
     &                  esemltph(iie  ,ii,iem)*(1. - wie) +
     &                  esemltph(iie+1,ii,iem)*      wie
                ex = ex + de*cos(alpha)
                ey = ey + de*sin(alpha)
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
c               --- quadrupole field 
                db = (msmmlt(iim  ,ii,imm)*(1.-wim)+ msmmlt(iim+1,ii,imm)*wim)*
     &               (mmltsc(im) + mmltsf(im))
                dbdxq = dbdxq + db
c                xdbdx = xdbdx + db*(x - mmltox(im))
c                ydbdy = ydbdy + db*(y - mmltoy(im))
                xdbdx = xdbdx + db*( (x - mmltox(im))*cos(2.*mmltph(im)) + 
     &                               (y - mmltoy(im))*sin(2.*mmltph(im)) )
                ydbdy = ydbdy + db*( (y - mmltoy(im))*cos(2.*mmltph(im)) - 
     &                               (x - mmltox(im))*sin(2.*mmltph(im)) )
                dbdxqc2p = dbdxqc2p + db*cos(2.*mmltph(im))
                dbdxqs2p = dbdxqs2p + db*sin(2.*mmltph(im))
              endif
              if ((nint(mmlt_n(ii)) == 1) .and. 
     &            (nint(mmlt_v(ii)) == 0)) then
c               --- dipole field
                db = (msmmlt(iim  ,ii,imm)*(1.-wim)+ msmmlt(iim+1,ii,imm)*wim)*
     &               (mmltsc(im) + mmltsf(im))
                alpha = mmltph(im) +
     &                  msmmltph(iim  ,ii,imm)*(1. - wim) +
     &                  msmmltph(iim+1,ii,imm)*      wim
                bx = bx + db*cos(alpha)
                by = by + db*sin(alpha)
              endif
              if ((nint(mmlt_n(ii)) == 0) .and. 
     &            (nint(mmlt_v(ii)) == 0)) then
c               db = 0.5*(msmmltp(iim  ,ii,imm)*(1.-wim) +
c    &                    msmmltp(iim+1,ii,imm)*wim)*
c    &                   (mmltsc(im) + mmltsf(im))
c               dbdrr = dbdrr + db
c               --- solenoidal field
                bz = bz + (msmmlt(iim  ,ii,imm)*(1.-wim) +
     &                     msmmlt(iim+1,ii,imm)*wim)*
     &                    (mmltsc(im) + mmltsf(im))
              endif
            enddo

          endif

        enddo

      endif

c     --- Extract any linear focusing field component from the gridded 
c     --- elements.
      
      if (egrds) then
         
c       --- find index of nearest egrd element
c       --- The cegrdid array calculated here is needed for the magnetic
c       --- elements described by E field data on a 3-D grid.
        do io=1,negrdol
          call getelemid(z,offset,negrd,egrdzs,egrdze,egrdol,
     &                   oegrdnn(io),oegrdoi(oegrdii(io)),oegrdio,
     &                   cegrdid(0,io),io)
c         --- load lattice info into comoving 1d array
          cegrdzs(0,io) = egrdzs(cegrdid(0,io)) + offset
          cegrdze(0,io) = egrdze(cegrdid(0,io)) + offset
c         --- Check whether this grid point is in the element.
c         --- The grid points intentionally overlap to avoid roundoff errors.
          if (cegrdzs(0,io) <= z .and. z <= cegrdze(0,io))
     &      linegrd(io) = .true.
        enddo
        
c       --- set linegrd(0) to true if any other elements are true
        linegrd(0) = any(linegrd)
        
c       --- Determine the x- and y-spacing used to calculate the E-field 
c       --- gradients.  This may not be the best way to find the right 
c       --- spacing if there are overlapping egrd elements, but it seems 
c       --- fairly robust.
        dxegrd = minval(egrddx)
        dyegrd = minval(egrddy)
        
c       --- Calculate electric fields and their gradients at a point (x,y,z)
c       --- using the interpolated field strengths at 
c       --- ( x +- dxegrd, y +- dyegrd, z )

c       Index point: 1,          2,      3,     4,          5
        xpg = (/ x - dxegrd,     x,      x,     x,      x + dxegrd /)
        ypg = (/     y,      y - dyegrd, y, y + dyegrd,     y      /)

c       --- Index points for finite differencing are arranged as follows:
c                      y
c                      ^
c                      |
c                  -   |       (4)
c           dyegrd |   |        
c                  -   |  (1)  (3)  (5)
c                      |        
c                      |       (2)
c                      |        
c                       ------------------> x
c                          |----|
c                          dxegrd

c       --- Initialize the electric fields on the 5-point grid illustrated 
c           above.  The suffix g stands for grid.  These must be set to 0, 
c           because applyegrd() adds the gridded electric field components 
c           to its last three arguments.
        exg = 0.
        eyg = 0.
        ezg = 0.
c       --- Calculate electric field strengths at index points
        call applyegrd(5,xpg,ypg,1,z,.true.,exg,eyg,ezg)
c       --- Extract total electric field components at central point (x,y)
        ex = ex + exg(3)
        ey = ey + eyg(3)
        ez = ez + ezg(3)
c       --- Calculate field gradients via centered finite differencing
c       --- Uniform field gradient
        dedrr = dedrr + 0.5*( (exg(5) - exg(1))/(2.*dxegrd) + 
     &                        (eyg(4) - eyg(2))/(2.*dyegrd) )
c       --- Quadrupole field gradient
c           WARNING: This should be the quadrupole field gradient in the 
c           rotated lens frame.  At present, this formula is only correct 
c           if the lens rotation angle is 0.
        dedxq = dedxq + 0.5*( (exg(5) - exg(1))/(2.*dxegrd) - 
     &                        (eyg(4) - eyg(2))/(2.*dyegrd) )
c       --- Take skew coupling into account
        dedxqc2p = dedxqc2p + 0.5*( (exg(5) - exg(1))/(2.*dxegrd) - 
     &                              (eyg(4) - eyg(2))/(2.*dyegrd) )
        dedxqs2p = dedxqs2p + 0.5*( (eyg(5) - eyg(1))/(2.*dxegrd) + 
     &                              (exg(4) - exg(2))/(2.*dyegrd) )
        
      endif


      if (bgrds) then
         
c       --- find index of nearest bgrd element
c       --- The cbgrdid array calculated here is needed for the magnetic
c       --- elements described by B field data on a 3-D grid.
        do io=1,nbgrdol
          call getelemid(z,offset,nbgrd,bgrdzs,bgrdze,bgrdol,
     &                   obgrdnn(io),obgrdoi(obgrdii(io)),obgrdio,
     &                   cbgrdid(0,io),io)
c         --- load lattice info into comoving 1d array
          cbgrdzs(0,io) = bgrdzs(cbgrdid(0,io)) + offset
          cbgrdze(0,io) = bgrdze(cbgrdid(0,io)) + offset
c         --- Check whether this grid point is in the element.
c         --- The grid points intentionally overlap to avoid roundoff errors.
          if (cbgrdzs(0,io) <= z .and. z <= cbgrdze(0,io))
     &      linbgrd(io) = .true.
        enddo
        
c       --- set linbgrd(0) to true if any other elements are true
        linbgrd(0) = any(linbgrd)
                 
c       --- Determine the x- and y-spacing used to calculate the B-field 
c       --- gradients.  This may not be the best way to find the right 
c       --- spacing if there are overlapping egrd elements, but it seems 
c       --- fairly robust.
        dxbgrd = minval(bgrddx)
        dybgrd = minval(bgrddy)
        
c       --- Calculate electric fields and their gradients at a point (x,y,z)
c       --- using the interpolated field strengths at 
c       --- ( x +- dxbgrd, y +- dybgrd, z )

c       Index point: 1,          2,      3,     4,          5
        xpg = (/ x - dxbgrd,     x,      x,     x,      x + dxbgrd /)
        ypg = (/     y,      y - dybgrd, y, y + dybgrd,     y      /)

c       --- Index points for finite differencing are arranged as follows:
c                      y
c                      ^
c                      |
c                  -   |       (4)
c           dybgrd |   |        
c                  -   |  (1)  (3)  (5)
c                      |        
c                      |       (2)
c                      |        
c                       ------------------> x
c                          |----|
c                          dxbgrd

c       --- Initialize the magnetic fields on the 5-point grid illustrated 
c           above.  The suffix g stands for grid.  These must be set to 0, 
c           because applybgrd() adds the gridded electric field components 
c           to its last three arguments.        
        bxg = 0.
        byg = 0.
        bzg = 0.
c       --- Calculate magnetic field strengths at index points
        call applybgrd(5,xpg,ypg,1,z,.true.,bxg,byg,bzg)
c       --- Extract magnetic field componenents at central point (x,y)
        bx = bx + bxg(3)
        by = by + byg(3)
        bz = bz + bzg(3)
c       --- Calculate field gradients via centered finite differencing
c           WARNING: This should be the quadrupole field gradient in the 
c           rotated lens frame.  At present, this formula is only correct 
c           if the lens rotation angle is 0.
        dbdxq = dbdxq + 0.5*( (byg(5) - byg(1))/(2.*dxbgrd) + 
     &                        (bxg(4) - bxg(2))/(2.*dybgrd) )
c       --- Take skew coupling into account
        dbdxqc2p = dbdxqc2p + 0.5*( (byg(5) - byg(1))/(2.*dxbgrd) + 
     &                              (bxg(4) - bxg(2))/(2.*dybgrd) )
        dbdxqs2p = dbdxqs2p + 0.5*( (byg(4) - byg(2))/(2.*dybgrd) - 
     &                              (bxg(5) - bxg(1))/(2.*dxbgrd) )
        
      endif


      if (pgrds) then
         
c       --- find index of nearest pgrd element
c       --- The cpgrdid array calculated here is needed for the electrostatic
c       --- elements described by potential data on a 3-D grid.
        do io=1,npgrdol
          call getelemid(z,offset,npgrd,pgrdzs,pgrdze,pgrdol,
     &                   opgrdnn(io),opgrdoi(opgrdii(io)),opgrdio,
     &                   cpgrdid(0,io),io)
c         --- load lattice info into comoving 1d array
          cpgrdzs(0,io) = pgrdzs(cpgrdid(0,io)) + offset
          cpgrdze(0,io) = pgrdze(cpgrdid(0,io)) + offset
c         --- Check whether this grid point is in the element.
c         --- The grid points intentionally overlap to avoid roundoff errors.
          if (cpgrdzs(0,io) <= z .and. z <= cpgrdze(0,io))
     &      linpgrd(io) = .true.
        enddo
        
c       --- set linpgrd(0) to true if any other elements are true
        linpgrd(0) = any(linpgrd)
        
c       --- Determine the x- and y-spacing used to calculate the E-field 
c       --- gradients.  This may not be the best way to find the right 
c       --- spacing if there are overlapping pgrd elements, but it seems 
c       --- fairly robust.
        dxpgrd = minval(pgrddx)
        dypgrd = minval(pgrddy)
        
c       --- Calculate electric fields and their gradients at a point (x,y,z)
c       --- using the interpolated field strengths at 
c       --- ( x +- dxpgrd, y +- dypgrd, z )

c       Index point: 1,          2,      3,     4,          5
        xpg = (/ x - dxpgrd,     x,      x,     x,      x + dxpgrd /)
        ypg = (/     y,      y - dypgrd, y, y + dypgrd,     y      /)
        
c       --- Index points for finite differencing are arranged as follows:
c                      y
c                      ^
c                      |
c                  -   |       (4)
c           dypgrd |   |        
c                  -   |  (1)  (3)  (5)
c                      |        
c                      |       (2)
c                      |        
c                       ------------------> x
c                          |----|
c                          dxpgrd

c       --- Initialize the electric fields on the 5-point grid illustrated 
c           above.  The suffix g stands for grid.  These must be set to 0, 
c           because applypgrd() adds the gridded electric field components 
c           to its last three arguments.        
        exg = 0.
        eyg = 0.
        ezg = 0.
c       --- Calculate electric field strengths at index points
        call applypgrd(5,xpg,ypg,1,z,.true.,exg,eyg,ezg)
c       --- Extract electric field components at central point (x,y)
        ex = ex + exg(3)
        ey = ey + eyg(3)
        ez = ez + ezg(3)
c       --- Calculate field gradients via centered finite differencing
c       --- Uniform field gradient
        dedrr = dedrr + 0.5*( (exg(5) - exg(1))/(2.*dxpgrd) + 
     &                        (eyg(4) - eyg(2))/(2.*dypgrd) )
c       --- Quadrupole field gradient
c           WARNING: This should be the quadrupole field gradient in the 
c           rotated lens frame.  At present, this formula is only correct 
c           if the lens rotation angle is 0.
        dedxq = dedxq + 0.5*( (exg(5) - exg(1))/(2.*dxpgrd) - 
     &                        (eyg(4) - eyg(2))/(2.*dypgrd) )
c       --- Take skew coupling into account
        dedxqc2p = dedxqc2p + 0.5*( (exg(5) - exg(1))/(2.*dxegrd) - 
     &                              (eyg(4) - eyg(2))/(2.*dyegrd) )
        dedxqs2p = dedxqs2p + 0.5*( (eyg(5) - eyg(1))/(2.*dxegrd) + 
     &                              (exg(4) - exg(2))/(2.*dyegrd) )
        
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

c     --- Apply scaling on the perveance to model the shorting of the self
c     --- fields near a conducting plate. This, for example, is used to correct
c     --- the envelope calculation as the beam impinges on the slit plate
c     --- of an emittance scanner.
      pscale = 1.
      if (n_STdiscs > 0) then
        do ii=1,n_STdiscs
c         --- Skip this disk if the fix is not to be applied.
          if (.not. lenvfix_STdisc(ii)) cycle
c         --- Calculate the z range over which the fix is to be applied
          rscale = sqrt(a*b)
          maxzscale = max(envfixzscale_STdisc(ii),
     &                    envfixsscale_STdisc(ii)*sqrt(a*b))
c         --- If the z is outside that range, then skip this disk
          if (abs(z - z_STdiscs(ii)) > maxzscale) cycle
c         --- Apply the fix. Note that fixes from multipe disks will be
c         --- multiplied.
c         --- The fix is scaled by the scale factor at the max distance
c         --- away from the plate that the fix is applied. This is done so
c         --- at the point, the scale factor converges to one.
          maxsscale = maxzscale/rscale
          maxpscale = maxsscale/sqrt(1.+maxsscale**2)
          sscale = abs(z - z_STdiscs(ii))/rscale
          pscale = pscale*(sscale/sqrt(1.+sscale**2))/maxpscale
        enddo
      endif

c     --- Calculate beam edge radii in the x- and y-planes, respectively
c         Note: these radii are in the tilted beam coordinates
      rx = 2.*sqrt(xx*cos(bph)**2 + yy*sin(bph)**2 + 2.*xy*sin(bph)*cos(bph))
      ry = 2.*sqrt(xx*sin(bph)**2 + yy*cos(bph)**2 - 2.*xy*sin(bph)*cos(bph))

c     --- Compute external uniform focusing (proportional to r and 
c     --- azimuthally symmetric) and quadrupole focusing strengths
c     --- Changed sign in front of dbdrr in fuy from - to +
c     --- Changed it back
      fux  = zion*echarge*(dedrr + dexdx - dbdrr*vz)/(aion*amu*gamma*vz**2)
      fuy  = zion*echarge*(dedrr + deydy - dbdrr*vz)/(aion*amu*gamma*vz**2)
      fq   = zion*echarge*(dedxq - dbdxq*vz)/(aion*amu*gamma*vz**2)
      fsol = (zion*echarge*bz/(2.*aion*amu*gamma*vz))**2
     
c     --- Calculate forces for advancing envelope, centroid, and test particles
c     --- NOTE: For the moment, forces calculated do not apply to mixed-type 
c         lattices, e.g. lattices with both quadrupole and solenoid focusing 
c         elements.  Will fix this later.
      feself = 2.*p*pscale/(a + b)
      fa = ( fq + fux - fsol)*a + feself + (emitnx_l*clight/(vz*gamma))**2/a**3
      fb = (-fq + fuy - fsol)*b + feself + (emitny_l*clight/(vz*gamma))**2/b**3
      fx = zion*echarge/(aion*amu*gamma*vz**2)*
     &      (xdedx - xdbdx*vz + yp*bz*vz + ex - by*vz)
      fy = zion*echarge/(aion*amu*gamma*vz**2)*
     &      (ydedy + ydbdy*vz - xp*bz*vz + ey + bx*vz)

c     These forces are used to advance the beam centroid coordinates and 
c     angles.  fxtrans and fytrans are due to transverse electric and magnetic 
c     field components, while flong stems from a longitudinal magnetic field.
      fxtrans = zion*echarge/(aion*amu*gamma*vz**2)*(xdedx-xdbdx*vz+ex-by*vz)
      fytrans = zion*echarge/(aion*amu*gamma*vz**2)*(ydedy+ydbdy*vz+ey+bx*vz)
      flong   = zion*echarge/(aion*amu*gamma*vz)*bz
      
      fexself = feself/a
c     Minus signs before the solenoid terms
      fx1  = (fq + fux - fsol + fexself)*x1 
      fx2  = (fq + fux - fsol + fexself)*x2
      fx01 = (fq + fux - fsol)*x01
      fx02 = (fq + fux - fsol)*x02
      
      feyself = feself/b
c     Minus signs before the solenoid and quadrupole terms
      fy1  = (-fq + fuy - fsol + feyself)*y1 
      fy2  = (-fq + fuy - fsol + feyself)*y2
      fy01 = (-fq + fuy - fsol)*y01
      fy02 = (-fq + fuy - fsol)*y02
      
      fqxplot =  fq
      fqyplot = -fq
      fuxplot =  fux
      fuyplot =  fuy 
      
c     These forces are useful in the second-order moment advance
c     --- Make sure that there is no divide by zero.
      if (rx == 0.) then
        fexxself = 0.
      else
        fexxself = 2.*p*pscale/((rx+ry)*rx)
      endif
      if (ry == 0.) then
        feyyself = 0.
      else
        feyyself = 2.*p*pscale/((rx+ry)*ry)
      endif
      
c     Define "k-function" force terms used in the second order moment 
c     advances.  Notation is similar to that used in Sven Chilton's MS thesis, 
c     "Implementation of an iterative matching scheme for the Kapchinskij-
c     Vladimirskij envelope equations in the WARP code."
c     "kf" terms stem from applied focusing effects, while "ks" terms stem 
c     from self-fields.
c     Caution: the "kf" terms here have opposite signs from those defined 
c     in the thesis.
      ksxx = fexxself*cos(bph)**2 + feyyself*sin(bph)**2
      ksxy = (fexxself-feyyself)*sin(bph)*cos(bph)
      ksyx = ksxy
      ksyy = fexxself*sin(bph)**2 + feyyself*cos(bph)**2
      kfxx = zion*echarge*(dedrr+dedxqc2p-dbdxqc2p*vz)/(aion*amu*gamma*vz**2)
      kfxy = zion*echarge*(dedxqs2p - dbdxqs2p*vz)/(aion*amu*gamma*vz**2)
      kfyx = kfxy
      kfyy = zion*echarge*(dedrr-dedxqc2p+dbdxqc2p*vz)/(aion*amu*gamma*vz**2)
      kxx  = ksxx + kfxx
      kxy  = ksxy + kfxy
      kyx  = ksyx + kfyx
      kyy  = ksyy + kfyy
      
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
      subroutine RK4HillSolve(karray,si,sf,xi,xpi,numsteps,xxparray)
c     Solves Hill's Equation x''(s) + k(s) = 0, x(si) = xi, x'(si) = xpi 
c     on interval [si,sf].
c     Inputs:
c       karray: array of k values with dimension 2*numsteps+1
c       si: initial s value
c       sf: final s value, usually equal to si+lperiod
c       xi: initial x value
c       xpi: initial x' value
c       numsteps: number of evenly spaced subintervals on [si,sf]
c       xxparray: dummy array to be modified
c     Output:
c       xxparray: array with dimensions (numsteps+1,2) containing x values
c       in column 1 and xp values in column 2
      integer(ISZ):: numsteps,i
      real(kind=8):: si,sf,xi,xpi,ds
      real(kind=8):: karray(2*numsteps+1),xxparray(numsteps+1,2)
      real(kind=8):: k1x,k1xp,k2x,k2xp,k3x,k3xp,k4x,k4xp
c     calculate step size 
      ds = (sf-si)/numsteps
c     set initial conditions
c     xxparray(1,1) = xi  = initial x  value
c     xxparray(1,2) = xpi = initial xp value
      xxparray(1,1) = xi
      xxparray(1,2) = xpi
      do i=1,numsteps
         k1x  = ds*xxparray(i,2)
         k1xp = ds*(-karray(2*i-1))*xxparray(i,1)
         k2x  = ds*(xxparray(i,2)+k1xp/2.)
         k2xp = ds*(-karray(2*i))*(xxparray(i,1)+k1x/2.)
         k3x  = ds*(xxparray(i,2)+k2xp/2.)
         k3xp = ds*(-karray(2*i))*(xxparray(i,1)+k2x/2.)
         k4x  = ds*(xxparray(i,2)+k3xp)
         k4xp = ds*(-karray(2*i+1))*(xxparray(i,1)+k3x)
         
         xxparray(i+1,1) = xxparray(i,1) + (k1x+ 2.*k2x+ 2.*k3x+ k4x)/6.
         xxparray(i+1,2) = xxparray(i,2) + (k1xp+2.*k2xp+2.*k3xp+k4xp)/6.
      end do
      return
      end
c========================================================================
      function kappax(z)
c     x-plane lattice focusing function
      use Picglb
      use Beam_acc
      real(kind=8):: z,dz
      real(kind=8):: fa,fb,fx,fy,fxtrans,fytrans,flong
      real(kind=8):: fx1,fx2,fx01,fx02,fy1,fy2,fy01,fy02,faccl
      real(kind=8):: fqxplot,fqyplot,fuxplot,fuyplot
      real(kind=8):: kfxx,kfxy,kfyx,kfyy,kxx,kxy,kyx,kyy
      real(kind=8):: kappax
      dz = SMALLPOS
      call envflatt(z,time,zion,aion,lrelativ,gammabar,vbeam,ekin,dz,
     & dedr,dexdx,deydy,dbdr,0.,
     & 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 
     & 0., 0., 0.,
     & fa, fb, fx, fy, fxtrans, fytrans, flong, 
     & fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     & fqxplot, fqyplot, fuxplot, fuyplot, 
     & kfxx, kfxy, kfyx, kfyy, kxx, kxy, kyx, kyy )
      kappax = -fa
      return
      end
c========================================================================
      function kappay(z)
c     y-plane lattice focusing function
      use Picglb
      use Beam_acc
      real(kind=8):: z,dz
      real(kind=8):: fa,fb,fx,fy,fxtrans,fytrans,flong
      real(kind=8):: fx1,fx2,fx01,fx02,fy1,fy2,fy01,fy02,faccl
      real(kind=8):: fqxplot,fqyplot,fuxplot,fuyplot
      real(kind=8):: kfxx,kfxy,kfyx,kfyy,kxx,kxy,kyx,kyy
      real(kind=8):: kappay
      dz = SMALLPOS
      call envflatt(z,time,zion,aion,lrelativ,gammabar,vbeam,ekin,dz,
     & dedr,dexdx,deydy,dbdr,0.,
     & 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 
     & 0., 0., 0., 
     & fa, fb, fx, fy, fxtrans, fytrans, flong, 
     & fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     & fqxplot, fqyplot, fuxplot, fuyplot, 
     & kfxx, kfxy, kfyx, kfyy, kxx, kxy, kyx, kyy )
      kappay = -fb
      return
      end
c========================================================================
      subroutine kappaxvec(zarray,karray,n)
c     Calculates kappax at each point in zarray and stores those values in 
c     karray
      use Picglb
      use Beam_acc
      real(kind=8):: z,dz
      real(kind=8):: fa,fb,fx,fy,fxtrans,fytrans,flong
      real(kind=8):: fx1,fx2,fx01,fx02,fy1,fy2,fy01,fy02,faccl
      real(kind=8):: fqxplot,fqyplot,fuxplot,fuyplot
      real(kind=8):: kfxx,kfxy,kfyx,kfyy,kxx,kxy,kyx,kyy
      integer(ISZ):: i,n
      real(kind=8):: zarray(n),karray(n)
      dz = SMALLPOS
      do i=1,n
        call envflatt(zarray(i),time,zion,aion,lrelativ,gammabar,vbeam,ekin,dz,
     & dedr,dexdx,deydy,dbdr,0.,
     & 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 
     & 0., 0., 0., 
     & fa, fb, fx, fy, fxtrans, fytrans, flong, 
     & fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     & fqxplot, fqyplot, fuxplot, fuyplot, 
     & kfxx, kfxy, kfyx, kfyy, kxx, kxy, kyx, kyy )
        karray(i) = -fa 
      end do
      return
      end
c========================================================================
      subroutine kappayvec(zarray,karray,n)
c     Calculates kappay at each point in zarray and stores those values in 
c     karray
      use Picglb
      use Beam_acc
      real(kind=8):: z,dz
      real(kind=8):: fa,fb,fx,fy,fxtrans,fytrans,flong
      real(kind=8):: fx1,fx2,fx01,fx02,fy1,fy2,fy01,fy02,faccl
      real(kind=8):: fqxplot,fqyplot,fuxplot,fuyplot
      real(kind=8):: kfxx,kfxy,kfyx,kfyy,kxx,kxy,kyx,kyy
      integer(ISZ):: i,n
      real(kind=8):: zarray(n),karray(n)
      dz = SMALLPOS
      do i=1,n
        call envflatt(zarray(i),time,zion,aion,lrelativ,gammabar,vbeam,ekin,dz,
     & dedr,dexdx,deydy,dbdr,0.,
     & 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 0., 
     & 0., 0., 0., 
     & fa, fb, fx, fy, fxtrans, fytrans, flong, 
     & fx1, fx2, fx01, fx02, fy1, fy2, fy01, fy02, faccl,
     & fqxplot, fqyplot, fuxplot, fuyplot, 
     & kfxx, kfxy, kfyx, kfyy, kxx, kxy, kyx, kyy )
        karray(i) = -fb 
      end do
      return
      end
c========================================================================
      subroutine ludcmp(a,n,np,indx,d)
      integer(ISZ):: n,np,indx(n)
      real(kind=8):: a(np,np),d

c     Perform LU decomposition on the matrix a
c     This is useful for determining the inverse of a

      integer(ISZ):: i,imax,j,k,nmax
      real(kind=8):: aamax,dum,sum,tiny
      parameter (nmax=500,tiny=SMALLPOS)
      real(kind=8):: vv(nmax)
      
      d = 1.
      do i=1,n
        aamax = 0.
        do j=1,n
          if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))
        end do
        if (aamax == 0.) then
          call kaboom("Singular matrix in ludcmp")
        end if
        vv(i) = 1./aamax
      end do
      do j=1,n
        do i=1,j-1
          sum = a(i,j)
          do k=1,i-1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
        end do
        aamax = 0.
        do i=j,n
          sum = a(i,j)
          do k=1,j-1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
          dum = vv(i)*abs(sum)
          if (dum >= aamax) then
            imax  = i
            aamax = dum
          end if
        end do
        if (j .ne. imax) then
          do k=1,n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          end do
          d = -d
          vv(imax) = vv(j)
        end if
        indx(j) = imax
        if (a(j,j) == 0.) a(j,j) = tiny
        if (j .ne. n) then
          dum = 1./a(j,j)
          do i=j+1,n
            a(i,j) = a(i,j)*dum
          end do
        end if
      end do
      return
      end
c========================================================================
      subroutine lubksb(a,n,np,indx,b)
      integer(ISZ):: n,np,indx(n)
      real(kind=8):: a(np,np),b(n)
      
c     Solves the set of n linear equations a \cdot x = b
c     Here, a is the LU decomposition of the original matrix a
c     as determined by the subroutine ludcmp.
c     Taken from Numerical Recipes
      
      integer(ISZ):: i,ii,j,ll
      real(kind=8):: sum
      
      ii = 0
      do i=1,n
        ll    = indx(i)
        sum   = b(ll)
        b(ll) = b(i)
        if (ii .ne. 0) then
          do j=ii,i-1
            sum = sum-a(i,j)*b(j)
          end do
        else if (sum .ne. 0.) then
          ii = 1
        end if
        b(i) = sum
      end do
      do i=n,1,-1
        sum = b(i)
        do j=i+1,n
          sum = sum-a(i,j)*b(j)
        end do
        b(i) = sum/a(i,i)
      end do
      return
      end
c========================================================================
      subroutine matinv(a,n,np,indx,d,y)
      integer(ISZ):: n,np,indx(n)
      real(kind=8):: a(np,np),d,y(np,np)
      
c     Find the inverse of the matrix a
c     Adapted from Numerical Recipes
      
      integer(ISZ):: i,j
      
      do i=1,n
        do j=1,n
          y(i,j) = 0.
        end do
        y(i,i) = 1.
      end do
      call ludcmp(a,n,np,indx,d)
      do j=1,n
        call lubksb(a,n,np,indx,y(:,j))
      end do
      return
      end
c========================================================================
      subroutine matprod(a,b,c,i,j,k)
      integer(ISZ):: i,j,k
      real(kind=8):: a(i,j),b(j,k),c(i,k)
c     Compute the product of two matrices a and b, 
c     with respective dimensions i by j and j by k
      integer(ISZ):: ii,jj,kk
      real(kind=8):: dum(j)
      
      do ii=1,i
        do kk=1,k
          do jj=1,j
            dum(jj) = a(ii,jj)*b(jj,kk)
          enddo
          c(ii,kk)  = sum(dum)
        enddo
      enddo
      return
      end
c========================================================================
      subroutine matvecprod(a,b,c,i,j)
      integer(ISZ):: i,j
      real(kind=8):: a(i,j),b(j),c(i)
c     Compute the matrix product of an i by j matrix a and a length j 
c     vector b
      integer(ISZ):: ii,jj
      real(kind=8):: dum(j)
            
      do ii=1,i
        do jj=1,j
          dum(jj) = a(ii,jj)*b(jj)
        enddo
        c(ii) = sum(dum)
      enddo
      
      return
      end
