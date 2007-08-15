#include "top.h"
c============================================================================
c@(#) File DTOP.F, version $Revision: 3.156 $, $Date: 2007/08/15 18:34:28 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is part of the main package TOP of code WARP
c  This file contains the diagnostic plotting routines.
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      subroutine minidiag (itt,tim,lspecial)
      use InDiag
      use Io
      use Z_arrays
      use Win_Moments
      use Z_Moments
      use Hist
      use Moments
      use Parallel,Only: my_index
      integer(ISZ):: itt
      real(kind=8):: tim
      logical(ISZ):: lspecial

c  Printout particle moments and field qty, writes into history buffers, etc.
c  Send one-liner to print file, and possibly to tty.
c  Makes call to hstall to utilize the hst package.


      logical(ISZ):: lgchange,checkzhistarray

      if (warpout > -1) then
        call oneliner (warpout,itt,tim,pz,ese,ek)
      endif

c     --- Call routine to print out diagnostics.  This can either be the
c     --- builting routine 'oneliner' or can be a Python function also
c     --- called 'oneliner'.
      if (loneliner) then
        call execuser("oneliner")
      else
        if (lspecial) call oneliner (STDOUT,itt,tim,pz,ese,ek)
      endif

      if (nhist > 0 .and. my_index == 0) then
c  Set builtin Hist arrays
        lgchange = .false.
c       --- increase the size of Hist arrays if necessary
        if (jhist == lenhist) then
          lenhist = lenhist + max(100,int(0.5*lenhist))
          lgchange = .true.
        endif

c       --- Make sure that the histories are set to have the same
c       --- number of species as the z moments.
c       --- This is not a good thing to do if nszmmnt was decreased since
c       --- some histories will be thrown away.
        if (nshist /= nszmmnt) then
          nshist = nszmmnt
          lgchange = .true.
        endif

c       --- Check whether any of the z histories have been turned on recently.
c       --- If so, space needs to be allocated for them.
        lgchange = checkzhistarray(lhlinechg,ihlinechg) .or. lgchange
        lgchange = checkzhistarray(lhvzofz,ihvzofz) .or. lgchange
        lgchange = checkzhistarray(lhcurrz,ihcurrz) .or. lgchange
        lgchange = checkzhistarray(lhpnumz,ihpnumz) .or. lgchange
        lgchange = checkzhistarray(lhepsxz,ihepsxz) .or. lgchange
        lgchange = checkzhistarray(lhepsyz,ihepsyz) .or. lgchange
        lgchange = checkzhistarray(lhepsnxz,ihepsnxz) .or. lgchange
        lgchange = checkzhistarray(lhepsnyz,ihepsnyz) .or. lgchange
        lgchange = checkzhistarray(lhepsgz,ihepsgz) .or. lgchange
        lgchange = checkzhistarray(lhepshz,ihepshz) .or. lgchange
        lgchange = checkzhistarray(lhepsngz,ihepsngz) .or. lgchange
        lgchange = checkzhistarray(lhepsnhz,ihepsnhz) .or. lgchange
        lgchange = checkzhistarray(lhxbarz,ihxbarz) .or. lgchange
        lgchange = checkzhistarray(lhybarz,ihybarz) .or. lgchange
        lgchange = checkzhistarray(lhxybarz,ihxybarz) .or. lgchange
        lgchange = checkzhistarray(lhxrmsz,ihxrmsz) .or. lgchange
        lgchange = checkzhistarray(lhyrmsz,ihyrmsz) .or. lgchange
        lgchange = checkzhistarray(lhrrmsz,ihrrmsz) .or. lgchange
        lgchange = checkzhistarray(lhxprmsz,ihxprmsz) .or. lgchange
        lgchange = checkzhistarray(lhyprmsz,ihyprmsz) .or. lgchange
        lgchange = checkzhistarray(lhxsqbarz,ihxsqbarz) .or. lgchange
        lgchange = checkzhistarray(lhysqbarz,ihysqbarz) .or. lgchange
        lgchange = checkzhistarray(lhvxbarz,ihvxbarz) .or. lgchange
        lgchange = checkzhistarray(lhvybarz,ihvybarz) .or. lgchange
        lgchange = checkzhistarray(lhvzbarz,ihvzbarz) .or. lgchange
        lgchange = checkzhistarray(lhxpbarz,ihxpbarz) .or. lgchange
        lgchange = checkzhistarray(lhypbarz,ihypbarz) .or. lgchange
        lgchange = checkzhistarray(lhvxrmsz,ihvxrmsz) .or. lgchange
        lgchange = checkzhistarray(lhvyrmsz,ihvyrmsz) .or. lgchange
        lgchange = checkzhistarray(lhvzrmsz,ihvzrmsz) .or. lgchange
        lgchange = checkzhistarray(lhxpsqbarz,ihxpsqbarz) .or. lgchange
        lgchange = checkzhistarray(lhypsqbarz,ihypsqbarz) .or. lgchange
        lgchange = checkzhistarray(lhxxpbarz,ihxxpbarz) .or. lgchange
        lgchange = checkzhistarray(lhyypbarz,ihyypbarz) .or. lgchange
        lgchange = checkzhistarray(lhxypbarz,ihxypbarz) .or. lgchange
        lgchange = checkzhistarray(lhyxpbarz,ihyxpbarz) .or. lgchange
        lgchange = checkzhistarray(lhxpypbarz,ihxpypbarz) .or. lgchange
        lgchange = checkzhistarray(lhxvzbarz,ihxvzbarz) .or. lgchange
        lgchange = checkzhistarray(lhyvzbarz,ihyvzbarz) .or. lgchange
        lgchange = checkzhistarray(lhvxvzbarz,ihvxvzbarz) .or. lgchange
        lgchange = checkzhistarray(lhvyvzbarz,ihvyvzbarz) .or. lgchange
        if (lgchange) call gchange("Hist",0)

c       --- save results of this time step in history buffers
c       --- must be done in seperate subroutine since the this subroutine
c       --- does not get effect of gchange
        if ( mod(itt,nhist) == 0) then
           call savehist(tim)
        endif

      endif

c  Make call to hst routine to utilize the hst package
      call hstall

      return
      end
c============================================================================
      logical function checkzhistarray(lhzarray,ihzarray)
      logical(ISZ):: lhzarray
      integer(ISZ):: ihzarray
c For histories of moments versus z, checks whether gchange should be called
c to allocate space. If the flag is turned on but the multiplier is zero,
c then space must be allocated for the history.  Also sets the multiplier
c to 1 if the flag has been set to collect the history.
      checkzhistarray = .false.
      if (lhzarray) then
        if (ihzarray < 0) checkzhistarray = .true.
        ihzarray = 1
      endif
      return
      end
c============================================================================
      subroutine savehist(tim)
      use Beam_acc
      use Picglb
      use Win_Moments
      use Moments
      use Z_Moments
      use Z_arrays
      use Hist
      real(kind=8):: tim
c  Save results of this time step in history buffers
      integer(ISZ):: iw,iz,js

      jhist = jhist + 1
      thist(jhist) = tim
      hzbeam(jhist) = zbeam
      hvbeam(jhist) = vbeam
      hefld(jhist) = ese

      do js=0,nshist
        hbmlen(jhist,js) = sqrt(max(SMALLPOS, zsqbar(0,js)-zbar(0,js)**2))
        hekzmbe(jhist,js) = ekzmbe(js)
        hekzbeam(jhist,js) = ekzbeam(js)
        hekperp(jhist,js) = ekperp(js)
        hxmaxp(jhist,js) = xmaxp(js)
        hxminp(jhist,js) = xminp(js)
        hymaxp(jhist,js) = ymaxp(js)
        hyminp(jhist,js) = yminp(js)
        hzmaxp(jhist,js) = zmaxp(js)
        hzminp(jhist,js) = zminp(js)
        hvxmaxp(jhist,js) = vxmaxp(js)
        hvxminp(jhist,js) = vxminp(js)
        hvymaxp(jhist,js) = vymaxp(js)
        hvyminp(jhist,js) = vyminp(js)
        hvzmaxp(jhist,js) = vzmaxp(js)
        hvzminp(jhist,js) = vzminp(js)

c       --- Hopefully, this comment no longer applies.
c       --- "This way" means using array notation, without the loops.
c       --- Write the code this way gets around an unexplained segmentation
c       --- fault that appears with the xlf compiler. It is almost certainly
c       --- a compiler bug since the seg fault goes away when using the -g
c       --- option, or when rewriting the code thusly.

        do iw=0,nzwind
          hepsx(iw,jhist,js)    = epsx(iw,js)
          hepsy(iw,jhist,js)    = epsy(iw,js)
          hepsz(iw,jhist,js)    = epsz(iw,js)
          hepsnx(iw,jhist,js)   = epsnx(iw,js)
          hepsny(iw,jhist,js)   = epsny(iw,js)
          hepsnz(iw,jhist,js)   = epsnz(iw,js)
          hepsg(iw,jhist,js)    = epsg(iw,js)
          hepsh(iw,jhist,js)    = epsh(iw,js)
          hepsng(iw,jhist,js)   = epsng(iw,js)
          hepsnh(iw,jhist,js)   = epsnh(iw,js)
          hpnum(iw,jhist,js)    = pnum(iw,js)
          hxbar(iw,jhist,js)    = xbar(iw,js)
          hybar(iw,jhist,js)    = ybar(iw,js)
          hzbar(iw,jhist,js)    = zbar(iw,js)
          hxrms(iw,jhist,js)    = xrms(iw,js)
          hyrms(iw,jhist,js)    = yrms(iw,js)
          hrrms(iw,jhist,js)    = rrms(iw,js)
          hzrms(iw,jhist,js)    = zrms(iw,js)
          hxprms(iw,jhist,js)   = xprms(iw,js)
          hyprms(iw,jhist,js)   = yprms(iw,js)
          hxsqbar(iw,jhist,js)  = xsqbar(iw,js)
          hysqbar(iw,jhist,js)  = ysqbar(iw,js)
          hxybar(iw,jhist,js)   = xybar(iw,js)
          hvxbar(iw,jhist,js)   = vxbar(iw,js)
          hvybar(iw,jhist,js)   = vybar(iw,js)
          hvzbar(iw,jhist,js)   = vzbar(iw,js)
          hxpbar(iw,jhist,js)   = xpbar(iw,js)
          hypbar(iw,jhist,js)   = ypbar(iw,js)
          hvxrms(iw,jhist,js)   = vxrms(iw,js)
          hvyrms(iw,jhist,js)   = vyrms(iw,js)
          hvzrms(iw,jhist,js)   = vzrms(iw,js)
          hxpsqbar(iw,jhist,js) = xpsqbar(iw,js) 
          hypsqbar(iw,jhist,js) = ypsqbar(iw,js)
          hxxpbar(iw,jhist,js)  = xxpbar(iw,js) 
          hyypbar(iw,jhist,js)  = yypbar(iw,js) 
          hxypbar(iw,jhist,js)  = xypbar(iw,js) 
          hyxpbar(iw,jhist,js)  = yxpbar(iw,js) 
          hxpypbar(iw,jhist,js) = xpypbar(iw,js) 
          hxvzbar(iw,jhist,js)  = xvzbar(iw,js) 
          hyvzbar(iw,jhist,js)  = yvzbar(iw,js) 
          hvxvzbar(iw,jhist,js) = vxvzbar(iw,js) 
          hvyvzbar(iw,jhist,js) = vyvzbar(iw,js) 
        enddo

c       --- Collect only those z histories specified by the user.
c       --- Making the hcurrz loop implicit gets around a problem
c       --- with the intel compiler.
        if (lhcurrz) hcurrz(:,jhist,js) = curr(:,js)
        do iz=0,nzmmnt
          if (lhpnumz)    hpnumz(iz,jhist,js)    = pnumz(iz,js)
          if (lhepsxz)    hepsxz(iz,jhist,js)    = epsxz(iz,js)
          if (lhepsyz)    hepsyz(iz,jhist,js)    = epsyz(iz,js)
          if (lhepsnxz)   hepsnxz(iz,jhist,js)   = epsnxz(iz,js)
          if (lhepsnyz)   hepsnyz(iz,jhist,js)   = epsnyz(iz,js)
          if (lhepsgz)    hepsgz(iz,jhist,js)    = epsgz(iz,js)
          if (lhepshz)    hepshz(iz,jhist,js)    = epshz(iz,js)
          if (lhepsngz)   hepsngz(iz,jhist,js)   = epsngz(iz,js)
          if (lhepsnhz)   hepsnhz(iz,jhist,js)   = epsnhz(iz,js)
          if (lhxbarz)    hxbarz(iz,jhist,js)    = xbarz(iz,js)
          if (lhybarz)    hybarz(iz,jhist,js)    = ybarz(iz,js)
          if (lhxybarz)   hxybarz(iz,jhist,js)   = xybarz(iz,js)
          if (lhxrmsz)    hxrmsz(iz,jhist,js)    = xrmsz(iz,js)
          if (lhyrmsz)    hyrmsz(iz,jhist,js)    = yrmsz(iz,js)
          if (lhrrmsz)    hrrmsz(iz,jhist,js)    = rrmsz(iz,js)
          if (lhxprmsz)   hxprmsz(iz,jhist,js)   = xprmsz(iz,js)
          if (lhyprmsz)   hyprmsz(iz,jhist,js)   = yprmsz(iz,js)
          if (lhxsqbarz)  hxsqbarz(iz,jhist,js)  = xsqbarz(iz,js)
          if (lhysqbarz)  hysqbarz(iz,jhist,js)  = ysqbarz(iz,js)
          if (lhvxbarz)   hvxbarz(iz,jhist,js)   = vxbarz(iz,js)
          if (lhvybarz)   hvybarz(iz,jhist,js)   = vybarz(iz,js)
          if (lhvzbarz)   hvzbarz(iz,jhist,js)   = vzbarz(iz,js)
          if (lhxpbarz)   hxpbarz(iz,jhist,js)   = xpbarz(iz,js)
          if (lhypbarz)   hypbarz(iz,jhist,js)   = ypbarz(iz,js)
          if (lhvxrmsz)   hvxrmsz(iz,jhist,js)   = vxrmsz(iz,js)
          if (lhvyrmsz)   hvyrmsz(iz,jhist,js)   = vyrmsz(iz,js)
          if (lhvzrmsz)   hvzrmsz(iz,jhist,js)   = vzrmsz(iz,js)
          if (lhxpsqbarz) hxpsqbarz(iz,jhist,js) = xpsqbarz(iz,js)
          if (lhypsqbarz) hypsqbarz(iz,jhist,js) = ypsqbarz(iz,js)
          if (lhxxpbarz)  hxxpbarz(iz,jhist,js)  = xxpbarz(iz,js)
          if (lhyypbarz)  hyypbarz(iz,jhist,js)  = yypbarz(iz,js)
          if (lhxypbarz)  hxypbarz(iz,jhist,js)  = xypbarz(iz,js)
          if (lhyxpbarz)  hyxpbarz(iz,jhist,js)  = yxpbarz(iz,js)
          if (lhxpypbarz) hxpypbarz(iz,jhist,js) = xpypbarz(iz,js)
          if (lhxvzbarz)  hxvzbarz(iz,jhist,js)  = xvzbarz(iz,js)
          if (lhyvzbarz)  hyvzbarz(iz,jhist,js)  = yvzbarz(iz,js)
          if (lhvxvzbarz) hvxvzbarz(iz,jhist,js) = vxvzbarz(iz,js)
          if (lhvyvzbarz) hvyvzbarz(iz,jhist,js) = vyvzbarz(iz,js)
        enddo
      enddo

c     --- These are the arrays which do not have the species index.
      do iw=0,nzwind
        hrhomid(iw,jhist)  = rhomid(iw)
        hrhomax(iw,jhist)  = rhomax(iw)
      enddo
      do iz=0,nzzarr
        if (lhlinechg)  hlinechg(iz,jhist) = linechg(iz)
        if (lhvzofz)    hvzofz(iz,jhist) = vzofz(iz)
      enddo

      return
      end
c=============================================================================
      subroutine oneliner (iunit,it,time,pz,ese,ek)
      use Io
      integer(ISZ):: iunit,it
      real(kind=8):: time,pz,ese,ek
      character*(240)::o_line


c  Summarize the state of the system in one line

c     --- return if line is not requested
      if (verbosity .lt. 2) return

      write (o_line,901) it,time,pz,ese,ek,ese+ek
 901  format("it =",i6," time =",1pe11.4," pz =",1pe11.4,
     & " ese =",1pe11.4," ek =",1pe11.4," et  =",1pe11.4)
      call remark(trim(o_line))

      return
      end
c=============================================================================
      subroutine prntpara(dxperp,dyperp,dz,lprntpara,pgroup)
      use ParticleGroupmodule
      use InGen
      use InPart
      use Io
      use Ch_var
      use Constant
      use Particles,Only: npmax
      use Picglb
      use OutParams
      use Beam_acc
      real(kind=8):: dxperp,dyperp,dz
      logical(ISZ):: lprntpara
      type(ParticleGroup):: pgroup

c  Calculates various parameters and optionally prints them out to a plot
c  frame and an output file or tty

      character(80):: textline

c  Calculate various quantities first

c  Current density at the center of the beam
      currdens = ibeam/dvnz(Pi*a0*b0)
c  Charge density at the center of the beam
      chrgdens = currdens/dvnz(vbeam)
c  Number density
      numbdens = chrgdens/dvnz(zion*echarge)
c  Plasma frequency
      omegap = sqrt(abs(numbdens*(zion*echarge)**2/dvnz(eps0)/dvnz(aion*amu)))
      omegapdt = omegap*dt
      omegaptq = omegap*(tunelen/dvnz(vbeam))
      taup = 2.*Pi/dvnz(omegap)
c  Transverse thermal velocities
      vthx = 0.5*vbeam*emitx/dvnz(a0)
      vthy = 0.5*vbeam*emity/dvnz(b0)
      vthxdt = vthx*dt
      vthydt = vthy*dt
      vthxdtodx = vthx*dt/dvnz(dxperp)
      vthydtody = vthy*dt/dvnz(dyperp)
c  Transverse Debye wavelength
      lamdebx = vthx/dvnz(omegap)
      lamdeby = vthy/dvnz(omegap)
      lamdebxodx = lamdebx/dvnz(dxperp)
      lamdebyody = lamdeby/dvnz(dyperp)
c  Longitudinal thermal velocity (rms)
      vthzdt = vthz*dt
      vthzdtodz = vthz*dt/dvnz(dz)
c  Longitudinal Debye wavelength
      lamdebz = vthz/dvnz(omegap)
      lamdebzodz = lamdebz/dvnz(dz)
c  Beam velocity
      vbeamoc = vbeam/dvnz(clight)
c  Totals
      npreal  = npmax*pgroup%sw(1)
      totmass = npmax*pgroup%sw(1)*(aion*amu)
      totchrg = npmax*pgroup%sw(1)*(zion*echarge)
c  Generalized perveance
      genperv = ibeam*(zion*echarge)/
     &          dvnz(2.*Pi*eps0*(aion*amu)*(vbeam*gammabar)**3)
c  Characteristic current
      charcurr = 4.*Pi*eps0*(aion*amu)*clight**3/dvnz((zion*echarge))
c  Budker parameter
      budker = ibeam/dvnz(charcurr*vbeam/clight)
c  Particle Betatron frequencies, both undepressed and depressed in the 
c  X- and Y- planes and in various units
      lambdab0x = tunelen*360./dvnz(sigma0x)
      lambdab0y = tunelen*360./dvnz(sigma0y)
      taub0x    = lambdab0x/dvnz(vbeam)
      taub0y    = lambdab0y/dvnz(vbeam)
      omegab0x  = 2.*Pi/dvnz(taub0x)
      omegab0y  = 2.*Pi/dvnz(taub0y)
      lambdabx  = tunelen*360./dvnz(sigmax)
      lambdaby  = tunelen*360./dvnz(sigmay)
      taubx     = lambdabx/dvnz(vbeam)
      tauby     = lambdaby/dvnz(vbeam)
      omegabx   = 2.*Pi/dvnz(taubx)
      omegaby   = 2.*Pi/dvnz(tauby)
c  Space charge wave velocity
      vwave = 0.5*sqrt(max(0.,gfactor))*omegap*sqrt(a0*b0)
c  Ratio of space charge and emittance forces
      femtxofscx = emitx**2/dvnz(2.*genperv)*(a0 + b0)/dvnz(a0**3)
      femtyofscy = emity**2/dvnz(2.*genperv)*(a0 + b0)/dvnz(b0**3)
c  Exit now if output parameters are not to be printed
      if (.not. lprntpara) return

c     --- Call script version of this routine.
      call execuser("printparameters")

c  Write to plot frame
   20 format(1x,a,1pe11.4,a)
   40 format(1x,a,1pe11.4,a,1pe11.4,a)
   30 format(1x,a,i8,a)
      write (textline,20) "Atomic number of ion = ",aion," "
      call remark(trim(textline))
      write (textline,20) "Charge state of ion  = ",zion," "
      call remark(trim(textline))
      write (textline,40) "Initial X, Y emittances = ",
     &      emitx,",  ",emity," m-rad"
      call remark(trim(textline))
      write (textline,40) "Initial X,Y envelope radii  = ",a0,",  ",b0," m"
      call remark(trim(textline))
      write (textline,40) "Initial X,Y envelope angles = ",ap0,",  ",bp0," rad"
      call remark(trim(textline))
      write (textline,20) "Input beam current = ",ibeam," amps"
      call remark(trim(textline))
      write (textline,20) "Current density = ",currdens," amps/m**2"
      call remark(trim(textline))
      write (textline,20) "Charge density = ",chrgdens," Coul/m**3"
      call remark(trim(textline))
      write (textline,20) "Number density = ",numbdens," "
      call remark(trim(textline))
      write (textline,20) "Plasma frequency     = ",omegap," 1/s"
      call remark(trim(textline))
      write (textline,20) "   times dt          = ",omegapdt," "
      call remark(trim(textline))
      write (textline,20) "   times quad period = ",omegaptq," "
      call remark(trim(textline))
      write (textline,20) "Plasma period        = ",taup," s"
      call remark(trim(textline))
      write (textline,40) "X-, Y-Thermal Velocities     = ",vthx,",  ",
     &                                                      vthy," m/s"
      call remark(trim(textline))
      write (textline,40) "   times dt                  = ",vthxdt,",  ", 
     &                                                      vthydt," m"
      call remark(trim(textline))
      write (textline,40) "   times dt/dx, dt/dy (X, Y) = ",vthxdtodx,",  ", 
     &                                                      vthydtody," "
      call remark(trim(textline))
      write (textline,40) "X-, Y-Debye Wavelengths  = ",lamdebx,",  ",
     &                                                  lamdeby," m"
      call remark(trim(textline))
      write (textline,40) "   over dx, dy (X and Y) = ",lamdebxodx,",  ", 
     &                                                  lamdebyody," "
      call remark(trim(textline))
      write (textline,20) "Longitudinal thermal velocity (rms) = ",vthz," m/s"
      call remark(trim(textline))
      write (textline,20) "   times dt                   = ",vthzdt," m"
      call remark(trim(textline))
      write (textline,20) "   times dt/dz                = ",vthzdtodz," "
      call remark(trim(textline))
      write (textline,20) "Longitudinal Debye wavelength = ",lamdebz," m"
      call remark(trim(textline))
      write (textline,20) "   over dz                    = ",lamdebzodz," "
      call remark(trim(textline))
c  Start a new frame since they all don't fit on one
      write (textline,20) "Beam velocity = ",vbeam," m/s"
      call remark(trim(textline))
      write (textline,20) "   over c     = ",vbeamoc," "
      call remark(trim(textline))
      write (textline,20) "Kinetic energy = ",ekin," eV"
      call remark(trim(textline))
      write (textline,20) "Weight of simulation particles = ",pgroup%sw(1)," "
      call remark(trim(textline))
      write (textline,30) "Number of simulation particles = ",npmax," "
      call remark(trim(textline))
      write (textline,20) "Number of real particles = ",npreal," "
      call remark(trim(textline))
      write (textline,20) "Total mass = ",totmass," kg"
      call remark(trim(textline))
      write (textline,20) "Total charge = ",totchrg," Coul"
      call remark(trim(textline))
      write (textline,20) "Generalized perveance = ",genperv," "
      call remark(trim(textline))
      write (textline,20) "Characteristic current = ",charcurr," amps"
      call remark(trim(textline))
      write (textline,20) "Budker parameter = ",budker," "
      call remark(trim(textline))
      write (textline,20) "Timestep size dt = ",dt," s"
      call remark(trim(textline))
      write (textline,20) "Tune length = ",tunelen," "
      call remark(trim(textline))
      write (textline,40) "Undep. X-, Y-Betatron frequencies  = ",
     &                    omegab0x,",  ",omegab0y," 1/s"
      call remark(trim(textline))
      write (textline,40) "Undep. X-, Y-Betatron periods      = ",
     &                    taub0x,",  ",taub0y," s"
      call remark(trim(textline))
      write (textline,40) "Undep. X-, Y-Betatron wavelengths  = ",
     &                    lambdab0x,",  ",lambdab0y," m"
      call remark(trim(textline))
      write (textline,40) "Dep.   X-, Y-Betatron frequencies  = ",
     &                    omegabx,",  ",omegaby," 1/s"
      call remark(trim(textline))
      write (textline,40) "Dep.   X-, Y-Betatron periods      = ",
     &                    taubx,",  ",tauby," s"
      call remark(trim(textline))
      write (textline,40) "Dep.   X-, Y-Betatron wavelengths  = ",
     &                    lambdabx,",  ",lambdaby," m"
      write (textline,40) "X-, Y-Tune Depressions (sigma/sigma0) = ",
     &                    sigmax/dvnz(sigma0x),",  ",sigmay/dvnz(sigma0y)," "
      call remark(trim(textline))
      write (textline,20) "Space charge wave velocity = ",vwave," m/s"
      call remark(trim(textline))
      write (textline,20) "Effective wall radius = ",rwall," m"
      call remark(trim(textline))
      write (textline,20) "Geometric factor = ",gfactor," "
      call remark(trim(textline))
      write (textline,40) "X-, Y-Emittance over Space charge forces = ",
     &                                  femtxofscx,",  ",femtyofscy," "
      call remark(trim(textline))

c  Write to text output file
      if (warpout > -1) then
        call edit (warpout, "runid")
        call edit (warpout, "it")
        call edit (warpout, "time")
        call edit (warpout, "aion")
        call edit (warpout, "zion")
        call edit (warpout, "emitx")
        call edit (warpout, "emity")
        call edit (warpout, "a0")
        call edit (warpout, "b0")
        call edit (warpout, "ap0")
        call edit (warpout, "bp0")
        call edit (warpout, "ibeam")
        call edit (warpout, "currdens")
        call edit (warpout, "chrgdens")
        call edit (warpout, "numbdens")
        call edit (warpout, "omegap")
        call edit (warpout, "omegapdt")
        call edit (warpout, "omegaptq")
        call edit (warpout, "taup")
        call edit (warpout, "vthx")
        call edit (warpout, "vthxdt")
        call edit (warpout, "vthxdtodx")
        call edit (warpout, "vthy")
        call edit (warpout, "vthydt")
        call edit (warpout, "vthydtody")
        call edit (warpout, "lamdebx")
        call edit (warpout, "lamdebxodx")
        call edit (warpout, "lamdeby")
        call edit (warpout, "lamdebyody")
        call edit (warpout, "vthz")
        call edit (warpout, "vthzdt")
        call edit (warpout, "vthzdtodz")
        call edit (warpout, "lamdebz")
        call edit (warpout, "lamdebzodz")
        call edit (warpout, "vbeam")
        call edit (warpout, "vbeamoc")
        call edit (warpout, "ekin")
        call edit (warpout, "sw")
        call edit (warpout, "np")
        call edit (warpout, "npreal")
        call edit (warpout, "totmass")
        call edit (warpout, "totchrg")
        call edit (warpout, "genperv")
        call edit (warpout, "charcurr")
        call edit (warpout, "budker")
        call edit (warpout, "dt")
        call edit (warpout, "tunelen")
        call edit (warpout, "omegab0x")
        call edit (warpout, "omegab0y")
        call edit (warpout, "taub0x")
        call edit (warpout, "taub0y")
        call edit (warpout, "lambdab0x")
        call edit (warpout, "lambdab0y")
        call edit (warpout, "omegabx")
        call edit (warpout, "omegaby")
        call edit (warpout, "taubx")
        call edit (warpout, "tauby")
        call edit (warpout, "lambdabx")
        call edit (warpout, "lambdaby")
        call edit (warpout, "vwave")
        call edit (warpout, "rwall")
        call edit (warpout, "gfactor")
        call edit (warpout, "femtxofscx")
        call edit (warpout, "femtyofscy")
      endif

      return
      end
c============================================================================
      subroutine psplots (freqflag)
      use InDiag
      use InPltCtl
      use Timers
      integer(ISZ):: freqflag

c  Phase space plots, both "frequent" ones and others


      real(kind=8):: timetemp,wtime
      integer(ISZ):: i

      timetemp = wtime()

      if (lpsplots) then
        if (freqflag == ALWAYS) call parsestr("psplotsalways()")
        if (freqflag == SELDOM) call parsestr("psplotsseldom()")
        plottime = plottime + (wtime() - timetemp)
        return
      endif

      plottime = plottime + (wtime() - timetemp)
      return
      end
c=============================================================================
      subroutine onedplts(freqflag)
      use InDiag
      use InPltCtl
      use Timers
      integer(ISZ):: freqflag
c  Plots all 1d (z dependent) quantities vs z
      real(kind=8):: timetemp,wtime
      timetemp = wtime()

      if (lonedplts) then
        if (freqflag == ALWAYS) call parsestr("onedpltsalways()")
        if (freqflag == SELDOM) call parsestr("onedpltsseldom()")
        return
      endif

      plottime = plottime + (wtime() - timetemp)
      return
      end
c=============================================================================
c=============================================================================
c  OTHER-THAN-PHYSICS-OR-GRAPHICS ROUTINES
c=============================================================================
      subroutine printpkgversion(iout,package,version)
      integer(ISZ):: iout
      character(*):: package
      character(19):: version
c Prints out the package name and version number. This routine is only
c used for the printing of the version number at the start up of the code.
      integer(ISZ):: i,lenv

c     --- Get length of version string. This is very specific to the cvs
c     --- versioning system. It assumes that the version string passed
c     --- in is of the format "$Revision: 3.156 $" where 'r' and 'v' are
c     --- the revision number of current version number.
      lenv = 12
      do i=12,19
        if (version(i:i) == "$") lenv = i
      enddo

c     --- Make the printout
      write (iout,'(" ******  ",a," version ",a)') package,version(12:lenv-2)

      return
      end
c=============================================================================
      subroutine parvers(iunit)
      use Code_version
      use Ch_var
      integer(ISZ):: iunit
      character(79):: version
      character(72):: author
      common /cverpar/ version, author

      version = "Particle simulation of heavy ion beams for inertial fusion"
      author =
     & "Alex Friedman, Dave Grote, Debbie Callahan, Bruce Langdon, Irv Haber"

      write(iunit,79) version, author, codeid, runid, runmaker
   79 format(a,/,a,/,a,/,a,/,a)

      write(iunit,89) pline3, pline2, pline1
   89 format(/,a,/,a,/,a,//)

      return
      end
c=============================================================================
      subroutine stepid (it, time, zbeam)
      use Ch_var
      integer(ISZ):: it
      real(kind=8):: time,zbeam

c  Loads character variable PLINE3 with time-dependent simulation params
c  associated with the current time step. These appear in the top line
c  of identifying information on each plot frame.

      integer(ISZ):: ilogt,ilogz,ilogtdigs,ilogzdigs
      character(3):: texp,zexp

c     --- The time and zbeam are printed out with an exponent that is a
c     --- multiple of 3 and the mantissa is between 1 and 999.
c     --- The extra if statement deals with the cases where the quantity
c     --- is very close to 10**(3i). The 5.e-5 is added since the f10.4 in the
c     --- write statements will round up numbers within 5.e-5 of 10**(3i).
      if (time == 0.) then
        ilogt = 0.
      else
        ilogt = int((log10(abs(time))-3.)/3.)*3
      endif
      if (time/10.**ilogt+5.e-5 >= 1000.) ilogt = ilogt + 3
      if (zbeam == 0.) then
        ilogz = 0.
      else
        ilogz = int((log10(abs(zbeam))-3.)/3.)*3
      endif
      if (zbeam/10.**ilogz+5.e-5 >= 1000.) ilogz = ilogz + 3

c     --- There is no left-justified format control, so it has to be dealt
c     --- explicitly. A normal write into a string is done and if the first
c     --- character is a space, it is removed. The ilogtdigs is used so that
c     --- the last character can be skipped.
c     --- The SP forces a + sign to be printed with ilogt > 0.
      write (texp,'(SP,i3)') ilogt
      ilogtdigs = 3
      if (texp(1:1) == " ") then
        texp(1:2) = texp(2:3)
        texp(3:3) = " "
        ilogtdigs = 2
      endif
      write (zexp,'(SP,i3)') ilogz
      ilogzdigs = 3
      if (zexp(1:1) == " ") then
        zexp(1:2) = zexp(2:3)
        zexp(3:3) = " "
        ilogzdigs = 2
      endif

      write ( pline3,
     & '("Step",i8,", T =",f10.4,"e",a," s, Zbeam =",f10.4,"e",a," m")' )
     & it,time/10.**ilogt,texp(1:ilogtdigs),zbeam/10.**ilogz,zexp(1:ilogzdigs)
      write ( frameti, '("Step",i7)' ) it

      return
      end
c=============================================================================
      logical function thisstep (it, icontrol, ncontrol)
      integer(ISZ):: it,ncontrol
      integer(ISZ):: icontrol(ncontrol)

      integer(ISZ):: istart,iend,iinc,i

      istart = icontrol(1)
      iend   = icontrol(2)
      iinc   = icontrol(3)
      thisstep = .false.

      if (iinc .ne. 0) then
         if ((it >= istart).and.(it <= iend)) then
            if (mod(it-istart,iinc) == 0) then
               thisstep = .true.
               return
            endif
         endif
      endif

      do i = 4,ncontrol
         if (it == icontrol(i)) thisstep = .true.
      enddo

      return
      end
c=============================================================================
      logical function thiszbeam(zbeaml,zbeamr,control,ncontrol)
      real(kind=8):: zbeaml,zbeamr
      integer(ISZ):: ncontrol
      real(kind=8):: control(ncontrol)

c Given a z range and a control array, this subroutine checks if a value
c of z given by the loop specified by the first three elements of the control
c array or a subsequent value in the control array are within that z range.
c This is used for the diagnostics which are done at specified z locations.
c
c The extra 'inc' added to the expressions to calculate il and ir is there
c so that for the case when zbeaml<start<zbeamr, a value of true is returned.
c Otherwise, the expression for il returns a negative number which is rounded
c up to zero and so is equal to ir.

      real(kind=8):: zstart,zend,zinc
      integer(ISZ):: i,il,ir

      zstart = control(1)
      zend   = control(2)
      zinc   = control(3)
      thiszbeam = .false.

      if (zinc .ne. 0.) then
        if ((zstart <= zbeamr) .and. (zbeaml <= zend)) then
          il = int((zbeaml - zstart + zinc)/dvnz(zinc))
          ir = int((zbeamr - zstart + zinc)/dvnz(zinc))
          if (ir > il) then
            thiszbeam = .true.
            return
          endif
        endif
      endif

      do i = 4,ncontrol
        if (zbeaml <= control(i) .and. control(i) < zbeamr)
     &    thiszbeam = .true.
      enddo

      return
      end
c=============================================================================
      logical function thistime(time1,time2,control,ncontrol)
      real(kind=8):: time1,time2
      integer(ISZ):: ncontrol
      real(kind=8):: control(ncontrol)

c Given a time range and a control array, this subroutine checks if a value
c of time given by the loop specified by the first three elements of the control
c array or a subsequent value in the control array are within that time range.
c This is used for the diagnostics which are done at specified times.
c
c The extra 'inc' added to the expressions to calculate il and ir is there
c so that for the case when time1<start<time2, a value of true is returned.
c Otherwise, the expression for il returns a negative number which is rounded
c up to zero and so is equal to ir.

      real(kind=8):: tstart,tend,tinc
      integer(ISZ):: i,il,ir

      tstart = control(1)
      tend   = control(2)
      tinc   = control(3)
      thistime = .false.

      if (tinc .ne. 0.) then
        if ((tstart <= time2) .and. (time1 <= tend)) then
          il = int((time1 - tstart + tinc)/dvnz(tinc))
          ir = int((time2 - tstart + tinc)/dvnz(tinc))
          if (ir > il) then
            thistime = .true.
            return
          endif
        endif
      endif

      do i = 4,ncontrol
        if (time1 <= control(i) .and. control(i) < time2)
     &    thistime = .true.
      enddo

      return
      end
c=============================================================================
      logical function dolabwn()
      use Picglb
      use Z_Moments
      use InDiag
      use Lab_Moments

c Check if the lab windows calculation will be done this step, checking
c each of the following conditions.
c   Are the lab windows turned on?
c   Have the correct number of timestep been passed?
c   Are there any lab windows in the beam frame?


      integer(ISZ):: il

c     --- Make sure the lab windows are initialized.
      call initlabwn(1)

      dolabwn = .false.

      if (nlabwn > 0) then
        if (mod(it,itlabwn) == 0) then

          do il=1,nlabwn
            if (zbeam+zmmntmin < zlw(il) .and.
     &                           zlw(il) < zbeam+zmmntmax) dolabwn = .true.
          enddo

        endif
      endif

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine tolabfrm(zcent,nn,xxx,zzz)
      use Lattice
      integer(ISZ):: nn
      real(kind=8):: zcent
      real(kind=8):: xxx(nn), zzz(nn)

c Converts data to lab frame using bend data in Lattice.
c The arrays xxx(nn) and zzz(nn) are converted in place.
c This routine works by first finding the bends over which the data extends
c and converts the bend locations into the lab frame.  The data is then
c converted relative to the nearby bends.
c
c The arrays contain the following data...
c   bendlats     Starts of bends in Lattice frame relative to zcent
c   bendlate     Ends of bends in Lattice frame relative to zcent
c   bend_a       Bend angle of each bend
c   bends_sum    Sum of bend angles from zcent up to bend start
c   cosbssum     Cosine of sums
c   sinbssum     Sine of sums
c   bende_sum    Sum of bend angles from zcent up to bend end
c   cosbesum     Cosine of sums
c   sinbesum     Sine of sums
c   labzs        Z location of bend starts in lab frame
c   labxs        X location of bend starts in lab frame
c   labze        Z location of bend ends in lab frame
c   labxe        X location of bend ends in lab frame


      integer(ISZ),allocatable:: ibend(:)
      real(kind=8),allocatable:: bendlats(:), bendlate(:)
      real(kind=8),allocatable:: bend_a(:)
      real(kind=8),allocatable:: bends_sum(:), cosbssum(:), sinbssum(:)
      real(kind=8),allocatable:: bende_sum(:), cosbesum(:), sinbesum(:)
      real(kind=8),allocatable:: labzs(:), labxs(:)
      real(kind=8),allocatable:: labze(:), labxe(:)
      real(kind=8):: offset,zzzmin,zzzmax,offsetmx,offsetmn,angle,aa,zz,tz
      integer(ISZ):: i,ib,numb,icentl,icentr,id,ibendl,ibendr,il,ir
      integer(ISZ):: ibendmn,ibendmx

c Return if there are no bends
      if (.not. bends) return

c Make sure that only bends within the range 0 to zlatperi are actually used.
      ibendmn = 0
      ibendmx = nbend
      do ib=0,nbend
        if (bendze(ib) < 0.) ibendmn = ib + 1
        if (bendzs(nbend-ib) > zlatperi .and. zlatperi > 0.) then
          ibendmx = nbend-ib - 1
        endif
      enddo

c Return if there are no usable bends
      if (ibendmx < ibendmn) return

c Calculate extent of data and offsets to get to lattice frame.
      zzzmin = zcent
      zzzmax = zcent
      do i=1,nn
        if (zzz(i) < zzzmin) zzzmin = zzz(i)
        if (zzz(i) > zzzmax) zzzmax = zzz(i)
      enddo
      if (zlatperi > 0.) then
        offsetmx = floor((zzzmax-zlatstrt)/zlatperi)*zlatperi
        offsetmn = floor((zzzmin-zlatstrt)/zlatperi)*zlatperi
        if (zzzmin <= bendzs(ibendmn)+offsetmn) offsetmn = offsetmn - zlatperi
        if (zzzmax >= bendze(ibendmx)+offsetmx) offsetmx = offsetmx + zlatperi
      else
        offsetmx = zlatstrt
        offsetmn = zlatstrt
      endif

c Find bend to the left of zzzmin and to the right of zzzmax.
      ibendl = 0
      ibendr = nbend
      do ib=ibendmn,ibendmx
        if (zzzmin > bendzs(ib)+offsetmn) ibendl = ib
        if (zzzmax < bendze(nbend-ib)+offsetmx) ibendr = nbend - ib
      enddo

c Make array to index lattice elements.
      numb = ibendr+1 - ibendl
      if (zlatperi > 0.) then
        numb = numb + nint((nbend+1)*(offsetmx - offsetmn)/zlatperi)
      endif

      allocate(ibend(numb),bendlats(numb),bendlate(numb),bend_a(numb))
      allocate(bends_sum(numb),cosbssum(numb),sinbssum(numb))
      allocate(bende_sum(numb),cosbesum(numb),sinbesum(numb))
      allocate(labzs(numb),labxs(numb),labze(numb),labxe(numb))

      ibend(1) = ibendl
      do i=2,numb
        ibend(i) = ibend(i-1) + 1
        if (ibend(i) > ibendmx) ibend(i) = ibendmn
      enddo

c Calculate bend starts and ends of bends in lattice frame.
c     --- They are constrained to be within 0 and zlatperi when zlatperi>0
      ib = ibendl
      offset = offsetmn
      do i=1,numb
        if (zlatperi == 0.) then
          bendlats(i) = offset + bendzs(ibend(i))
          bendlate(i) = offset + bendze(ibend(i))
        else
          bendlats(i) = offset + max(0.,bendzs(ibend(i)))
          bendlate(i) = offset + min(zlatperi,bendze(ibend(i)))
        endif
        ib = ib + 1
        if (ib > ibendmx) then
          ib = ibendmn
          offset = offset + zlatperi
        endif
      enddo

c Calculate bend angles.
      do i=1,numb
        ib = ibend(i)
        bend_a(i) = (bendze(ib) - bendzs(ib))/bendrc(ib)
      enddo

c Calculate index of bend left and right of frame center.
      icentl = 0
      icentr = numb + 1
      do i=1,numb
        if (zcent > bendlats(i)) icentl = i
        if (zcent < bendlate(numb-i+1)) icentr = numb-i+1
      enddo

c Calculate sum of bend angles from bends surrounding zcent up to the
c bend previous to the current bend.
      if (icentl > 0) then
        i = icentl
        zz = min((bendlate(i)-bendlats(i)),(zcent - bendlats(i)))
        bends_sum(i) = zz/bendrc(ibend(i))
        bende_sum(i) = max(0.,(bendlate(i) - zcent))/bendrc(ibend(i))
        cosbssum(i) = cos(bends_sum(i))
        sinbssum(i) = sin(bends_sum(i))
        cosbesum(i) = cos(bende_sum(i))
        sinbesum(i) = sin(bende_sum(i))
        do i = icentl-1,1,-1
          bende_sum(i) = bends_sum(i+1)
          bends_sum(i) = bende_sum(i) + bend_a(i)
          cosbssum(i) = cos(bends_sum(i))
          sinbssum(i) = sin(bends_sum(i))
          cosbesum(i) = cos(bende_sum(i))
          sinbesum(i) = sin(bende_sum(i))
        enddo
      endif

      if (icentr <= numb) then
        i = icentr
        zz = min((bendlate(i)-bendlats(i)),(bendlate(i) - zcent))
        bends_sum(i) = max(0.,(zcent - bendlats(i)))/bendrc(ibend(i))
        bende_sum(i) = zz/bendrc(ibend(i))
        cosbssum(i) = cos(bends_sum(i))
        sinbssum(i) = sin(bends_sum(i))
        cosbesum(i) = cos(bende_sum(i))
        sinbesum(i) = sin(bende_sum(i))
        do i=icentr+1,numb
          bends_sum(i) = bende_sum(i-1)
          bende_sum(i) = bends_sum(i) + bend_a(i)
          cosbssum(i) = cos(bends_sum(i))
          sinbssum(i) = sin(bends_sum(i))
          cosbesum(i) = cos(bende_sum(i))
          sinbesum(i) = sin(bende_sum(i))
        enddo
      endif

c Calculate location of starts and ends of bends in lab frame.
c --- assumes that zcent is not in a bend
c --- left of zcent
      if (icentl > 0) then
        i = icentl
        if (zcent >= bendlate(i)) then
          labze(i) = bendlate(i)
          labxe(i) = 0.
          labzs(i) = labze(i) - bendrc(ibend(i))*sin(bend_a(i))
          labxs(i) = labxe(i) - bendrc(ibend(i))*(1. - cos(bend_a(i)))
          angle = bend_a(i)
        else 
          aa = (bendlate(i)-zcent)/bendrc(ibend(i))
          labze(i) = zcent + bendrc(ibend(i))*sin(aa)
          labxe(i) = 0. - bendrc(ibend(i))*(1. - cos(aa))
          aa = (zcent - bendlats(i))/bendrc(ibend(i))
          labzs(i) = zcent - bendrc(ibend(i))*sin(aa)
          labxs(i) = 0. - bendrc(ibend(i))*(1. - cos(aa))
          angle = aa
        endif
        do i=icentl-1,1,-1
          labze(i) = labzs(i+1) - (bendlats(i+1) - bendlate(i))*cos(angle)
          labxe(i) = labxs(i+1) - (bendlats(i+1) - bendlate(i))*sin(angle)
          labzs(i) = labze(i) - bendrc(ibend(i))*(sin(angle+bend_a(i)) -
     &               sin(angle))
          labxs(i) = labxe(i) - bendrc(ibend(i))*(cos(angle) -
     &               cos(angle+bend_a(i)))
          angle = angle + bend_a(i)
        enddo
      endif

c --- right of zcent
      if (icentr <= numb) then
        i = icentr
        if (zcent <= bendlats(i)) then
          labzs(i) = bendlats(i)
          labxs(i) = 0.
          labze(i) = labzs(i) + bendrc(ibend(i))*sin(bend_a(i))
          labxe(i) = labxs(i) - bendrc(ibend(i))*(1. - cos(bend_a(i)))
          angle = bend_a(i)
        else 
          aa = (zcent - bendlats(i))/bendrc(ibend(i))
          labzs(i) = zcent - bendrc(ibend(i))*sin(aa)
          labxs(i) = 0. - bendrc(ibend(i))*(1. - cos(aa))
          aa = (bendlate(i)-zcent)/bendrc(ibend(i))
          labze(i) = zcent + bendrc(ibend(i))*sin(aa)
          labxe(i) = 0. - bendrc(ibend(i))*(1. - cos(aa))
          angle = aa
        endif
        do i=icentr+1,numb
          labzs(i) = labze(i-1) + (bendlats(i) - bendlate(i-1))*cos(angle)
          labxs(i) = labxe(i-1) - (bendlats(i) - bendlate(i-1))*sin(angle)
          labze(i) = labzs(i) + bendrc(ibend(i))*(sin(angle+bend_a(i)) -
     &               sin(angle))
          labxe(i) = labxs(i) - bendrc(ibend(i))*(cos(angle) -
     &               cos(angle+bend_a(i)))
          angle = angle + bend_a(i)
        enddo
      endif

c Convert data to lab frame.
      do id=1,nn

c       --- Find bend starting to the left and bend ending to the right of zzz.
        il = 0
        ir = 0
        do i=1,numb
          if (zzz(id) > bendlats(i)) il = i
          if (zzz(id) < bendlate(numb-i+1)) ir = numb-i+1
        enddo

c       --- Do work for cases where zzz < end of bend left of zcent.
c       --- Only two cases: zzz is in bend, so il=ir; or zzz is in drift so
c       --- ir = il+1.
        if (zzz(id) < bendlate(icentl)) then
          if (il == ir) then
            aa = bends_sum(il) - (zzz(id) - bendlats(il))/bendrc(ibend(il))
            tz = (labzs(il) - bendrc(ibend(il))*(sin(aa) - sinbssum(il)) -
     &            xxx(id)*sin(aa))
            xxx(id) = (labxs(il) - bendrc(ibend(il))*(cosbssum(il) - cos(aa)) +
     &                 xxx(id)*cos(aa))
            zzz(id) = tz
          else
            tz = (labzs(ir) - (bendlats(ir) - zzz(id))*cosbssum(ir) -
     &            xxx(id)*sinbssum(ir))
            xxx(id) = (labxs(ir) - (bendlats(ir)-zzz(id))*sinbssum(ir) +
     &                 xxx(id)*cosbssum(ir))
            zzz(id) = tz
          endif
c       --- Do work for cases where zzz > start of bend right of zcent.
c       --- Only two cases: zzz is in bend, so il=ir; or zzz is in drift so
c       --- ir = il+1.
        elseif (zzz(id) > bendlats(icentr)) then
          if (il == ir) then
            aa = bende_sum(il) - (bendlate(il) - zzz(id))/bendrc(ibend(il))
            tz = (labze(il) + bendrc(ibend(il))*(sin(aa) - sinbesum(il)) +
     &            xxx(id)*sin(aa))
            xxx(id) = (labxe(il) - bendrc(ibend(il))*(cosbesum(il) - cos(aa)) +
     &                 xxx(id)*cos(aa))
            zzz(id) = tz
          else
            tz = (labze(il) + (zzz(id) - bendlate(il))*cosbesum(il) +
     &            xxx(id)*sinbesum(il))
            xxx(id) = (labxe(il) - (zzz(id) - bendlate(il))*sinbesum(il) +
     &                 xxx(id)*cosbesum(il))
            zzz(id) = tz
          endif
        endif
      enddo

      deallocate(ibend,bendlats,bendlate,bend_a)
      deallocate(bends_sum,cosbssum,sinbssum)
      deallocate(bende_sum,cosbesum,sinbesum)
      deallocate(labzs,labxs,labze,labxe)

      return
      end
c=============================================================================
c=============================================================================
      subroutine getpsgrd(np,xp,uxp,nw,nh,psgrd,
     &                    wmin,wmax,hmin,hmax,zl,zr,zp,uzp,slope)
      integer(ISZ):: nw,nh,np
      real(kind=8):: xp(np), zp(np), uxp(np), uzp(np)
      real(kind=8):: psgrd(0:nw,0:nh)
      real(kind=8):: slope,wmax,wmin,hmax,hmin,zl,zr

c lays particles down onto slanted mesh in phase space.  Slanted by slope.
c Inputs are position, xp, velocity, uxp, axial position and velocity, zp,
c and uzp, slope, and grid size and limits, nw, nh, wmin, wmax, hmin, hmax.
c output is psgrd.

      integer(ISZ):: ip,iw,ih
      real(kind=8):: gw,gh
      real(kind=8):: dwi,dhi,ww,wh

c set grid cell inverse sizes
      dwi = nw/(wmax - wmin)
      dhi = nh/(hmax - hmin)

c loop over particles
      do ip=1,np

c       --- if within z grid cell and not lost, include
        if (zl < zp(ip) .and. zp(ip) < zr .and. uzp(ip) /= 0.) then

c         --- find location on grid
          gw = (xp(ip) - wmin)*dwi
          gh = (uxp(ip)/uzp(ip) - xp(ip)*slope - hmin)*dhi

c         --- if within grid, accumulate
          if (gw >= 0. .and. gw <= nw .and.
     &        gh >= 0. .and. gh <= nh) then
            iw = int(gw)
            ih = int(gh)
            ww = gw - iw
            wh = gh - ih
            psgrd(iw  ,ih  ) = psgrd(iw  ,ih  ) + (1. - ww)*(1. - wh)
            psgrd(iw+1,ih  ) = psgrd(iw+1,ih  ) +       ww *(1. - wh)
            psgrd(iw  ,ih+1) = psgrd(iw  ,ih+1) + (1. - ww)*      wh
            psgrd(iw+1,ih+1) = psgrd(iw+1,ih+1) +       ww *      wh
          endif
        endif
      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine setgrid1d(np,x,nx,grid,xmin,xmax)
      integer(ISZ):: nx,np
      real(kind=8):: x(np)
      real(kind=8):: grid(0:nx)
      real(kind=8):: xmin,xmax

c Lays particles down onto 1-D mesh.

      integer(ISZ):: ip,ix
      real(kind=8):: gx,dxi,wx

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx) then
          ix = int(gx)
          wx = gx - ix
          grid(ix  ) = grid(ix  ) + (1. - wx)
          grid(ix+1) = grid(ix+1) +       wx
        endif

      enddo

      return
      end
c=============================================================================
      subroutine deposgrid1d(itask,np,x,z,nx,grid,gridcount,xmin,xmax)
      integer(ISZ):: itask,nx,np
      real(kind=8):: x(np), z(np)
      real(kind=8):: grid(0:nx), gridcount(0:nx)
      real(kind=8):: xmin,xmax

c Deposits particle data onto 1-D mesh.
c When itask is zero, grid is reset to zero, when not zero, grid is not reset.
c Also, if itask is zero, then the data is divided by the number density
c of particles.
c This allows data to be accumulated to obtain better statistics.

      integer(ISZ):: ip,ix
      real(kind=8):: gx,dxi,wx

c     --- Reset grids if requested
      if (itask == 0) then
        grid = 0.
        gridcount = 0.
      endif

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx) then
          ix = int(gx)
          wx = gx - ix
          grid(ix  ) = grid(ix  ) + z(ip)*(1. - wx)
          grid(ix+1) = grid(ix+1) + z(ip)*      wx
          gridcount(ix  ) = gridcount(ix  ) + (1. - wx)
          gridcount(ix+1) = gridcount(ix+1) +       wx
        endif

      enddo

c     --- Divide out the number of particles at each grid point, only if
c     --- itask == 0.
      if (itask == 0) then
        do ix=0,nx
          if (gridcount(ix) > 0.) then
            grid(ix) = grid(ix)/gridcount(ix)
          endif
        enddo
      endif

      return
      end
c=============================================================================
      subroutine getgrid1d(np,x,z,nx,grid,xmin,xmax)
      integer(ISZ):: nx,np
      real(kind=8):: x(np), z(np)
      real(kind=8):: grid(0:nx)
      real(kind=8):: xmin,xmax

c Gathers data from a 1-D mesh onto particle positions.

      integer(ISZ):: ip,ix
      real(kind=8):: gx,dxi,wx

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi

c       --- if within grid, gather
        if (0. <= gx .and. gx <= nx) then
          ix = int(gx)
          wx = gx - ix
          z(ip) = grid(ix  )*(1. - wx) + grid(ix+1)*wx
        endif

      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine setgrid2d(np,x,y,nx,ny,grid,xmin,xmax,ymin,ymax)
      integer(ISZ):: nx,ny,np
      real(kind=8):: x(np), y(np)
      real(kind=8):: grid(0:nx,0:ny)
      real(kind=8):: xmin,xmax,ymin,ymax

c Lays particles down onto 2-D mesh.

      integer(ISZ):: ip,ix,iy
      real(kind=8):: gx,gy,dxi,dyi,wx,wy

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gy = (y(ip) - ymin)*dyi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny) then
          ix = int(gx)
          iy = int(gy)
          wx = gx - ix
          wy = gy - iy
          grid(ix  ,iy  ) = grid(ix  ,iy  ) + (1. - wx)*(1. - wy)
          grid(ix+1,iy  ) = grid(ix+1,iy  ) +       wx *(1. - wy)
          grid(ix  ,iy+1) = grid(ix  ,iy+1) + (1. - wx)*      wy
          grid(ix+1,iy+1) = grid(ix+1,iy+1) +       wx *      wy
        endif

      enddo

      return
      end
c=============================================================================
      subroutine setgrid2dw(np,x,y,w,nx,ny,grid,xmin,xmax,ymin,ymax)
      integer(ISZ):: nx,ny,np
      real(kind=8):: x(np), y(np), w(np)
      real(kind=8):: grid(0:nx,0:ny)
      real(kind=8):: xmin,xmax,ymin,ymax

c Lays particles with variable weights down onto 2-D mesh.

      integer(ISZ):: ip,ix,iy
      real(kind=8):: gx,gy,dxi,dyi,wx,wy

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gy = (y(ip) - ymin)*dyi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny) then
          ix = int(gx)
          iy = int(gy)
          wx = gx - ix
          wy = gy - iy
          grid(ix  ,iy  ) = grid(ix  ,iy  ) + (1. - wx)*(1. - wy) * w(ip)
          grid(ix+1,iy  ) = grid(ix+1,iy  ) +       wx *(1. - wy) * w(ip)
          grid(ix  ,iy+1) = grid(ix  ,iy+1) + (1. - wx)*      wy  * w(ip)
          grid(ix+1,iy+1) = grid(ix+1,iy+1) +       wx *      wy  * w(ip)
        endif

      enddo

      return
      end
c=============================================================================
      subroutine deposgrid2d(itask,np,x,y,z,nx,ny,grid,gridcount,
     &                       xmin,xmax,ymin,ymax)
      integer(ISZ):: itask,nx,ny,np
      real(kind=8):: x(np), y(np), z(np)
      real(kind=8):: grid(0:nx,0:ny), gridcount(0:nx,0:ny)
      real(kind=8):: xmin,xmax,ymin,ymax

c Deposits particle data onto 2-D mesh.
c When itask is zero, grid is reset to zero, when not zero, grid is not reset.
c Also, if itask is zero, then the data is divided by the number density
c of particles.
c This allows data to be accumulated to obtain better statistics.

      integer(ISZ):: ip,ix,iy
      real(kind=8):: gx,gy,dxi,dyi,wx,wy

c     --- Reset grids if requested
      if (itask == 0) then
        grid = 0.
        gridcount = 0.
      endif

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gy = (y(ip) - ymin)*dyi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny) then
          ix = int(gx)
          iy = int(gy)
          wx = gx - ix
          wy = gy - iy
          grid(ix  ,iy  ) = grid(ix  ,iy  ) + z(ip)*(1. - wx)*(1. - wy)
          grid(ix+1,iy  ) = grid(ix+1,iy  ) + z(ip)*      wx *(1. - wy)
          grid(ix  ,iy+1) = grid(ix  ,iy+1) + z(ip)*(1. - wx)*      wy
          grid(ix+1,iy+1) = grid(ix+1,iy+1) + z(ip)*      wx *      wy
          gridcount(ix  ,iy  ) = gridcount(ix  ,iy  ) + (1. - wx)*(1. - wy)
          gridcount(ix+1,iy  ) = gridcount(ix+1,iy  ) +       wx *(1. - wy)
          gridcount(ix  ,iy+1) = gridcount(ix  ,iy+1) + (1. - wx)*      wy
          gridcount(ix+1,iy+1) = gridcount(ix+1,iy+1) +       wx *      wy
        endif

      enddo

c     --- Divide out the number of particles at each grid point, only if
c     --- itask == 0.
      if (itask == 0) then
        do iy=0,ny
          do ix=0,nx
            if (gridcount(ix,iy) > 0.) then
              grid(ix,iy) = grid(ix,iy)/gridcount(ix,iy)
            endif
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine deposgrid2dw(itask,np,x,y,z,w,nx,ny,grid,gridcount,
     &                       xmin,xmax,ymin,ymax)
      integer(ISZ):: itask,nx,ny,np
      real(kind=8):: x(np), y(np), z(np), w(np)
      real(kind=8):: grid(0:nx,0:ny), gridcount(0:nx,0:ny)
      real(kind=8):: xmin,xmax,ymin,ymax

c Deposits particle data onto 2-D mesh.
c When itask is zero, grid is reset to zero, when not zero, grid is not reset.
c Also, if itask is zero, then the data is divided by the number density
c of particles.
c This allows data to be accumulated to obtain better statistics.

      integer(ISZ):: ip,ix,iy
      real(kind=8):: gx,gy,dxi,dyi,wx,wy

c     --- Reset grids if requested
      if (itask == 0) then
        grid = 0.
        gridcount = 0.
      endif

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gy = (y(ip) - ymin)*dyi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny) then
          ix = int(gx)
          iy = int(gy)
          wx = gx - ix
          wy = gy - iy
          grid(ix  ,iy  ) = grid(ix  ,iy  ) + z(ip)*(1. - wx)*(1. - wy) * w(ip)
          grid(ix+1,iy  ) = grid(ix+1,iy  ) + z(ip)*      wx *(1. - wy) * w(ip)
          grid(ix  ,iy+1) = grid(ix  ,iy+1) + z(ip)*(1. - wx)*      wy  * w(ip)
          grid(ix+1,iy+1) = grid(ix+1,iy+1) + z(ip)*      wx *      wy  * w(ip)
          gridcount(ix  ,iy  ) = gridcount(ix  ,iy  ) + (1. - wx)*(1. - wy) * w(ip)
          gridcount(ix+1,iy  ) = gridcount(ix+1,iy  ) +       wx *(1. - wy) * w(ip)
          gridcount(ix  ,iy+1) = gridcount(ix  ,iy+1) + (1. - wx)*      wy  * w(ip)
          gridcount(ix+1,iy+1) = gridcount(ix+1,iy+1) +       wx *      wy  * w(ip)
        endif

      enddo

c     --- Divide out the number of particles at each grid point, only if
c     --- itask == 0.
      if (itask == 0) then
        do iy=0,ny
          do ix=0,nx
            if (gridcount(ix,iy) > 0.) then
              grid(ix,iy) = grid(ix,iy)/gridcount(ix,iy)
            endif
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine deposgridrzvect(itask,np,x,y,z,vx,vy,vz,w,nx,ny,grid,gridcount,
     &                       xmin,xmax,ymin,ymax)
      integer(ISZ):: itask,nx,ny,np
      real(kind=8):: x(np), y(np), z(np), vx(np), vy(np), vz(np), w(np)
      real(kind=8):: grid(0:nx,0:ny,3), gridcount(0:nx,0:ny)
      real(kind=8):: xmin,xmax,ymin,ymax,r,vr,vt

c Deposits particle data onto 2-D mesh.
c When itask is zero, grid is reset to zero, when not zero, grid is not reset.
c Also, if itask is zero, then the data is divided by the number density
c of particles.
c This allows data to be accumulated to obtain better statistics.

      integer(ISZ):: ip,ix,iy
      real(kind=8):: gx,gy,dxi,dyi,wx,wy

c     --- Reset grids if requested
      if (itask == 0) then
        grid = 0.
        gridcount = 0.
      endif

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        r = sqrt(x(ip)*x(ip)+y(ip)*y(ip))
        if (r*dxi>1.e-10) then
          vr =  vx(ip)*x(ip)/r+vy(ip)*y(ip)/r
          vt = -vx(ip)*y(ip)/r+vy(ip)*x(ip)/r
        else
          vr = sqrt(vx(ip)*vx(ip)+vy(ip)*vy(ip))
          vt = 0.
        end if
        gx = (r - xmin)*dxi
        gy = (z(ip) - ymin)*dyi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny) then
          ix = int(gx)
          iy = int(gy)
          wx = gx - ix
          wy = gy - iy
          grid(ix  ,iy  ,1) = grid(ix  ,iy  ,1) + vr*(1. - wx)*(1. - wy) * w(ip)
          grid(ix+1,iy  ,1) = grid(ix+1,iy  ,1) + vr*      wx *(1. - wy) * w(ip)
          grid(ix  ,iy+1,1) = grid(ix  ,iy+1,1) + vr*(1. - wx)*      wy  * w(ip)
          grid(ix+1,iy+1,1) = grid(ix+1,iy+1,1) + vr*      wx *      wy  * w(ip)
          grid(ix  ,iy  ,2) = grid(ix  ,iy  ,2) + vt*(1. - wx)*(1. - wy) * w(ip)
          grid(ix+1,iy  ,2) = grid(ix+1,iy  ,2) + vt*      wx *(1. - wy) * w(ip)
          grid(ix  ,iy+1,2) = grid(ix  ,iy+1,2) + vt*(1. - wx)*      wy  * w(ip)
          grid(ix+1,iy+1,2) = grid(ix+1,iy+1,2) + vt*      wx *      wy  * w(ip)
          grid(ix  ,iy  ,3) = grid(ix  ,iy  ,3) + vz(ip)*(1. - wx)*(1. - wy) * w(ip)
          grid(ix+1,iy  ,3) = grid(ix+1,iy  ,3) + vz(ip)*      wx *(1. - wy) * w(ip)
          grid(ix  ,iy+1,3) = grid(ix  ,iy+1,3) + vz(ip)*(1. - wx)*      wy  * w(ip)
          grid(ix+1,iy+1,3) = grid(ix+1,iy+1,3) + vz(ip)*      wx *      wy  * w(ip)
          gridcount(ix  ,iy  ) = gridcount(ix  ,iy  ) + (1. - wx)*(1. - wy) * w(ip)
          gridcount(ix+1,iy  ) = gridcount(ix+1,iy  ) +       wx *(1. - wy) * w(ip)
          gridcount(ix  ,iy+1) = gridcount(ix  ,iy+1) + (1. - wx)*      wy  * w(ip)
          gridcount(ix+1,iy+1) = gridcount(ix+1,iy+1) +       wx *      wy  * w(ip)
        endif

      enddo

c     --- Divide out the number of particles at each grid point, only if
c     --- itask == 0.
      if (itask == 0) then
        do iy=0,ny
          do ix=0,nx
            if (gridcount(ix,iy) > 0.) then
              grid(ix,iy,:) = grid(ix,iy,:)/gridcount(ix,iy)
            endif
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine getgrid2d(np,x,y,z,nx,ny,grid,xmin,xmax,ymin,ymax)
      integer(ISZ):: nx,ny,np
      real(kind=8):: x(np), y(np), z(np)
      real(kind=8):: grid(0:nx,0:ny)
      real(kind=8):: xmin,xmax,ymin,ymax

c Gathers data from a 2-D mesh onto particle positions.

      integer(ISZ):: ip,ix,iy
      real(kind=8):: gx,gy,dxi,dyi,wx,wy

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gy = (y(ip) - ymin)*dyi

c       --- if within grid, gather
        if (0. <= gx .and. gx <= nx .and.
     &      0. <= gy .and. gy <= ny) then
          ix = int(gx)
          iy = int(gy)
          wx = gx - ix
          wy = gy - iy
          z(ip) = grid(ix  ,iy  )*(1. - wx)*(1. - wy) +
     &            grid(ix+1,iy  )*      wx *(1. - wy) +
     &            grid(ix  ,iy+1)*(1. - wx)*      wy  +
     &            grid(ix+1,iy+1)*      wx *      wy
        endif

      enddo

      return
      end
c=============================================================================
      subroutine getgridngp2d(np,x,y,z,nx,ny,grid,xmin,xmax,ymin,ymax)
      integer(ISZ):: nx,ny,np
      real(kind=8):: x(np), y(np), z(np)
      real(kind=8):: grid(0:nx,0:ny)
      real(kind=8):: xmin,xmax,ymin,ymax

c Gathers data from a 2-D mesh onto particle positions.

      integer(ISZ):: ip,ix,iy
      real(kind=8):: gx,gy,dxi,dyi

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gy = (y(ip) - ymin)*dyi

c       --- if within grid, gather
        if (0. <= gx .and. gx <= nx .and.
     &      0. <= gy .and. gy <= ny) then
          ix = nint(gx)
          iy = nint(gy)
          z(ip) = grid(ix,iy)
        endif

      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine setgrid3d(np,x,y,z,nx,ny,nz,grid,xmin,xmax,ymin,ymax,zmin,zmax)
      integer(ISZ):: nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np)
      real(kind=8):: grid(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax

c Lays particles down onto 3-D mesh.

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz,dxi,dyi,dzi,wx,wy,wz

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gy = (y(ip) - ymin)*dyi
        gz = (z(ip) - zmin)*dzi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny .and.
     &      0. <= gz .and. gz < nz) then
          ix = int(gx)
          iy = int(gy)
          iz = int(gz)
          wx = gx - ix
          wy = gy - iy
          wz = gz - iz
          grid(ix  ,iy  ,iz  ) = grid(ix  ,iy  ,iz  ) + (1.-wx)*(1.-wy)*(1.-wz)
          grid(ix+1,iy  ,iz  ) = grid(ix+1,iy  ,iz  ) +     wx *(1.-wy)*(1.-wz)
          grid(ix  ,iy+1,iz  ) = grid(ix  ,iy+1,iz  ) + (1.-wx)*    wy *(1.-wz)
          grid(ix+1,iy+1,iz  ) = grid(ix+1,iy+1,iz  ) +     wx *    wy *(1.-wz)
          grid(ix  ,iy  ,iz+1) = grid(ix  ,iy  ,iz+1) + (1.-wx)*(1.-wy)*    wz
          grid(ix+1,iy  ,iz+1) = grid(ix+1,iy  ,iz+1) +     wx *(1.-wy)*    wz
          grid(ix  ,iy+1,iz+1) = grid(ix  ,iy+1,iz+1) + (1.-wx)*    wy *    wz
          grid(ix+1,iy+1,iz+1) = grid(ix+1,iy+1,iz+1) +     wx *    wy *    wz
        endif

      enddo

      return
      end
c=============================================================================
      subroutine deposgrid3d(itask,np,x,y,z,q,nx,ny,nz,grid,gridcount,
     &                       xmin,xmax,ymin,ymax,zmin,zmax)
      integer(ISZ):: itask,nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np), q(np)
      real(kind=8):: grid(0:nx,0:ny,0:nz), gridcount(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax

c Deposits particle data onto 3-D mesh.
c When itask is zero, grid is reset to zero, when not zero, grid is not reset.
c Also, if itask is zero, then the data is divided by the number density
c of particles.
c This allows data to be accumulated to obtain better statistics.

      integer(ISZ):: ip,ix,iy,iz,ixp1,iyp1,izp1
      integer(ISZ):: nxtemp,nytemp,nztemp
      real(kind=8):: gx,gy,gz,dxi,dyi,dzi,wx,wy,wz

c     --- Reset grids if requested
      if (itask == 0) then
        grid = 0.
        gridcount = 0.
      endif

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      if (ymax == ymin) then
        dyi = 1.
      else
        dyi = ny/(ymax - ymin)
      endif
      dzi = nz/(zmax - zmin)

      nxtemp = nx
      nytemp = ny
      nztemp = nz
      if (nx == 0) nxtemp = 1
      if (ny == 0) nytemp = 1
      if (nz == 0) nztemp = 1

      gx = 0.
      gy = 0.
      gz = 0.
      ix = 0
      iy = 0
      iz = 0
      ixp1 = 0
      iyp1 = 0
      izp1 = 0
      wx = 0.
      wy = 0.
      wz = 0.

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        if (nx > 0) gx = (x(ip) - xmin)*dxi
        if (ny > 0) gy = (y(ip) - ymin)*dyi
        if (nz > 0) gz = (z(ip) - zmin)*dzi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nxtemp .and.
     &      0. <= gy .and. gy < nytemp .and.
     &      0. <= gz .and. gz < nztemp) then
          if (nx > 0) then
            ix = int(gx)
            ixp1 = ix + 1
            wx = gx - ix
          endif
          if (ny > 0) then
            iy = int(gy)
            iyp1 = iy + 1
            wy = gy - iy
          endif
          if (nz > 0) then
            izp1 = iz + 1
            iz = int(gz)
            wz = gz - iz
          endif
          grid(ix  ,iy  ,iz)   = grid(ix  ,iy  ,iz)   + q(ip)*(1. - wx)*(1. - wy)*(1. - wz)
          grid(ixp1,iy  ,iz)   = grid(ixp1,iy  ,iz)   + q(ip)*      wx *(1. - wy)*(1. - wz)
          grid(ix  ,iyp1,iz)   = grid(ix  ,iyp1,iz)   + q(ip)*(1. - wx)*      wy *(1. - wz)
          grid(ixp1,iyp1,iz)   = grid(ixp1,iyp1,iz)   + q(ip)*      wx *      wy *(1. - wz)
          grid(ix  ,iy  ,izp1) = grid(ix  ,iy  ,izp1) + q(ip)*(1. - wx)*(1. - wy)*      wz
          grid(ixp1,iy  ,izp1) = grid(ixp1,iy  ,izp1) + q(ip)*      wx *(1. - wy)*      wz
          grid(ix  ,iyp1,izp1) = grid(ix  ,iyp1,izp1) + q(ip)*(1. - wx)*      wy *      wz
          grid(ixp1,iyp1,izp1) = grid(ixp1,iyp1,izp1) + q(ip)*      wx *      wy *      wz 
          gridcount(ix  ,iy  ,iz)   = gridcount(ix  ,iy  ,iz)   + (1. - wx)*(1. - wy)*(1. - wz)
          gridcount(ixp1,iy  ,iz)   = gridcount(ixp1,iy  ,iz)   +       wx *(1. - wy)*(1. - wz)
          gridcount(ix  ,iyp1,iz)   = gridcount(ix  ,iyp1,iz)   + (1. - wx)*      wy *(1. - wz)
          gridcount(ixp1,iyp1,iz)   = gridcount(ixp1,iyp1,iz)   +       wx *      wy *(1. - wz)
          gridcount(ix  ,iy  ,izp1) = gridcount(ix  ,iy  ,izp1) + (1. - wx)*(1. - wy)*      wz
          gridcount(ixp1,iy  ,izp1) = gridcount(ixp1,iy  ,izp1) +       wx *(1. - wy)*      wz
          gridcount(ix  ,iyp1,izp1) = gridcount(ix  ,iyp1,izp1) + (1. - wx)*      wy *      wz
          gridcount(ixp1,iyp1,izp1) = gridcount(ixp1,iyp1,izp1) +       wx *      wy *      wz
        endif

      enddo

c     --- Divide out the number of particles at each grid point, only if
c     --- itask == 0.
      if (itask == 0) then
        do iz=0,nz
          do iy=0,ny
            do ix=0,nx
              if (gridcount(ix,iy,iz) > 0.) then
                grid(ix,iy,iz) = grid(ix,iy,iz)/gridcount(ix,iy,iz)
              endif
            enddo
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine deposgrid3dvect(itask,np,x,y,z,vx,vy,vz,w,nx,ny,nz,grid,gridcount,
     &                       xmin,xmax,ymin,ymax,zmin,zmax)
      integer(ISZ):: itask,nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np), vx(np), vy(np), vz(np), w(np)
      real(kind=8):: grid(0:nx,0:ny,0:nz,3), gridcount(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax
      logical(ISZ):: l_radial

c Deposits particle velocities onto 3-D mesh, returning averages (or fluid) velocities on grid.
c When itask is zero, grid is reset to zero, when not zero, grid is not reset.
c Also, if itask is zero, then the data is divided by the number density
c of particles.
c This allows data to be accumulated to obtain better statistics.

      integer(ISZ):: ip,ix,iy,iz,ixp1,iyp1,izp1
      integer(ISZ):: nxtemp,nytemp,nztemp
      real(kind=8):: gx,gy,gz,dxi,dyi,dzi,wx,wy,wz

c     --- Reset grids if requested
      if (itask == 0) then
        grid = 0.
        gridcount = 0.
      endif

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      if (ymax == ymin) then
        dyi = 1.
      else
        dyi = ny/(ymax - ymin)
      endif
      dzi = nz/(zmax - zmin)

      nxtemp = nx
      nytemp = ny
      nztemp = nz
      if (nx == 0) nxtemp = 1
      if (ny == 0) nytemp = 1
      if (nz == 0) nztemp = 1

      gx = 0.
      gy = 0.
      gz = 0.
      ix = 0
      iy = 0
      iz = 0
      ixp1 = 0
      iyp1 = 0
      izp1 = 0
      wx = 0.
      wy = 0.
      wz = 0.

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        if (nx > 0) gx = (x(ip) - xmin)*dxi
        if (ny > 0) gy = (y(ip) - ymin)*dyi
        if (nz > 0) gz = (z(ip) - zmin)*dzi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nxtemp .and.
     &      0. <= gy .and. gy < nytemp .and.
     &      0. <= gz .and. gz < nztemp) then
          if (nx > 0) then
            ix = int(gx)
            ixp1 = ix + 1
            wx = gx - ix
          endif
          if (ny > 0) then
            iy = int(gy)
            iyp1 = iy + 1
            wy = gy - iy
          endif
          if (nz > 0) then
            izp1 = iz + 1
            iz = int(gz)
            wz = gz - iz
          endif
          grid(ix  ,iy  ,iz  , 1) = grid(ix  ,iy  ,iz  , 1) + w(ip)*vx(ip)*(1. - wx)*(1. - wy)*(1. - wz)
          grid(ixp1,iy  ,iz  , 1) = grid(ixp1,iy  ,iz  , 1) + w(ip)*vx(ip)*      wx *(1. - wy)*(1. - wz)
          grid(ix  ,iyp1,iz  , 1) = grid(ix  ,iyp1,iz  , 1) + w(ip)*vx(ip)*(1. - wx)*      wy *(1. - wz)
          grid(ixp1,iyp1,iz  , 1) = grid(ixp1,iyp1,iz  , 1) + w(ip)*vx(ip)*      wx *      wy *(1. - wz)
          grid(ix  ,iy  ,izp1, 1) = grid(ix  ,iy  ,izp1, 1) + w(ip)*vx(ip)*(1. - wx)*(1. - wy)*      wz
          grid(ixp1,iy  ,izp1, 1) = grid(ixp1,iy  ,izp1, 1) + w(ip)*vx(ip)*      wx *(1. - wy)*      wz
          grid(ix  ,iyp1,izp1, 1) = grid(ix  ,iyp1,izp1, 1) + w(ip)*vx(ip)*(1. - wx)*      wy *      wz
          grid(ixp1,iyp1,izp1, 1) = grid(ixp1,iyp1,izp1, 1) + w(ip)*vx(ip)*      wx *      wy *      wz 
          grid(ix  ,iy  ,iz  , 2) = grid(ix  ,iy  ,iz  , 2) + w(ip)*vy(ip)*(1. - wx)*(1. - wy)*(1. - wz)
          grid(ixp1,iy  ,iz  , 2) = grid(ixp1,iy  ,iz  , 2) + w(ip)*vy(ip)*      wx *(1. - wy)*(1. - wz)
          grid(ix  ,iyp1,iz  , 2) = grid(ix  ,iyp1,iz  , 2) + w(ip)*vy(ip)*(1. - wx)*      wy *(1. - wz)
          grid(ixp1,iyp1,iz  , 2) = grid(ixp1,iyp1,iz  , 2) + w(ip)*vy(ip)*      wx *      wy *(1. - wz)
          grid(ix  ,iy  ,izp1, 2) = grid(ix  ,iy  ,izp1, 2) + w(ip)*vy(ip)*(1. - wx)*(1. - wy)*      wz
          grid(ixp1,iy  ,izp1, 2) = grid(ixp1,iy  ,izp1, 2) + w(ip)*vy(ip)*      wx *(1. - wy)*      wz
          grid(ix  ,iyp1,izp1, 2) = grid(ix  ,iyp1,izp1, 2) + w(ip)*vy(ip)*(1. - wx)*      wy *      wz
          grid(ixp1,iyp1,izp1, 2) = grid(ixp1,iyp1,izp1, 2) + w(ip)*vy(ip)*      wx *      wy *      wz 
          grid(ix  ,iy  ,iz  , 3) = grid(ix  ,iy  ,iz  , 3) + w(ip)*vz(ip)*(1. - wx)*(1. - wy)*(1. - wz)
          grid(ixp1,iy  ,iz  , 3) = grid(ixp1,iy  ,iz  , 3) + w(ip)*vz(ip)*      wx *(1. - wy)*(1. - wz)
          grid(ix  ,iyp1,iz  , 3) = grid(ix  ,iyp1,iz  , 3) + w(ip)*vz(ip)*(1. - wx)*      wy *(1. - wz)
          grid(ixp1,iyp1,iz  , 3) = grid(ixp1,iyp1,iz  , 3) + w(ip)*vz(ip)*      wx *      wy *(1. - wz)
          grid(ix  ,iy  ,izp1, 3) = grid(ix  ,iy  ,izp1, 3) + w(ip)*vz(ip)*(1. - wx)*(1. - wy)*      wz
          grid(ixp1,iy  ,izp1, 3) = grid(ixp1,iy  ,izp1, 3) + w(ip)*vz(ip)*      wx *(1. - wy)*      wz
          grid(ix  ,iyp1,izp1, 3) = grid(ix  ,iyp1,izp1, 3) + w(ip)*vz(ip)*(1. - wx)*      wy *      wz
          grid(ixp1,iyp1,izp1, 3) = grid(ixp1,iyp1,izp1, 3) + w(ip)*vz(ip)*      wx *      wy *      wz 
          gridcount(ix  ,iy  ,iz  ) = gridcount(ix  ,iy  ,iz  ) + w(ip)*(1. - wx)*(1. - wy)*(1. - wz)
          gridcount(ixp1,iy  ,iz  ) = gridcount(ixp1,iy  ,iz  ) + w(ip)*      wx *(1. - wy)*(1. - wz)
          gridcount(ix  ,iyp1,iz  ) = gridcount(ix  ,iyp1,iz  ) + w(ip)*(1. - wx)*      wy *(1. - wz)
          gridcount(ixp1,iyp1,iz  ) = gridcount(ixp1,iyp1,iz  ) + w(ip)*      wx *      wy *(1. - wz)
          gridcount(ix  ,iy  ,izp1) = gridcount(ix  ,iy  ,izp1) + w(ip)*(1. - wx)*(1. - wy)*      wz
          gridcount(ixp1,iy  ,izp1) = gridcount(ixp1,iy  ,izp1) + w(ip)*      wx *(1. - wy)*      wz
          gridcount(ix  ,iyp1,izp1) = gridcount(ix  ,iyp1,izp1) + w(ip)*(1. - wx)*      wy *      wz
          gridcount(ixp1,iyp1,izp1) = gridcount(ixp1,iyp1,izp1) + w(ip)*      wx *      wy *      wz
        endif

      enddo

c     --- Divide out the number of particles at each grid point, only if
c     --- itask == 0.
      if (itask == 0) then
        do iz=0,nz
          do iy=0,ny
            do ix=0,nx
              if (gridcount(ix,iy,iz) > 0.) then
                grid(ix,iy,iz,:) = grid(ix,iy,iz,:)/gridcount(ix,iy,iz)
              endif
            enddo
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine getgrid3d(np,x,y,z,f,nx,ny,nz,grid,
     &                     xmin,xmax,ymin,ymax,zmin,zmax,l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np), f(np)
      real(kind=8):: grid(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax
      logical(ISZ):: l2symtry,l4symtry

c Gathers data from a 3-D mesh onto particle positions.

      integer(ISZ):: ip,ix,iy,iz,ixp1,iyp1,izp1
      real(kind=8):: xx,yy,gx,gy,gz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz
      real(kind=8):: xsymmetryplane,ysymmetryplane

c     --- These now default to zero, but could be input quantities.
      xsymmetryplane = 0.
      ysymmetryplane = 0.

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      if (ymax == ymin) then
        dyi = 1.
      else
        dyi = ny/(ymax - ymin)
      endif
      dzi = nz/(zmax - zmin)

      gx = 0.
      gy = 0.
      gz = 0.
      ix = 0
      iy = 0
      iz = 0
      ixp1 = 0
      iyp1 = 0
      izp1 = 0
      wx = 0.
      wy = 0.
      wz = 0.

c     --- loop over particles
      do ip=1,np

c       --- find location on grid, taking into account any symmetries.
        if (l4symtry) then
          xx = abs(x(ip) - xsymmetryplane) + xsymmetryplane
          yy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else if (l2symtry) then
          xx = x(ip)
          yy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else
          xx = x(ip)
          yy = y(ip)
        endif
        if (nx > 0) gx = (xx    - xmin)*dxi
        if (ny > 0) gy = (yy    - ymin)*dyi
        if (nz > 0) gz = (z(ip) - zmin)*dzi

c       --- if within grid, gather
        if (0. <= gx .and. gx <= nx .and.
     &      0. <= gy .and. gy <= ny .and.
     &      0. <= gz .and. gz <= nz) then
          if (nx > 0) then
            ix = int(gx)
            ixp1 = ix + 1
            wx = gx - ix
          endif
          if (ny > 0) then
            iy = int(gy)
            iyp1 = iy + 1
            wy = gy - iy
          endif
          if (nz > 0) then
            izp1 = iz + 1
            iz = int(gz)
            wz = gz - iz
          endif
          f(ip) = grid(ix  ,iy  ,iz  )*(1. - wx)*(1. - wy)*(1. - wz) +
     &            grid(ixp1,iy  ,iz  )*      wx *(1. - wy)*(1. - wz) +
     &            grid(ix  ,iyp1,iz  )*(1. - wx)*      wy *(1. - wz) +
     &            grid(ixp1,iyp1,iz  )*      wx *      wy *(1. - wz) +
     &            grid(ix  ,iy  ,izp1)*(1. - wx)*(1. - wy)*      wz  +
     &            grid(ixp1,iy  ,izp1)*      wx *(1. - wy)*      wz  +
     &            grid(ix  ,iyp1,izp1)*(1. - wx)*      wy *      wz  +
     &            grid(ixp1,iyp1,izp1)*      wx *      wy *      wz
        endif

      enddo

      return
      end
c=============================================================================
      subroutine getgridngp3d(np,x,y,z,f,nx,ny,nz,grid,
     &                        xmin,xmax,ymin,ymax,zmin,zmax,zgrid,
     &                        l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np), f(np)
      real(kind=8):: grid(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax,zgrid
      logical(ISZ):: l2symtry,l4symtry

c Gathers data from a 3-D mesh onto particle positions.

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz
      real(kind=8):: xsymmetryplane,ysymmetryplane

c     --- These now default to zero, but could be input quantities.
      xsymmetryplane = 0.
      ysymmetryplane = 0.

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid, taking into account any symmetries.
        if (l4symtry) then
          gx = abs(x(ip) - xsymmetryplane) + xsymmetryplane
          gy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else if (l2symtry) then
          gx = x(ip)
          gy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else
          gx = x(ip)
          gy = y(ip)
        endif
        gx = (gx    - xmin)*dxi
        gy = (gy    - ymin)*dyi
        gz = (z(ip) - zmin - zgrid)*dzi

c       --- if within grid, gather
        if (0. <= gx .and. gx <= nx .and.
     &      0. <= gy .and. gy <= ny .and.
     &      0. <= gz .and. gz <= nz) then
          ix = nint(gx)
          iy = nint(gy)
          iz = nint(gz)
          f(ip) = grid(ix,iy,iz)
        endif

      enddo

      return
      end
c=============================================================================
      subroutine getgridngp3di(np,x,y,z,f,nx,ny,nz,grid,
     &                         xmin,xmax,ymin,ymax,zmin,zmax,zgrid,
     &                         l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np)
      integer(ISZ):: grid(0:nx,0:ny,0:nz), f(np)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax,zgrid
      logical(ISZ):: l2symtry,l4symtry

c Gathers data from a 3-D mesh onto particle positions.

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz
      real(kind=8):: xsymmetryplane,ysymmetryplane

c     --- These now default to zero, but could be input quantities.
      xsymmetryplane = 0.
      ysymmetryplane = 0.

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid, taking into account any symmetries.
        if (l4symtry) then
          gx = abs(x(ip) - xsymmetryplane) + xsymmetryplane
          gy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else if (l2symtry) then
          gx = x(ip)
          gy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else
          gx = x(ip)
          gy = y(ip)
        endif
        gx = (gx    - xmin)*dxi
        gy = (gy    - ymin)*dyi
        gz = (z(ip) - zmin - zgrid)*dzi

c       --- if within grid, gather
        if (0. <= gx .and. gx <= nx .and.
     &      0. <= gy .and. gy <= ny .and.
     &      0. <= gz .and. gz <= nz) then
          ix = nint(gx)
          iy = nint(gy)
          iz = nint(gz)
          f(ip) = grid(ix,iy,iz)
        endif

      enddo

      return
      end
c=============================================================================
      subroutine gridtogrid3d(nxin,nyin,nzin,
     &                        xminin,xmaxin,yminin,ymaxin,zminin,zmaxin,gridin,
     &                        nxout,nyout,nzout,
     &                        xminout,xmaxout,yminout,ymaxout,zminout,zmaxout,
     &                        gridout)
      integer(ISZ):: nxin,nyin,nzin
      real(kind=8):: xminin,xmaxin,yminin,ymaxin,zminin,zmaxin
      real(kind=8):: gridin(0:nxin,0:nyin,0:nzin)
      integer(ISZ):: nxout,nyout,nzout
      real(kind=8):: xminout,xmaxout,yminout,ymaxout,zminout,zmaxout
      real(kind=8):: gridout(0:nxout,0:nyout,0:nzout)

c Does linear interpolation from one grid to another. Note that 1d and 2d grids
c can be passed in with two or one of the n's zero.

      real(kind=8):: xin,yin,zin
      integer(ISZ):: ixin,iyin,izin
      integer(ISZ):: ixout,iyout,izout
      real(kind=8):: dxiin,dyiin,dziin
      real(kind=8):: dxout,dyout,dzout
      real(kind=8):: wx1,wy1,wz1
      real(kind=8):: wx0,wy0,wz0
      integer(ISZ):: ixinp1,iyinp1,izinp1

c     --- If the min and max of the input grid are the same, then set the
c     --- inverse grid cell size to zero so that the index will always
c     --- be zero (independent of nin)
      if (xmaxin == xminin) then
        dxiin = 0.
      else
        dxiin = nxin/(xmaxin - xminin)
      endif
      if (ymaxin == yminin) then
        dyiin = 0.
      else
        dyiin = nyin/(ymaxin - yminin)
      endif
      if (zmaxin == zminin) then
        dziin = 0.
      else
        dziin = nzin/(zmaxin - zminin)
      endif

c     --- If nout is zero, then the output grid cell size doesn't
c     --- matter since it will only be multiplied by zero.
      if (nxout > 0) then
        dxout = (xmaxout - xminout)/nxout
      else
        dxout = 0.
      endif
      if (nyout > 0) then
        dyout = (ymaxout - yminout)/nyout
      else
        dyout = 0.
      endif
      if (nzout > 0) then
        dzout = (zmaxout - zminout)/nzout
      else
        dzout = 0.
      endif

!$OMP DO
      do izout=0,nzout

c       --- Get z position of output grid cell relative to the input mesh.
        zin = ((zminout + izout*dzout) - zminin)*dziin

c       --- If it is outside of the input mesh, then force this plane to
c       --- zero and go to the next plane.
        if (zin < 0. .or. zin > nzin) then
          gridout(:,:,izout) = 0.
          cycle
        endif

c       --- Get the grid cell number and weights
        izin = int(zin)
        izinp1 = izin + 1
        wz1 =  zin - izin
        wz0 =  1. - wz1

c       --- If the grid cell is at the upper boundary of the input mesh,
c       --- then set izinp1 to be within the mesh, but force the weight to
c       --- zero (it should actually already be zero anyway though).
        if (izin == nzin) then
          izinp1 = izin
          wz1 = 0.
        endif

        do iyout=0,nyout

c         --- See comments above for z.
          yin = ((yminout + iyout*dyout) - yminin)*dyiin

          if (yin < 0. .or. yin > nyin) then
            gridout(:,iyout,izout) = 0.
            cycle
          endif

          iyin = int(yin)
          iyinp1 = iyin + 1
          wy1 =  yin - iyin
          wy0 =  1. - wy1

          if (iyin == nyin) then
            iyinp1 = iyin
            wy1 = 0.
          endif

          do ixout=0,nxout

c           --- See comments above for z.
            xin = ((xminout + ixout*dxout) - xminin)*dxiin

            if (xin < 0. .or. xin > nxin) then
              gridout(ixout,iyout,izout) = 0.
              cycle
            endif

            ixin = int(xin)
            ixinp1 = ixin + 1
            wx1 =  xin - ixin
            wx0 =  1. - wx1

            if (ixin == nxin) then
              ixinp1 = ixin
              wx1 = 0.
            endif

            gridout(ixout,iyout,izout) = 
     &             wx0*wy0*wz0*gridin(ixin  ,iyin  ,izin  ) +
     &             wx1*wy0*wz0*gridin(ixinp1,iyin  ,izin  ) +
     &             wx0*wy1*wz0*gridin(ixin  ,iyinp1,izin  ) +
     &             wx1*wy1*wz0*gridin(ixinp1,iyinp1,izin  ) +
     &             wx0*wy0*wz1*gridin(ixin  ,iyin  ,izinp1) +
     &             wx1*wy0*wz1*gridin(ixinp1,iyin  ,izinp1) +
     &             wx0*wy1*wz1*gridin(ixin  ,iyinp1,izinp1) +
     &             wx1*wy1*wz1*gridin(ixinp1,iyinp1,izinp1)

          enddo
        enddo
      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine take2dint(a,n1,n2,i,j,n,b)
      integer(ISZ):: n1,n2,n
      integer(ISZ):: a(0:n1-1,0:n2-1)
      integer(ISZ):: i(n),j(n)
      integer(ISZ):: b(n)

      integer(ISZ):: k
      do k=1,n
        b(k) = a(i(k),j(k))
      enddo

      return
      end
c============================================================================
c============================================================================
      subroutine emitthresh(pgroup,n,threshold,is,iw,ngridw,ngridh,tepsx,tepsy)
      use ParticleGroupmodule
      use InPart
      use InDiag
      use Picglb
      use Win_Moments
      type(ParticleGroup):: pgroup
      integer(ISZ):: n
      real(kind=8):: threshold(n)
      integer(ISZ):: is,iw,ngridw,ngridh
      real(kind=8):: tepsx(n),tepsy(n)

c --- Calculates the emittance with thesholding. Particles in regions where
c --- the density is less than threshold are not included. The thresholding is
c --- done in a way which mimics experiment - the particles are binned onto a
c --- grid which has the slope subtracted (to minimize empty space).
c --- Input:
c ---   - n is number of thresholds
c ---   - threshold are threshold values
c ---   - is species of particles to include
c ---   - iw z-window to select particles
c ---   - ngridw number of grid points position is binned into
c ---   - ngridh number of grid points velocity is binned into
c --- Output:
c ---   - tepsx,tepsy

      real(kind=8):: xxpmesh(0:ngridw,0:ngridh)
      real(kind=8):: yypmesh(0:ngridw,0:ngridh)
      integer(ISZ):: i1,i2,ip,i,j,jw,jh,js
      real(kind=8):: zw1,zw2,vms,ww,wh,density
      real(kind=8):: slopex,mmaxx,dwx,dhx
      real(kind=8):: slopey,mmaxy,dwy,dhy
      real(kind=8):: maxes(4),mines(4)
      integer(ISZ):: nn(2,n)
      real(kind=8):: aves(2,5,n),te(2,n)
      real(kind=8):: uzpi
      real(kind=8):: dwxi,dhxi,dwyi,dhyi

c     --- Set some temporaries
      i1 = pgroup%ins(is)
      i2 = pgroup%ins(is)+pgroup%nps(is)-1
      zw1 = zwindows(1,iw)+zbeam
      zw2 = zwindows(2,iw)+zbeam
      js = min(is,nswind)
      slopex = (xxpbar(iw,js)-xbar(iw,js)*xpbar(iw,js))/dvnz(xrms(iw,js)**2)
      slopey = (yypbar(iw,js)-ybar(iw,js)*ypbar(iw,js))/dvnz(yrms(iw,js)**2)

c     --- Find min and max of position and velocity
      mines = +LARGEPOS
      maxes = -LARGEPOS
      do ip=i1,i2
        if (pgroup%uzp(ip) /= 0. .and.
     &      zw1 <= pgroup%zp(ip) .and. pgroup%zp(ip) <= zw2) then
          uzpi = 1./pgroup%uzp(ip)
          mines(1) = min(mines(1),pgroup%xp(ip))
          maxes(1) = max(maxes(1),pgroup%xp(ip))
          mines(2) = min(mines(2),pgroup%yp(ip))
          maxes(2) = max(maxes(2),pgroup%yp(ip))
          vms = pgroup%uxp(ip)*uzpi - pgroup%xp(ip)*slopex
          mines(3) = min(mines(3),vms)
          maxes(3) = max(maxes(3),vms)
          vms = pgroup%uyp(ip)*uzpi - pgroup%yp(ip)*slopey
          mines(4) = min(mines(4),vms)
          maxes(4) = max(maxes(4),vms)
        endif
      enddo

#ifdef MPIPARALLEL
      call parallelmaxrealarray(maxes,4)
      call parallelminrealarray(mines,4)
#endif

c     --- Extend the maximum a tiny bit to make sure that all particles are
c     --- within range (min,max]. Calculate grid cell sizes.
      maxes = maxes + (maxes - mines)*1.e-6
      dwx = (maxes(1) - mines(1))/ngridw
      dhx = (maxes(3) - mines(3))/ngridh
      dwxi = 1./dwx
      dhxi = 1./dhx
      dwy = (maxes(2) - mines(2))/ngridw
      dhy = (maxes(4) - mines(4))/ngridh
      dwyi = 1./dwy
      dhyi = 1./dhy

c     --- Bin up the data onto a 2-D grid
      xxpmesh = 0.
      yypmesh = 0.
      i1 = pgroup%ins(is)
      i2 = pgroup%ins(is) + pgroup%nps(is) - 1
      call getpsgrd(pgroup%nps(is),pgroup%xp(i1:i2),pgroup%uxp(i1:i2),
     &              ngridw,ngridh,xxpmesh,
     &              mines(1),maxes(1),mines(3),maxes(3),
     &              zw1,zw2,pgroup%zp(i1:i2),pgroup%uzp(i1:i2),slopex)
      call getpsgrd(pgroup%nps(is),pgroup%yp(i1:i2),pgroup%uyp(i1:i2),
     &              ngridw,ngridh,yypmesh,
     &              mines(2),maxes(2),mines(4),maxes(4),
     &              zw1,zw2,pgroup%zp(i1:i2),pgroup%uzp(i1:i2),slopey)

#ifdef MPIPARALLEL
      call parallelsumrealarray(xxpmesh,(1+ngridw)*(1+ngridh))
      call parallelsumrealarray(yypmesh,(1+ngridw)*(1+ngridh))
#endif
c     --- Find max of data and rescale it so the max is 1.
      mmaxx = maxval(xxpmesh)
      mmaxy = maxval(yypmesh)
      xxpmesh = xxpmesh/mmaxx
      yypmesh = yypmesh/mmaxy

c     --- Calculate average of particles in regions meeting threshold.
      nn = 0
      aves = 0.
      do ip=i1,i2
        if (pgroup%uzp(ip) /= 0. .and.
     &      zw1 <= pgroup%zp(ip) .and. pgroup%zp(ip) <= zw2) then
          uzpi = 1./pgroup%uzp(ip)
          jw = int((pgroup%xp(ip) - mines(1))*dwxi)
          ww =    ((pgroup%xp(ip) - mines(1))*dwxi) - jw
          vms = pgroup%uxp(ip)*uzpi - slopex*pgroup%xp(ip)
          jh = int((vms - mines(3))*dhxi)
          wh =    ((vms - mines(3))*dhxi) - jh
          density = xxpmesh(jw  ,jh  )*(1.-ww)*(1.-wh) +
     &              xxpmesh(jw+1,jh  )*(   ww)*(1.-wh) +
     &              xxpmesh(jw  ,jh+1)*(1.-ww)*(   wh) +
     &              xxpmesh(jw+1,jh+1)*(   ww)*(   wh)
          where (density > threshold)
            nn(1,:) = nn(1,:) + 1
            aves(1,1,:) = aves(1,1,:) + pgroup%xp(ip)
            aves(1,2,:) = aves(1,2,:) + vms
            aves(1,3,:) = aves(1,3,:) + pgroup%xp(ip)**2
            aves(1,4,:) = aves(1,4,:) + vms**2
            aves(1,5,:) = aves(1,5,:) + pgroup%xp(ip)*vms
          end where
          jw = int((pgroup%yp(ip) - mines(2))*dwyi)
          ww =    ((pgroup%yp(ip) - mines(2))*dwyi) - jw
          vms = pgroup%uyp(ip)*uzpi - slopey*pgroup%yp(ip)
          jh = int((vms - mines(4))*dhyi)
          wh =    ((vms - mines(4))*dhyi) - jh
          density = yypmesh(jw  ,jh  )*(1.-ww)*(1.-wh) +
     &              yypmesh(jw+1,jh  )*(   ww)*(1.-wh) +
     &              yypmesh(jw  ,jh+1)*(1.-ww)*(   wh) +
     &              yypmesh(jw+1,jh+1)*(   ww)*(   wh)
          where (density > threshold)
            nn(2,:) = nn(2,:) + 1
            aves(2,1,:) = aves(2,1,:) + pgroup%yp(ip)
            aves(2,2,:) = aves(2,2,:) + vms
            aves(2,3,:) = aves(2,3,:) + pgroup%yp(ip)**2
            aves(2,4,:) = aves(2,4,:) + vms**2
            aves(2,5,:) = aves(2,5,:) + pgroup%yp(ip)*vms
          end where
        endif
      enddo
#ifdef MPIPARALLEL
      call parallelsumintegerarray(nn,2*n)
      call parallelsumrealarray(aves,2*5*n)
#endif
      where (nn > 0)
        aves(:,1,:) = aves(:,1,:)/nn
        aves(:,2,:) = aves(:,2,:)/nn
        aves(:,3,:) = aves(:,3,:)/nn
        aves(:,4,:) = aves(:,4,:)/nn
        aves(:,5,:) = aves(:,5,:)/nn
      end where

c     --- Calculate emittance from the remaining particles.
      te = ((aves(:,3,:) - aves(:,1,:)**2)*(aves(:,4,:) - aves(:,2,:)**2) -
     &      (aves(:,5,:) - aves(:,1,:)*aves(:,2,:))**2)
      tepsx = 4.*sqrt(te(1,:))
      tepsy = 4.*sqrt(te(2,:))

      return
      end
c=============================================================================
c=============================================================================
c The following two subroutines were modified by S.M. Lund from 
c C.M. Celata's emitellipse and unshear for improved numerical efficiency, 
c elliminate lost particles from the accumulations, to 
c better conformity with WARP program structures, and to allow fixed 
c emittance bin specifications.     
c=============================================================================
      subroutine emitfrac(xp,uxp,uzp,np,
     &                    xbar,xpbar,xsqbar,xxpbar,xpsqbar,
     &                    fracbin,emitbin,npts,emitbinmax,
     &                    tx,txp,emitp,rwork,iwork) 
      use Constant
      integer(ISZ):: np,npts
      real(kind=8):: xp(np),uxp(np),uzp(np)
      real(kind=8):: xbar,xsqbar,xxpbar,xpsqbar,xpbar
      real(kind=8):: fracbin(0:npts),emitbin(0:npts)
      real(kind=8):: emitbinmax 
      real(kind=8):: tx(np),txp(np),emitp(np)
      real(kind=8):: rwork(3,npts)
      integer(ISZ):: iwork(npts)

c This subroutine calculates the points for a graph of emittance vs.
c the fraction of the beam particles in phase space enclosed by
c nested ellipses (with the rms equivalent beam) in phase space.
c Also, the single particle emittances are returned for each
c particle.  These emittances are calculated assuming that the instantaneous
c emittance ellipse of each particle is nested with that of the rms
c equivalent beam phase-space ellipse.
c The subroutine prin is used to rotate and translate the phase space
c ellipse first, so that it is upright in phase space.
c The emittance calculated is the unnormalized rms edge emittance
c = 4*Sqrt(<x**2><x'**2>-<x*x'>**2).  The emittance bins
c are set uniformly (npts bins) relative to emitbinmax as set by
c input.  If emitbinmax is zero, the maximum single particle emittance
c (live particle only) will be employed for the emittance bins.  Note
c that this will result in varying bin sizes as the maximum particle
c excursion in phase-space evolves.

c Variables are:
c Input:
c    xp(np) = x coordinates of particles
c   uxp(np) = gamma*vx      of particles
c   uzp(np) = gamma*vz      of particles
c   np  = number of particles
c
c   xbar    = <x>        of input particles
c   xsqbar  = <x**2>     of input particles
c   xxpbar  = <x*x'>     of input particles
c   xpsqbar = <x'**2>    of input particles
c   xpbar   = <x'>       of input particles
c
c   vbeam    = mean axial beam velocity
c   gammabar = mean axial beam gamma (paraxial)
c
c   emitbinmax = max value of emittance bins ().  If set, npts bins will be
c                generated uniformly from 0 to emitbinmax.  If emitbinmax
c                is zero, the maximum single particle emittance (live
c                particles only) will be employed to set the bins.
c
c   npts = number of nested ellipses into which to divide phase space
c
c Output:
c   fracbin(0:npts) = fraction beam within an ellipse with emittance emit(npts)
c                     nested with beam rms emittance ellipse
c   emitbin(0:npts) = Emittance of particles within an ellipse
c
c   tx(np)    = transformed x  coordinates of the particles (mixed units)
c   txp(np)   = transformed x' coordinates of the particles (mixed units)
c                   Note: This transfomation employs a translation and
c                         rotation to remove the <x>, <x'>, and <x*x'>
c                         moments of the particles. It (for simplicity)
c                         employs a mixed unit symplectic transform that
c                         perserves area.  However, direct physical
c                         interpretation of the coordinates must be taken
c                         with care.
c
c   emitp(np) = unormalized emittance of particle in x-x' phase space
c                 assuming that the phase space ellipse is nested with that
c                 of the rms equivalent beam.
c
c Work Arrays:
c   rwork(3,npts) = real(kind=8):: work array.  Used to store
c                   rwork(1,npts) = sum of tx*tx   in elliptical shell
c                   rwork(2,npts) = sum of txp*txp in elliptical shell
c                   rwork(3,npts) = sum of tx*txp  in elliptical shell
c   iwork(npts)   = integer work array. Used to store
c                   iwork(npts) = sum of particles in elliptical shell
c

      integer(ISZ):: ip,ipt,aindex
      integer(ISZ):: tnptot,nplive 
      real(kind=8):: emitmax,eps,epsinvsq
      real(kind=8):: txsqbar,txpsqbar,txxpbar,txrms,txprms
      real(kind=8):: txp2,txxp,tx2
      real(kind=8):: txp2tot,txxptot,tx2tot

c Rotate and translate the particles coordinates and moments to
c principal axes.  Note that the rotated coordinates and moments
c are expressed in mixed units.

      call prin(xp,uxp,uzp,np,xsqbar,xpsqbar,xxpbar,xbar,xpbar,
     &          tx,txp,txsqbar,txpsqbar,txxpbar)

      txrms  = sqrt(txsqbar-xbar**2)
      txprms = sqrt(txpsqbar-xpbar**2)

c     --- ellipticity of rms equivalent beam phase space ellipse
      eps = txprms/txrms

c For each particle, calculate the area/pi = single particle emittance
c within an ellipse in tx-txp phase space with the same
c ellipticity as the ellipse with axes 2*txrms and 2*txprms which passes
c through the particle's coordinates.  This information will be used to
c place the particle in the set of nested ellipses.  Note that since the
c transformation is symplectic, this mixed unit transform yields the
c correct physical area measure in x-x' phase space for the emittance.

      emitmax = 0.
      epsinvsq = 1./(eps*eps) 

      do ip = 1,np
c       --- check that particle is live
        if (uzp(ip) /= 0.) then 
          emitp(ip) = eps*( tx(ip)**2 + epsinvsq*(txp(ip)**2) )
          emitmax   = max(emitmax,emitp(ip))
        endif       
      enddo

#ifdef MPIPARALLEL
      call parallelmaxrealarray(emitmax,1)
#endif

      if (emitbinmax .gt. 0.) emitmax = emitbinmax

c Find the nested ellipses, and locate the particles within them.  Add
c each particle's moments to the moments for the elliptical shell it is
c in.  These will be summed later to find the emittance within each ellipse.

      nplive = 0 

c     --- zero moment arrays
      do ipt = 1,npts
         iwork(ipt) = 0
         rwork(1,ipt) = 0.
         rwork(2,ipt) = 0.
         rwork(3,ipt) = 0.
      enddo

c     --- accumulate moment sums
      do ip = 1,np
c       --- sample live particles only
        if (uzp(ip) /= 0.) then 
c         --- count live particles in set analyzed
          nplive = nplive + 1
c         --- find area index
          aindex = aint(npts*emitp(ip)/dvnz(emitmax)) + 1 
c         --- accumulate moments
          if (aindex .le. npts) then  
            iwork(aindex) =  iwork(aindex) + 1
            rwork(1,aindex) = rwork(1,aindex) + tx(ip)**2
            rwork(2,aindex) = rwork(2,aindex) + txp(ip)**2
            rwork(3,aindex) = rwork(3,aindex) + tx(ip)*txp(ip)
          endif 
        endif 
      enddo

#ifdef MPIPARALLEL
      call parallelsumintegerarray(nplive,1)
      call parallelsumrealarray(rwork,npts*3)
      call parallelsumintegerarray(iwork,npts)
#endif

c Find the number of particles within each nested ellipse and their emittance

       tx2tot = 0.
      txp2tot = 0.
      txxptot = 0.
       tnptot = 0.

      tx2  = 0.
      txp2 = 0.
      txxp = 0.

      fracbin(0) = 0.
      emitbin(0) = 0.
      
      do ipt = 1,npts

         tnptot =  tnptot + iwork(ipt)
         tx2tot =  tx2tot + rwork(1,ipt)
        txp2tot = txp2tot + rwork(2,ipt)
        txxptot = txxptot + rwork(3,ipt)

        if (tnptot .ne. 0) then 
          tx2  = tx2tot/real(tnptot)
          txp2 = txp2tot/real(tnptot)
          txxp = txxptot/real(tnptot)
        endif
c       --- express emittance as a fraction of *live* particles
        fracbin(ipt) = real(tnptot)/real(nplive)
        emitbin(ipt) = 4.*sqrt(tx2*txp2-txxp**2)

      enddo

      return
      end
c=============================================================================
      subroutine prin(xp,uxp,uzp,np,
     &                xsqbar,xpsqbar,xxpbar,xbar,xpbar,
     &                tx,txp,
     &                txsqbar,txpsqbar,txxpbar)
      integer(ISZ):: np
      real(kind=8):: xp(np),uxp(np),uzp(np)
      real(kind=8):: xsqbar,xpsqbar,xxpbar,xbar,xpbar
      real(kind=8):: tx(np),txp(np)
      real(kind=8):: txsqbar,txpsqbar,txxpbar

c This subroutine transforms a phase ellipse which has non-zero <xx'>,
c <x> and <x'> and rotates and translates it to produce a distribution
c with zero <x*x'>, <x>, and <x'> in the transformed (principal axis)
c coordinate system.  Note that this invloves a rotation and translation
c -- the extent of the ellipse in x is not the same as before
c the transformation.  Also, the coordinate transformation employs mixed
c units but is still symplectic (area preserving).
c
c Input variables are:
c
c   xp(np)  = coordinates of particles
c   uxp(np) = gamma*vx    of particles
c   uzp(np) = gamma*vz    of particles
c   np      = number of particles
c
c   xsqbar  = <x**2>  of particle distribution input
c   xpsqbar = <x'**2> of particle distribution input
c   xxpbar  = <x*x'>  of particle distribution input
c   xbar    = <x>     of particle distribution input
c   xpbar   = <x'>    of particle distribution input
c
c Output variables are:
c
c   tx(np)  = "coordinates" of particles rotated to principal axes
c             (mixed units)
c   txp(np) = "angles"      of particles rotated to principal axes
c             (mixed units)
c   txsqbar  = <x**2>  of rotated particle distribution (mixed units)
c   txpsqbar = <x'**2> of rotated particle distribution (mixed units)
c   txxpbar  = <x*x'>  of rotated particle distribution (mixed units)

      real(kind=8):: theta,sint,cost,sintcost,sint2,cost2
      real(kind=8):: txsqoff,txpsqoff,txxpoff
      integer(ISZ):: ip

c Find rotation angle to make ellipse upright

      theta = 0.5*atan2(2.*(xxpbar-xbar*xpbar),
     &                  dvnz(xsqbar-xpsqbar-xbar**2+xpbar**2))
c
      cost = cos(theta)
      sint = sin(theta)

c Rotate coordinates of all live particles

      do ip = 1,np
c       --- check that particle is live
        if (uzp(ip) /= 0.) then 
          tx(ip)  =  (xp(ip) - xbar)*cost + (uxp(ip)/uzp(ip) - xpbar)*sint
          txp(ip) = -(xp(ip) - xbar)*sint + (uxp(ip)/uzp(ip) - xpbar)*cost
        endif 
      enddo

c Rotate moments

      cost2 = cost**2
      sint2 = sint**2
      sintcost = sint*cost
c
      txsqoff  =  xsqbar - xbar**2
      txpsqoff = xpsqbar - xpbar**2
      txxpoff  =  xxpbar - xbar*xpbar
c
      txsqbar  = cost2*txsqoff + sint2*txpsqoff + 2*sintcost*txxpoff + xbar**2
      txpsqbar = sint2*txsqoff + cost2*txpsqoff - 2*sintcost*txxpoff + xpbar**2
      txxpbar  = sintcost*(-txsqoff+txpsqoff) + 
     &            (cost2-sint2)*txxpoff + xbar*xpbar

      return
      end
c=============================================================================
c=============================================================================
c The following two subroutines were created by C. M. Celata
c=============================================================================
      subroutine emitellipse(xpshear,vxgam,npart,xbar,xsqbar,xxpbar,
     &                       xpsqbar,xpbar,npts,vbeam,gamma,percent,emitt,
     &                       xp,xprime,area,x2sum,xp2sum,xxpsum,xsum,xpsum,
     &                       upercent,uemitt,nshell)
      integer(ISZ):: npart,npts
      real(kind=8):: xpshear(npart),vxgam(npart)
      real(kind=8):: xbar,xsqbar,xxpbar,xpsqbar,xpbar,vbeam,gamma
      real(kind=8):: percent(0:npts),emitt(0:npts)
      real(kind=8):: xp(npart),xprime(npart),area(npart)
      real(kind=8):: x2sum(npts),xp2sum(npts),xxpsum(npts),xsum(npts),xpsum(npts)
      real(kind=8):: upercent(0:npts),uemitt(0:npts)
      integer(ISZ):: nshell(npts)

c Comments:  This subroutine calculates the points for a graph of
c emittance vs. percent of the beam current enclosed by nested ellipses in
c phase space. It uses subroutine Unshear to rotate the phase space
c ellipse first, so that it is upright in phase space.  Unshear leaves the
c center of the ellipse where it is.  Emittance is assumed to be
c normalized edge rms, or 4*gamma*beta*Sqrt(<x**2><x'**2>-<xx'>**2).

c Variables are:
c Input:
c   xpshear = array of dim = npart containing the original x coordinates of the
c             particles, before the ellipse is rotated
c   vxgam = array of dimension=npart containing the gamma*vx of the particles
c   npart = number of particles
c   xbar = <x>
c   xsqbar = <x**2>
c   xxpbar = <x*x'>
c   xpsqbar =  <x'**2> of all the particles
c   xpbar = <x'>
c   npts = number of nested ellipses into which to divide phase space
c   vbeam = vz
c   gamma = relativistic factor
c Output:
c   percent = array of dim = npts containing percent beam within an ellipse
c   emitt = array of dim = npts containing the emittance of particles within
c           an ellipse
c Work Arrays:
c   xp = array of dimension=npart containing the x coordinate of the
c        particles after ellipse rotated so that it is upright
c   xprime = array of dimension=npart containing the x' of the particles
c   area
c   x2sum = sum of (x-xbar)**2 for each particle within the elliptical shell
c   xp2sum = sum of (x'-xpbar)**2 for each particle within the elliptical shell
c   xxpsum = sum of (x-xbar)*(x'-xpbar) for each particle within the elliptical
c            shell
c   xsum
c   xpsum
c   upercent
c   uemitt
c   nshell
c Local:
c   xrms = sqrt(<x**2>-<x>**2)

      integer(ISZ):: i,ipart,ipt,areaindex,ntotal
      real(kind=8):: areamax,beta,ellipticity
      real(kind=8):: xp2,xxp,x2,xrms,xprms
      real(kind=8):: x2total,xp2total,xxptotal,xsumtotal,xpsumtotal

C Unshear the ellipse

      call unshear(xpshear,vxgam,npart,xsqbar,xpsqbar,xxpbar,xbar,xpbar,
     &             xp,xprime,vbeam,gamma)

Cdo i=1,npart
Cwrite (21,700) xp(i),xprime(i)
C700  format(2(1pe13.5,1x))
Cenddo

      xrms=sqrt(xsqbar-xbar**2)
      xprms=sqrt(xpsqbar-xpbar**2)
      ellipticity=xprms/xrms
      beta=vbeam/2.997925e+8

C    write (19,890)
C890 format (" In Emittellip"/)
C       write (19,900) xrms,xprms,npts
C900    format (1h ,"xrms=",1pe13.5,1x,"xprms=",1pe13.5,1x,"npts=",i4)
C       write (19,910) vbeam,gamma
C910    format ("vbeam=",1pe13.5,1x,"gamma=",1pe13.5)
C    write (19,960)ellipticity
C960 format("ellipticity=",1pe13.5)

c For each particle, calculate the area within an ellipse with the same
c ellipticity as the ellipse with axes 2*xrms and 2*xprms which passes
c through the particle's coordinates.  This information will be used to
c place the particle in the set of nested ellipses.  Note:  what is
c actually calculated is the area divided by pi times the ellipticity.

      areamax=0.0e0

      do ipart=1,npart
      
        area(ipart)=((xp(ipart)-xbar)*ellipticity)**2+(xprime(ipart)-xpbar)**2
        areamax=max(areamax,area(ipart))
      
      end do

C Debugging
C
Cdo ipart=1,20
C
C              write (19,920) ipart, xp(ipart), xprime(ipart), area(ipart)
C920 format ("i=",i4,1x,"x=",1pe13.5,1x,"xprime=",1pe13.5,"area=",1pe13.5)
C
Cend do
C              write (19,930) areamax
C930    format ("area max=",1pe13.5)

C Zero moment arrays, percent array, emittance array, number array

      percent(0)=0.0e0
      emitt(0)=0.0e0
      do i=1,npts
        x2sum(i)=0.0e0
        xp2sum(i)=0.0e0
        xxpsum(i)=0.0e0
        xpsum(i)=0.0e0
        xsum(i)=0.0e0
        nshell(i)=0
        percent(i)=0.0e0
        emitt(i)=0.0e0
      enddo
      x2total=0.e0
      xp2total=0.e0
      xxptotal=0.0e0
      ntotal=0.0e0
      xsumtotal=0.0e0
      xpsumtotal=0.0e0

C  Zero moments

      x2=0.0e0
      xp2=0.0e0
      xxp=0.0e0

c Find the nested ellipses, and locate the particles within them.  Add
c each particle's moments to the moments for the elliptical shell it is
c in.  These will be summed later to find the emittance within each ellipse.  

      do ipart=1,npart
      
        areaindex=aint(npts*area(ipart)/areamax)+1
        if (areaindex==npts+1) areaindex=npts
        nshell(areaindex)=nshell(areaindex)+1
        xsum(areaindex)=xsum(areaindex)+xp(ipart)
        xpsum(areaindex)=xpsum(areaindex)+xprime(ipart)
        x2sum(areaindex)=x2sum(areaindex)+(xp(ipart)-xbar)**2
        xp2sum(areaindex)=xp2sum(areaindex)+(xprime(ipart)-xpbar)**2
        xxpsum(areaindex)=xxpsum(areaindex)+(xp(ipart)-xbar)*(xprime(ipart)-xpbar)

C Debugging

C       write (19,940) ipart, areaindex, nshell(areaindex)
C940 format ("i=",i4,1x,"index=", i4,1x,"number in shell=",i4)
Cwrite (19,950) xsum(areaindex),xpsum(areaindex)
C950 format ("for shell, centroid=",1pe13.5,1x,"v centroid=",1pe13.5)
Cwrite (19,960) x2sum(areaindex),xp2sum(areaindex),xxpsum(areaindex)
C960    format ("x2sum=",1pe13.5,1x,"xp2sum=",1pe13.5,1x,"xxpsum=",1pe13.5)

      end do

C Find the number of particles within each nested ellipse and their emittance

      percent(0)=0.0e0
      emitt(0)=0.0e0
      
      do ipt=1,npts
      
        x2total=x2total+x2sum(ipt)
        xp2total=xp2total+xp2sum(ipt)
        xxptotal=xxptotal+xxpsum(ipt)
        xsumtotal=xsumtotal+xsum(ipt)
        xpsumtotal=xpsumtotal+xpsum(ipt)
        ntotal=ntotal+nshell(ipt)
C       write (19,500) ipt,ntotal
C500 format (1h ,"ipt=",i3,1x,"ntotal=",i4)
        if (ntotal.ne.0) then

                x2=x2total/ntotal
                xp2=(xp2total/ntotal)
                xxp=xxptotal/ntotal
      
        end if
      
        percent(ipt)=real(ntotal)/real(npart)
        emitt(ipt)=4*gamma*beta*sqrt(x2*xp2-xxp**2)
C       write (19,550) nshell(ipt),percent(ipt),x2,xp2,emitt(ipt)
C550 format (1h ,"nshell,percent, x2,xp2,emitt=",i4,1x,4(1pe13.5,1x))
      
      end do

C Debugging

Cdo ipt=1,npts
C    write (19,799)
C799 format (" Percent and emittance=")
C    write (19,800) percent(ipt),emitt(ipt)
C800 format(1h ,1pe13.5,1x,1pe13.5)
Cend do

c Calculate the curve of emittance vs. percent for a uniform beam with
c the same rms radii

      upercent(0)=0.0e0
      uemitt(0)=0.0e0
      
      do ipt=1,npts
      
        upercent(ipt)=real(ipt)/real(npts)
        uemitt(ipt)=ipt*beta*4.0e0*xrms*xprms/npts
      
      end do

C Debugging
C      write (19,650)
C650   format ("Percent and Emittance for uniform beam")
C      do ipt=1,npts
C         write (19,655) upercent(ipt),uemitt(ipt)
C655   format(2(1pe13.5,1x))
C      end do

      return
      end
c=============================================================================
      subroutine unshear(xp,xprime,npart,xsqbar,xpsqbar,xxpbar,xbar,xpbar,xun,
     &                   xpun,vbeam,gamma)
      integer(ISZ):: npart
      real(kind=8):: xp(npart),xprime(npart),xun(npart),xpun(npart)
      real(kind=8):: xsqbar,xpsqbar,xxpbar,xbar,xpbar,vbeam,gamma

c This subroutine takes a phase ellipse which has non-zero <xx'> and
c rotates it to produce an ellipse with its axes along the coordinate axes
c and center at the same position as the original ellipse.  Note that this
c is a rotation--the extent of the ellipse in x is not the same as before
c it was rotated.  

c Input variables are:
c   xp = array of dimension = npart, containing the coordinates of all the
c        particles
c   xprime = array of dimension = npart, containing gamma*vx for all the
c            particles
c   npart = number of particles
c   xsqbar = <x**2> for the original distribution
c   xpsqbar = <xprime**2> for the original distribution
c   xxpbar = <x*xprime> for the original distribution
c   xbar = <x>
c   xpbar=<xprime>

c Output variables are:
c   xun = array of dimension = npart, containing final distribution particle
c         coordinates
c   xpun = array of dimension = npart, containing final distribution v/vz values

      real(kind=8):: numerator,denominator,theta,sint,cost,vfac
      real(kind=8):: x2ck,xp2ck,xxpck,sintcost,sint2,cost2,xsqoff,xpsqoff,xxpoff
      integer(ISZ):: i
      
      vfac=1.0e0/(gamma*vbeam)

C Debugging - input variables

C    write (19,200)
C200 format ("In Unshear")
C    write (19,220)npart,xsqbar,xpsqbar,xxpbar
C220 format (1h ,"npart=",i4,1x,"xsqbar=",1pe13.5,1x,"xpsqbar=",1pe13.5,1x,"xxpbar=",1pe13.5)
C    write (19,230)xbar,xpbar
C230 format (1h ,"xbar=",1pe13.5,1x,"xpbar=",1pe13.5)

C Check moments - debugging

      x2ck=0.0e0
      xp2ck=0.0e0
      xxpck=0.0e0
      
      do i=1,npart
        x2ck=x2ck+xp(i)**2
        xp2ck=xp2ck+xprime(i)**2
        xxpck=xxpck+xp(i)*xprime(i)
      end do
        x2ck=x2ck/npart
        xp2ck=xp2ck*vfac**2/npart
        xxpck=xxpck*vfac/npart

C    write (19,440)x2ck,xp2ck,xxpck

C Calculate the rotation angle for the ellipse

      numerator=2.0e0*(xxpbar-xbar*xpbar)
      denominator=xsqbar-xpsqbar-xbar**2+xpbar**2
    
C    write (19,260)numerator,denominator
C260 format(1h ,"numerator=",1pe13.5,1x,"denom=",1pe13.5)

C  If the ellipse is upright, don't rotate

      if (abs(denominator).lt.10**6*abs(numerator)) then

        theta=atan(numerator/denominator)/2.0e0
        cost=cos(theta)
        sint=sin(theta)

C       write (19,10)theta
C10     format (1h ,"theta=",1pe13.5)

C Rotate coordinates of all particles

      do i=1,npart
      
        xun(i)=(xp(i)-xbar)*cost+(vfac*xprime(i)-xpbar)*sint+xbar
        xpun(i)=-(xp(i)-xbar)*sint+(vfac*xprime(i)-xpbar)*cost+xpbar
      
      end do

C Rotate moments

      cost2=cost**2
      sint2=sint**2
      sintcost=sint*cost
      xsqoff=xsqbar-xbar**2
      xpsqoff=xpsqbar-xpbar**2
      xxpoff=xxpbar-xbar*xpbar
      xsqbar=cost2*xsqoff+sint2*xpsqoff+2*sintcost*xxpoff+xbar**2
      xpsqbar=sint2*xsqoff+cost2*xpsqoff-2*sintcost*xxpoff+xpbar**2
      xxpbar=sintcost*(-xsqoff+xpsqoff)+(cost2-sint2)*xxpoff+xbar*xpbar

C    write (19,400) xsqbar,xpsqbar,xxpbar
C400 format (1h ,"rotated moments x2, xp2, xxp=",3(1pe13.5,1x))

C Check moments - debugging

      x2ck=0.0e0
      xp2ck=0.0e0
      xxpck=0.0e0
      
      do i=1,npart
        x2ck=x2ck+xun(i)**2
        xp2ck=xp2ck+xpun(i)**2
        xxpck=xxpck+xun(i)*xpun(i)
      end do
        x2ck=x2ck/npart
        xp2ck=xp2ck/npart
        xxpck=xxpck/npart

C    write (19,440)x2ck,xp2ck,xxpck
C440 format(1h ,"xsq,xpsq,xxp actual moments are:",3(1pe13.5,1x))

      else
      
        do i=1,npart
          xun(i)=xp(i)
          xpun(i)=vfac*xprime(i)
        end do
      
      end if
      
      return
      end
c=============================================================================
      subroutine grid2grid(unew, nxnew, nynew, xminnew, xmaxnew, yminnew, ymaxnew,
     &                     uold, nxold, nyold, xminold, xmaxold, yminold, ymaxold)
c project field from one grid to another

      implicit none
      INTEGER(ISZ), INTENT(IN) :: nxnew, nynew, nxold, nyold
      REAL(8), INTENT(IN) :: uold(0:nxold,0:nyold)
      REAL(8), INTENT(OUT) :: unew(0:nxnew,0:nynew)
      REAL(8), INTENT(IN) :: xminold, xmaxold, yminold, ymaxold,
     &                       xminnew, xmaxnew, yminnew, ymaxnew

      INTEGER(ISZ) :: jnew, knew, j, k
      REAL(8) :: x, y, xx, yy, dxold, dyold, dxnew, dynew, invdxold, invdyold, delx, dely
      REAL(8) :: ddx(0:nxnew), ddy(0:nynew)
      INTEGER(ISZ) :: jold(0:nxnew), kold(0:nynew)

c --- computes mesh size for both grids
      dxold = (xmaxold-xminold) / nxold
      dyold = (ymaxold-yminold) / nyold
      dxnew = (xmaxnew-xminnew) / nxnew
      dynew = (ymaxnew-yminnew) / nynew

c --- check if new grid boundaries enclosed into old grid
c --- if not, issue error message
      if(xminnew<xminold .and. abs(xminnew-xminold)/dxnew>1.e-6) then
        call kaboom("Error in grid2grid: xminnew < xminold")
        return
      end if
      if(xmaxnew>xmaxold .and. abs(xmaxnew-xmaxold)/dxnew>1.e-6) then
        call kaboom("Error in grid2grid: xmaxnew > xmaxold")
        return
      end if
      if(yminnew<yminold .and. abs(yminnew-yminold)/dynew>1.e-6) then
        call kaboom("Error in grid2grid: yminnew < yminold")
        return
      end if
      if(ymaxnew>ymaxold .and. abs(ymaxnew-ymaxold)/dynew>1.e-6) then
        call kaboom("Error in grid2grid: ymaxnew > ymaxold")
        return
      end if

c --- computes temporaries

      invdxold = 1./dxold
      invdyold = 1./dyold

      do knew = 0, nynew
        y = yminnew+knew*dynew
        yy = (y-yminold) * invdyold
        kold(knew) = MIN(nyold-1,INT(yy))
        ddy(knew) = yy-real(kold(knew))
      END do
      do jnew = 0, nxnew
        x = xminnew+jnew*dxnew
        xx = (x-xminold) * invdxold
        jold(jnew) = MIN(nxold-1,INT(xx))
        ddx(jnew) = xx-real(jold(jnew))
      END do

c --- interpolate field on new grid node using linear interpolation from coarse grid

      do knew = 0, nynew
        k = kold(knew)
        dely = ddy(knew)
        do jnew = 0, nxnew
          j = jold(jnew)
          delx = ddx(jnew)
          unew(jnew,knew) = uold(j,  k)   * (1.-delx) * (1.-dely)
     &                    + uold(j+1,k)   * delx      * (1.-dely)
     &                    + uold(j,  k+1) * (1.-delx) * dely
     &                    + uold(j+1,k+1) * delx      * dely
        end do
      END do

      return
      END subroutine grid2grid
c=============================================================================
      subroutine sum_neighbors3d(fin,fout,nx,ny,nz)
      INTEGER(ISZ), INTENT(IN) :: nx, ny, nz
      INTEGER(ISZ) :: fin(nx+1,ny+1,nz+1)
      INTEGER(ISZ) :: fout(nx+1,ny+1,nz+1)

      INTEGER(ISZ) :: j,k,l,jmin,jmax,kmin,kmax,lmin,lmax

      do l = 1, nz+1
        lmin = MAX(1,l-1)
        lmax = MIN(nz+1,l+1)
        do k = 1, ny+1
          kmin = MAX(1,k-1)
          kmax = MIN(ny+1,k+1)
          do j = 1, nx+1
            jmin = MAX(1,j-1)
            jmax = MIN(nx+1,j+1)
            fout(j,k,l) = SUM(fin(jmin:jmax,kmin:kmax,lmin:lmax))
          end do
        end do
      end do

      end subroutine sum_neighbors3d
c=============================================================================
c=============================================================================
      subroutine reduceisinsidegrid(isinside,reducedisinside,nx,ny,nz)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: isinside(0:nx,0:ny,0:nz)
      real(kind=8):: reducedisinside(0:nx,0:ny,0:nz)

c reducedisinside is a copy of isinside but will be modified to remove
c redundant information. This provides an optimization of the routines
c which find intersections with conductors. Normally, a particle is
c compared against the conductors that the grid point surrounding it
c are in. If more than one of those grid points are in the same
c conductor, the particle will be checked against that conductor
c multiple times. This is a waste of CPU time. The reducing routine
c checks if a grid point is between two grid points that are in the
c same conductor as itself. If so, then the fact that the grid point
c is inside that conductor can be ignored, since particles nearby
c will get a reference to the conductor from the neighboring grid
c points. Note that the routine never ignores grid points that have
c nx,ny,nz all even.

      integer(ISZ):: ix,iy,iz

      do iz=0,nz
        do iy=0,ny
          do ix=0,nx
c           --- Don't ignore the data when nx,ny,nz are all even.
c           --- These grid points are needed to preserve the minimal
c           --- information.
            if (mod(ix,2) == 0 .and. mod(iy,2) == 0 .and. mod(iz,2) == 0) cycle

c           --- Check along each dimension if the two neighbors are inside
c           --- of the same conductor. If so, the the fact that the center
c           --- point is inside if a conductor can be ignored.
            if (ix-1 >= 0 .and. ix+1 <= nx) then
              if (isinside(ix-1,iy,iz) == isinside(ix,iy,iz) .and.
     &            isinside(ix+1,iy,iz) == isinside(ix,iy,iz)) then
                reducedisinside(ix,iy,iz) = 0.
              endif
            endif

            if (iy-1 >= 0 .and. iy+1 <= ny) then
              if (isinside(ix,iy-1,iz) == isinside(ix,iy,iz) .and.
     &            isinside(ix,iy+1,iz) == isinside(ix,iy,iz)) then
                reducedisinside(ix,iy,iz) = 0.
              endif
            endif

            if (iz-1 >= 0 .and. iz+1 <= nz) then
              if (isinside(ix,iy,iz-1) == isinside(ix,iy,iz) .and.
     &            isinside(ix,iy,iz+1) == isinside(ix,iy,iz)) then
                reducedisinside(ix,iy,iz) = 0.
              endif
            endif

          enddo
        enddo
      enddo

      return
      end subroutine reduceisinsidegrid
c=============================================================================


