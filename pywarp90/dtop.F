#include "top.h"
c============================================================================
c@(#) File DTOP.M, version $Revision: 3.98 $, $Date: 2003/07/18 20:50:53 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is part of the main package TOP of code WARP
c  This file contains the diagnostic plotting routines.
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      subroutine minidiag (itt,tim,lmoments)
      use InDiag
      use Io
      use Z_arrays
      use Win_Moments
      use Z_Moments
      use Hist
      use Moments
      integer(ISZ):: itt
      real(kind=8):: tim
      logical(ISZ):: lmoments

c  Printout particle moments and field qty, writes into history buffers, etc.
c  Send one-liner to print file, and possibly to tty.
c  Makes call to hstall to utilize the hst package.


      logical(ISZ):: lgchange,checkzhistarray

      if (warpout > -1) then
        call oneliner (warpout,itt,tim,pz,ese,ek)
      endif

c     --- Call routine to print out diagnostics.  This can either be the
c     --- builting routine 'oneliner' or can be a basis function also
c     --- called 'oneliner'.
      if (loneliner) then
        call execuser("oneliner")
      else
        if (lmoments) call oneliner (STDOUT,itt,tim,pz,ese,ek)
      endif

      if (nhist .gt. 0) then
c  Set builtin Hist arrays
        lgchange = .false.
c       --- increase the size of Hist arrays if necessary
        if (jhist == lenhist) then
          lenhist = lenhist + 100
          lgchange = .true.
        endif

c       --- Check whether any of the z histories have been turned on recently.
c       --- If so, space needs to be allocated for them.
        lgchange = checkzhistarray(lhlinechg,ihlinechg) .or. lgchange
        lgchange = checkzhistarray(lhvzofz,ihvzofz) .or. lgchange
        lgchange = checkzhistarray(lhcurrz,ihcurrz) .or. lgchange
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
        if (ihzarray == 0) checkzhistarray = .true.
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
      integer(ISZ):: iw,iz

      jhist = jhist + 1
      thist(jhist) = tim
      hzbeam(jhist) = zbeam
      hvbeam(jhist) = vbeam
      hefld(jhist) = ese
      hbmlen(jhist) = sqrt(max(SMALLPOS, zsqbar(0)-zbar(0)**2))
      hekzmbe(jhist) = ekzmbe
      hekzbeam(jhist) = ekzbeam
      hekperp(jhist) = ekperp
      hxmaxp(jhist) = xmaxp
      hxminp(jhist) = xminp
      hymaxp(jhist) = ymaxp
      hyminp(jhist) = yminp
      hzmaxp(jhist) = zmaxp
      hzminp(jhist) = zminp
      hvxmaxp(jhist) = vxmaxp
      hvxminp(jhist) = vxminp
      hvymaxp(jhist) = vymaxp
      hvyminp(jhist) = vyminp
      hvzmaxp(jhist) = vzmaxp
      hvzminp(jhist) = vzminp
      do iw=0,nzwind
        hepsx(iw,jhist)    = epsx(iw)
        hepsy(iw,jhist)    = epsy(iw)
        hepsz(iw,jhist)    = epsz(iw)
        hepsnx(iw,jhist)   = epsnx(iw)
        hepsny(iw,jhist)   = epsny(iw)
        hepsnz(iw,jhist)   = epsnz(iw)
        hepsg(iw,jhist)    = epsg(iw)
        hepsh(iw,jhist)    = epsh(iw)
        hepsng(iw,jhist)   = epsng(iw)
        hepsnh(iw,jhist)   = epsnh(iw)
        hpnum(iw,jhist)    = pnum(iw)
        hxbar(iw,jhist)    = xbar(iw)
        hybar(iw,jhist)    = ybar(iw)
        hxrms(iw,jhist)    = xrms(iw)
        hyrms(iw,jhist)    = yrms(iw)
        hrrms(iw,jhist)    = rrms(iw)
        hxprms(iw,jhist)   = xprms(iw)
        hyprms(iw,jhist)   = yprms(iw)
        hxsqbar(iw,jhist)  = xsqbar(iw)
        hysqbar(iw,jhist)  = ysqbar(iw)
        hxybar(iw,jhist)   = xybar(iw)
        hvxbar(iw,jhist)   = vxbar(iw)
        hvybar(iw,jhist)   = vybar(iw)
        hvzbar(iw,jhist)   = vzbar(iw)
        hxpbar(iw,jhist)   = xpbar(iw)
        hypbar(iw,jhist)   = ypbar(iw)
        hvxrms(iw,jhist)   = vxrms(iw)
        hvyrms(iw,jhist)   = vyrms(iw)
        hvzrms(iw,jhist)   = vzrms(iw)
        hxpsqbar(iw,jhist) = xpsqbar(iw) 
        hypsqbar(iw,jhist) = ypsqbar(iw)
        hxxpbar(iw,jhist)  = xxpbar(iw) 
        hyypbar(iw,jhist)  = yypbar(iw) 
        hxypbar(iw,jhist)  = xypbar(iw) 
        hyxpbar(iw,jhist)  = yxpbar(iw) 
        hxpypbar(iw,jhist) = xpypbar(iw) 
        hxvzbar(iw,jhist)  = xvzbar(iw) 
        hyvzbar(iw,jhist)  = yvzbar(iw) 
        hvxvzbar(iw,jhist) = vxvzbar(iw) 
        hvyvzbar(iw,jhist) = vyvzbar(iw) 
        hrhomid(iw,jhist)  = rhomid(iw)
        hrhomax(iw,jhist)  = rhomax(iw)
      enddo

c     --- Collect only those z histories specified by the user.
      do iz=0,nzzarr
        if (lhlinechg)  hlinechg(iz,jhist) = linechg(iz)
        if (lhvzofz)    hvzofz(iz,jhist) = vzofz(iz)
        if (lhcurrz)    hcurrz(iz,jhist) = curr(iz)
      enddo
      do iz=0,nzmmnt
        if (lhepsxz)    hepsxz(iz,jhist)    = epsxz(iz)
        if (lhepsyz)    hepsyz(iz,jhist)    = epsyz(iz)
        if (lhepsnxz)   hepsnxz(iz,jhist)   = epsnxz(iz)
        if (lhepsnyz)   hepsnyz(iz,jhist)   = epsnyz(iz)
        if (lhepsgz)    hepsgz(iz,jhist)    = epsgz(iz)
        if (lhepshz)    hepshz(iz,jhist)    = epshz(iz)
        if (lhepsngz)   hepsngz(iz,jhist)   = epsngz(iz)
        if (lhepsnhz)   hepsnhz(iz,jhist)   = epsnhz(iz)
        if (lhxbarz)    hxbarz(iz,jhist)    = xbarz(iz)
        if (lhybarz)    hybarz(iz,jhist)    = ybarz(iz)
        if (lhxybarz)   hxybarz(iz,jhist)   = xybarz(iz)
        if (lhxrmsz)    hxrmsz(iz,jhist)    = xrmsz(iz)
        if (lhyrmsz)    hyrmsz(iz,jhist)    = yrmsz(iz)
        if (lhrrmsz)    hrrmsz(iz,jhist)    = rrmsz(iz)
        if (lhxprmsz)   hxprmsz(iz,jhist)   = xprmsz(iz)
        if (lhyprmsz)   hyprmsz(iz,jhist)   = yprmsz(iz)
        if (lhxsqbarz)  hxsqbarz(iz,jhist)  = xsqbarz(iz)
        if (lhysqbarz)  hysqbarz(iz,jhist)  = ysqbarz(iz)
        if (lhvxbarz)   hvxbarz(iz,jhist)   = vxbarz(iz)
        if (lhvybarz)   hvybarz(iz,jhist)   = vybarz(iz)
        if (lhvzbarz)   hvzbarz(iz,jhist)   = vzbarz(iz)
        if (lhxpbarz)   hxpbarz(iz,jhist)   = xpbarz(iz)
        if (lhypbarz)   hypbarz(iz,jhist)   = ypbarz(iz)
        if (lhvxrmsz)   hvxrmsz(iz,jhist)   = vxrmsz(iz)
        if (lhvyrmsz)   hvyrmsz(iz,jhist)   = vyrmsz(iz)
        if (lhvzrmsz)   hvzrmsz(iz,jhist)   = vzrmsz(iz)
        if (lhxpsqbarz) hxpsqbarz(iz,jhist) = xpsqbarz(iz)
        if (lhypsqbarz) hypsqbarz(iz,jhist) = ypsqbarz(iz)
        if (lhxxpbarz)  hxxpbarz(iz,jhist)  = xxpbarz(iz)
        if (lhyypbarz)  hyypbarz(iz,jhist)  = yypbarz(iz)
        if (lhxypbarz)  hxypbarz(iz,jhist)  = xypbarz(iz)
        if (lhyxpbarz)  hyxpbarz(iz,jhist)  = yxpbarz(iz)
        if (lhxpypbarz) hxpypbarz(iz,jhist) = xpypbarz(iz)
        if (lhxvzbarz)  hxvzbarz(iz,jhist)  = xvzbarz(iz)
        if (lhyvzbarz)  hyvzbarz(iz,jhist)  = yvzbarz(iz)
        if (lhvxvzbarz) hvxvzbarz(iz,jhist) = vxvzbarz(iz)
        if (lhvyvzbarz) hvyvzbarz(iz,jhist) = vyvzbarz(iz)
      enddo

      return
      end
c=============================================================================
      subroutine oneliner (iunit,it,time,pz,ese,ek)
      use Io
      integer(ISZ):: iunit,it
      real(kind=8):: time,pz,ese,ek


c  Summarize the state of the system in one line

c     --- return if line is not requested
      if (verbosity .lt. 2) return

      write (iunit,901) it,time,pz,ese,ek,ese+ek
 901  format("it =",i6," time =",1pe11.4," pz =",1pe11.4,
     & " ese =",1pe11.4," ek =",1pe11.4," et  =",1pe11.4)

      return
      end
c=============================================================================
      subroutine prntpara(dxperp,dyperp,dz,lprntpara)
      use InGen
      use InPart
      use Io
      use Ch_var
      use Constant
      use Particles
      use Picglb
      use OutParams
      use Beam_acc
      real(kind=8):: dxperp,dyperp,dz
      logical(ISZ):: lprntpara

c  Calculates various parameters and optionally prints them out to a plot
c  frame and an output file or tty

      character(80):: textline

c  Calculate various quantities first

c  Current density at the center of the beam
      currdens = ibeam/dvnz(Pi*a0*b0)
c  Charge density at the center of the beam
      chrgdens = currdens/dvnz(vbeam)
c  Number density
      numbdens = chrgdens/dvnz(sq(1))
c  Plasma frequency
      omegap = sqrt(numbdens*sq(1)**2/dvnz(eps0)/dvnz(sm(1)))
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
      npreal  = np*sw(1)
      totmass = np*sw(1)*sm(1)
      totchrg = np*sw(1)*sq(1)
c  Generalized perveance
      genperv = ibeam*sq(1)/dvnz(2.*Pi*eps0*sm(1)*(vbeam*gammabar)**3)
c  Characteristic current
      charcurr = 4.*Pi*eps0*sm(1)*clight**3/dvnz(sq(1))
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
      vwave = 0.5*sqrt(gfactor)*omegap*sqrt(a0*b0)
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
      call remark  (textline)
      write (textline,20) "Charge state of ion  = ",zion," "
      call remark  (textline)
      write (textline,40) "Initial X, Y emittances = ",
     &      emitx,",  ",emity," m-rad"
      call remark  (textline)
      write (textline,40) "Initial X,Y envelope radii  = ",a0,",  ",b0," m"
      call remark  (textline)
      write (textline,40) "Initial X,Y envelope angles = ",ap0,",  ",bp0," rad"
      call remark  (textline)
      write (textline,20) "Input beam current = ",ibeam," amps"
      call remark  (textline)
      write (textline,20) "Current density = ",currdens," amps/m**2"
      call remark  (textline)
      write (textline,20) "Charge density = ",chrgdens," Coul/m**3"
      call remark  (textline)
      write (textline,20) "Number density = ",numbdens," "
      call remark  (textline)
      write (textline,20) "Plasma frequency     = ",omegap," 1/s"
      call remark  (textline)
      write (textline,20) "   times dt          = ",omegapdt," "
      call remark  (textline)
      write (textline,20) "   times quad period = ",omegaptq," "
      call remark  (textline)
      write (textline,20) "Plasma period        = ",taup," s"
      call remark  (textline)
      write (textline,40) "X-, Y-Thermal Velocities     = ",vthx,",  ",
     &                                                      vthy," m/s"
      call remark  (textline)
      write (textline,40) "   times dt                  = ",vthxdt,",  ", 
     &                                                      vthydt," m"
      call remark  (textline)
      write (textline,40) "   times dt/dx, dt/dy (X, Y) = ",vthxdtodx,",  ", 
     &                                                      vthydtody," "
      call remark  (textline)
      write (textline,40) "X-, Y-Debye Wavelengths  = ",lamdebx,",  ",
     &                                                  lamdeby," m"
      call remark  (textline)
      write (textline,40) "   over dx, dy (X and Y) = ",lamdebxodx,",  ", 
     &                                                  lamdebyody," "
      call remark  (textline)
      write (textline,20) "Longitudinal thermal velocity (rms) = ",vthz," m/s"
      call remark  (textline)
      write (textline,20) "   times dt                   = ",vthzdt," m"
      call remark  (textline)
      write (textline,20) "   times dt/dz                = ",vthzdtodz," "
      call remark  (textline)
      write (textline,20) "Longitudinal Debye wavelength = ",lamdebz," m"
      call remark  (textline)
      write (textline,20) "   over dz                    = ",lamdebzodz," "
      call remark  (textline)
c  Start a new frame since they all don't fit on one
      write (textline,20) "Beam velocity = ",vbeam," m/s"
      call remark  (textline)
      write (textline,20) "   over c     = ",vbeamoc," "
      call remark  (textline)
      write (textline,20) "Kinetic energy = ",ekin," eV"
      call remark  (textline)
      write (textline,20) "Weight of simulation particles = ",sw(1)," "
      call remark  (textline)
      write (textline,30) "Number of simulation particles = ",np," "
      call remark  (textline)
      write (textline,20) "Number of real particles = ",npreal," "
      call remark  (textline)
      write (textline,20) "Total mass = ",totmass," kg"
      call remark  (textline)
      write (textline,20) "Total charge = ",totchrg," Coul"
      call remark  (textline)
      write (textline,20) "Generalized perveance = ",genperv," "
      call remark  (textline)
      write (textline,20) "Characteristic current = ",charcurr," amps"
      call remark  (textline)
      write (textline,20) "Budker parameter = ",budker," "
      call remark  (textline)
      write (textline,20) "Timestep size dt = ",dt," s"
      call remark  (textline)
      write (textline,20) "Tune length = ",tunelen," "
      call remark  (textline)
      write (textline,40) "Undep. X-, Y-Betatron frequencies  = ",
     &                    omegab0x,",  ",omegab0y," 1/s"
      call remark  (textline)
      write (textline,40) "Undep. X-, Y-Betatron periods      = ",
     &                    taub0x,",  ",taub0y," s"
      call remark  (textline)
      write (textline,40) "Undep. X-, Y-Betatron wavelengths  = ",
     &                    lambdab0x,",  ",lambdab0y," m"
      call remark  (textline)
      write (textline,40) "Dep.   X-, Y-Betatron frequencies  = ",
     &                    omegabx,",  ",omegaby," 1/s"
      call remark  (textline)
      write (textline,40) "Dep.   X-, Y-Betatron periods      = ",
     &                    taubx,",  ",tauby," s"
      call remark  (textline)
      write (textline,40) "Dep.   X-, Y-Betatron wavelengths  = ",
     &                    lambdabx,",  ",lambdaby," m"
      write (textline,40) "X-, Y-Tune Depressions (sigma/sigma0) = ",
     &                    sigmax/dvnz(sigma0x),",  ",sigmay/dvnz(sigma0y)," "
      call remark  (textline)
      write (textline,20) "Space charge wave velocity = ",vwave," m/s"
      call remark  (textline)
      write (textline,20) "Effective wall radius = ",rwall," m"
      call remark  (textline)
      write (textline,20) "Geometric factor = ",gfactor," "
      call remark  (textline)
      write (textline,40) "X-, Y-Emittance over Space charge forces = ",
     &                                  femtxofscx,",  ",femtyofscy," "
      call remark  (textline)

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
      subroutine psubsets
      use InDiag
      use InPart
      use Pspwork
      use Particles

c  Establishes particle subsets used for phase space plots.
c  The same subsets are used for all time, since we don't want
c  movies to flicker as subsets change.


      integer(ISZ):: iset,ip,ndone,nnew
      real(kind=8):: fracplot,wranf

      do iset = 1, NSUBSETS
         if (npplot(iset).gt.0 .and. inclump(iset).gt.0) then
            if (inclump(iset)*(npsplt/inclump(iset)) .ne. npsplt) then
               call remark("PSUBSETS Error: npsplt not divisible by inclump.")
               call kaboom (0)
            endif
            fracplot = min( 1., 1.*npplot(iset)/max(1.*np_s(1), 1.) )
            ntopick(iset) = fracplot * inclump(iset) + 0.5
            if (ntopick(iset) .lt. 1) then
               call remark ("PSUBSETS Warning: ntopick < 1, set to 1")
               call remark ("  so plots are not devoid of particles.")
               call remark ("  Try a bigger clump size (inclump).")
               ntopick(iset) = 1
            endif
            do ip=1,npsplt
              isubset(ip,iset) = 0
            enddo
            do ip = 1, npsplt, inclump(iset)
               ndone = 0
               do while (ndone .lt. ntopick(iset))
                  nnew = inclump(iset) * wranf()
                  if (isubset(ip+nnew,iset) == 0) then
                     ndone = ndone + 1
                     isubset(ip+nnew,iset) = 1
                  endif
               end do
            enddo
         endif
      enddo

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
c     --- in is of the format "$Revision: 3.98 $" where 'r' and 'v' are
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

c  This is where "Basis version ..... <some person>....." gets printed.
c  We pervert it for other info.


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
c  associated with the current time step.  These appear in the top line
c  of identifying information on each plot frame.


      if (zbeam < 2000.) then
        write ( pline3,
     &   '("Step",i7,", T =",f9.4,"x10|S2|-6  s, Zbeam =",f9.4," m")' )
     &   it, time*1.e6, zbeam
      else
        write ( pline3,
     &   '("Step",i7,", T =",f9.4,"x10|S2|-6  s, Zbeam =",f9.4," km")' )
     &   it, time*1.e6, zbeam/1000.
      endif
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
         if ((it.ge.istart).and.(it.le.iend)) then
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
        if ((zstart .le. zbeamr) .and. (zbeaml .le. zend)) then
          il = int((zbeaml - zstart + zinc)/dvnz(zinc))
          ir = int((zbeamr - zstart + zinc)/dvnz(zinc))
          if (ir .gt. il) then
            thiszbeam = .true.
            return
          endif
        endif
      endif

      do i = 4,ncontrol
        if (zbeaml .le. control(i) .and. control(i) .lt. zbeamr)
     &    thiszbeam = .true.
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
c The arrays xxx(nn) and zzz(nn) are converted in place.  Currently, data
c limited to extended over at most 100 bends.
c This routine works by first finding the bends over which the data extends
c and converts the bend locations into the lab frame.  The data is then
c converted relative to the nearby bends.
c
c The arrays contain the following data...
c   bendlats     Starts of bends in Lattice frame relative to zcent
c   bendlate     Ends of bends in Lattice frame relative to zcent
c   bend_a       Bend angle of each bend
c   bend_sum     Sum of bend angles from zcent up to previous bend
c   cosbsum      Cosine of sum
c   sinbsum      Sine of sum
c   labzs        Z location of bend starts in lab frame
c   labxs        X location of bend starts in lab frame
c   labze        Z location of bend ends in lab frame
c   labxe        X location of bend ends in lab frame


      integer(ISZ):: maxnumb
      parameter (maxnumb=100)
      integer(ISZ):: ibend(maxnumb)
      real(kind=8):: bendlats(maxnumb), bendlate(maxnumb)
      real(kind=8):: bend_a(maxnumb), bend_sum(0:maxnumb+1)
      real(kind=8):: cosbsum(0:maxnumb+1), sinbsum(0:maxnumb+1)
      real(kind=8):: labzs(maxnumb), labxs(maxnumb)
      real(kind=8):: labze(maxnumb), labxe(maxnumb)
      real(kind=8):: offset,zzzmin,zzzmax,offsetmx,offsetmn,angle,aa,tz
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
        offsetmx = (int(10.+(zzzmax-zlatstrt)/zlatperi)-10)*zlatperi + zlatstrt
        offsetmn = (int(10.+(zzzmin-zlatstrt)/zlatperi)-10)*zlatperi + zlatstrt
        if (zzzmin < bendzs(ibendmn)+offsetmn) offsetmn = offsetmn - zlatperi
        if (zzzmax > bendze(ibendmx)+offsetmx) offsetmx = offsetmx + zlatperi
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
      if (zlatperi > 0.) then
        numb = ibendr + int(nbend*(offsetmx - offsetmn)/zlatperi + 0.5) -
     &         ibendl + 1
      else
        numb = ibendr - ibendl + 1
      endif
      ibend(1) = ibendl
      do i=2,numb
        ibend(i) = ibend(i-1) + 1
        if (ibend(i) > ibendmx) ibend(i) = ibendmn
      enddo

c Calculate bend starts and ends of bends in lattice frame.
      ib = ibendl
      offset = offsetmn
      do i=1,numb
        bendlats(i) = offset + bendzs(ibend(i))
        bendlate(i) = offset + bendze(ibend(i))
        ib = ib + 1
        if (ib .gt. ibendmx) then
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
      bend_sum(icentl) = 0.
      cosbsum(icentl) = cos(bend_sum(icentl))
      sinbsum(icentl) = sin(bend_sum(icentl))
      do i=icentl-1,0,-1
        bend_sum(i) = bend_sum(i+1) + bend_a(i+1)
        cosbsum(i) = cos(bend_sum(i))
        sinbsum(i) = sin(bend_sum(i))
      enddo
      bend_sum(icentr) = 0.
      cosbsum(icentr) = cos(bend_sum(icentr))
      sinbsum(icentr) = sin(bend_sum(icentr))
      do i=icentr+1,numb+1
        bend_sum(i) = bend_sum(i-1) + bend_a(i-1)
        cosbsum(i) = cos(bend_sum(i))
        sinbsum(i) = sin(bend_sum(i))
      enddo

c Calculate location of starts and ends of bends in lab frame.
c --- assumes that zcent is not in a bend
c --- left of zcent
      if (icentl > 0) then
        i = icentl
        labze(i) = bendlate(i)
        labxe(i) = 0.
        labzs(i) = labze(i) - bendrc(ibend(i))*sin(bend_a(i))
        labxs(i) = labxe(i) - bendrc(ibend(i))*(1. - cos(bend_a(i)))
        angle = bend_a(i)
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
        labzs(i) = bendlats(i)
        labxs(i) = 0.
        labze(i) = labzs(i) + bendrc(ibend(i))*sin(bend_a(i))
        labxe(i) = labxs(i) - bendrc(ibend(i))*(1. - cos(bend_a(i)))
        angle = bend_a(i)
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
            aa = bend_sum(il) + (bendlate(il) - zzz(id))/bendrc(ibend(il))
            tz = labze(il) - bendrc(ibend(il))*(sin(aa) - sinbsum(il)) -
     &           xxx(id)*sin(aa)
            xxx(id) = labxe(il) - bendrc(ibend(il))*(cosbsum(il) - cos(aa)) +
     &                xxx(id)*cos(aa)
            zzz(id) = tz
          else
            tz = labzs(ir) - (bendlats(ir) - zzz(id))*cosbsum(ir-1) -
     &           xxx(id)*sinbsum(ir-1)
            xxx(id) = labxs(ir) - (bendlats(ir)-zzz(id))*sinbsum(ir-1) +
     &                xxx(id)*cosbsum(ir-1)
            zzz(id) = tz
          endif
c       --- Do work for cases where zzz > start of bend right of zcent.
c       --- Only two cases: zzz is in bend, so il=ir; or zzz is in drift so
c       --- ir = il+1.
        elseif (zzz(id) > bendlats(icentr)) then
          if (il == ir) then
            aa = bend_sum(il) + (zzz(id) - bendlats(il))/bendrc(ibend(il))
            tz = labzs(il) + bendrc(ibend(il))*(sin(aa) - sinbsum(il)) +
     &           xxx(id)*sin(aa)
            xxx(id) = labxs(il) - bendrc(ibend(il))*(cosbsum(il) - cos(aa)) +
     &                xxx(id)*cos(aa)
            zzz(id) = tz
          else
            tz = labze(il) + (zzz(id) - bendlate(il))*cosbsum(il+1) +
     &           xxx(id)*sinbsum(il+1)
            xxx(id) = labxe(il) - (zzz(id) - bendlate(il))*sinbsum(il+1) +
     &                xxx(id)*cosbsum(il+1)
            zzz(id) = tz
          endif
        endif
      enddo

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
          if (gw >= 0. .and. gw < nw .and.
     &        gh >= 0. .and. gh < nh) then
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
        if (0. <= gx .and. gx < nx) then
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
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny) then
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
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny) then
          ix = int(gx)
          iy = int(gy)
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

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz,dxi,dyi,dzi,wx,wy,wz

c     --- Reset grids if requested
      if (itask == 0) then
        grid = 0.
        gridcount = 0.
      endif

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
          grid(ix  ,iy  ,iz)   = grid(ix  ,iy  ,iz)   + q(ip)*(1. - wx)*(1. - wy)*(1. - wz)
          grid(ix+1,iy  ,iz)   = grid(ix+1,iy  ,iz)   + q(ip)*      wx *(1. - wy)*(1. - wz)
          grid(ix  ,iy+1,iz)   = grid(ix  ,iy+1,iz)   + q(ip)*(1. - wx)*      wy *(1. - wz)
          grid(ix+1,iy+1,iz)   = grid(ix+1,iy+1,iz)   + q(ip)*      wx *      wy *(1. - wz)
          grid(ix  ,iy  ,iz+1) = grid(ix  ,iy  ,iz+1) + q(ip)*(1. - wx)*(1. - wy)*      wz
          grid(ix+1,iy  ,iz+1) = grid(ix+1,iy  ,iz+1) + q(ip)*      wx *(1. - wy)*      wz
          grid(ix  ,iy+1,iz+1) = grid(ix  ,iy+1,iz+1) + q(ip)*(1. - wx)*      wy *      wz
          grid(ix+1,iy+1,iz+1) = grid(ix+1,iy+1,iz+1) + q(ip)*      wx *      wy *      wz 
          gridcount(ix  ,iy  ,iz)   = gridcount(ix  ,iy  ,iz)   + (1. - wx)*(1. - wy)*(1. - wz)
          gridcount(ix+1,iy  ,iz)   = gridcount(ix+1,iy  ,iz)   +       wx *(1. - wy)*(1. - wz)
          gridcount(ix  ,iy+1,iz)   = gridcount(ix  ,iy+1,iz)   + (1. - wx)*      wy *(1. - wz)
          gridcount(ix+1,iy+1,iz)   = gridcount(ix+1,iy+1,iz)   +       wx *      wy *(1. - wz)
          gridcount(ix  ,iy  ,iz+1) = gridcount(ix  ,iy  ,iz+1) + (1. - wx)*(1. - wy)*      wz
          gridcount(ix+1,iy  ,iz+1) = gridcount(ix+1,iy  ,iz+1) +       wx *(1. - wy)*      wz
          gridcount(ix  ,iy+1,iz+1) = gridcount(ix  ,iy+1,iz+1) + (1. - wx)*      wy *      wz
          gridcount(ix+1,iy+1,iz+1) = gridcount(ix+1,iy+1,iz+1) +       wx *      wy *      wz
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
      subroutine getgrid3d(np,x,y,z,f,nx,ny,nz,grid,
     &                     xmin,xmax,ymin,ymax,zmin,zmax,l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np), f(np)
      real(kind=8):: grid(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax
      logical(ISZ):: l2symtry,l4symtry

c Gathers data from a 3-D mesh onto particle positions.

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        if (l4symtry) then
          gx = abs((x(ip) - xmin)*dxi)
          gy = abs((y(ip) - ymin)*dyi)
          gz =     (z(ip) - zmin)*dzi
        else if (l2symtry) then
          gx =     (x(ip) - xmin)*dxi
          gy = abs((y(ip) - ymin)*dyi)
          gz =     (z(ip) - zmin)*dzi
        else
          gx =     (x(ip) - xmin)*dxi
          gy =     (y(ip) - ymin)*dyi
          gz =     (z(ip) - zmin)*dzi
        endif

c       --- if within grid, gather
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny .and.
     &      0. <= gz .and. gz < nz) then
          ix = int(gx)
          iy = int(gy)
          iz = int(gz)
          wx = gx - ix
          wy = gy - iy
          wz = gz - iz
          f(ip) = grid(ix  ,iy  ,iz  )*(1. - wx)*(1. - wy)*(1. - wz) +
     &            grid(ix+1,iy  ,iz  )*      wx *(1. - wy)*(1. - wz) +
     &            grid(ix  ,iy+1,iz  )*(1. - wx)*      wy *(1. - wz) +
     &            grid(ix+1,iy+1,iz  )*      wx *      wy *(1. - wz) +
     &            grid(ix  ,iy  ,iz+1)*(1. - wx)*(1. - wy)*      wz  +
     &            grid(ix+1,iy  ,iz+1)*      wx *(1. - wy)*      wz  +
     &            grid(ix  ,iy+1,iz+1)*(1. - wx)*      wy *      wz  +
     &            grid(ix+1,iy+1,iz+1)*      wx *      wy *      wz
        endif

      enddo

      return
      end
c=============================================================================
      subroutine getgridngp3d(np,x,y,z,f,nx,ny,nz,grid,
     &                        xmin,xmax,ymin,ymax,zmin,zmax,l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np), f(np)
      real(kind=8):: grid(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax
      logical(ISZ):: l2symtry,l4symtry

c Gathers data from a 3-D mesh onto particle positions.

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dyi = ny/(ymax - ymin)
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        if (l4symtry) then
          gx = abs((x(ip) - xmin)*dxi)
          gy = abs((y(ip) - ymin)*dyi)
          gz =     (z(ip) - zmin)*dzi
        else if (l2symtry) then
          gx =     (x(ip) - xmin)*dxi
          gy = abs((y(ip) - ymin)*dyi)
          gz =     (z(ip) - zmin)*dzi
        else
          gx =     (x(ip) - xmin)*dxi
          gy =     (y(ip) - ymin)*dyi
          gz =     (z(ip) - zmin)*dzi
        endif

c       --- if within grid, gather
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny .and.
     &      0. <= gz .and. gz < nz) then
          ix = int(gx)
          iy = int(gy)
          iz = int(gz)
          f(ip) = grid(ix,iy,iz)
        endif

      enddo

      return
      end
c=============================================================================
c===========================================================================
      subroutine emitthresh(n,threshold,js,iw,ngridw,ngridh,tepsx,tepsy)
      use InPart
      use InDiag
      use Picglb
      use Win_Moments
      use Particles
      integer(ISZ):: n
      real(kind=8):: threshold(n)
      integer(ISZ):: js,iw,ngridw,ngridh
      real(kind=8):: tepsx(n),tepsy(n)

c --- Calculates the emittance with thesholding. Particles in regions where
c --- the density is less than threshold are not included. The thresholding is
c --- done in a way which mimics experiment - the particles are binned onto a
c --- grid which has the slope subtracted (to minimize empty space).
c --- Input:
c ---   - n is number of thresholds
c ---   - threshold are threshold values
c ---   - js species of particles to include
c ---   - iw z-window to select particles
c ---   - ngridw number of grid points position is binned into
c ---   - ngridh number of grid points velocity is binned into
c --- Output:
c ---   - tepsx,tepsy

      real(kind=8):: xxpmesh(0:ngridw,0:ngridh)
      real(kind=8):: yypmesh(0:ngridw,0:ngridh)
      integer(ISZ):: i1,i2,ip,i,j,jw,jh
      real(kind=8):: zw1,zw2,vms,ww,wh,density
      real(kind=8):: slopex,mmaxx,dwx,dhx
      real(kind=8):: slopey,mmaxy,dwy,dhy
      real(kind=8):: maxes(4),mines(4)
      integer(ISZ):: nn(2,n)
      real(kind=8):: aves(2,5,n),te(2,n)
      real(kind=8):: uzpi
      real(kind=8):: dwxi,dhxi,dwyi,dhyi

c     --- Set some temporaries
      i1 = ins(js)
      i2 = ins(js)+nps(js)-1
      zw1 = zwindows(1,iw)+zbeam
      zw2 = zwindows(2,iw)+zbeam
      slopex = (xxpbar(iw)-xbar(iw)*xpbar(iw))/dvnz(xrms(iw)**2)
      slopey = (yypbar(iw)-ybar(iw)*ypbar(iw))/dvnz(yrms(iw)**2)

c     --- Find min and max of position and velocity
      mines = +LARGEPOS
      maxes = -LARGEPOS
      do ip=i1,i2
        if (uzp(ip) /= 0. .and. zw1 <= zp(ip) .and. zp(ip) <= zw2) then
          uzpi = 1./uzp(ip)
          mines(1) = min(mines(1),xp(ip))
          maxes(1) = max(maxes(1),xp(ip))
          mines(2) = min(mines(2),yp(ip))
          maxes(2) = max(maxes(2),yp(ip))
          vms = uxp(ip)*uzpi - xp(ip)*slopex
          mines(3) = min(mines(3),vms)
          maxes(3) = max(maxes(3),vms)
          vms = uyp(ip)*uzpi - yp(ip)*slopey
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
      call getpsgrd(nps(js),xp(ins(js)),uxp(ins(js)),ngridw,ngridh,xxpmesh,
     &              mines(1),maxes(1),mines(3),maxes(3),
     &              zw1,zw2,zp(ins(js)),uzp(ins(js)),slopex)
      call getpsgrd(nps(js),yp(ins(js)),uyp(ins(js)),ngridw,ngridh,yypmesh,
     &              mines(2),maxes(2),mines(4),maxes(4),
     &              zw1,zw2,zp(ins(js)),uzp(ins(js)),slopey)

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
        if (uzp(ip) /= 0. .and. zw1 <= zp(ip) .and. zp(ip) <= zw2) then
          uzpi = 1./uzp(ip)
          jw = int((xp(ip) - mines(1))*dwxi)
          ww =    ((xp(ip) - mines(1))*dwxi) - jw
          vms = uxp(ip)*uzpi - slopex*xp(ip)
          jh = int((vms - mines(3))*dhxi)
          wh =    ((vms - mines(3))*dhxi) - jh
          density = xxpmesh(jw  ,jh  )*(1.-ww)*(1.-wh) +
     &              xxpmesh(jw+1,jh  )*(   ww)*(1.-wh) +
     &              xxpmesh(jw  ,jh+1)*(1.-ww)*(   wh) +
     &              xxpmesh(jw+1,jh+1)*(   ww)*(   wh)
          where (density > threshold)
            nn(1,:) = nn(1,:) + 1
            aves(1,1,:) = aves(1,1,:) + xp(ip)
            aves(1,2,:) = aves(1,2,:) + vms
            aves(1,3,:) = aves(1,3,:) + xp(ip)**2
            aves(1,4,:) = aves(1,4,:) + vms**2
            aves(1,5,:) = aves(1,5,:) + xp(ip)*vms
          end where
          jw = int((yp(ip) - mines(2))*dwyi)
          ww =    ((yp(ip) - mines(2))*dwyi) - jw
          vms = uyp(ip)*uzpi - slopey*yp(ip)
          jh = int((vms - mines(4))*dhyi)
          wh =    ((vms - mines(4))*dhyi) - jh
          density = yypmesh(jw  ,jh  )*(1.-ww)*(1.-wh) +
     &              yypmesh(jw+1,jh  )*(   ww)*(1.-wh) +
     &              yypmesh(jw  ,jh+1)*(1.-ww)*(   wh) +
     &              yypmesh(jw+1,jh+1)*(   ww)*(   wh)
          where (density > threshold)
            nn(2,:) = nn(2,:) + 1
            aves(2,1,:) = aves(2,1,:) + yp(ip)
            aves(2,2,:) = aves(2,2,:) + vms
            aves(2,3,:) = aves(2,3,:) + yp(ip)**2
            aves(2,4,:) = aves(2,4,:) + vms**2
            aves(2,5,:) = aves(2,5,:) + yp(ip)*vms
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
c --- if not, issue error message and stop
      if(xminnew<xminold .and. abs(xminnew-xminold)/dxnew>1.e-6) then
        write(0,*) 'Error in grid2grid: xminnew < xminold'
        stop
      end if
      if(xmaxnew>xmaxold .and. abs(xmaxnew-xmaxold)/dxnew>1.e-6) then
        write(0,*) 'Error in grid2grid: xmaxnew > xmaxold'
        stop
      end if
      if(yminnew<yminold .and. abs(yminnew-yminold)/dynew>1.e-6) then
        write(0,*) 'Error in grid2grid: yminnew < yminold'
        stop
      end if
      if(ymaxnew>ymaxold .and. abs(ymaxnew-ymaxold)/dynew>1.e-6) then
        write(0,*) 'Error in grid2grid: ymaxnew > ymaxold'
        stop
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
