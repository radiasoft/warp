#include "top.h"
c=============================================================================
c@(#) File TOP.F, version $Revision: 3.204 $, $Date: 2007/02/02 01:14:58 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  Main source file for the package TOP of the PIC code WARP.
c  It handles global version control, and loads computed global variables.
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      subroutine topinit
      use Io
      use Ch_var
      use InGen

c  Called at first reference to package.

c  Echo the package version

      call topvers (STDOUT)

c  Initialize text output file, print versions of all packages
      if (warpout > -1) then
        call outfile (warpout, "Warp PIC output files")
        call topvers (warpout)
c       call w3dvers (warpout)
c       call wrzvers (warpout)
c       call envvers (warpout)
c       call f3dvers (warpout)
c       call frzvers (warpout)
      endif

      return
      end
c=============================================================================
      subroutine topvers (iout)
      use TOPversion
c  Echoes code version,etc. to output files when they're created
      integer(ISZ):: iout
      call printpkgversion(iout,"Main package TOP",verstop)
      return
      end
c=============================================================================
      subroutine cigar(np,zunifrm,zpunifrm,z,zp,perpscal,straight,
     &                 szp,phi,zpm,zpms)
      use Constant
      integer(ISZ):: np
      real(kind=8):: straight
      real(kind=8):: zunifrm(np),zpunifrm(np),z(np),zp(np),perpscal(np)
      real(kind=8):: szp(np),phi(np),zpm(np),zpms(np)

c     This is a self-contained subroutine which calculates the normalized
c     velocities and positions for a "neuffer-type" distribution[1],
c     consisting of a uniform central region and quadratic end caps.
c     The line charge in the end caps is assumed symmetric.
c     The output transverse pseudo-radius can be used with your favorite
c     algorithm to generate the transverse distribution.  If i/epsilon
c     scaling is assumed for the current in the end caps, a single
c     solution of the envelope equations can be used to de-scale the
c     transverse variables.
c
c                                              I. Haber
c
c     Cubic solver of  x**3+a2*x**2+a1*x+a0=0 from Abromiwitz and Stegun p. 17
c     greatly simplified for the case needed here.
c     This subroutine is also vectorized?
c                                             D.P.Grote
c
c  input arguments:
c
c      zuniform   a floating point number on the interval (0.,1.) which is
c                 mapped into the longitudinal position on a beam bunch in
c                 interval (0.,1.)
c      zpunifrm   a floating point number in the interval (-.5,.5) which is
c                 mapped into a z velocity appropriate to the z position on
c                 the bunch.
c      straight   the fraction of the beam length occupied by a uniform
c                 line charge section at the bunch center.
c
c  output arguments:
c
c      z          longitudinal position of a particle in z relative to a bunch
c                 on the interval (0.,1.).
c      zp         longitudinal velocity at the position z, normalized to the
c                 maximum velocity at the bunch center.
c      perpscal   the sqrt of line charge at position z which can be used to
c                 scale transverse positions and velocities.
c
c  others:
c
c      szp, phi, zpm, zpms
c                 all scratch arrays of size np
c
c  reference:
c
c     David Neuffer, "longitudinal motion in high current ion beams - a
c     self-consistent phase space distribution with an envelope equation,"
c     IEEE Trans. Nucl. Sci. vol ns-26, June 1979, p. 3031.
c
c
c
c     initialization section.  the logic only works if local
c     storage is non-volatile.
c
      real(kind=8):: pion2,oneonpi,zs,zend,zendr,zmap,zl,zendi
      real(kind=8):: third,root3
      integer(ISZ):: i

      pion2 = pi*0.5
      oneonpi = 1./pi
      zs = straight
      zend = (1.-zs)*0.5
      zendr = 1.-zend
      zmap = zs+(4./3.)*zend
      zl = zend/3.
      zendi = 1./dvnz(zend)

c initialize stuff for cubic solver
      third = 1./3.
      root3 = sqrt(3.)
c     a1 = -3.
c     a2 = 0.
c     q = a1*third-a2**2/9.0

      do i=1,np
c
c     map the center of the beam into the uniform region
c
          z(i) = zl + zunifrm(i)*zmap
          zpm(i) = 0.
c
c     solve a cubic to map the two parabolic end caps
c
          if(z(i) < zend) then
c  Solve for the third cubic root
c            t = q**3+r**2<0 always, so three real(kind=8):: roots
             phi(i) = atan2(sqrt(1.-(1.5*(z(i)-zend)*zendi)**2),
     &                     -1.5*(z(i)-zend)*zendi)*third
             zpm(i) = -cos(phi(i))+root3*sin(phi(i))
             z(i) = (1.+zpm(i))*zend
          elseif(z(i) > zendr) then
c  Solve for the third cubic root
c            t = q**3+r**2<0 always, so three real(kind=8):: roots
             phi(i) = atan2(sqrt(1.-(1.5*(z(i)-zendr)*zendi)**2),
     &                      -1.5*(z(i)-zendr)*zendi)*third
             zpm(i) = -cos(phi(i))+root3*sin(phi(i))
             z(i) = zpm(i)*zend+zendr
          endif
c
c      newton-raphson iteration to solve for velocity mapping of the
c      interval -0.5 < zpunifrm < +0.5 .
c
          zpms(i) = 1.0-zpm(i)*zpm(i)
          perpscal(i) = sqrt(zpms(i))
c
c     initial guess (phi is a temp for zp to allow in place calculation)
c
          phi(i) = pion2*perpscal(i)*zpunifrm(i)
c
c      The four iterations of the loop are expanded for vectorization.
c      Original expression before reducing number of divides.
c         phi(i) = phi(i) - ((phi(i)*szp(i)+zpms(i)*asin(phi(i)/perpscal(i)))/
c    &            (zpms(i)*pi) - zpunifrm(i))/(twoonpi/zpms(i)*szp(i))

          szp(i) = sqrt(zpms(i)-phi(i)**2)
          phi(i) = phi(i) - ((phi(i)*szp(i)+zpms(i)*asin(phi(i)/perpscal(i)))*
     &             oneonpi - zpunifrm(i)*zpms(i))*pion2/szp(i)

          szp(i) = sqrt(zpms(i)-phi(i)**2)
          phi(i) = phi(i) - ((phi(i)*szp(i)+zpms(i)*asin(phi(i)/perpscal(i)))*
     &             oneonpi - zpunifrm(i)*zpms(i))*pion2/szp(i)

          szp(i) = sqrt(zpms(i)-phi(i)**2)
          phi(i) = phi(i) - ((phi(i)*szp(i)+zpms(i)*asin(phi(i)/perpscal(i)))*
     &             oneonpi - zpunifrm(i)*zpms(i))*pion2/szp(i)

          szp(i) = sqrt(zpms(i)-phi(i)**2)
          phi(i) = phi(i) - ((phi(i)*szp(i)+zpms(i)*asin(phi(i)/perpscal(i)))*
     &             oneonpi - zpunifrm(i)*zpms(i))*pion2/szp(i)

          zp(i) = phi(i)

      enddo

      return
      end
c=============================================================================
      subroutine species
      use InPart
      use Constant
      use Beam_acc

c Sets up species related quantities which were not set by the user, based
c on the rest of the user's input.    


      integer(ISZ):: is

c     First, set parameters for whole beam (for envelope solver).  If
c     they are unset, set to parameters of species 1.

c     --- Beam size in X
      if (a0 == 0. .and. a0_s(1) /= LARGEPOS) a0 = a0_s(1)
c     --- Beam divergence in X
      if (ap0 == 0. .and. ap0_s(1) /= LARGEPOS) ap0 = ap0_s(1)
c     --- Beam size in Y
      if (b0 == 0. .and. b0_s(1) /= LARGEPOS) b0 = b0_s(1)
c     --- Beam divergence in Y
      if (bp0 == 0. .and. bp0_s(1) /= LARGEPOS) bp0 = bp0_s(1)
c     --- Beam centroid in X
      if (x0 == 0. .and. x0_s(1) /= LARGEPOS) x0 = x0_s(1)
c     --- Beam centroid angle in X
      if (xp0 == 0. .and. xp0_s(1) /= LARGEPOS) xp0 = xp0_s(1)
c     --- Beam centroid in Y
      if (y0 == 0. .and. y0_s(1) /= LARGEPOS) y0 = y0_s(1)
c     --- Beam centroid angle in Y
      if (yp0 == 0. .and. yp0_s(1) /= LARGEPOS) yp0 = yp0_s(1)
c     --- Atomic number
      if (aion == 0. .and. aion_s(1) /= LARGEPOS) aion = aion_s(1)
c     --- Particle energy
      if (ekin == 0. .and. ekin_s(1) /= LARGEPOS) ekin = ekin_s(1)
c     --- Emittances
      if (emit == 0. .and. emit_s(1) /= LARGEPOS) emit = emit_s(1)
      if (emitx == 0. .and. emitx_s(1) /= LARGEPOS) emitx = emitx_s(1)
      if (emity == 0. .and. emity_s(1) /= LARGEPOS) emity = emity_s(1)
c     --- Normalized emittances
      if (emitn == 0. .and. emitn_s(1) /= LARGEPOS) emitn = emitn_s(1)
      if (emitnx == 0. .and. emitnx_s(1) /= LARGEPOS) emitnx = emitnx_s(1)
      if (emitny == 0. .and. emitny_s(1) /= LARGEPOS) emitny = emitny_s(1)
c     --- Current
      if (ibeam == 0. .and. ibeam_s(1) /= LARGEPOS) ibeam = ibeam_s(1)
c     --- Charge state
      if (zion == 0. .and. zion_s(1) /= LARGEPOS) zion = zion_s(1)
c     --- Axial velocity
      if (vbeam == 0. .and. vbeam_s(1) /= LARGEPOS) vbeam = vbeam_s(1)
c     --- Axial velocity tilt
      if (vtilt == 0. .and. vtilt_s(1) /= LARGEPOS) vtilt = vtilt_s(1)
c     --- Transverse thermal spread
      if (vthperp == 0. .and. vthperp_s(1) /= LARGEPOS) vthperp = vthperp_s(1)
c     --- Axial thermal spread
      if (vthz == 0. .and. vthz_s(1) /= LARGEPOS) vthz = vthz_s(1)
c     --- Minimum initial z of beam
      if (zimin == 0. .and. zimin_s(1) /= LARGEPOS) zimin = zimin_s(1)
c     --- Maximum initial z of beam
      if (zimax == 0. .and. zimax_s(1) /= LARGEPOS) zimax = zimax_s(1)
c     --- Fractional length of uniform part of beam
      if (straight == 0. .and. straight_s(1) /= LARGEPOS) straight = straight_s(1)

c     Now, set parameters for species.  If they are unset, set to the whole
c     beam parameters. Also, always reset the values of species number one
c     since the single species input parameters are kept consistent with
c     species one parameters. This alleviates many of the problems of having a
c     seperate variable for input and for internal use.

c     --- Make sure that space has been allocated for the arrays.
      if (ns > 1) then
        call gchange("InPart",0)
      endif

      do is=1,ns
c       --- Beam size in X
        if (a0_s(is) == LARGEPOS .or. is == 1) then
          a0_s(is) = a0
        endif
c       --- Beam divergence in X
        if (ap0_s(is) == LARGEPOS .or. is == 1) then
          ap0_s(is) = ap0
        endif
c       --- Beam size in Y
        if (b0_s(is) == LARGEPOS .or. is == 1) then
          b0_s(is) = b0
        endif
c       --- Beam divergence in Y
        if (bp0_s(is) == LARGEPOS .or. is == 1) then
          bp0_s(is) = bp0
        endif
c       --- Beam centroid in X
        if (x0_s(is) == LARGEPOS .or. is == 1) then
          x0_s(is) = x0
        endif
c       --- Beam centroid angle in X
        if (xp0_s(is) == LARGEPOS .or. is == 1) then
          xp0_s(is) = xp0
        endif
c       --- Beam centroid in Y
        if (y0_s(is) == LARGEPOS .or. is == 1) then
          y0_s(is) = y0
        endif
c       --- Beam centroid angle in Y
        if (yp0_s(is) == LARGEPOS .or. is == 1) then
          yp0_s(is) = yp0
        endif
c       --- Atomic number
        if (aion_s(is) == LARGEPOS .or. is == 1) then
          aion_s(is) = aion
        endif
c       --- Particle energy
        if (ekin_s(is) == LARGEPOS .or. is == 1) then
          ekin_s(is) = ekin
        endif
c       --- Emittance
        if (emit_s(is) == LARGEPOS .or. is == 1) then
          emit_s(is) = emit
        endif
        if (emitx_s(is) == LARGEPOS .or. is == 1) then
          emitx_s(is) = emitx
        endif
        if (emity_s(is) == LARGEPOS .or. is == 1) then
          emity_s(is) = emity
        endif
c       --- Normalized emittance
        if (emitn_s(is) == LARGEPOS .or. is == 1) then
          emitn_s(is) = emitn
        endif
        if (emitnx_s(is) == LARGEPOS .or. is == 1) then
          emitnx_s(is) = emitnx
        endif
        if (emitny_s(is) == LARGEPOS .or. is == 1) then
          emitny_s(is) = emitny
        endif
c       --- Current
        if (ibeam_s(is) == LARGEPOS .or. is == 1) then
          ibeam_s(is) = ibeam
        endif
c       --- Charge state
        if (zion_s(is) == LARGEPOS .or. is == 1) then
          zion_s(is) = zion
        endif
c       --- Axial velocity
        if (vbeam_s(is) == LARGEPOS .or. is == 1) then
          vbeam_s(is) = vbeam
        endif
c       --- Axial velocity tilt
        if (vtilt_s(is) == LARGEPOS .or. is == 1) then
          vtilt_s(is) = vtilt
        endif
c       --- Transverse thermal spread
        if (vthperp_s(is) == LARGEPOS .or. is == 1) then
          vthperp_s(is) = vthperp
        endif
c       --- Axial thermal spread
        if (vthz_s(is) == LARGEPOS .or. is == 1) then
          vthz_s(is) = vthz
        endif
c       --- Minimum initial z of beam
        if (zimin_s(is) == LARGEPOS .or. is == 1) then
          zimin_s(is) = zimin
        endif
c       --- Maximum initial z of beam
        if (zimax_s(is) == LARGEPOS .or. is == 1) then
          zimax_s(is) = zimax
        endif
c       --- Fractional length of uniform part of beam
        if (straight_s(is) == LARGEPOS .or. is == 1) then
          straight_s(is) = straight
        endif
      enddo

      return
      end
c=============================================================================
      subroutine derivqty
      use InPart
      use Constant
      use Beam_acc

c  Sets globally derived quantities that need computation.


      integer(ISZ):: is
      real(kind=8):: ke,u

c     Set species data 
      call species  

c     Set constants that are derived from one another
c     --- Magnetic constant = 4*pi*1.e-7
      mu0 = 4.*pi*1.e-7
c     --- Conversion factor from joules to eV is just echarge
      jperev = echarge
c     --- Epsilon_0 calculated from speed of light and mu_0
      eps0 = 1./(mu0*clight*clight)

c     --- Compute gammabar and vbeam or ekin (beam kinetic energy) from
c     --- whichever of the two quantities, vbeam or ekin, the user has
c     --- chosen to set.  (If ekin is set, then vbeam = 0, and vice-versa.)
c     --- Include both relativistic and nonrelativistic cases.
c     --- The do loop covers each species and what follows covers the
c     --- beam in general.
      if (lrelativ) then

c       --- loop over species
        do is=1,ns
          if (aion_s(is) .ne. 0.) then
            if (vbeam_s(is) == 0.) then

c             --- Beam energy in units of mc**2
              ke = jperev*ekin_s(is)/dvnz(aion_s(is)*amu*clight**2)
              gammabar = 1. + ke

c             --- The expression for vbeam was rewritten to avoid the
c             --- problem of taking the difference of nearly equal numbers
c             --- (i.e.  (1-1/gammabar) where gammabar ~ 1).
c             --- vbeam_s(is) = clight * sqrt(1.-1./gammbar**2)
              vbeam_s(is) = clight * sqrt((2*ke+ke**2)/gammabar**2)

            elseif(ekin_s(is) == 0.) then

              u = (vbeam_s(is)/clight)**2
              gammabar = 1. / sqrt (1. - u)

c             --- The expression for ekin was rewritten to avoid the
c             --- problem of taking the difference of nearly equal numbers.
c             --- ekin = (aion_s(is)*amu*clight**2)*(gammabar - 1.)/jperev
              ekin = (aion_s(is)*amu*clight**2)*(u/(sqrt(1.-u)+1.-u))/jperev

            endif
          endif

        enddo

c       --- for the beam in general
        if (vbeam == 0.) then
          ke = jperev * ekin / dvnz(aion * amu * clight**2)
          gammabar = 1. + ke
          vbeam = clight * sqrt((2*ke+ke**2)/gammabar**2)
        elseif(ekin == 0.) then
          u = (vbeam/clight)**2
          gammabar = 1. / sqrt (1. - u)
c         --- ekin = (aion * amu * clight**2) * (gammabar - 1.) / jperev
          ekin = (aion * amu * clight**2)*(u/(sqrt(1.-u)+1.-u))/jperev
        endif

      else

c       --- non-relativistic
c       --- Note that in the expression for vbeam, amu is outside of the dvnz
c       --- macro since it causes a loss of accuracy since amu is so small.
c       --- This makes the assumption that amu would never be set to zero
c       --- (it should never even be changed and should in fact be truly a
c       --- constant).

        gammabar = 1.

c       --- loop over species
        do is=1,ns
          if (aion_s(is) .ne. 0.) then
            if (vbeam_s(is) == 0.) then
              vbeam_s(is) = sqrt(2.*ekin_s(is)*jperev/dvnz(aion_s(is))/amu)
            elseif(ekin_s(is) == 0.) then
              ekin_s(is) = 0.5*(aion_s(is)*amu*vbeam_s(is)**2)/jperev
            endif
          endif
        enddo

c       --- for the beam in general
        if (vbeam == 0.) then
          vbeam = sqrt(2.*ekin*jperev/dvnz(aion)/amu)
        elseif(ekin == 0.) then
          ekin = 0.5*(aion*amu*vbeam**2)/jperev
        endif

      endif

c     --- Convert between emittance and normalized emittance
c     --- and set the x- and y- plane emittances.  Here a somewhat
c     --- confusing case structure is employed.  Note that ouside of
c     --- envelope matching routines, the code only employs emitx and emity.

c     --- Set x- and y-plane emittances from whole beam emittances if
c     --- either is set, or from the x- and y- plane normalized emittances.
c     --- This is done this way so that the user can later change emit for
c     --- example and have emitx and emity changed appropriately when
c     --- this routine is called again.
      if (emit .ne. 0.) then
        emitx = emit
        emity = emit
      elseif (emitn .ne. 0) then
        emitx = emitn*clight/dvnz(vbeam*gammabar)
        emity = emitn*clight/dvnz(vbeam*gammabar)
      elseif (emitnx .ne. 0. .or. emitny .ne. 0.) then
        emitx = emitnx*clight/dvnz(vbeam*gammabar)
        emity = emitny*clight/dvnz(vbeam*gammabar)
      endif

c     --- Loop over species with the same emittance construction above.
      do is=1,ns

c       --- Set relativistic gamma factor for this species.
        if (lrelativ) gammabar = 1./sqrt(1.-(vbeam_s(is)/clight)**2) 

c       --- (See comment above.)
        if (emit_s(is) .ne. 0.) then
          emitx_s(is) = emit_s(is)
          emity_s(is) = emit_s(is)
        elseif (emitn_s(is) .ne. 0) then
          emitx_s(is) = emitn_s(is)*clight/dvnz(vbeam_s(is)*gammabar)
          emity_s(is) = emitn_s(is)*clight/dvnz(vbeam_s(is)*gammabar)
        elseif (emitnx_s(is) .ne. 0. .or. emitny_s(is) .ne. 0.) then
          emitx_s(is) = emitnx_s(is)*clight/dvnz(vbeam_s(is)*gammabar)
          emity_s(is) = emitny_s(is)*clight/dvnz(vbeam_s(is)*gammabar)
        endif

      enddo

c     --- Reset relativistic gamma factor for the whole beam. For the
c     --- non-relativistic case, it has already been set to 1.
      if (lrelativ) gammabar = 1./sqrt(1.-(vbeam/clight)**2)

c     --- compute geometric factor needed for wave speed and (possibly) 
c     --- ears calculation 
      if (rwall > 0. .and. gfactor == 0.) then
         gfactor = log(rwall**2/dvnz(a0*b0))
      endif

      return
      end
c=============================================================================
      subroutine getzmmnt(np,xp,yp,zp,uxp,uyp,uzp,gaminv,q,m,w,dt,dtscale,
     &                    itask,nplive,
     &                    uxpo,uypo,uzpo,is,isid,ismax,maxp,minp,zmmnts0,zmmnts)
      use Constant
      use Beam_acc
      use InDiag
      use Z_Moments
      use Win_Moments
      use Moments
      use Picglb
      use ExtPart
      use Timers
      use Particles, only: wpid
      integer(ISZ):: np,itask,nplive,is,isid,ismax
      real(kind=8):: q,m,w,dt,dtscale
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np)
      real(kind=8):: maxp(6,0:nszmmnt),minp(6,0:nszmmnt)
      real(kind=8):: zmmnts0(NUMZMMNT,0:nszmmnt)
      real(kind=8):: zmmnts(0:nzmmnt,NUMZMMNT,0:nszmmnt)

c  Sets moments for species 1 in as a function of z.
c  Interpolation is done to grid centers
c  Note that lost particles are to the left of the z grid.
c  Note: "window zero" includes all particles, and no extrapolation is done
c  Note: zbar and zsqbar calculations are meaningless except for 
c        window zero (WE SHOULD DELETE ?)
c  Note: Vectorized over moments.
c        particle moment calculations are vectorized over particles
c        moments are summed into a 2-D array, vectorized over moments
c        then they are dumped into the seperate 1-D arrays in itask=3
c  Three parts
c  When itask=1  zeros out all moments
c       itask=2  sums moments from particles
c       itask=3  divides by number of particles, calculates emittances and rms

c  Note that the moments are scaled by the dtscale factor. This accounts for
c  the fact that particles with larger dtscale would normally contribute more
c  if nominal dt would be used since they would spend proportionately more
c  time in any location.

      real(kind=8),allocatable:: pp(:,:)
      integer(ISZ),allocatable:: iz(:)
      integer(ISZ):: im,izm,iw,ip,icell,i,js
      real(kind=8):: dti,dtip1,wz1,wz0,zwin
      real(kind=8):: oneondt,clighti,vzi
      real(kind=8):: pnumi
      real(kind=8):: delxsq,delxxp,delxpsq
      real(kind=8):: delysq,delyyp,delypsq
      real(kind=8):: delxy,delxyp
      real(kind=8):: delyxp,delxpyp
      real(kind=8):: delzsq,delvzsq,delzvz
      real(kind=8):: delvxsq,delvysq
      real(kind=8):: tg,th
      real(kind=8):: gamma
      real(kind=8):: timetemp,wtime
      real(kind=8):: zmmntscopy(NUMZMMNT)
#ifdef J90
      real(kind=8):: pmmnts(NUMZMMNT,np)
#endif

      timetemp = wtime()

      if (ifzmmnt == 0) return

      oneondt = 1./dvnz(dt)

c     --- If lspeciesmoments is true, then calculate the moments for each
c     --- species separately. The moments for all species combined will be
c     --- done afterwards. Otherwise, only calculate the combined moments.
      if (lspeciesmoments) then
c       --- Note that species index for these arrays is zero based.
        js = isid - 1
      else
        js = 0
      endif

      if (itask == 1 .and. .not. laccumulate_zmoments) then

         if (lspeciesmoments) then
c          --- Check if the moments are to be calculated separately for
c          --- each species. If so, check if nszmmnt already has been set
c          --- appropriately. If not, set it and allocate the arrays.
c          --- If only one species, then don't have separate species data.
           if (nszmmnt < ismax .and. ismax > 1) then
             nszmmnt = ismax
             call gchange("Z_Moments",0)
           endif
         else
           if (nszmmnt /= 0) then
             nszmmnt = 0
             call gchange("Z_Moments",0)
           endif
         endif

c        --- This is a very kludgy thing to do. This assumes that this
c        --- getzmmnt routine will always be called with the temp arrays
c        --- passed in when itask==1. The only reason the work arrays are
c        --- passed in anyway is because of the OpenMP code which requires
c        --- the work arrays to be thread private when called with itask==2.
         tempmaxp = -LARGEPOS
         tempminp = +LARGEPOS
         tempzmmnts0 = 0.
         tempzmmnts = 0.
         if (nepwin > 0) nep = 0
         if(wpid>0) nplive = 0
         zmomentscalculated = .false.
         zmomentscalculated(nszmmnt) = .true.
      endif

      if (itask == 2) then

        allocate(pp(0:17,np),iz(np))

c       --- Save the mass and weight
        zmmntsq(js) = q
        zmmntsm(js) = m
        zmmntsw(js) = w

c       --- Set flag that moments for this species has been calculated.
        zmomentscalculated(js) = .true.

c       --- Set maximum and minimum particle coordinate variables
c       --- set multiplier so lost particles don't contribute
        do ip = 1, np
          if (uzp(ip) /= 0.) then
            maxp(1,js) = max(maxp(1,js), xp(ip))
            minp(1,js) = min(minp(1,js), xp(ip))
            maxp(2,js) = max(maxp(2,js), yp(ip))
            minp(2,js) = min(minp(2,js), yp(ip))
            maxp(3,js) = max(maxp(3,js), zp(ip))
            minp(3,js) = min(minp(3,js), zp(ip))
            maxp(4,js) = max(maxp(4,js), gaminv(ip)*uxp(ip))
            minp(4,js) = min(minp(4,js), gaminv(ip)*uxp(ip))
            maxp(5,js) = max(maxp(5,js), gaminv(ip)*uyp(ip))
            minp(5,js) = min(minp(5,js), gaminv(ip)*uyp(ip))
            maxp(6,js) = max(maxp(6,js), gaminv(ip)*uzp(ip))
            minp(6,js) = min(minp(6,js), gaminv(ip)*uzp(ip))

            pp(0,ip) = dtscale
            if (ifzmmnt == 1) then
              if ((zp(ip) - zbeam) < zmmntmin .or.
     &            (zp(ip) - zbeam) > zmmntmax) pp(0,ip) = 0.
            else
              if ((zp(ip) - zbeam) <  zmmntmin .or.
     &            (zp(ip) - zbeam) >= zmmntmax) pp(0,ip) = 0.
            endif

          else
            pp(0,ip) = 0.
          endif
        enddo

c       --- set multiplier so lost particles don't contribute
c       ---   not added to total moments
c       ---   weights set to zero for z moments
c       if (ifzmmnt == 1) then
c         do ip=1,np
c           pp(0,ip) = 1.
c           if (uzp(ip) == 0. .or.
c    &          (zp(ip) - zbeam) < zmmntmin .or.
c    &          (zp(ip) - zbeam) > zmmntmax) pp(0,ip) = 0.
c         enddo
c       else
c         do ip=1,np
c           pp(0,ip) = 1.
c           if (uzp(ip) == 0. .or.
c    &          (zp(ip) - zbeam) < zmmntmin .or.
c    &          (zp(ip) - zbeam) >= zmmntmax) pp(0,ip) = 0.
c         enddo
c       endif

c       --- Accumate global moments first.

c       --- Used to add in dtscale factor below.
        if (dtscale /= 1.) zmmntscopy = zmmnts0(:,js)

c-------------------------------------------------------------------------
#ifdef J90
c       --- This method is faster for vectorized machines, primarily the J90.
c       --- put particle moments into pmmnts array for total moments
        do ip=1,np
          if (pp(0,ip) /= 0.) then
            vzi = 1./(uzp(ip)+SMALLPOS)
            pmmnts( 1,ip) = 1.
            pmmnts( 2,ip) = xp(ip)
            pmmnts( 3,ip) = yp(ip)
            pmmnts( 4,ip) = zp(ip)
            pmmnts( 7,ip) = uxp(ip)*gaminv(ip)
            pmmnts( 8,ip) = uyp(ip)*gaminv(ip)
            pmmnts( 9,ip) = uzp(ip)*gaminv(ip)
            pmmnts(10,ip) = xp(ip)*yp(ip)
            pmmnts(12,ip) = xp(ip)**2
            pmmnts(13,ip) = yp(ip)**2
            pmmnts(14,ip) = zp(ip)**2
            pmmnts(17,ip) = (gaminv(ip)*uxp(ip))**2
            pmmnts(18,ip) = (gaminv(ip)*uyp(ip))**2
            pmmnts(19,ip) = (gaminv(ip)*uzp(ip))**2
            pmmnts(22,ip) = zp(ip)*gaminv(ip)*uzp(ip)
            pmmnts(25,ip) = xp(ip)*uzp(ip)*gaminv(ip)
            pmmnts(26,ip) = yp(ip)*uzp(ip)*gaminv(ip)
            pmmnts(27,ip) = uxp(ip)*uzp(ip)*gaminv(ip)**2
            pmmnts(28,ip) = uyp(ip)*uzp(ip)*gaminv(ip)**2
            pmmnts( 5,ip) = uxp(ip)*vzi
            pmmnts( 6,ip) = uyp(ip)*vzi
            pmmnts(11,ip) = uxp(ip)*uyp(ip)*vzi**2
            pmmnts(15,ip) = (uxp(ip)*vzi)**2
            pmmnts(16,ip) = (uyp(ip)*vzi)**2
            pmmnts(20,ip) = xp(ip)*uxp(ip)*vzi
            pmmnts(21,ip) = yp(ip)*uyp(ip)*vzi
            pmmnts(23,ip) = xp(ip)*uyp(ip)*vzi
            pmmnts(24,ip) = yp(ip)*uxp(ip)*vzi
          endif
        enddo

c       --- sum total moments
c       --- Timings show that this loop runs much faster if the inner loop is
c       --- completely unrolled.  It is faster on both HP workstations and
c       --- CRAY C90.  WARNING: if NUMZMMNT changes, this loop needs to be
c       --- adjusted to reflect that change.  Note that the similar loops
c       --- below are still faster in their present form.
        do ip=1,np
          if (pp(0,ip) /= 0.) then
c           do im=1,NUMZMMNT
c             zmmnts0(im) = zmmnts0(im) + pmmnts(ip,im)
c           enddo
            zmmnts0( 1,js) = zmmnts0( 1,js) + pmmnts( 1,ip)
            zmmnts0( 2,js) = zmmnts0( 2,js) + pmmnts( 2,ip)
            zmmnts0( 3,js) = zmmnts0( 3,js) + pmmnts( 3,ip)
            zmmnts0( 4,js) = zmmnts0( 4,js) + pmmnts( 4,ip)
            zmmnts0( 5,js) = zmmnts0( 5,js) + pmmnts( 5,ip)
            zmmnts0( 6,js) = zmmnts0( 6,js) + pmmnts( 6,ip)
            zmmnts0( 7,js) = zmmnts0( 7,js) + pmmnts( 7,ip)
            zmmnts0( 8,js) = zmmnts0( 8,js) + pmmnts( 8,ip)
            zmmnts0( 9,js) = zmmnts0( 9,js) + pmmnts( 9,ip)
            zmmnts0(10,js) = zmmnts0(10,js) + pmmnts(10,ip)
            zmmnts0(11,js) = zmmnts0(11,js) + pmmnts(11,ip)
            zmmnts0(12,js) = zmmnts0(12,js) + pmmnts(12,ip)
            zmmnts0(13,js) = zmmnts0(13,js) + pmmnts(13,ip)
            zmmnts0(14,js) = zmmnts0(14,js) + pmmnts(14,ip)
            zmmnts0(15,js) = zmmnts0(15,js) + pmmnts(15,ip)
            zmmnts0(16,js) = zmmnts0(16,js) + pmmnts(16,ip)
            zmmnts0(17,js) = zmmnts0(17,js) + pmmnts(17,ip)
            zmmnts0(18,js) = zmmnts0(18,js) + pmmnts(18,ip)
            zmmnts0(19,js) = zmmnts0(19,js) + pmmnts(19,ip)
            zmmnts0(20,js) = zmmnts0(20,js) + pmmnts(20,ip)
            zmmnts0(21,js) = zmmnts0(21,js) + pmmnts(21,ip)
            zmmnts0(22,js) = zmmnts0(22,js) + pmmnts(22,ip)
            zmmnts0(23,js) = zmmnts0(23,js) + pmmnts(23,ip) 
            zmmnts0(24,js) = zmmnts0(24,js) + pmmnts(24,ip) 
            zmmnts0(25,js) = zmmnts0(25,js) + pmmnts(25,ip) 
            zmmnts0(26,js) = zmmnts0(26,js) + pmmnts(26,ip) 
            zmmnts0(27,js) = zmmnts0(27,js) + pmmnts(27,ip) 
            zmmnts0(28,js) = zmmnts0(28,js) + pmmnts(28,ip) 
          endif
        enddo

#else

c-------------------------------------------------------------------------
c       --- This method is faster for RISC and pentium machines.
c       --- sum particle moments into total moments
        do ip=1,np
          if (pp(0,ip) /= 0.) then
            vzi = 1./(uzp(ip)+SMALLPOS)
            zmmnts0( 1,js) = zmmnts0( 1,js) + 1.
            zmmnts0( 2,js) = zmmnts0( 2,js) + xp(ip)
            zmmnts0( 3,js) = zmmnts0( 3,js) + yp(ip)
            zmmnts0( 4,js) = zmmnts0( 4,js) + zp(ip)
            zmmnts0( 7,js) = zmmnts0( 7,js) + uxp(ip)*gaminv(ip)
            zmmnts0( 8,js) = zmmnts0( 8,js) + uyp(ip)*gaminv(ip)
            zmmnts0( 9,js) = zmmnts0( 9,js) + uzp(ip)*gaminv(ip)
            zmmnts0(10,js) = zmmnts0(10,js) + xp(ip)*yp(ip)
            zmmnts0(12,js) = zmmnts0(12,js) + xp(ip)**2
            zmmnts0(13,js) = zmmnts0(13,js) + yp(ip)**2
            zmmnts0(14,js) = zmmnts0(14,js) + zp(ip)**2
            zmmnts0(17,js) = zmmnts0(17,js) + (gaminv(ip)*uxp(ip))**2
            zmmnts0(18,js) = zmmnts0(18,js) + (gaminv(ip)*uyp(ip))**2
            zmmnts0(19,js) = zmmnts0(19,js) + (gaminv(ip)*uzp(ip))**2
            zmmnts0(22,js) = zmmnts0(22,js) + zp(ip)*gaminv(ip)*uzp(ip)
            zmmnts0(25,js) = zmmnts0(25,js) + xp(ip)*uzp(ip)*gaminv(ip)
            zmmnts0(26,js) = zmmnts0(26,js) + yp(ip)*uzp(ip)*gaminv(ip)
            zmmnts0(27,js) = zmmnts0(27,js) + uxp(ip)*uzp(ip)*gaminv(ip)**2
            zmmnts0(28,js) = zmmnts0(28,js) + uyp(ip)*uzp(ip)*gaminv(ip)**2
            zmmnts0( 5,js) = zmmnts0( 5,js) + uxp(ip)*vzi
            zmmnts0( 6,js) = zmmnts0( 6,js) + uyp(ip)*vzi
            zmmnts0(11,js) = zmmnts0(11,js) + uxp(ip)*uyp(ip)*vzi**2
            zmmnts0(15,js) = zmmnts0(15,js) + (uxp(ip)*vzi)**2
            zmmnts0(16,js) = zmmnts0(16,js) + (uyp(ip)*vzi)**2
            zmmnts0(20,js) = zmmnts0(20,js) + xp(ip)*uxp(ip)*vzi
            zmmnts0(21,js) = zmmnts0(21,js) + yp(ip)*uyp(ip)*vzi
            zmmnts0(23,js) = zmmnts0(23,js) + xp(ip)*uyp(ip)*vzi
            zmmnts0(24,js) = zmmnts0(24,js) + yp(ip)*uxp(ip)*vzi
          endif
        enddo

c-------------------------------------------------------------------------
#endif

c       --- Adds in dtscale factor efficiently. Somewhat kludgy though.
        if (dtscale /= 1.)
     &    zmmnts0(:,js) = zmmntscopy + (zmmnts0(:,js) - zmmntscopy)*dtscale


c       --- Now, accumulate moments onto Z grid.
        if (ifzmmnt == 2) then

c         --- Extrapolate quantities onto grid points so all particles
c         ---   contributing to a moment are at same z.
c         --- dti and dtip1 calculated by assuming constant Vz and 
c         ---   computing times at which z crosses the two nearest grid points.
c         --- V extrapolated using current and previous values of V
c         ---   (this avoids having to use E and B to extrapolate V).
c         --- X extrapolated assuming constant V (uses V at particle position).
c         --- After extrapolation, data is linearly weighted onto grid by
c         ---   particle's z location.  This has the effect that the less
c         ---   accurately extrapolated particles, which are far from the
c         ---   grid point, are weighted less, hopefully helping smoothness.

c         --- extrapolate to two nearest cells
          do ip=1,np
            if (pp(0,ip) /= 0.) then
              vzi = 1./(uzp(ip)*gaminv(ip)+SMALLPOS)
              iz(ip)  = (zp(ip) - zbeam - zmmntmin)*dzmi
              dti = (iz(ip)*dzm-zp(ip)+zbeam+zmmntmin)*vzi
              dti = max(-zmmntdtextmax*dt,min(zmmntdtextmax*dt,dti))
              pp(3,ip) = uxp(ip)*(1. + dti*oneondt) - uxpo(ip)*dti*oneondt
              pp(7,ip) = uyp(ip)*(1. + dti*oneondt) - uypo(ip)*dti*oneondt
              pp(11,ip) = uzp(ip)*(1. + dti*oneondt) - uzpo(ip)*dti*oneondt
              pp(1,ip) = xp(ip) + uxp(ip)*dti*gaminv(ip)
              pp(5,ip) = yp(ip) + uyp(ip)*dti*gaminv(ip)
              pp(9,ip) = zp(ip) + uzp(ip)*dti*gaminv(ip)
c             pp(9,ip) = iz(ip)*dzm + zbeam + zmmntmin

              dtip1=((iz(ip)+1)*dzm-zp(ip)+zbeam+zmmntmin)*vzi
              dtip1 = max(-zmmntdtextmax*dt,min(zmmntdtextmax*dt,dtip1))
              pp(4,ip) = uxp(ip)*(1. + dtip1*oneondt) - uxpo(ip)*dtip1*oneondt
              pp(8,ip) = uyp(ip)*(1. + dtip1*oneondt) - uypo(ip)*dtip1*oneondt
              pp(12,ip) = uzp(ip)*(1. + dtip1*oneondt) - uzpo(ip)*dtip1*oneondt
              pp(2,ip) = xp(ip) + uxp(ip)*dtip1*gaminv(ip)
              pp(6,ip) = yp(ip) + uyp(ip)*dtip1*gaminv(ip)
              pp(10,ip) = zp(ip) + uzp(ip)*dtip1*gaminv(ip)
c             pp(10,ip) = (iz(ip)+1)*dzm + zbeam + zmmntmin

              pp(13,ip) = gaminv(ip)
              pp(14,ip) = (1. - (zp(ip)-zbeam-zmmntmin)*dzmi + iz(ip))*pp(0,ip)
              pp(15,ip) = ((zp(ip) - zbeam - zmmntmin)*dzmi - iz(ip))*pp(0,ip)
              pp(16,ip) = dti
              pp(17,ip) = dtip1
            else
              iz(ip) = 0
              pp(14,ip) = 0.
              pp(15,ip) = 0.
            endif
          enddo

c         --- Add particles to ExtPart arrays (skips if nepwin = 0).
c         --- This allows access to the particle data when it is extrapolated
c         --- to the moments grid cell centers.
c         do iw=1,nepwin
c           if (izepwin(iw) < 0 .or. izepwin(iw) > nzmmnt) cycle

c            --- Make sure that there is enough room - assuming that all
c            --- particles could be saved.
c            if (nep(iw,is)+np > nepmax) then
c              nepmax = nep(iw,is)+np
c              call gchange("ExtPart",0)
c            endif
c            if (npidepmax .ne. npid) then
c              call gchange("ExtPart",0)
c            endif

c           do ip=1,np
c             if (iz(ip) == izepwin(iw) .and. nep(iw,is) < nepmax .and.
c    &            pp(0,ip) /= 0.) then
c$OMP ATOMIC
c               nep(iw,is) = nep(iw,is) + 1
c               tep(nep(iw,is),iw,is) = pp(16,ip) + time
c               xep(nep(iw,is),iw,is) = pp(1,ip)
c               yep(nep(iw,is),iw,is) = pp(5,ip)
c               uxep(nep(iw,is),iw,is) = pp(3,ip)
c               uyep(nep(iw,is),iw,is) = pp(7,ip)
c               uzep(nep(iw,is),iw,is) = pp(11,ip)
c               if (npid > 0) pidep(nep(iw,is),:,iw,is) = pid(ip,:)
c             else if (iz(ip)+1 == izepwin(iw) .and.
c    &                 nep(iw,is) < nepmax .and.
c    &                 pp(0,ip) /= 0.) then
c$OMP ATOMIC
c               nep(iw,is) = nep(iw,is) + 1
c               tep(nep(iw,is),iw,is) = pp(17,ip) + time
c               xep(nep(iw,is),iw,is) = pp(2,ip)
c               yep(nep(iw,is),iw,is) = pp(6,ip)
c               uxep(nep(iw,is),iw,is) = pp(4,ip)
c               uyep(nep(iw,is),iw,is) = pp(8,ip)
c               uzep(nep(iw,is),iw,is) = pp(12,ip)
c               if (npid > 0) pidep(nep(iw,is),:,iw,is) = pid(ip,:)
c             endif
c           enddo
c         enddo

c-------------------------------------------------------------------------
#ifdef J90
c       --- This method is faster for vectorized machines, primarily the J90.
c         --- put particle moments into pmmnts array for zmmnts(iz,)
          do ip=1,np
            vzi = 1./(pp(11,ip)+SMALLPOS)
            wz0 = pp(14,ip)
            pmmnts( 1,ip) = wz0
            pmmnts( 2,ip) = pp(1,ip)*wz0
            pmmnts( 3,ip) = pp(5,ip)*wz0
            pmmnts( 4,ip) = pp(9,ip)*wz0
            pmmnts( 7,ip) = pp(3,ip)*pp(13,ip)*wz0
            pmmnts( 8,ip) = pp(7,ip)*pp(13,ip)*wz0
            pmmnts( 9,ip) = pp(11,ip)*pp(13,ip)*wz0
            pmmnts(10,ip) = pp(1,ip)*pp(5,ip)*wz0
            pmmnts(12,ip) = pp(1,ip)**2*wz0
            pmmnts(13,ip) = pp(5,ip)**2*wz0
            pmmnts(14,ip) = pp(9,ip)**2*wz0
            pmmnts(17,ip) = (pp(13,ip)*pp(3,ip))**2*wz0
            pmmnts(18,ip) = (pp(13,ip)*pp(7,ip))**2*wz0
            pmmnts(19,ip) = (pp(13,ip)*pp(11,ip))**2*wz0
            pmmnts(22,ip) = pp(9,ip)*pp(13,ip)*pp(11,ip)*wz0
            pmmnts(25,ip) = pp(1,ip)*pp(11,ip)*pp(13,ip)*wz0
            pmmnts(26,ip) = pp(5,ip)*pp(11,ip)*pp(13,ip)*wz0
            pmmnts(27,ip) = pp(3,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            pmmnts(28,ip) = pp(7,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            pmmnts( 5,ip) = pp(3,ip)*vzi*wz0
            pmmnts( 6,ip) = pp(7,ip)*vzi*wz0
            pmmnts(11,ip) = pp(3,ip)*pp(7,ip)*vzi**2*wz0
            pmmnts(15,ip) = (pp(3,ip)*vzi)**2*wz0
            pmmnts(16,ip) = (pp(7,ip)*vzi)**2*wz0
            pmmnts(20,ip) = pp(1,ip)*pp(3,ip)*vzi*wz0
            pmmnts(21,ip) = pp(5,ip)*pp(7,ip)*vzi*wz0
            pmmnts(23,ip) = pp(1,ip)*pp(7,ip)*vzi*wz0
            pmmnts(24,ip) = pp(5,ip)*pp(3,ip)*vzi*wz0
          enddo

c         --- deposit moments into array
          do ip=1,np
            do im=1,NUMZMMNT
              zmmnts(iz(ip),im,js) = zmmnts(iz(ip),im,js) + pmmnts(im,ip)
            enddo
          enddo

c         --- put particle moments into pmmnts array for zmmnts(iz+1,)
          do ip=1,np
            vzi = 1./(pp(12,ip)+SMALLPOS)
            wz1 = pp(15,ip)
            pmmnts( 1,ip) = wz1
            pmmnts( 2,ip) = pp(2,ip)*wz1
            pmmnts( 3,ip) = pp(6,ip)*wz1
            pmmnts( 4,ip) = pp(10,ip)*wz1
            pmmnts( 7,ip) = pp(4,ip)*pp(13,ip)*wz1
            pmmnts( 8,ip) = pp(8,ip)*pp(13,ip)*wz1
            pmmnts( 9,ip) = pp(12,ip)*pp(13,ip)*wz1
            pmmnts(10,ip) = pp(2,ip)*pp(6,ip)*wz1
            pmmnts(12,ip) = pp(2,ip)**2*wz1
            pmmnts(13,ip) = pp(6,ip)**2*wz1
            pmmnts(14,ip) = pp(10,ip)**2*wz1
            pmmnts(17,ip) = (pp(13,ip)*pp(4,ip))**2*wz1
            pmmnts(18,ip) = (pp(13,ip)*pp(8,ip))**2*wz1
            pmmnts(19,ip) = (pp(13,ip)*pp(12,ip))**2*wz1
            pmmnts(22,ip) = pp(10,ip)*pp(13,ip)*pp(12,ip)*wz1
            pmmnts(25,ip) = pp(2,ip)*pp(12,ip)*pp(13,ip)*wz1
            pmmnts(26,ip) = pp(6,ip)*pp(12,ip)*pp(13,ip)*wz1
            pmmnts(27,ip) = pp(4,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            pmmnts(28,ip) = pp(8,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            pmmnts( 5,ip) = pp(4,ip)*vzi*wz1
            pmmnts( 6,ip) = pp(8,ip)*vzi*wz1
            pmmnts(11,ip) = pp(4,ip)*pp(8,ip)*vzi**2*wz1
            pmmnts(15,ip) = (pp(4,ip)*vzi)**2*wz1
            pmmnts(16,ip) = (pp(8,ip)*vzi)**2*wz1
            pmmnts(20,ip) = pp(2,ip)*pp(4,ip)*vzi*wz1
            pmmnts(21,ip) = pp(6,ip)*pp(8,ip)*vzi*wz1
            pmmnts(23,ip) = pp(2,ip)*pp(8,ip)*vzi*wz1
            pmmnts(24,ip) = pp(6,ip)*pp(4,ip)*vzi*wz1
          enddo

c         --- deposit moments into array
          do ip=1,np
            do im=1,NUMZMMNT
              zmmnts(iz(ip)+1,im,js) = zmmnts(iz(ip)+1,im,js) + pmmnts(im,ip)
            enddo
          enddo

#else

c-------------------------------------------------------------------------
c       --- This method is faster for RISC and pentium machines.

          do ip=1,np
c           --- sum particle moments into zmmnts(iz,)
            vzi = 1./(pp(11,ip)+SMALLPOS)
            wz0 = pp(14,ip)
            i = iz(ip)
            zmmnts(i, 1,js) = zmmnts(i, 1,js) + wz0
            zmmnts(i, 2,js) = zmmnts(i, 2,js) + pp(1,ip)*wz0
            zmmnts(i, 3,js) = zmmnts(i, 3,js) + pp(5,ip)*wz0
            zmmnts(i, 4,js) = zmmnts(i, 4,js) + pp(9,ip)*wz0
            zmmnts(i, 7,js) = zmmnts(i, 7,js) + pp(3,ip)*pp(13,ip)*wz0
            zmmnts(i, 8,js) = zmmnts(i, 8,js) + pp(7,ip)*pp(13,ip)*wz0
            zmmnts(i, 9,js) = zmmnts(i, 9,js) + pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,10,js) = zmmnts(i,10,js) + pp(1,ip)*pp(5,ip)*wz0
            zmmnts(i,12,js) = zmmnts(i,12,js) + pp(1,ip)**2*wz0
            zmmnts(i,13,js) = zmmnts(i,13,js) + pp(5,ip)**2*wz0
            zmmnts(i,14,js) = zmmnts(i,14,js) + pp(9,ip)**2*wz0
            zmmnts(i,17,js) = zmmnts(i,17,js) + (pp(13,ip)*pp(3,ip))**2*wz0
            zmmnts(i,18,js) = zmmnts(i,18,js) + (pp(13,ip)*pp(7,ip))**2*wz0
            zmmnts(i,19,js) = zmmnts(i,19,js) + (pp(13,ip)*pp(11,ip))**2*wz0
            zmmnts(i,22,js) = zmmnts(i,22,js) + pp(9,ip)*pp(13,ip)*pp(11,ip)*wz0
            zmmnts(i,25,js) = zmmnts(i,25,js) + pp(1,ip)*pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,26,js) = zmmnts(i,26,js) + pp(5,ip)*pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,27,js) = zmmnts(i,27,js) + pp(3,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            zmmnts(i,28,js) = zmmnts(i,28,js) + pp(7,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            zmmnts(i, 5,js) = zmmnts(i, 5,js) + pp(3,ip)*vzi*wz0
            zmmnts(i, 6,js) = zmmnts(i, 6,js) + pp(7,ip)*vzi*wz0
            zmmnts(i,11,js) = zmmnts(i,11,js) + pp(3,ip)*pp(7,ip)*vzi**2*wz0
            zmmnts(i,15,js) = zmmnts(i,15,js) + (pp(3,ip)*vzi)**2*wz0
            zmmnts(i,16,js) = zmmnts(i,16,js) + (pp(7,ip)*vzi)**2*wz0
            zmmnts(i,20,js) = zmmnts(i,20,js) + pp(1,ip)*pp(3,ip)*vzi*wz0
            zmmnts(i,21,js) = zmmnts(i,21,js) + pp(5,ip)*pp(7,ip)*vzi*wz0
            zmmnts(i,23,js) = zmmnts(i,23,js) + pp(1,ip)*pp(7,ip)*vzi*wz0
            zmmnts(i,24,js) = zmmnts(i,24,js) + pp(5,ip)*pp(3,ip)*vzi*wz0

c           --- sum particle moments into zmmnts(iz+1,,js)
            vzi = 1./(pp(12,ip)+SMALLPOS)
            wz1 = pp(15,ip)
            i = iz(ip) + 1
            zmmnts(i, 1,js) = zmmnts(i, 1,js) + wz1
            zmmnts(i, 2,js) = zmmnts(i, 2,js) + pp(2,ip)*wz1
            zmmnts(i, 3,js) = zmmnts(i, 3,js) + pp(6,ip)*wz1
            zmmnts(i, 4,js) = zmmnts(i, 4,js) + pp(10,ip)*wz1
            zmmnts(i, 7,js) = zmmnts(i, 7,js) + pp(4,ip)*pp(13,ip)*wz1
            zmmnts(i, 8,js) = zmmnts(i, 8,js) + pp(8,ip)*pp(13,ip)*wz1
            zmmnts(i, 9,js) = zmmnts(i, 9,js) + pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,10,js) = zmmnts(i,10,js) + pp(2,ip)*pp(6,ip)*wz1
            zmmnts(i,12,js) = zmmnts(i,12,js) + pp(2,ip)**2*wz1
            zmmnts(i,13,js) = zmmnts(i,13,js) + pp(6,ip)**2*wz1
            zmmnts(i,14,js) = zmmnts(i,14,js) + pp(10,ip)**2*wz1
            zmmnts(i,17,js) = zmmnts(i,17,js) + (pp(13,ip)*pp(4,ip))**2*wz1
            zmmnts(i,18,js) = zmmnts(i,18,js) + (pp(13,ip)*pp(8,ip))**2*wz1
            zmmnts(i,19,js) = zmmnts(i,19,js) + (pp(13,ip)*pp(12,ip))**2*wz1
            zmmnts(i,22,js) = zmmnts(i,22,js) + pp(10,ip)*pp(13,ip)*pp(12,ip)*wz1
            zmmnts(i,25,js) = zmmnts(i,25,js) + pp(2,ip)*pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,26,js) = zmmnts(i,26,js) + pp(6,ip)*pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,27,js) = zmmnts(i,27,js) + pp(4,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            zmmnts(i,28,js) = zmmnts(i,28,js) + pp(8,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            zmmnts(i, 5,js) = zmmnts(i, 5,js) + pp(4,ip)*vzi*wz1
            zmmnts(i, 6,js) = zmmnts(i, 6,js) + pp(8,ip)*vzi*wz1
            zmmnts(i,11,js) = zmmnts(i,11,js) + pp(4,ip)*pp(8,ip)*vzi**2*wz1
            zmmnts(i,15,js) = zmmnts(i,15,js) + (pp(4,ip)*vzi)**2*wz1
            zmmnts(i,16,js) = zmmnts(i,16,js) + (pp(8,ip)*vzi)**2*wz1
            zmmnts(i,20,js) = zmmnts(i,20,js) + pp(2,ip)*pp(4,ip)*vzi*wz1
            zmmnts(i,21,js) = zmmnts(i,21,js) + pp(6,ip)*pp(8,ip)*vzi*wz1
            zmmnts(i,23,js) = zmmnts(i,23,js) + pp(2,ip)*pp(8,ip)*vzi*wz1
            zmmnts(i,24,js) = zmmnts(i,24,js) + pp(6,ip)*pp(4,ip)*vzi*wz1

          enddo

c-------------------------------------------------------------------------
#endif

c         --- End of Z moments accumulation
        endif

        deallocate(pp,iz)

      endif

      if (itask == 3 .and. .not. laccumulate_zmoments) then

c       --- If the moments for each species was calculated separately, sum
c       --- them to get the combined moments.
        if (nszmmnt > 0) then
          maxp(:,nszmmnt) = maxval(maxp(:,0:nszmmnt-1),2)
          minp(:,nszmmnt) = minval(minp(:,0:nszmmnt-1),2)
          zmmnts0(:,nszmmnt) = sum(zmmnts0(:,0:nszmmnt-1),2)
          zmmnts(:,:,nszmmnt) = sum(zmmnts(:,:,0:nszmmnt-1),3)
        endif


        clighti = 1./clight
        gamma = 1./sqrt(1. - (vbeam*clighti)**2)

c       --- Make sure that the moments are set to have the same
c       --- number of species as the z moments.
        if (nsmmnt /= nszmmnt) then
          nsmmnt = nszmmnt
          call gchange("Moments",0)
        endif

c       --- Make sure that the window moments are set to have the same
c       --- number of species as the z moments.
        if (nswind /= nszmmnt) then
          nswind = nszmmnt
          call gchange("Win_Moments",0)
        endif

        xmaxp  = maxp(1,:)
        xminp  = minp(1,:)
        ymaxp  = maxp(2,:)
        yminp  = minp(2,:)
        zmaxp  = maxp(3,:)
        zminp  = minp(3,:)
        vxmaxp = maxp(4,:)
        vxminp = minp(4,:)
        vymaxp = maxp(5,:)
        vyminp = minp(5,:)
        vzmaxp = maxp(6,:)
        vzminp = minp(6,:)

c       --- For slave, call routine which sums moments over processors.
#ifdef MPIPARALLEL
        call parallel_sum_mmnts(zmmnts0,zmmnts)
#endif

c       --- Complete the calculation of moments: divide by particle number

c       --- Global moments first
        do js=0,nsmmnt
          pnumi = 1./(zmmnts0( 1,js)+SMALLPOS)
          pnum(0,js)    = zmmnts0( 1,js)
          xbar(0,js)    = zmmnts0( 2,js)*pnumi
          ybar(0,js)    = zmmnts0( 3,js)*pnumi
          zbar(0,js)    = zmmnts0( 4,js)*pnumi
          xpbar(0,js)   = zmmnts0( 5,js)*pnumi
          ypbar(0,js)   = zmmnts0( 6,js)*pnumi
          vxbar(0,js)   = zmmnts0( 7,js)*pnumi
          vybar(0,js)   = zmmnts0( 8,js)*pnumi
          vzbar(0,js)   = zmmnts0( 9,js)*pnumi
          xybar(0,js)   = zmmnts0(10,js)*pnumi
          xpypbar(0,js) = zmmnts0(11,js)*pnumi
          xsqbar(0,js)  = zmmnts0(12,js)*pnumi
          ysqbar(0,js)  = zmmnts0(13,js)*pnumi
          zsqbar(0,js)  = zmmnts0(14,js)*pnumi
          xpsqbar(0,js) = zmmnts0(15,js)*pnumi
          ypsqbar(0,js) = zmmnts0(16,js)*pnumi
          vxsqbar(0,js) = zmmnts0(17,js)*pnumi
          vysqbar(0,js) = zmmnts0(18,js)*pnumi
          vzsqbar(0,js) = zmmnts0(19,js)*pnumi
          xxpbar(0,js)  = zmmnts0(20,js)*pnumi
          yypbar(0,js)  = zmmnts0(21,js)*pnumi
          zvzbar(0,js)  = zmmnts0(22,js)*pnumi
          xypbar(0,js)  = zmmnts0(23,js)*pnumi 
          yxpbar(0,js)  = zmmnts0(24,js)*pnumi 
          xvzbar(0,js)  = zmmnts0(25,js)*pnumi
          yvzbar(0,js)  = zmmnts0(26,js)*pnumi
          vxvzbar(0,js) = zmmnts0(27,js)*pnumi
          vyvzbar(0,js) = zmmnts0(28,js)*pnumi

c         --- Compute second order moments with averages subtracted
          delxsq = xsqbar(0,js) - xbar(0,js)*xbar(0,js)
          delxxp = xxpbar(0,js) - xbar(0,js)*xpbar(0,js)
          delxpsq = xpsqbar(0,js) - xpbar(0,js)*xpbar(0,js)
          delysq = ysqbar(0,js) - ybar(0,js)*ybar(0,js)
          delyyp = yypbar(0,js) - ybar(0,js)*ypbar(0,js)
          delypsq = ypsqbar(0,js) - ypbar(0,js)*ypbar(0,js)
          delxy = xybar(0,js) - xbar(0,js)*ybar(0,js)
          delxyp = xypbar(0,js) - xbar(0,js)*ypbar(0,js)
          delyxp = yxpbar(0,js) - ybar(0,js)*xpbar(0,js)
          delxpyp = xpypbar(0,js) - xpbar(0,js)*ypbar(0,js)

          delvxsq = vxsqbar(0,js) - vxbar(0,js)**2
          delvysq = vysqbar(0,js) - vybar(0,js)**2

          delzsq = zsqbar(0,js) - zbar(0,js)**2
          delvzsq = vzsqbar(0,js) - vzbar(0,js)**2
          delzvz = zvzbar(0,js) - zbar(0,js)*vzbar(0,js)

c         --- Compute overall true RMS positions
          xrms(0,js) = sqrt(max(0.,delxsq))
          yrms(0,js) = sqrt(max(0.,delysq))
          zrms(0,js) = sqrt(max(0.,delzsq))
          rrms(0,js) = sqrt(xrms(0,js)**2 + yrms(0,js)**2)

c         --- Compute overall true RMS velocities
          xprms(0,js) = sqrt(max(0.,delxpsq))
          yprms(0,js) = sqrt(max(0.,delypsq))
          vxrms(0,js) = sqrt(max(0.,delvxsq))
          vyrms(0,js) = sqrt(max(0.,delvysq))
          vzrms(0,js) = sqrt(max(0.,delvzsq))

c         --- Compute overall emittances
          epsx(0,js) = 4.*sqrt(max(0.,delxsq*delxpsq - delxxp**2))
          epsy(0,js) = 4.*sqrt(max(0.,delysq*delypsq - delyyp**2))
          epsz(0,js) = 4.*sqrt(max(0.,delzsq*delvzsq - delzvz**2))
          epsnx(0,js) = epsx(0,js)*gamma*vzbar(0,js)*clighti*1.e6
          epsny(0,js) = epsy(0,js)*gamma*vzbar(0,js)*clighti*1.e6
          epsnz(0,js) = epsz(0,js)*gamma*vzbar(0,js)*clighti*1.e6

c         --- Compute generalized emittances
          tg = 0.5*(epsx(0,js)**2 + epsy(0,js)**2) +
     &             16*(delxy*delxpyp - delxyp*delyxp)
          epsg(0,js) = sqrt(max(0.,tg))
          th = epsx(0,js)**2*epsy(0,js)**2 +
     &       256*((delxy*delxpyp)**2 + (delxyp*delyxp)**2 -
     &       delxsq*delysq*(delxpyp)**2 - delxsq*delypsq*(delyxp)**2 -
     &       delxpsq*delysq*(delxyp)**2 - delxpsq*delypsq*(delxy)**2 -
     &       2*delxy*delxyp*delyxp*delxpyp + 2*delxxp*delypsq*delxy*delyxp -
     &       2*delxxp*delyyp*delxy*delxpyp - 2*delxxp*delyyp*delxyp*delyxp +
     &       2*delxpsq*delyyp*delxy*delxyp + 2*delxsq*delyyp*delyxp*delxpyp +
     &       2*delxxp*delysq*delxpyp*delxyp)
          epsh(0,js) = sqrt(sqrt(max(0.,th)))
          epsng(0,js) = epsg(0,js)*gamma*vzbar(0,js)*clighti*1.e6
          epsnh(0,js) = epsh(0,js)*gamma*vzbar(0,js)*clighti*1.e6

c         ---  Axial z momentum in vbeam frame (classical)
          pz(js) = zmmntsm(js)*zmmntsw(js)*pnum(0,js) * (vzbar(0,js)-vbeam)

c         ---  Total Z kinetic energy in lab frame minus beam energy
c         ---  (classical) using <vz^2>-vbeam^2;
          ekzmbe(js) = 0.5*zmmntsm(js)*zmmntsw(js)*pnum(0,js)*(vzsqbar(0,js) -
     &                 vbeam**2)

c         ---  Z kinetic energy in beam frame (classical), using <(vz-vbeam)^2>;
          ekzbeam(js) = 0.5*zmmntsm(js)*zmmntsw(js)*pnum(0,js)*
     &                  (vzsqbar(0,js) - 2.*vzbar(0,js)*vbeam + vbeam**2)

c         ---  Perp kinetic energy 
          ekperp(js) =  0.5*zmmntsm(js)*zmmntsw(js)*pnum(0,js)*(vxsqbar(0,js)+
     &                  vysqbar(0,js))

c         --- Total kinetic energy in beam frame
          ek(js) = ekzbeam(js) + ekperp(js)

        enddo

c       ---  Number of live particles (integer)
        if(wpid==0) nplive = pnum(0,nswind) + 0.5


c       --- Now scale the full Z moments by particle number
        if (ifzmmnt == 2) then
          do js=0,nszmmnt
            do izm = 0, nzmmnt
              pnumi = 1./(zmmnts(izm, 1,js)+SMALLPOS)
              pnumz(izm,js)    = zmmnts(izm, 1,js)
              xbarz(izm,js)    = zmmnts(izm, 2,js)*pnumi
              ybarz(izm,js)    = zmmnts(izm, 3,js)*pnumi
              zbarz(izm,js)    = zmmnts(izm, 4,js)*pnumi
              xpbarz(izm,js)   = zmmnts(izm, 5,js)*pnumi
              ypbarz(izm,js)   = zmmnts(izm, 6,js)*pnumi
              vxbarz(izm,js)   = zmmnts(izm, 7,js)*pnumi
              vybarz(izm,js)   = zmmnts(izm, 8,js)*pnumi
              vzbarz(izm,js)   = zmmnts(izm, 9,js)*pnumi
              xybarz(izm,js)   = zmmnts(izm,10,js)*pnumi
              xpypbarz(izm,js) = zmmnts(izm,11,js)*pnumi
              xsqbarz(izm,js)  = zmmnts(izm,12,js)*pnumi
              ysqbarz(izm,js)  = zmmnts(izm,13,js)*pnumi
              zsqbarz(izm,js)  = zmmnts(izm,14,js)*pnumi
              xpsqbarz(izm,js) = zmmnts(izm,15,js)*pnumi
              ypsqbarz(izm,js) = zmmnts(izm,16,js)*pnumi
              vxsqbarz(izm,js) = zmmnts(izm,17,js)*pnumi
              vysqbarz(izm,js) = zmmnts(izm,18,js)*pnumi
              vzsqbarz(izm,js) = zmmnts(izm,19,js)*pnumi
              xxpbarz(izm,js)  = zmmnts(izm,20,js)*pnumi
              yypbarz(izm,js)  = zmmnts(izm,21,js)*pnumi
              zvzbarz(izm,js)  = zmmnts(izm,22,js)*pnumi
              xypbarz(izm,js)  = zmmnts(izm,23,js)*pnumi 
              yxpbarz(izm,js)  = zmmnts(izm,24,js)*pnumi 
              xvzbarz(izm,js)  = zmmnts(izm,25,js)*pnumi
              yvzbarz(izm,js)  = zmmnts(izm,26,js)*pnumi
              vxvzbarz(izm,js) = zmmnts(izm,27,js)*pnumi
              vyvzbarz(izm,js) = zmmnts(izm,28,js)*pnumi

c             --- Compute second order moments with averages subtracted
              delxsq = xsqbarz(izm,js) - xbarz(izm,js)*xbarz(izm,js)
              delxxp = xxpbarz(izm,js) - xbarz(izm,js)*xpbarz(izm,js)
              delxpsq = xpsqbarz(izm,js) - xpbarz(izm,js)*xpbarz(izm,js)
              delysq = ysqbarz(izm,js) - ybarz(izm,js)*ybarz(izm,js)
              delyyp = yypbarz(izm,js) - ybarz(izm,js)*ypbarz(izm,js)
              delypsq = ypsqbarz(izm,js) - ypbarz(izm,js)*ypbarz(izm,js)
              delxy = xybarz(izm,js) - xbarz(izm,js)*ybarz(izm,js)
              delxyp = xypbarz(izm,js) - xbarz(izm,js)*ypbarz(izm,js)
              delyxp = yxpbarz(izm,js) - ybarz(izm,js)*xpbarz(izm,js)
              delxpyp = xpypbarz(izm,js) - xpbarz(izm,js)*ypbarz(izm,js)

              delvxsq = vxsqbarz(izm,js) - vxbarz(izm,js)**2
              delvysq = vysqbarz(izm,js) - vybarz(izm,js)**2

              delzsq = zsqbarz(izm,js) - zbarz(izm,js)**2
              delvzsq = vzsqbarz(izm,js) - vzbarz(izm,js)**2
              delzvz = zvzbarz(izm,js) - zbarz(izm,js)*vzbarz(izm,js)

c             --- Compute overall true RMS positions
              xrmsz(izm,js) = sqrt(max(0.,delxsq))
              yrmsz(izm,js) = sqrt(max(0.,delysq))
              zrmsz(izm,js) = sqrt(max(0.,delzsq))
              rrmsz(izm,js) = sqrt(xrmsz(izm,js)**2 + yrmsz(izm,js)**2)

c             --- Compute overall true RMS velocities
              xprmsz(izm,js) = sqrt(max(0.,delxpsq))
              yprmsz(izm,js) = sqrt(max(0.,delypsq))
              vxrmsz(izm,js) = sqrt(max(0.,delvxsq))
              vyrmsz(izm,js) = sqrt(max(0.,delvysq))
              vzrmsz(izm,js) = sqrt(max(0.,delvzsq))

c             --- Compute overall emittances
              epsxz(izm,js) = 4.*sqrt(max(0.,delxsq*delxpsq - delxxp**2))
              epsyz(izm,js) = 4.*sqrt(max(0.,delysq*delypsq - delyyp**2))
              epszz(izm,js) = 4.*sqrt(max(0.,delzsq*delvzsq - delzvz**2))
              epsnxz(izm,js) = epsxz(izm,js)*gamma*vzbarz(izm,js)*clighti*1.e6
              epsnyz(izm,js) = epsyz(izm,js)*gamma*vzbarz(izm,js)*clighti*1.e6
              epsnzz(izm,js) = epszz(izm,js)*gamma*vzbarz(izm,js)*clighti*1.e6

c             --- Compute generalized emittances
              tg = 0.5*(epsxz(izm,js)**2 + epsyz(izm,js)**2) +
     &           16*(delxy*delxpyp - delxyp*delyxp)
              epsgz(izm,js) = sqrt(max(0.,tg))
              th = epsxz(izm,js)**2*epsyz(izm,js)**2 +
     &          256*((delxy*delxpyp)**2 + (delxyp*delyxp)**2 -
     &          delxsq*delysq*(delxpyp)**2 - delxsq*delypsq*(delyxp)**2 -
     &          delxpsq*delysq*(delxyp)**2 - delxpsq*delypsq*(delxy)**2 -
     &          2*delxy*delxyp*delyxp*delxpyp + 2*delxxp*delypsq*delxy*delyxp -
     &          2*delxxp*delyyp*delxy*delxpyp - 2*delxxp*delyyp*delxyp*delyxp +
     &          2*delxpsq*delyyp*delxy*delxyp + 2*delxsq*delyyp*delyxp*delxpyp +
     &          2*delxxp*delysq*delxpyp*delxyp)
              epshz(izm,js) = sqrt(sqrt(max(0.,th)))
              epsngz(izm,js) = epsgz(izm,js)*gamma*vzbarz(izm,js)*clighti*1.e6
              epsnhz(izm,js) = epshz(izm,js)*gamma*vzbarz(izm,js)*clighti*1.e6

            enddo

c           --- Save window moments.  Linearly interpolate the Z moments data
c           --- to the Z window center.
c           --- If the zwindow is outside the zmoments range, don't set
c           --- anything.
            do iw = 1,nzwind
              zwin = 0.5*(zwindows(1,iw) + zwindows(2,iw))
              if (zmmntmin <= zwin .and. zwin < zmmntmax) then
                icell = int((zwin - zmmntmin)*dzmi)
                wz0 = (zwin - zmmntmin)*dzmi - icell
                wz1 = 1. - wz0
                pnum(iw,js)   =pnumz(icell,js)*wz1    + pnumz(icell+1,js)*wz0
                xbar(iw,js)   =xbarz(icell,js)*wz1    + xbarz(icell+1,js)*wz0
                ybar(iw,js)   =ybarz(icell,js)*wz1    + ybarz(icell+1,js)*wz0
                zbar(iw,js)   =zbarz(icell,js)*wz1    + zbarz(icell+1,js)*wz0
                xpbar(iw,js)  =xpbarz(icell,js)*wz1   + xpbarz(icell+1,js)*wz0
                ypbar(iw,js)  =ypbarz(icell,js)*wz1   + ypbarz(icell+1,js)*wz0
                vxbar(iw,js)  =vxbarz(icell,js)*wz1   + vxbarz(icell+1,js)*wz0
                vybar(iw,js)  =vybarz(icell,js)*wz1   + vybarz(icell+1,js)*wz0
                vzbar(iw,js)  =vzbarz(icell,js)*wz1   + vzbarz(icell+1,js)*wz0
                xybar(iw,js)  =xybarz(icell,js)*wz1   + xybarz(icell+1,js)*wz0
                xypbar(iw,js) =xypbarz(icell,js)*wz1  + xypbarz(icell+1,js)*wz0
                yxpbar(iw,js) =yxpbarz(icell,js)*wz1  + yxpbarz(icell+1,js)*wz0
                xpypbar(iw,js)=xpypbarz(icell,js)*wz1 + xpypbarz(icell+1,js)*wz0
                xsqbar(iw,js) =xsqbarz(icell,js)*wz1  + xsqbarz(icell+1,js)*wz0
                ysqbar(iw,js) =ysqbarz(icell,js)*wz1  + ysqbarz(icell+1,js)*wz0
                zsqbar(iw,js) =zsqbarz(icell,js)*wz1  + zsqbarz(icell+1,js)*wz0
                xpsqbar(iw,js)=xpsqbarz(icell,js)*wz1 + xpsqbarz(icell+1,js)*wz0
                ypsqbar(iw,js)=ypsqbarz(icell,js)*wz1 + ypsqbarz(icell+1,js)*wz0
                vxsqbar(iw,js)=vxsqbarz(icell,js)*wz1 + vxsqbarz(icell+1,js)*wz0
                vysqbar(iw,js)=vysqbarz(icell,js)*wz1 + vysqbarz(icell+1,js)*wz0
                vzsqbar(iw,js)=vzsqbarz(icell,js)*wz1 + vzsqbarz(icell+1,js)*wz0
                xxpbar(iw,js) =xxpbarz(icell,js)*wz1  + xxpbarz(icell+1,js)*wz0
                yypbar(iw,js) =yypbarz(icell,js)*wz1  + yypbarz(icell+1,js)*wz0
                xvzbar(iw,js) =xvzbarz(icell,js)*wz1  + xvzbarz(icell+1,js)*wz0
                yvzbar(iw,js) =yvzbarz(icell,js)*wz1  + yvzbarz(icell+1,js)*wz0
                zvzbar(iw,js) =zvzbarz(icell,js)*wz1  + zvzbarz(icell+1,js)*wz0
                vxvzbar(iw,js)=vxvzbarz(icell,js)*wz1 + vxvzbarz(icell+1,js)*wz0
                vyvzbar(iw,js)=vyvzbarz(icell,js)*wz1 + vyvzbarz(icell+1,js)*wz0
                xrms(iw,js)   =xrmsz(icell,js)*wz1    + xrmsz(icell+1,js)*wz0
                yrms(iw,js)   =yrmsz(icell,js)*wz1    + yrmsz(icell+1,js)*wz0
                zrms(iw,js)   =zrmsz(icell,js)*wz1    + zrmsz(icell+1,js)*wz0
                rrms(iw,js)   =rrmsz(icell,js)*wz1    + rrmsz(icell+1,js)*wz0
                xprms(iw,js)  =xprmsz(icell,js)*wz1   + xprmsz(icell+1,js)*wz0
                yprms(iw,js)  =yprmsz(icell,js)*wz1   + yprmsz(icell+1,js)*wz0
                epsx(iw,js)   =epsxz(icell,js)*wz1    + epsxz(icell+1,js)*wz0
                epsy(iw,js)   =epsyz(icell,js)*wz1    + epsyz(icell+1,js)*wz0
                epsz(iw,js)   =epszz(icell,js)*wz1    + epszz(icell+1,js)*wz0
                epsnx(iw,js)  =epsnxz(icell,js)*wz1   + epsnxz(icell+1,js)*wz0
                epsny(iw,js)  =epsnyz(icell,js)*wz1   + epsnyz(icell+1,js)*wz0
                epsnz(iw,js)  =epsnzz(icell,js)*wz1   + epsnzz(icell+1,js)*wz0
                epsg(iw,js)   =epsgz(icell,js)*wz1    + epsgz(icell+1,js)*wz0
                epsh(iw,js)   =epshz(icell,js)*wz1    + epshz(icell+1,js)*wz0
                epsng(iw,js)  =epsngz(icell,js)*wz1   + epsngz(icell+1,js)*wz0
                epsnh(iw,js)  =epsnhz(icell,js)*wz1   + epsnhz(icell+1,js)*wz0
                vxrms(iw,js)  =vxrmsz(icell,js)*wz1   + vxrmsz(icell+1,js)*wz0
                vyrms(iw,js)  =vyrmsz(icell,js)*wz1   + vyrmsz(icell+1,js)*wz0
                vzrms(iw,js)  =vzrmsz(icell,js)*wz1   + vzrmsz(icell+1,js)*wz0
              endif
            enddo
          enddo

        endif

      endif

!$OMP ATOMIC
      momentstime = momentstime + (wtime() - timetemp)

      return
      end
c=============================================================================
      subroutine getzmmnt_weights(np,xp,yp,zp,uxp,uyp,uzp,gaminv,wp,q,m,w,dt,
     &                            dtscale,itask,nplive,
     &                    uxpo,uypo,uzpo,is,isid,ismax,maxp,minp,zmmnts0,zmmnts)
      use Constant
      use Beam_acc
      use InDiag
      use Z_Moments
      use Win_Moments
      use Moments
      use Picglb
      use ExtPart
      use Timers
      integer(ISZ):: np,itask,nplive,is,isid,ismax
      real(kind=8):: q,m,w,dt,dtscale
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np), wp(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np)
      real(kind=8):: maxp(6,0:nszmmnt),minp(6,0:nszmmnt)
      real(kind=8):: zmmnts0(NUMZMMNT,0:nszmmnt)
      real(kind=8):: zmmnts(0:nzmmnt,NUMZMMNT,0:nszmmnt)

c  Sets moments for species 1 in as a function of z.
c  Interpolation is done to grid centers
c  Note that lost particles are to the left of the z grid.
c  Note: "window zero" includes all particles, and no extrapolation is done
c  Note: zbar and zsqbar calculations are meaningless except for 
c        window zero (WE SHOULD DELETE ?)
c  Note: Vectorized over moments.
c        particle moment calculations are vectorized over particles
c        moments are summed into a 2-D array, vectorized over moments
c        then they are dumped into the seperate 1-D arrays in itask=3
c  Three parts
c  When itask=1  zeros out all moments
c       itask=2  sums moments from particles
c       itask=3  divides by number of particles, calculates emittances and rms

c  Note that the moments are scaled by the dtscale factor. This accounts for
c  the fact that particles with larger dtscale would normally contribute more
c  if nominal dt would be used since they would spend proportionately more
c  time in any location.

      real(kind=8):: pp(1:17,np)
      integer(ISZ):: iz(np)
      integer(ISZ):: im,izm,iw,ip,icell,i,js
      real(kind=8):: dti,dtip1,wz1,wz0,zwin
      real(kind=8):: oneondt,clighti,vzi
      real(kind=8):: pnumi(0:nszmmnt)
      real(kind=8):: delxsq(0:nszmmnt),delxxp(0:nszmmnt),delxpsq(0:nszmmnt)
      real(kind=8):: delysq(0:nszmmnt),delyyp(0:nszmmnt),delypsq(0:nszmmnt)
      real(kind=8):: delxy(0:nszmmnt),delxyp(0:nszmmnt)
      real(kind=8):: delyxp(0:nszmmnt),delxpyp(0:nszmmnt)
      real(kind=8):: delzsq(0:nszmmnt),delvzsq(0:nszmmnt),delzvz(0:nszmmnt)
      real(kind=8):: delvxsq(0:nszmmnt),delvysq(0:nszmmnt)
      real(kind=8):: tg(0:nszmmnt),th(0:nszmmnt)
      real(kind=8):: gamma
      real(kind=8):: timetemp,wtime
      real(kind=8):: zmmntscopy(NUMZMMNT)
#ifdef J90
      real(kind=8):: pmmnts(NUMZMMNT,np)
#endif

      timetemp = wtime()

      if (ifzmmnt == 0) return

      oneondt = 1./dvnz(dt)

c     --- If lspeciesmoments is true, then calculate the moments for each
c     --- species separately. The moments for all species combined will be
c     --- done afterwards. Otherwise, only calculate the combined moments.
      if (lspeciesmoments) then
        js = isid - 1
      else
        js = 0
      endif

      if (itask == 1 .and. .not. laccumulate_zmoments) then
        call getzmmnt(np,xp,yp,zp,uxp,uyp,uzp,gaminv,q,m,w,dt,dtscale,
     &                itask,nplive,
     &                uxpo,uypo,uzpo,is,isid,ismax,maxp,minp,zmmnts0,zmmnts)
      endif

      if (itask == 2) then

c       --- Save the mass and weight
c       --- Note that the charge is not needed and should be removed from the
c       --- arugment list.
        zmmntsq(js) = q
        zmmntsm(js) = m
        zmmntsw(js) = w

c       --- Set flag that moments for this species has been calculated.
        zmomentscalculated(js) = .true.

c       --- Set maximum and minimum particle coordinate variables
        do ip = 1, np
          if (uzp(ip) /= 0.) then
            maxp(1,js) = max(maxp(1,js), xp(ip))
            minp(1,js) = min(minp(1,js), xp(ip))
            maxp(2,js) = max(maxp(2,js), yp(ip))
            minp(2,js) = min(minp(2,js), yp(ip))
            maxp(3,js) = max(maxp(3,js), zp(ip))
            minp(3,js) = min(minp(3,js), zp(ip))
            maxp(4,js) = max(maxp(4,js), gaminv(ip)*uxp(ip))
            minp(4,js) = min(minp(4,js), gaminv(ip)*uxp(ip))
            maxp(5,js) = max(maxp(5,js), gaminv(ip)*uyp(ip))
            minp(5,js) = min(minp(5,js), gaminv(ip)*uyp(ip))
            maxp(6,js) = max(maxp(6,js), gaminv(ip)*uzp(ip))
            minp(6,js) = min(minp(6,js), gaminv(ip)*uzp(ip))
            nplive = nplive + 1
          endif
        enddo

c       --- Accumate global moments first.

c       --- Used to add in dtscale factor below.
        if (dtscale /= 1.) zmmntscopy = zmmnts0(:,js)

c-------------------------------------------------------------------------
#ifdef J90
c       --- This method is faster for vectorized machines, primarily the J90.
c       --- put particle moments into pmmnts array for total moments
        do ip=1,np
          vzi = 1./(uzp(ip)+SMALLPOS)
          if (wp(ip) > 0.) then
            pmmnts( 1,ip) = wp(ip)*wp(ip)
            pmmnts( 2,ip) = xp(ip)*wp(ip)
            pmmnts( 3,ip) = yp(ip)*wp(ip)
            pmmnts( 4,ip) = zp(ip)*wp(ip)
            pmmnts( 7,ip) = uxp(ip)*gaminv(ip)*wp(ip)
            pmmnts( 8,ip) = uyp(ip)*gaminv(ip)*wp(ip)
            pmmnts( 9,ip) = uzp(ip)*gaminv(ip)*wp(ip)
            pmmnts(10,ip) = xp(ip)*yp(ip)*wp(ip)
            pmmnts(12,ip) = xp(ip)**2*wp(ip)
            pmmnts(13,ip) = yp(ip)**2*wp(ip)
            pmmnts(14,ip) = zp(ip)**2*wp(ip)
            pmmnts(17,ip) = (gaminv(ip)*uxp(ip))**2*wp(ip)
            pmmnts(18,ip) = (gaminv(ip)*uyp(ip))**2*wp(ip)
            pmmnts(19,ip) = (gaminv(ip)*uzp(ip))**2*wp(ip)
            pmmnts(22,ip) = zp(ip)*gaminv(ip)*uzp(ip)*wp(ip)
            pmmnts(25,ip) = xp(ip)*uzp(ip)*gaminv(ip)*wp(ip)
            pmmnts(26,ip) = yp(ip)*uzp(ip)*gaminv(ip)*wp(ip)
            pmmnts(27,ip) = uxp(ip)*uzp(ip)*gaminv(ip)**2*wp(ip)
            pmmnts(28,ip) = uyp(ip)*uzp(ip)*gaminv(ip)**2*wp(ip)
            pmmnts( 5,ip) = uxp(ip)*vzi*wp(ip)
            pmmnts( 6,ip) = uyp(ip)*vzi*wp(ip)
            pmmnts(11,ip) = uxp(ip)*uyp(ip)*vzi**2*wp(ip)
            pmmnts(15,ip) = (uxp(ip)*vzi)**2*wp(ip)
            pmmnts(16,ip) = (uyp(ip)*vzi)**2*wp(ip)
            pmmnts(20,ip) = xp(ip)*uxp(ip)*vzi*wp(ip)
            pmmnts(21,ip) = yp(ip)*uyp(ip)*vzi*wp(ip)
            pmmnts(23,ip) = xp(ip)*uyp(ip)*vzi*wp(ip)
            pmmnts(24,ip) = yp(ip)*uxp(ip)*vzi*wp(ip)
          endif
        enddo

c       --- sum total moments
c       --- Timings show that this loop runs much faster if the inner loop is
c       --- completely unrolled.  It is faster on both HP workstations and
c       --- CRAY C90.  WARNING: if NUMZMMNT changes, this loop needs to be
c       --- adjusted to reflect that change.  Note that the similar loops
c       --- below are still faster in their present form.
        do ip=1,np
          if (wp(ip) > 0.) then
c           do im=1,NUMZMMNT
c             zmmnts0(im,js) = zmmnts0(im,js) + pmmnts(ip,im)
c           enddo
            zmmnts0( 1,js) = zmmnts0( 1,js) + pmmnts( 1,ip)
            zmmnts0( 2,js) = zmmnts0( 2,js) + pmmnts( 2,ip)
            zmmnts0( 3,js) = zmmnts0( 3,js) + pmmnts( 3,ip)
            zmmnts0( 4,js) = zmmnts0( 4,js) + pmmnts( 4,ip)
            zmmnts0( 5,js) = zmmnts0( 5,js) + pmmnts( 5,ip)
            zmmnts0( 6,js) = zmmnts0( 6,js) + pmmnts( 6,ip)
            zmmnts0( 7,js) = zmmnts0( 7,js) + pmmnts( 7,ip)
            zmmnts0( 8,js) = zmmnts0( 8,js) + pmmnts( 8,ip)
            zmmnts0( 9,js) = zmmnts0( 9,js) + pmmnts( 9,ip)
            zmmnts0(10,js) = zmmnts0(10,js) + pmmnts(10,ip)
            zmmnts0(11,js) = zmmnts0(11,js) + pmmnts(11,ip)
            zmmnts0(12,js) = zmmnts0(12,js) + pmmnts(12,ip)
            zmmnts0(13,js) = zmmnts0(13,js) + pmmnts(13,ip)
            zmmnts0(14,js) = zmmnts0(14,js) + pmmnts(14,ip)
            zmmnts0(15,js) = zmmnts0(15,js) + pmmnts(15,ip)
            zmmnts0(16,js) = zmmnts0(16,js) + pmmnts(16,ip)
            zmmnts0(17,js) = zmmnts0(17,js) + pmmnts(17,ip)
            zmmnts0(18,js) = zmmnts0(18,js) + pmmnts(18,ip)
            zmmnts0(19,js) = zmmnts0(19,js) + pmmnts(19,ip)
            zmmnts0(20,js) = zmmnts0(20,js) + pmmnts(20,ip)
            zmmnts0(21,js) = zmmnts0(21,js) + pmmnts(21,ip)
            zmmnts0(22,js) = zmmnts0(22,js) + pmmnts(22,ip)
            zmmnts0(23,js) = zmmnts0(23,js) + pmmnts(23,ip) 
            zmmnts0(24,js) = zmmnts0(24,js) + pmmnts(24,ip) 
            zmmnts0(25,js) = zmmnts0(25,js) + pmmnts(25,ip) 
            zmmnts0(26,js) = zmmnts0(26,js) + pmmnts(26,ip) 
            zmmnts0(27,js) = zmmnts0(27,js) + pmmnts(27,ip) 
            zmmnts0(28,js) = zmmnts0(28,js) + pmmnts(28,ip) 
          endif
        enddo

#else

c-------------------------------------------------------------------------
c       --- This method is faster for RISC and pentium machines.
c       --- sum particle moments into total moments
        do ip=1,np
          vzi = 1./(uzp(ip)+SMALLPOS)
          if (wp(ip) > 0.) then
            zmmnts0( 1,js) = zmmnts0( 1,js) + wp(ip)*wp(ip)
            zmmnts0( 2,js) = zmmnts0( 2,js) + xp(ip)*wp(ip)
            zmmnts0( 3,js) = zmmnts0( 3,js) + yp(ip)*wp(ip)
            zmmnts0( 4,js) = zmmnts0( 4,js) + zp(ip)*wp(ip)
            zmmnts0( 7,js) = zmmnts0( 7,js) + uxp(ip)*gaminv(ip)*wp(ip)
            zmmnts0( 8,js) = zmmnts0( 8,js) + uyp(ip)*gaminv(ip)*wp(ip)
            zmmnts0( 9,js) = zmmnts0( 9,js) + uzp(ip)*gaminv(ip)*wp(ip)
            zmmnts0(10,js) = zmmnts0(10,js) + xp(ip)*yp(ip)*wp(ip)
            zmmnts0(12,js) = zmmnts0(12,js) + xp(ip)**2*wp(ip)
            zmmnts0(13,js) = zmmnts0(13,js) + yp(ip)**2*wp(ip)
            zmmnts0(14,js) = zmmnts0(14,js) + zp(ip)**2*wp(ip)
            zmmnts0(17,js) = zmmnts0(17,js) + (gaminv(ip)*uxp(ip))**2*wp(ip)
            zmmnts0(18,js) = zmmnts0(18,js) + (gaminv(ip)*uyp(ip))**2*wp(ip)
            zmmnts0(19,js) = zmmnts0(19,js) + (gaminv(ip)*uzp(ip))**2*wp(ip)
            zmmnts0(22,js) = zmmnts0(22,js) + zp(ip)*gaminv(ip)*uzp(ip)*wp(ip)
            zmmnts0(25,js) = zmmnts0(25,js) + xp(ip)*uzp(ip)*gaminv(ip)*wp(ip)
            zmmnts0(26,js) = zmmnts0(26,js) + yp(ip)*uzp(ip)*gaminv(ip)*wp(ip)
            zmmnts0(27,js) = zmmnts0(27,js) + uxp(ip)*uzp(ip)*gaminv(ip)**2*wp(ip)
            zmmnts0(28,js) = zmmnts0(28,js) + uyp(ip)*uzp(ip)*gaminv(ip)**2*wp(ip)
            zmmnts0( 5,js) = zmmnts0( 5,js) + uxp(ip)*vzi*wp(ip)
            zmmnts0( 6,js) = zmmnts0( 6,js) + uyp(ip)*vzi*wp(ip)
            zmmnts0(11,js) = zmmnts0(11,js) + uxp(ip)*uyp(ip)*vzi**2*wp(ip)
            zmmnts0(15,js) = zmmnts0(15,js) + (uxp(ip)*vzi)**2*wp(ip)
            zmmnts0(16,js) = zmmnts0(16,js) + (uyp(ip)*vzi)**2*wp(ip)
            zmmnts0(20,js) = zmmnts0(20,js) + xp(ip)*uxp(ip)*vzi*wp(ip)
            zmmnts0(21,js) = zmmnts0(21,js) + yp(ip)*uyp(ip)*vzi*wp(ip)
            zmmnts0(23,js) = zmmnts0(23,js) + xp(ip)*uyp(ip)*vzi*wp(ip)
            zmmnts0(24,js) = zmmnts0(24,js) + yp(ip)*uxp(ip)*vzi*wp(ip)
          endif
        enddo

c-------------------------------------------------------------------------
#endif

c       --- Adds in dtscale factor efficiently. Somewhat kludgy though.
        if (dtscale /= 1.)
     &    zmmnts0(:,js) = zmmntscopy + (zmmnts0(:,js) - zmmntscopy)*dtscale


c       --- Now, accumulate moments onto Z grid.
        if (ifzmmnt == 2) then

c         --- Extrapolate quantities onto grid points so all particles
c         ---   contributing to a moment are at same z.
c         --- dti and dtip1 calculated by assuming constant Vz and 
c         ---   computing times at which z crosses the two nearest grid points.
c         --- V extrapolated using current and previous values of V
c         ---   (this avoids having to use E and B to extrapolate V).
c         --- X extrapolated assuming constant V (uses V at particle position).
c         --- After extrapolation, data is linearly weighted onto grid by
c         ---   particle's z location.  This has the effect that the less
c         ---   accurately extrapolated particles, which are far from the
c         ---   grid point, are weighted less, hopefully helping smoothness.

c         --- extrapolate to two nearest cells
          do ip=1,np
            vzi = 1./(uzp(ip)*gaminv(ip)+SMALLPOS)
            iz(ip)  = (zp(ip) - zbeam - zmmntmin)*dzmi
            dti = (iz(ip)*dzm-zp(ip)+zbeam+zmmntmin)*vzi
            pp(3,ip) = uxp(ip)*(1. + dti*oneondt) - uxpo(ip)*dti*oneondt
            pp(7,ip) = uyp(ip)*(1. + dti*oneondt) - uypo(ip)*dti*oneondt
            pp(11,ip) = uzp(ip)*(1. + dti*oneondt) - uzpo(ip)*dti*oneondt
            pp(1,ip) = xp(ip) + uxp(ip)*dti*gaminv(ip)
            pp(5,ip) = yp(ip) + uyp(ip)*dti*gaminv(ip)
            pp(9,ip) = zp(ip) + uzp(ip)*dti*gaminv(ip)
c           pp(9,ip) = iz(ip)*dzm + zbeam + zmmntmin

            dtip1=((iz(ip)+1)*dzm-zp(ip)+zbeam+zmmntmin)*vzi
            pp(4,ip) = uxp(ip)*(1. + dtip1*oneondt) - uxpo(ip)*dtip1*oneondt
            pp(8,ip) = uyp(ip)*(1. + dtip1*oneondt) - uypo(ip)*dtip1*oneondt
            pp(12,ip) = uzp(ip)*(1. + dtip1*oneondt) - uzpo(ip)*dtip1*oneondt
            pp(2,ip) = xp(ip) + uxp(ip)*dtip1*gaminv(ip)
            pp(6,ip) = yp(ip) + uyp(ip)*dtip1*gaminv(ip)
            pp(10,ip) = zp(ip) + uzp(ip)*dtip1*gaminv(ip)
c           pp(10,ip) = (iz(ip)+1)*dzm + zbeam + zmmntmin

            pp(13,ip) = gaminv(ip)
            pp(14,ip) = (1. - (zp(ip)-zbeam-zmmntmin)*dzmi + iz(ip))*wp(ip)*
     &                  dtscale
            pp(15,ip) = ((zp(ip) - zbeam - zmmntmin)*dzmi - iz(ip))*wp(ip)*
     &                  dtscale
            pp(16,ip) = dti
            pp(17,ip) = dtip1
          enddo

c         --- Add particles to ExtPart arrays (skips if nepwin = 0).
c         --- This allows access to the particle data when it is extrapolated
c         --- to the moments grid cell centers.
c         do iw=1,nepwin
c           if (izepwin(iw) < 0 .or. izepwin(iw) > nzmmnt) cycle

c            --- Make sure that there is enough room - assuming that all
c            --- particles could be saved.
c            if (nep(iw,is)+np > nepmax) then
c              nepmax = nep(iw,is)+np
c              call gchange("ExtPart",0)
c            endif
c            if (npidepmax .ne. npid) then
c              call gchange("ExtPart",0)
c            endif

c           do ip=1,np
c             if (iz(ip) == izepwin(iw) .and. nep(iw,is) < nepmax .and.
c    &            wp(ip) > 0.) then
c$OMP ATOMIC
c               nep(iw,is) = nep(iw,is) + 1
c               tep(nep(iw,is),iw,is) = pp(16,ip) + time
c               xep(nep(iw,is),iw,is) = pp(1,ip)
c               yep(nep(iw,is),iw,is) = pp(5,ip)
c               uxep(nep(iw,is),iw,is) = pp(3,ip)
c               uyep(nep(iw,is),iw,is) = pp(7,ip)
c               uzep(nep(iw,is),iw,is) = pp(11,ip)
c               if (npid > 0) pidep(nep(iw,is),:,iw,is) = pid(ip,:)
c             else if (iz(ip)+1 == izepwin(iw) .and.
c    &                 nep(iw,is) < nepmax .and.
c    &                 wp(ip) > 0.) then
c$OMP ATOMIC
c               nep(iw,is) = nep(iw,is) + 1
c               tep(nep(iw,is),iw,is) = pp(17,ip) + time
c               xep(nep(iw,is),iw,is) = pp(2,ip)
c               yep(nep(iw,is),iw,is) = pp(6,ip)
c               uxep(nep(iw,is),iw,is) = pp(4,ip)
c               uyep(nep(iw,is),iw,is) = pp(8,ip)
c               uzep(nep(iw,is),iw,is) = pp(12,ip)
c               if (npid > 0) pidep(nep(iw,is),:,iw,is) = pid(ip,:)
c             endif
c           enddo
c         enddo

c-------------------------------------------------------------------------
#ifdef J90
c       --- This method is faster for vectorized machines, primarily the J90.
c         --- put particle moments into pmmnts array for zmmnts(iz,)
          do ip=1,np
            vzi = 1./(pp(11,ip)+SMALLPOS)
            wz0 = pp(14,ip)
            pmmnts( 1,ip) = wz0
            pmmnts( 2,ip) = pp(1,ip)*wz0
            pmmnts( 3,ip) = pp(5,ip)*wz0
            pmmnts( 4,ip) = pp(9,ip)*wz0
            pmmnts( 7,ip) = pp(3,ip)*pp(13,ip)*wz0
            pmmnts( 8,ip) = pp(7,ip)*pp(13,ip)*wz0
            pmmnts( 9,ip) = pp(11,ip)*pp(13,ip)*wz0
            pmmnts(10,ip) = pp(1,ip)*pp(5,ip)*wz0
            pmmnts(12,ip) = pp(1,ip)**2*wz0
            pmmnts(13,ip) = pp(5,ip)**2*wz0
            pmmnts(14,ip) = pp(9,ip)**2*wz0
            pmmnts(17,ip) = (pp(13,ip)*pp(3,ip))**2*wz0
            pmmnts(18,ip) = (pp(13,ip)*pp(7,ip))**2*wz0
            pmmnts(19,ip) = (pp(13,ip)*pp(11,ip))**2*wz0
            pmmnts(22,ip) = pp(9,ip)*pp(13,ip)*pp(11,ip)*wz0
            pmmnts(25,ip) = pp(1,ip)*pp(11,ip)*pp(13,ip)*wz0
            pmmnts(26,ip) = pp(5,ip)*pp(11,ip)*pp(13,ip)*wz0
            pmmnts(27,ip) = pp(3,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            pmmnts(28,ip) = pp(7,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            pmmnts( 5,ip) = pp(3,ip)*vzi*wz0
            pmmnts( 6,ip) = pp(7,ip)*vzi*wz0
            pmmnts(11,ip) = pp(3,ip)*pp(7,ip)*vzi**2*wz0
            pmmnts(15,ip) = (pp(3,ip)*vzi)**2*wz0
            pmmnts(16,ip) = (pp(7,ip)*vzi)**2*wz0
            pmmnts(20,ip) = pp(1,ip)*pp(3,ip)*vzi*wz0
            pmmnts(21,ip) = pp(5,ip)*pp(7,ip)*vzi*wz0
            pmmnts(23,ip) = pp(1,ip)*pp(7,ip)*vzi*wz0
            pmmnts(24,ip) = pp(5,ip)*pp(3,ip)*vzi*wz0
          enddo

c         --- deposit moments into array
          do ip=1,np
            do im=1,NUMZMMNT
              zmmnts(iz(ip),im,js) = zmmnts(iz(ip),im,js) + pmmnts(im,ip)
            enddo
          enddo

c         --- put particle moments into pmmnts array for zmmnts(iz+1,)
          do ip=1,np
            vzi = 1./(pp(12,ip)+SMALLPOS)
            wz1 = pp(15,ip)
            pmmnts( 1,ip) = wz1
            pmmnts( 2,ip) = pp(2,ip)*wz1
            pmmnts( 3,ip) = pp(6,ip)*wz1
            pmmnts( 4,ip) = pp(10,ip)*wz1
            pmmnts( 7,ip) = pp(4,ip)*pp(13,ip)*wz1
            pmmnts( 8,ip) = pp(8,ip)*pp(13,ip)*wz1
            pmmnts( 9,ip) = pp(12,ip)*pp(13,ip)*wz1
            pmmnts(10,ip) = pp(2,ip)*pp(6,ip)*wz1
            pmmnts(12,ip) = pp(2,ip)**2*wz1
            pmmnts(13,ip) = pp(6,ip)**2*wz1
            pmmnts(14,ip) = pp(10,ip)**2*wz1
            pmmnts(17,ip) = (pp(13,ip)*pp(4,ip))**2*wz1
            pmmnts(18,ip) = (pp(13,ip)*pp(8,ip))**2*wz1
            pmmnts(19,ip) = (pp(13,ip)*pp(12,ip))**2*wz1
            pmmnts(22,ip) = pp(10,ip)*pp(13,ip)*pp(12,ip)*wz1
            pmmnts(25,ip) = pp(2,ip)*pp(12,ip)*pp(13,ip)*wz1
            pmmnts(26,ip) = pp(6,ip)*pp(12,ip)*pp(13,ip)*wz1
            pmmnts(27,ip) = pp(4,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            pmmnts(28,ip) = pp(8,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            pmmnts( 5,ip) = pp(4,ip)*vzi*wz1
            pmmnts( 6,ip) = pp(8,ip)*vzi*wz1
            pmmnts(11,ip) = pp(4,ip)*pp(8,ip)*vzi**2*wz1
            pmmnts(15,ip) = (pp(4,ip)*vzi)**2*wz1
            pmmnts(16,ip) = (pp(8,ip)*vzi)**2*wz1
            pmmnts(20,ip) = pp(2,ip)*pp(4,ip)*vzi*wz1
            pmmnts(21,ip) = pp(6,ip)*pp(8,ip)*vzi*wz1
            pmmnts(23,ip) = pp(2,ip)*pp(8,ip)*vzi*wz1
            pmmnts(24,ip) = pp(6,ip)*pp(4,ip)*vzi*wz1
          enddo

c         --- deposit moments into array
          do ip=1,np
            do im=1,NUMZMMNT
              zmmnts(iz(ip)+1,im,js) = zmmnts(iz(ip)+1,im,js) + pmmnts(im,ip)
            enddo
          enddo

#else

c-------------------------------------------------------------------------
c       --- This method is faster for RISC and pentium machines.

          do ip=1,np
c           --- sum particle moments into zmmnts(iz,,js)
            vzi = 1./(pp(11,ip)+SMALLPOS)
            wz0 = pp(14,ip)
            i = iz(ip)
            zmmnts(i, 1,js) = zmmnts(i, 1,js) + wz0
            zmmnts(i, 2,js) = zmmnts(i, 2,js) + pp(1,ip)*wz0
            zmmnts(i, 3,js) = zmmnts(i, 3,js) + pp(5,ip)*wz0
            zmmnts(i, 4,js) = zmmnts(i, 4,js) + pp(9,ip)*wz0
            zmmnts(i, 7,js) = zmmnts(i, 7,js) + pp(3,ip)*pp(13,ip)*wz0
            zmmnts(i, 8,js) = zmmnts(i, 8,js) + pp(7,ip)*pp(13,ip)*wz0
            zmmnts(i, 9,js) = zmmnts(i, 9,js) + pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,10,js) = zmmnts(i,10,js) + pp(1,ip)*pp(5,ip)*wz0
            zmmnts(i,12,js) = zmmnts(i,12,js) + pp(1,ip)**2*wz0
            zmmnts(i,13,js) = zmmnts(i,13,js) + pp(5,ip)**2*wz0
            zmmnts(i,14,js) = zmmnts(i,14,js) + pp(9,ip)**2*wz0
            zmmnts(i,17,js) = zmmnts(i,17,js) + (pp(13,ip)*pp(3,ip))**2*wz0
            zmmnts(i,18,js) = zmmnts(i,18,js) + (pp(13,ip)*pp(7,ip))**2*wz0
            zmmnts(i,19,js) = zmmnts(i,19,js) + (pp(13,ip)*pp(11,ip))**2*wz0
            zmmnts(i,22,js) = zmmnts(i,22,js) + pp(9,ip)*pp(13,ip)*pp(11,ip)*wz0
            zmmnts(i,25,js) = zmmnts(i,25,js) + pp(1,ip)*pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,26,js) = zmmnts(i,26,js) + pp(5,ip)*pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,27,js) = zmmnts(i,27,js) + pp(3,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            zmmnts(i,28,js) = zmmnts(i,28,js) + pp(7,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            zmmnts(i, 5,js) = zmmnts(i, 5,js) + pp(3,ip)*vzi*wz0
            zmmnts(i, 6,js) = zmmnts(i, 6,js) + pp(7,ip)*vzi*wz0
            zmmnts(i,11,js) = zmmnts(i,11,js) + pp(3,ip)*pp(7,ip)*vzi**2*wz0
            zmmnts(i,15,js) = zmmnts(i,15,js) + (pp(3,ip)*vzi)**2*wz0
            zmmnts(i,16,js) = zmmnts(i,16,js) + (pp(7,ip)*vzi)**2*wz0
            zmmnts(i,20,js) = zmmnts(i,20,js) + pp(1,ip)*pp(3,ip)*vzi*wz0
            zmmnts(i,21,js) = zmmnts(i,21,js) + pp(5,ip)*pp(7,ip)*vzi*wz0
            zmmnts(i,23,js) = zmmnts(i,23,js) + pp(1,ip)*pp(7,ip)*vzi*wz0
            zmmnts(i,24,js) = zmmnts(i,24,js) + pp(5,ip)*pp(3,ip)*vzi*wz0

c           --- sum particle moments into zmmnts(iz+1,,js)
            vzi = 1./(pp(12,ip)+SMALLPOS)
            wz1 = pp(15,ip)
            i = iz(ip) + 1
            zmmnts(i, 1,js) = zmmnts(i, 1,js) + wz1
            zmmnts(i, 2,js) = zmmnts(i, 2,js) + pp(2,ip)*wz1
            zmmnts(i, 3,js) = zmmnts(i, 3,js) + pp(6,ip)*wz1
            zmmnts(i, 4,js) = zmmnts(i, 4,js) + pp(10,ip)*wz1
            zmmnts(i, 7,js) = zmmnts(i, 7,js) + pp(4,ip)*pp(13,ip)*wz1
            zmmnts(i, 8,js) = zmmnts(i, 8,js) + pp(8,ip)*pp(13,ip)*wz1
            zmmnts(i, 9,js) = zmmnts(i, 9,js) + pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,10,js) = zmmnts(i,10,js) + pp(2,ip)*pp(6,ip)*wz1
            zmmnts(i,12,js) = zmmnts(i,12,js) + pp(2,ip)**2*wz1
            zmmnts(i,13,js) = zmmnts(i,13,js) + pp(6,ip)**2*wz1
            zmmnts(i,14,js) = zmmnts(i,14,js) + pp(10,ip)**2*wz1
            zmmnts(i,17,js) = zmmnts(i,17,js) + (pp(13,ip)*pp(4,ip))**2*wz1
            zmmnts(i,18,js) = zmmnts(i,18,js) + (pp(13,ip)*pp(8,ip))**2*wz1
            zmmnts(i,19,js) = zmmnts(i,19,js) + (pp(13,ip)*pp(12,ip))**2*wz1
            zmmnts(i,22,js) = zmmnts(i,22,js) + pp(10,ip)*pp(13,ip)*pp(12,ip)*wz1
            zmmnts(i,25,js) = zmmnts(i,25,js) + pp(2,ip)*pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,26,js) = zmmnts(i,26,js) + pp(6,ip)*pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,27,js) = zmmnts(i,27,js) + pp(4,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            zmmnts(i,28,js) = zmmnts(i,28,js) + pp(8,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            zmmnts(i, 5,js) = zmmnts(i, 5,js) + pp(4,ip)*vzi*wz1
            zmmnts(i, 6,js) = zmmnts(i, 6,js) + pp(8,ip)*vzi*wz1
            zmmnts(i,11,js) = zmmnts(i,11,js) + pp(4,ip)*pp(8,ip)*vzi**2*wz1
            zmmnts(i,15,js) = zmmnts(i,15,js) + (pp(4,ip)*vzi)**2*wz1
            zmmnts(i,16,js) = zmmnts(i,16,js) + (pp(8,ip)*vzi)**2*wz1
            zmmnts(i,20,js) = zmmnts(i,20,js) + pp(2,ip)*pp(4,ip)*vzi*wz1
            zmmnts(i,21,js) = zmmnts(i,21,js) + pp(6,ip)*pp(8,ip)*vzi*wz1
            zmmnts(i,23,js) = zmmnts(i,23,js) + pp(2,ip)*pp(8,ip)*vzi*wz1
            zmmnts(i,24,js) = zmmnts(i,24,js) + pp(6,ip)*pp(4,ip)*vzi*wz1

          enddo

c-------------------------------------------------------------------------
#endif

c         --- End of Z moments accumulation
        endif

      endif

      if (itask == 3 .and. .not. laccumulate_zmoments) then
        call getzmmnt(np,xp,yp,zp,uxp,uyp,uzp,gaminv,q,m,w,dt,dtscale,
     &                itask,nplive,
     &                uxpo,uypo,uzpo,is,isid,ismax,maxp,minp,zmmnts0,zmmnts)
      endif

!$OMP ATOMIC
      momentstime = momentstime + (wtime() - timetemp)

      return
      end
c=============================================================================
      subroutine getextrapolatedparticles(ipmin,np,pgroup,npid,
     &                                    dt,uxpo,uypo,uzpo,isid,time,zbeam)
      use ParticleGroupmodule
      use Z_Moments
      use ExtPart
      use Timers
      integer(ISZ):: ipmin,np,npid,isid
      type(ParticleGroup):: pgroup
      real(kind=8):: dt,time,zbeam
      real(kind=8):: uxpo(np),uypo(np),uzpo(np)

c Given current and previous velocities, extrapolates particles to the
c positions given (in the lab frame) by zzepwin.
c The extrapolation is approximate. A better version might use the full
c advance capabilities (from padvnc3d) but that would probably be overkill.
c --- az calculated from current and previous Vz and dt
c --- dti calculated by assuming constant az and 
c ---   computing times at which z crosses lab window
c --- V extrapolated using current V and estimated A
c ---   (this avoids having to use E and B to extrapolate V).
c --- X extrapolated using constant V and estimated A

      real(kind=8):: zzep,wzep
      real(kind=8):: oneondt,dtstar,delz,vx,vy,vz,ax,ay,az,dp,dm
      real(kind=8):: vznp1,zznp1,delznp1,dtstarnp1
      integer(ISZ):: iw,ip,nn

      real(kind=8):: timetemp,wtime

c     --- Check if the pgroup had been initialized.
      if (isid == 0) return

      oneondt = 1./dt

      if (npidepmax .ne. npid) then
        npidepmax = npid
        call gchange("ExtPart",0)
      endif

c     --- Loop over windows
      do iw=1,nepwin

c       --- Get the location of the extrapolation plane.
        if (izepwin(iw) >= 0) then
          zzep = zmmntmin + izepwin(iw)*dzm + zbeam
          wzep = dzm
        else
          zzep = zzepwin(iw)
          wzep = wzepwin(iw)
        endif

c       --- Skip windows that are outside of the grid frame.
        if (zzep+wzep <= zmmntmin+zbeam .or.
     &      zzep-wzep >= zmmntmax+zbeam) cycle

c       --- Make sure that there is enough room - assuming that all particles
c       --- could be saved.
        if (nep(iw,isid)+np > nepmax) then
          nepmax = nep(iw,isid)+np
          call gchange("ExtPart",0)
        endif

c       --- extrapolate to window if near it
        do ip=ipmin,ipmin+np-1
          if (pgroup%zp(ip) < zzep-wzep .or.
     &        pgroup%zp(ip) > zzep+wzep) cycle

c         --- Check if particle has already been saved in this window.
c         --- This only works if multiple windows don't overlap.
          if (lepsaveonce) then
            if (pgroup%pid(ip,epflagpid) == iw) cycle
          endif

          vx = pgroup%uxp(ip)*pgroup%gaminv(ip)
          vy = pgroup%uyp(ip)*pgroup%gaminv(ip)
          vz = pgroup%uzp(ip)*pgroup%gaminv(ip)
          ax = (pgroup%uxp(ip) - uxpo(ip-ipmin+1))*oneondt
          ay = (pgroup%uyp(ip) - uypo(ip-ipmin+1))*oneondt
          az = (pgroup%uzp(ip) - uzpo(ip-ipmin+1))*oneondt

          delz = (pgroup%zp(ip) - zzep)
          if (az == 0.) then
            dtstar = -delz/vz
          else
            if ((vz**2 - 2.*az*delz) > 0.) then
              dm = (-vz - sqrt(vz**2 - 2.*az*delz))/az
              dp = (-vz + sqrt(vz**2 - 2.*az*delz))/az
              if (abs(dp) < abs(dm)) then
                dtstar = dp
              else
                dtstar = dm
              endif
            else
c             --- In this case, because of az, the particle will never
c             --- cross the zzep plane, so effectively decrease az
c             --- so that the particle would just touch the plane.
c             --- This is very kludgy and ad-hoc, but the calculation
c             --- of az is already equally kludgy.
              dtstar = -vz/az
            endif
          endif

c         --- Make an estimate of where the particle will be on the next step.
c         --- This assumes that az will be constant.
c         --- If it will be much closer to the z plane (smaller dtstar),
c         --- then wait until then to save the particle.
          if (lepsaveonce) then
            vznp1 = vz + 0.5*az*dt
            zznp1 = pgroup%zp(ip) + vznp1*dt
            vznp1 = vznp1 + 0.5*az*dt
            delznp1 = (zznp1 - zzep)
            dtstarnp1 = 2.*delznp1/(sqrt(vznp1**2 - 2.*az*delznp1) + vznp1)
            if (az == 0.) then
              dtstarnp1 = -delznp1/vznp1
            else
              if ((vznp1**2 - 2.*az*delznp1) > 0.) then
                dm = (-vznp1 - sqrt(vznp1**2 - 2.*az*delznp1))/az
                dp = (-vznp1 + sqrt(vznp1**2 - 2.*az*delznp1))/az
                if (abs(dp) < abs(dm)) then
                  dtstarnp1 = dp
                else
                  dtstarnp1 = dm
                endif
              else
                dtstarnp1 = -vznp1/az
              endif
            endif
c           --- The factor of 0.75 is an adhoc parameter, an attempt to
c           --- prevent particles from being missed if the estimate of
c           --- dtstarnp1 is way off.
            if (abs(dtstarnp1) < abs(dtstar)*epclosenessfactor) cycle
          endif

          nep(iw,isid) = nep(iw,isid) + 1
          nn = nep(iw,isid)

          if (lepsaveonce) then
            pgroup%pid(ip,epflagpid) = iw
          endif

          tep(nn,iw,isid) = time + dtstar
          xep(nn,iw,isid) = pgroup%xp(ip) + vx*dtstar + 0.5*ax*dtstar**2
          yep(nn,iw,isid) = pgroup%yp(ip) + vy*dtstar + 0.5*ay*dtstar**2
          uxep(nn,iw,isid) = vx + ax*dtstar
          uyep(nn,iw,isid) = vy + ay*dtstar
          uzep(nn,iw,isid) = vz + az*dtstar
          if (npid > 0) pidep(nn,:,iw,isid) = pgroup%pid(ip,:)

        enddo
      enddo

      return
      end
c=============================================================================
      subroutine initlabwn(ntlabwnestimate)
      use InDiag
      use Hist
      use Lab_Moments
      integer(ISZ):: ntlabwnestimate

c Initialize the lab window moments

      integer(ISZ):: iwin,nl

      nl = nlabwn
      nlabwn = 0
      do iwin = 1, nl
         if (zlw(iwin) < LARGEPOS) nlabwn = iwin
      enddo

      if (nlabwn > 0) then
c       --- Note that the original code used abs(nhist), but it was giving
c       --- the HPUX f90 compiler fits. The negative nhist feature was
c       --- made obsolete anyway.
        if (itlabwn == 0) itlabwn = max(1,nhist)
c       --- ntlabwnestimate is an estimate of the number of time steps across
c       --- the system, i.e. the number of times that data will be saved
c       --- assuming nhist=1. That then needs to be scale by how often the
c       --- data is collected.
        if (ntlabwn == 0) ntlabwn = max(1,ntlabwnestimate/itlabwn)
      endif
      call gchange("Lab_Moments", 0)

      return
      end
c=============================================================================
      subroutine getlabwn()
      use Beam_acc
      use InDiag
      use Picglb
      use Z_Moments
      use Z_arrays
      use Lab_Moments
 
c  Sets moments for species 1 in all axial windows in the lab frame
c  Also saves appropriate moments in time history arrays.
c  Get info from z moments
c  The data is only saved every itlabwn steps, even if the moments
c  are actually calculated more or less often.

      real(kind=8):: w0,w1,w0z,w1z
      integer(ISZ):: i,ilw,iz,izz,js

      if (nlabwn == 0) return
      if (mod(it,itlabwn) > 0) return

c     --- Make sure that the lab window moments are set to have the same
c     --- number of species as the z moments.
      if (nslabwn /= nszmmnt) then
        nslabwn = nszmmnt
        call gchange("Lab_Moments",0)
      endif

      do i = 1, nlabwn
        if (zlw(i) > (zbeam+zmmntmin) .and.
     &      zlw(i) < (zbeam+zmmntmax)) then
          if (maxval(ilabwn(i,:)) == ntlabwn) then
            ntlabwn = ntlabwn + max(10,int(ntlabwn/10))
            call gchange("Lab_Moments",0)
          endif
          iz = (zlw(i) - zbeam - zmmntmin)*dzmi
          w0 = (zlw(i) - zbeam - zmmntmin)*dzmi - iz
          w1 = 1. - w0
          izz = (zlw(i) - zbeam - zzmin)*dzzi
          w0z = (zlw(i) - zbeam - zzmin)*dzzi - izz
          w1z = 1. - w0z
          do js=0,nslabwn
            if (zmomentscalculated(js)) then
              ilabwn(i,js) = ilabwn(i,js) + 1
              ilw = ilabwn(i,js)
              timelw(ilw,i,js) = time
              pnumlw(ilw,i,js)   = w1*pnumz(iz,js)   + w0*pnumz(iz+1,js)
              xbarlw(ilw,i,js)   = w1*xbarz(iz,js)   + w0*xbarz(iz+1,js)
              ybarlw(ilw,i,js)   = w1*ybarz(iz,js)   + w0*ybarz(iz+1,js)
              xrmslw(ilw,i,js)   = w1*xrmsz(iz,js)   + w0*xrmsz(iz+1,js)
              yrmslw(ilw,i,js)   = w1*yrmsz(iz,js)   + w0*yrmsz(iz+1,js)
              rrmslw(ilw,i,js)   = w1*rrmsz(iz,js)   + w0*rrmsz(iz+1,js)
              xxpbarlw(ilw,i,js) = w1*xxpbarz(iz,js) + w0*xxpbarz(iz+1,js)
              yypbarlw(ilw,i,js) = w1*yypbarz(iz,js) + w0*yypbarz(iz+1,js)
              vzbarlw(ilw,i,js)  = w1*vzbarz(iz,js)  + w0*vzbarz(iz+1,js)
              epsxlw(ilw,i,js)   = w1*epsxz(iz,js)   + w0*epsxz(iz+1,js)
              epsylw(ilw,i,js)   = w1*epsyz(iz,js)   + w0*epsyz(iz+1,js)
              epszlw(ilw,i,js)   = w1*epszz(iz,js)   + w0*epszz(iz+1,js)
              vxrmslw(ilw,i,js)  = w1*vxrmsz(iz,js)  + w0*vxrmsz(iz+1,js)
              vyrmslw(ilw,i,js)  = w1*vyrmsz(iz,js)  + w0*vyrmsz(iz+1,js)
              vzrmslw(ilw,i,js)  = w1*vzrmsz(iz,js)  + w0*vzrmsz(iz+1,js)
              currlw(ilw,i,js)  = w1z*curr(izz,js) + w0z*curr(izz+1,js)
              lostparslw(ilw,i,js)=w1z*lostpars(izz,js)+w0z*lostpars(izz+1,js)
              if (js < nslabwn) then
c               --- This special coding is needed since the linecharge
c               --- is not calculated for each species.
                linechglw(ilw,i,js)  =
     &                     w1z*pnumz(izz  ,js)*dzmi*zmmntsw(js)*zmmntsq(js) +
     &                     w0z*pnumz(izz+1,js)*dzmi*zmmntsw(js)*zmmntsq(js)
              else
                linechglw(ilw,i,js) = w1z*linechg(izz  ) + w0z*linechg(izz+1)
              endif
            endif
          enddo
        endif
      enddo

      return
      end
c=============================================================================
      subroutine gammaadv(np,gaminv,uxp,uyp,uzp,gamadv,lrelativ)
      use Constant
      integer(ISZ):: np
      real(kind=8):: gaminv(np),uxp(np),uyp(np),uzp(np)
      character(8):: gamadv
      logical(ISZ):: lrelativ

c  Advance Gamma Inverse from it-3/2 to it-1/2 using method given by gamadv


      real(kind=8):: clghtisq,usq
      integer(ISZ):: ip

c  If not using relativity then don't do anything
      if (lrelativ) then

      clghtisq = 1./clight**2

c  Standard gamma advance
      if (gamadv == "stndrd") then
        do ip=1,np
          usq = (uxp(ip)**2 + uyp(ip)**2 + uzp(ip)**2)*clghtisq
          gaminv(ip) = 1./sqrt(1. + usq)
        enddo
      endif

c  Fast gamma advance, version 1
      if (gamadv == "fast 1") then
        do ip=1,np
          usq = (uxp(ip)**2 + uyp(ip)**2 + uzp(ip)**2)*clghtisq
          gaminv(ip) = gaminv(ip)*1.5 - (1.+usq)/gaminv(ip)*0.5
        enddo
      endif

c  Fast gamma advance, version 2, fastest, but least accurate
      if (gamadv == "fast 2") then
        do ip=1,np
          gaminv(ip)  = 1. - 0.5*(uxp(ip)**2+uyp(ip)**2+uzp(ip)**2)*clghtisq
        enddo
      endif

      endif

      return
      end
c=============================================================================
      subroutine checkz_arrays(ns)
      use Z_arrays
      use InDiag,Only: lspeciesmoments
      integer(ISZ):: ns

c     --- Set nszarr and update array allocation if needed
      if (lspeciesmoments) then
c       --- Check if the moments are to be calculated separately for
c       --- each species. If so, check if nszarr already has been set
c       --- appropriately. If not, set it and allocate the arrays.
c       --- If only one species, then don't have separate species data.
        if (nszarr < ns .and. ns > 1) then
          nszarr = ns
          call gchange("Z_arrays",0)
        endif
      else
        if (nszarr /= 0) then
          nszarr = 0
          call gchange("Z_arrays",0)
        endif
      endif

      return
      end
c=============================================================================
      subroutine setcurr(pgroup,zbeam,ns,wpid,lspeciesmoments,lzero,bound0)
      use ParticleGroupmodule
      use Z_arrays
      type(ParticleGroup):: pgroup
      integer(ISZ):: ns,wpid,bound0
      real(kind=8):: zbeam
      logical(ISZ):: lspeciesmoments,lzero

c  Sets 1d beam current as a function of z, directly from particle data

      integer(ISZ):: kk
      real(kind=8):: g,gw,w1,w0
      integer(ISZ):: is,js,ip,isid

      call checkz_arrays(ns)

      if (lzero) call zeroarry(curr,(nzzarr+1)*(1+nszarr))

      dzzi = 1. / dzz

!$OMP PARALLEL PRIVATE(ip,js,g,gw,kk,w1,w0)
      do is=1,pgroup%ns
        if (pgroup%sid(is-1) == -1) cycle

        if (lspeciesmoments) then
          js = pgroup%sid(is-1)
        else
          js = 0
        endif

        g = pgroup%sw(is)*pgroup%dtscale(is)*pgroup%sq(is)*dzzi
        gw = g
!$OMP DO
        do ip = pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
          kk = (pgroup%zp(ip) - zbeam - zzmin)*dzzi
          if (pgroup%uzp(ip) == 0 .or. kk < 0 .or. kk+1 > nzzarr) cycle
          w1 = (pgroup%zp(ip) - zbeam - zzmin)*dzzi - kk
          w0 = 1. - w1
          if (wpid > 0) gw = g*pgroup%pid(ip,wpid)
          curr(kk  ,js) = curr(kk  ,js)+w0*gw*pgroup%gaminv(ip)*pgroup%uzp(ip)
          curr(kk+1,js) = curr(kk+1,js)+w1*gw*pgroup%gaminv(ip)*pgroup%uzp(ip)
        enddo

!$OMP END DO
      enddo
!$OMP END PARALLEL

      if (lzero) call percurr(curr,nzzarr,nszarr,bound0)

      return
      end
c=============================================================================
      subroutine percurr(curr,nzzarr,nszarr,bound0)
      use GlobalVars
      integer(ISZ):: nzzarr,nszarr,bound0
      real(kind=8):: curr(0:nzzarr,0:nszarr)

c Make current periodic.

#ifdef MPIPARALLEL
c Parallel version sums curr over all processors.
      call parallelsumrealarray(curr,nzzarr+1)
#endif

      if (nszarr > 0) then
c       --- Sum up the current from all of the speices.
        curr(:,nszarr) = sum(curr(:,0:nszarr-1),2)
      endif

      if (bound0==periodic) then
        curr(0,:) = curr(0,:) + curr(nzzarr,:)
        curr(nzzarr,:) = curr(0,:)
      else
        curr(0,:) = 2.*curr(0,:)
        curr(nzzarr,:) = 2.*curr(nzzarr,:)
      endif

      return
      end
c=============================================================================
      subroutine getvzofz
      use Z_arrays

c  Computes the mean beam z velocity from the current and line charge density

      where (linechg /= 0.)
        vzofz = curr(:,nszarr)/linechg
      elsewhere
        vzofz = 0.
      end where

      return
      end
c=============================================================================
      subroutine setegap
      use Beam_acc
      use InGaps
      use Z_arrays

c  Computes the "gap" electric field as a function of z, using the "smeared gap"
c  approximation.  We sweep back from head of beam, since that's the direction
c  in which information moves.

      integer(ISZ):: i
      real(kind=8):: d

      if (ifgap) then

        do i = nzzarr-1, 1, -1
           d = dvnz(2.*rgap*dzz)
           d = dvnz(cgap*vbeam**2/dzz**2 + vbeam/d)
           egap(i-1) = (- cgap * vbeam**2 * (-2.*egap(i)+egap(i+1)) / dzz**2
     &                  + vbeam * egap(i+1) / dvnz(2.*rgap*dzz)
     &                  - egap(i) / dvnz(lgap)
     &                  + vbeam*(curr(i+1,nszarr) - curr(i-1,nszarr))/(2.*dzz)
     &                 ) / d
        enddo

      endif

      return
      end
c=============================================================================
      subroutine gapfield (np,zp,ez,zbeam,zzmin,egap,dzz)
      use InGaps
      integer(ISZ):: np
      real(kind=8):: zbeam,zzmin,dzz
      real(kind=8):: zp(np),ez(np),egap(0:*)

c  Interpolates the spatially averaged gap electric field onto the particle
c  positions, then uses it to augment the individual particle Ez's.

      real(kind=8):: dzzi,w1,w0
      integer(ISZ):: ip,kk

      if (ifgap) then

        dzzi = 1./dzz
        do ip=1,np
           kk = (zp(ip) - zbeam - zzmin) * dzzi
           w1 = (zp(ip) - zbeam - zzmin) * dzzi - kk
           w0 = 1. - w1
c          if (kk < 0 .or. kk+1 > nzzarr) cycle
           ez(ip) = ez(ip) + w0*egap(kk) + w1*egap(kk+1)
        enddo

      endif

      return
      end
c=============================================================================
      subroutine setgamma(pgroup,lrelativ)
      use ParticleGroupmodule
      use Constant
      use InPart
      type(ParticleGroup):: pgroup
      logical(ISZ):: lrelativ
c  Converts v to u, sets inverse gamma factor for all particles


      real(kind=8):: clghtisq,gamma
      integer(ISZ):: is,ip

      if (lrelativ) then
        clghtisq = 1./clight**2
        do is=1,pgroup%ns
          do ip = pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
             gamma = 1. / sqrt ( 1. -
     &    (pgroup%uxp(ip)**2 + pgroup%uyp(ip)**2 + pgroup%uzp(ip)**2)*clghtisq)
             pgroup%uxp(ip) = gamma*pgroup%uxp(ip)
             pgroup%uyp(ip) = gamma*pgroup%uyp(ip)
             pgroup%uzp(ip) = gamma*pgroup%uzp(ip)
             pgroup%gaminv(ip) = 1./gamma
          enddo
        enddo
      else
        do is=1,pgroup%ns
          do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
             pgroup%gaminv(ip) = 1.0
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine zpartbnd(pgroup,zmmax,zmmin,dz)
      use ParticleGroupmodule
      use Picglb, only: it
      use GlobalVars
      use InGen
      use InPart
      use Subcycling,Only: ndtstorho,zgridndts
      type(ParticleGroup):: pgroup
      real(kind=8):: zmmax,zmmin,dz

c  Impose boundary conditions on zp

      real(kind=8):: syslen,zgrid
      integer(ISZ):: is,ip,indts

c Parallel version makes call to special routine zpartbnd_slave
#ifdef MPIPARALLEL
      call zpartbnd_slave(pgroup,zmmax,zmmin,dz)
#else

      if (pbound0==periodic) then
c     --- periodic boundary condition
        syslen = zmmax - zmmin
        do is=1,pgroup%ns
          if (.not. pgroup%ldts(is-1)) cycle
          indts = ndtstorho(pgroup%ndts(is-1))
          zgrid = zgridndts(indts)
          do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
            if (pgroup%zp(ip) < zgrid + zmmin) then
              pgroup%zp(ip) = pgroup%zp(ip) +
     &                      syslen*int((zgrid + zmmax - pgroup%zp(ip))/syslen)
            endif
            if (pgroup%zp(ip) > zgrid + zmmax) then
              pgroup%zp(ip) = pgroup%zp(ip) -
     &                      syslen*int((pgroup%zp(ip) - zgrid - zmmin)/syslen)
            endif
          enddo
        enddo
      else
c     --- other boundary condition (particles absorbed or reflected)
        do is=1,pgroup%ns
          if (.not. pgroup%ldts(is-1)) cycle
          indts = ndtstorho(pgroup%ndts(is-1))
          zgrid = zgridndts(indts)
          do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
            if ((pgroup%zp(ip)-zgrid) >= zmmax) then
              if(pboundnz==absorb) then
                pgroup%gaminv(ip) = 0.
              elseif(pboundnz==reflect) then
                pgroup%zp(ip) = 2.*zgrid + 2.*zmmax - pgroup%zp(ip)
                pgroup%uzp(ip) = -pgroup%uzp(ip)
              end if
            elseif ((pgroup%zp(ip)-zgrid) <= zmmin) then
              if(pbound0==absorb) then
                pgroup%gaminv(ip) = 0.
              elseif(pbound0==reflect) then
                pgroup%zp(ip) = 2.*zgrid + 2.*zmmin - pgroup%zp(ip)
                pgroup%uzp(ip) = -pgroup%uzp(ip)
              end if
            endif
          enddo
        enddo
      endif

#endif
c --- end of ifelse above for paralle version.

      return
      end
c=============================================================================
      subroutine zpartbndwithdata(n,z,uz,gaminv,zmmax,zmmin,dz,zgrid)
      use GlobalVars
      use InGen
      use InPart
      integer(ISZ):: n
      real(kind=8):: z(n),uz(n),gaminv(n)
      real(kind=8):: zmmax,zmmin,dz,zgrid

c  Impose boundary conditions on z

      real(kind=8):: syslen
      integer(ISZ):: is,ip

      if (pbound0==periodic) then
c     --- periodic boundary condition
        syslen = zmmax - zmmin
        do ip=1,n
          if (z(ip) < zgrid + zmmin)
     &      z(ip) = z(ip) + syslen*int((zgrid + zmmax - z(ip))/syslen)
          if (z(ip) > zgrid + zmmax)
     &      z(ip) = z(ip) - syslen*int((z(ip) - zgrid - zmmin)/syslen)
        enddo
      else
c     --- other boundary condition (particles absorbed or reflected)
        do ip=1,n
          if ((z(ip)-zgrid) >= zmmax) then
            if(pboundnz==absorb) then
              gaminv(ip) = 0.
            elseif(pboundnz==reflect) then
              z(ip) = 2.*zgrid + 2.*zmmax - z(ip)
              uz(ip) = -uz(ip)
            end if
          elseif ((z(ip)-zgrid) <= zmmin) then
            if(pbound0==absorb) then
              gaminv(ip) = 0.
            elseif(pbound0==reflect) then
              z(ip) = 2.*zgrid + 2.*zmmin - z(ip)
              uz(ip) = -uz(ip)
            end if
          endif
        enddo
      endif

      return
      end
c=============================================================================
      subroutine reorgparticles(pgroup)
      use ParticleGroupmodule
      type(ParticleGroup):: pgroup
#ifdef MPIPARALLEL
      call reorgparticles_parallel(pgroup)
#endif
      return
      end
c=============================================================================
      subroutine semitransparent_disc(pgroup,dz)
      use ParticleGroupmodule
      use Picglb, only: it
      use InGen
      use Moments
      use InPart
      use SemiTransparentDisc
      type(ParticleGroup):: pgroup
      real(kind=8), intent(in) :: dz 

c Random absorption of particles passing through a semitranparent disc.
c This assumes that no particle travels more than dz in one time step.
c dz is used only for efficiency purpose: only particles within dz 
c of semi-transparent disc location are checked for disc crossing.

      integer(ISZ) :: is, ip, izd, idisc
      real(kind=8) :: zpold, dtcross, xpatdisc, ypatdisc, r, wranf
     &              , zmind(n_STdiscs),zmaxd(n_STdiscs)

      if(n_STdiscs==0) return
      if(clearlostpart==0) clearlostpart = 1

      do idisc = 1, n_STdiscs
        zmind(idisc) = z_STdiscs(idisc)-dz
        zmaxd(idisc) = z_STdiscs(idisc)+dz
      end do

c     --- loop on species
      do is=1,pgroup%ns
c       --- do it only if species was advanced on this time step
        if (.not. pgroup%ldts(is-1)) cycle
c       --- loop on particles
        do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
c         --- loop on number of discs
          do idisc = 1, n_STdiscs
c           --- do only the work for particles close to disc location
            if(pgroup%zp(ip)>=zmind(idisc) .and.
     &         pgroup%zp(ip)<=zmaxd(idisc)) then
              zpold = pgroup%zp(ip)-dt*pgroup%uzp(ip)*pgroup%gaminv(ip)
c             --- test if particle crossed disc location
              if(zpold<z_STdiscs(idisc) .and.
     &           pgroup%zp(ip)>=z_STdiscs(idisc)) then
                dtcross = dt*(pgroup%zp(ip) - z_STdiscs(idisc))/
     &                       (pgroup%zp(ip) - zpold)
                xpatdisc = pgroup%xp(ip) -
     &                     dtcross*pgroup%uxp(ip)*pgroup%gaminv(ip)
                ypatdisc = pgroup%yp(ip) -
     &                     dtcross*pgroup%uyp(ip)*pgroup%gaminv(ip)
                r = sqrt(xpatdisc**2+ypatdisc**2)
c               --- test if particle passed through disc
                if(r<=r_STdiscs(idisc)) then
c                  --- randomly absorb particles according to transparency
c                  --- coefficient
                   if(wranf()>t_STdiscs(idisc)) then
                     pgroup%gaminv(ip) = 0.
                   end if
                end if
              end if
            end if
          end do
        end do
      end do
      
      end
c=============================================================================
c=============================================================================
c The following two subroutines are made obsolete by zpartbnd above, but
c are maintained for use by the wrz and wxy packages.
c=============================================================================
      subroutine periz(np,zp,zgrid,zmmax,zmmin)
      integer(ISZ):: np
      real(kind=8):: zgrid,zmmax,zmmin
      real(kind=8):: zp(np)

c  Imposes periodic boundary conditions on zp

      real(kind=8):: syslen,sysleni
      integer(ISZ):: ip

      syslen = zmmax - zmmin
      sysleni = 1./syslen
      do ip=1,np
        if (zp(ip) < zgrid + zmmin)
     &    zp(ip) = zp(ip) + syslen*int((zgrid + zmmax - zp(ip))*sysleni)
        if (zp(ip) > zgrid + zmmax)
     &    zp(ip) = zp(ip) - syslen*int((zp(ip) - zgrid - zmmin)*sysleni)
      enddo

      return
      end
c=============================================================================
      subroutine stckyz(np,zp,zmmax,zmmin,dz,uxp,uyp,uzp,gaminv,zgrid)
      integer(ISZ):: np
      real(kind=8):: zmmax,zmmin,dz,zgrid
      real(kind=8):: zp(np)
      real(kind=8):: uxp(np),uyp(np),uzp(np),gaminv(np)

c  Enforces sticky b.c.'s on the z walls.
c  Particles stick at the extreme edges of the grid.

      integer(ISZ):: ip

      do ip=1,np
        if ((zp(ip)-zgrid) >= zmmax) then
          gaminv(ip) = 0.
        elseif ((zp(ip)-zgrid) <= zmmin) then
          gaminv(ip) = 0.
        endif
      enddo

      return
      end
c=============================================================================
c=============================================================================
      real(kind=8) function getbeamcom(pgroup)
      use ParticleGroupmodule
      use InPart
      type(ParticleGroup):: pgroup

c Calculate the center of mass of the beam along the z axis from the particles

      integer(ISZ):: is,ip
      real(kind=8):: totz,totmass

      totz = 0.
      totmass = 0.
      do is=1,pgroup%ns
        do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
          if (pgroup%uzp(ip) /= 0.) then
            totz = totz + pgroup%zp(ip)*pgroup%sm(is)
            totmass = totmass + pgroup%sm(is)
          endif
        enddo
      enddo

      getbeamcom = totz/dvnz(totmass)

      return
      end
c=============================================================================
c=============================================================================
      subroutine setuppgroup(pgroup)
      use ParticleGroupmodule
      use Particles,Only: npid
      type(ParticleGroup):: pgroup

      pgroup%npid = npid
c     call gchange("Particles",0)
      call ParticleGroupchange(pgroup)

      return
      end
c=============================================================================
      subroutine alotpart(pgroup)
      use ParticleGroupmodule
      use Subtimerstop
      use InPart
      use Particles,Only: npmax,npid
      type(ParticleGroup):: pgroup

c Allocate space for particles and set up associated arrays.

      integer(ISZ):: is,sid(0:ns)
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Either npmax or np_s can be used as input quantities. If one of
c     --- them is zero, set it from the other. If both are set, the code
c     --- uses np_s, ignoring npmax.
      do is=1,ns
        if (np_s(is) == 0) np_s(is) = npmax*sp_fract(is)
      enddo
      if (pgroup%npmax == 0) pgroup%npmax = sum(np_s)

c     --- Force pgroup%npmax to be large enough to hold all of the particles
c     --- to be loaded.
      pgroup%npmax = sum(np_s)

c     --- Allocate space. A gchange is done primarily so that any data
c     --- in ipmax_s is preserved.
      pgroup%ns = ns
      call setuppgroup(pgroup)

c     --- setup ipmax_s array
      pgroup%ipmax_s(0) = 0
      do is=1,pgroup%ns
        pgroup%ipmax_s(is) = pgroup%ipmax_s(is-1) + np_s(is)
        pgroup%ins(is) = pgroup%ipmax_s(is-1) + 1
        pgroup%nps(is) = 0
        if (pgroup%sid(is-1) == -1) pgroup%sid(is-1) = is-1
      enddo

!$OMP MASTER
      if (ltoptimesubs) timealotpart = timealotpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine chckpart(pgroup,is,nlower,nhigher,lfullshft)
      use ParticleGroupmodule
      use Subtimerstop
      use InPart
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,nlower,nhigher
      logical(ISZ):: lfullshft

c Make sure that there is enough space in the particle arrays for nlower
c new particles below and nhigher above the live particles.  Returns if
c there is already enough space above and below.  If there is enough total
c space but not enough room above or below, the particles are shifted
c appropriately.  If lfullshft is true, the particles are shifted all the
c way to the end of the species block minus nhigher.  If there is not
c enough space, add more to the arrays. Particle data is shifted
c appropriately.

      integer(ISZ):: nadd,i,ishift
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Call this here in case chckpart has been called before alotpart.
      call setuppgroup(pgroup)

c     --- Make sure that ins has a sensible value
      if (pgroup%ins(is) == 0 .and. pgroup%ipmax_s(is) == 0) pgroup%ins(is) = 1

c     --- If there is already enough space, then return.
      if (pgroup%ins(is)-pgroup%ipmax_s(is-1)-1 >= nlower .and.
     &    pgroup%ipmax_s(is)-pgroup%ins(is)-pgroup%nps(is)+1 >= nhigher) return

c     --- If there is not enough space, then allocate new space.
      if (pgroup%ins(is)-pgroup%ipmax_s(is-1)-1 +
     &    pgroup%ipmax_s(is)-pgroup%ins(is)-pgroup%nps(is)+1 <
     &      nlower + nhigher) then

c       --- Amount of new space to add: this is calculated to be a sizable
c       --- chunk, but not too big.  By default, add 10 times the space needed.
c       --- Limit that to 10000, except, if the space needed is greater than
c       --- 10000, add the amount of space needed.
        nadd = nlower + nhigher -
     &   (pgroup%ins(is)-pgroup%ipmax_s(is-1)-1 +
     &    pgroup%ipmax_s(is)-pgroup%ins(is)-pgroup%nps(is)+1)
        nadd = min(10*nadd, max(10000, nadd))

c       --- change the particle array dimension variables and allot space
        pgroup%npmax = pgroup%npmax + nadd
        call ParticleGroupchange(pgroup)

c       --- Loop over species species above 'is', shifting
c       --- them up to make space below for species is.
        do i=pgroup%ns,is+1,-1
          pgroup%ipmax_s(i) = pgroup%ipmax_s(i) + nadd
          call shftpart(pgroup,i,nadd)
        enddo

c       --- Increase size of the block for species 'is'.
        pgroup%ipmax_s(is) = pgroup%ipmax_s(is) + nadd

      endif

c     --- Now, there is enough room to shift the particles appropriately.
c     --- Even though there may be enough space needed for the new particles,
c     --- the shift should still be done for the case when lfullshft is true.
c     --- Note that is some cases no shift will be needed.  In those cases,
c     --- ishift comes out to be zero and shftpart returns immediately when
c     --- ishift is zero.
      if (lfullshft) then
        ishift = pgroup%ipmax_s(is) - pgroup%ins(is) - pgroup%nps(is) + 1 - nhigher
      else
        ishift = max(0,nlower - (pgroup%ins(is)-pgroup%ipmax_s(is-1)-1)) +
     &           min(0,(pgroup%ipmax_s(is)-pgroup%ins(is)-pgroup%nps(is)+1) -
     &                 nhigher)
      endif
      call shftpart(pgroup,is,ishift)

!$OMP MASTER
      if (ltoptimesubs) timechckpart = timechckpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine shftpart(pgroup,is,ishift)
      use ParticleGroupmodule
      use InPart
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,ishift
      if (pgroup%npid > 0) then
        call shftpartwork(is,ishift,pgroup%npmax,pgroup%ns,
     &                    pgroup%ipmax_s,pgroup%ins,pgroup%nps,pgroup%npid,
     &                    pgroup%xp,pgroup%yp,pgroup%zp,
     &                    pgroup%uxp,pgroup%uyp,pgroup%uzp,pgroup%gaminv,
     &                    .true.,
     &                    pgroup%ex,pgroup%ey,pgroup%ez,
     &                    pgroup%bx,pgroup%by,pgroup%bz,
     &                    pgroup%pid)
      else
        call shftpartwork(is,ishift,pgroup%npmax,pgroup%ns,
     &                    pgroup%ipmax_s,pgroup%ins,pgroup%nps,pgroup%npid,
     &                    pgroup%xp,pgroup%yp,pgroup%zp,
     &                    pgroup%uxp,pgroup%uyp,pgroup%uzp,pgroup%gaminv,
     &                    .true.,
     &                    pgroup%ex,pgroup%ey,pgroup%ez,
     &                    pgroup%bx,pgroup%by,pgroup%bz,
     &                    0.)
      endif
      return
      end
c=============================================================================
      subroutine shftpartwork(is,ishift,npmax,ns,ipmax_s,ins,nps,
     &                        npid,xp,yp,zp,uxp,uyp,uzp,gaminv,
     &                        lshiftfields,ex,ey,ez,bx,by,bz,pid)
      use Subtimerstop
      integer(ISZ):: is,ishift
      integer(ISZ):: npmax,ns,ipmax_s(0:ns),ins(ns),nps(ns),npid
      real(kind=8):: xp(npmax),yp(npmax),zp(npmax)
      real(kind=8):: uxp(npmax),uyp(npmax),uzp(npmax)
      real(kind=8):: gaminv(npmax),pid(npmax,npid)
      logical(ISZ):: lshiftfields
      real(kind=8):: ex(npmax),ey(npmax),ez(npmax)
      real(kind=8):: bx(npmax),by(npmax),bz(npmax)

c Shift particles by 'ishift'.  All of the data from the low end of the
c species block to the top of the live particles are moved.  The area
c below the live particles may include newly injected particles which
c have not been included into the simulation yet.
c Assumes that the is enought space in the particle arrays to make the
c shift without clobbering other particles or going past the ends of
c the arrays.

      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- return if shift is zero
      if (ishift == 0) return

c     --- Positive and negative shifts are treated differently since
c     --- a negative shift may result in lost data.
      if (ishift > 0) then

c       --- Loop from top of live particles to bottom of species block.
        do ip=ins(is)+nps(is)-1, ipmax_s(is-1)+1, -1
          xp(ip+ishift) = xp(ip)
          yp(ip+ishift) = yp(ip)
          zp(ip+ishift) = zp(ip)
          uxp(ip+ishift) = uxp(ip)
          uyp(ip+ishift) = uyp(ip)
          uzp(ip+ishift) = uzp(ip)
          gaminv(ip+ishift) = gaminv(ip)
          if (lshiftfields) then
            ex(ip+ishift) = ex(ip)
            ey(ip+ishift) = ey(ip)
            ez(ip+ishift) = ez(ip)
            bx(ip+ishift) = bx(ip)
            by(ip+ishift) = by(ip)
            bz(ip+ishift) = bz(ip)
          endif
          if (npid > 0) then
            pid(ip+ishift,:) = pid(ip,:)
          endif
        enddo
c       --- Move this into a seperate loop since some compilers were screwing
c       --- up. This will also be slightly faster too, maybe.
        do ip=ipmax_s(is-1)+1,ipmax_s(is-1)+1+ishift-1
          uzp(ip) = 0.
          if (npid > 0) pid(ip,:) = 0
        enddo

c       --- change lower bound of live particles
        ins(is) = ins(is) + ishift

      else
c       --- Loop from -ishift above bottom of species block to top of live
c       --- particles.  Note that particles will be clobbered if
c       --- (ins(is) - ipmax_s(is-1)) < -ishift
        if (ins(is) - ipmax_s(is-1) < -ishift) then
          print*,"WARNING: shftpart: particle data will be lost since there is"
          print*,"                   not enough space in the arrays."
          stop
        endif

        do ip=ipmax_s(is-1)+1-ishift, ins(is)+nps(is)-1
          xp(ip+ishift) = xp(ip)
          yp(ip+ishift) = yp(ip)
          zp(ip+ishift) = zp(ip)
          uxp(ip+ishift) = uxp(ip)
          uyp(ip+ishift) = uyp(ip)
          uzp(ip+ishift) = uzp(ip)
          gaminv(ip+ishift) = gaminv(ip)
          if (lshiftfields) then
            ex(ip+ishift) = ex(ip)
            ey(ip+ishift) = ey(ip)
            ez(ip+ishift) = ez(ip)
            bx(ip+ishift) = bx(ip)
            by(ip+ishift) = by(ip)
            bz(ip+ishift) = bz(ip)
          endif
          if (npid > 0) then
            pid(ip+ishift,:) = pid(ip,:)
          endif
        enddo
c       --- Move this into a seperate loop since some compilers were screwing
c       --- up. This will also be slightly faster too, maybe.
        do ip=ins(is)+nps(is)-1,ins(is)+nps(is)-1+ishift+1,-1
          uzp(ip) = 0.
          if (npid > 0) pid(ip,:) = 0
        enddo

c       --- change number of particles if data is lost
        if (ins(is) - ipmax_s(is-1) < -ishift) then
          nps(is) = nps(is) - (-ishift - (ins(is) - ipmax_s(is-1)))
          ins(is) = ipmax_s(is-1) + 1
        else
c         --- change lower bound only
          ins(is) = ins(is) + ishift
        endif

      endif

!$OMP MASTER
      if (ltoptimesubs) timeshftpartwork = timeshftpartwork + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine copypart(pgroup,it,nn,ii,istart)
      use ParticleGroupmodule
      use Subtimerstop
      use InPart
      type(ParticleGroup):: pgroup
      integer(ISZ):: it,nn,istart
      integer(ISZ):: ii(0:nn-1)

c Copies particle data from locations given by ii to the locations starting
c at 'it', or if istart > 0, copy sequential particles from locations
c starting at 'istart' to locations starting at 'it'.

      integer(ISZ):: i,in,ij
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (istart > 0) then
        do i=0,nn-1
          ij = istart + i
          in = it + i
          pgroup%xp(in) = pgroup%xp(ij)
          pgroup%yp(in) = pgroup%yp(ij)
          pgroup%zp(in) = pgroup%zp(ij)
          pgroup%uxp(in) = pgroup%uxp(ij)
          pgroup%uyp(in) = pgroup%uyp(ij)
          pgroup%uzp(in) = pgroup%uzp(ij)
          pgroup%gaminv(in) = pgroup%gaminv(ij)
          pgroup%ex(in) = pgroup%ex(ij)
          pgroup%ey(in) = pgroup%ey(ij)
          pgroup%ez(in) = pgroup%ez(ij)
          pgroup%bx(in) = pgroup%bx(ij)
          pgroup%by(in) = pgroup%by(ij)
          pgroup%bz(in) = pgroup%bz(ij)
          if (pgroup%npid > 0) then
            pgroup%pid(in,:) = pgroup%pid(ij,:)
          endif
        enddo
      else
        do i=0,nn-1
          ij = ii(i)
          in = it + i
          pgroup%xp(in) = pgroup%xp(ij)
          pgroup%yp(in) = pgroup%yp(ij)
          pgroup%zp(in) = pgroup%zp(ij)
          pgroup%uxp(in) = pgroup%uxp(ij)
          pgroup%uyp(in) = pgroup%uyp(ij)
          pgroup%uzp(in) = pgroup%uzp(ij)
          pgroup%gaminv(in) = pgroup%gaminv(ij)
          pgroup%ex(in) = pgroup%ex(ij)
          pgroup%ey(in) = pgroup%ey(ij)
          pgroup%ez(in) = pgroup%ez(ij)
          pgroup%bx(in) = pgroup%bx(ij)
          pgroup%by(in) = pgroup%by(ij)
          pgroup%bz(in) = pgroup%bz(ij)
          if (pgroup%npid > 0) then
            pgroup%pid(in,:) = pgroup%pid(ij,:)
          endif
        enddo
      endif

!$OMP MASTER
      if (ltoptimesubs) timecopypart = timecopypart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine addpart(pgroup,nn,nnpid,x,y,z,vx,vy,vz,gi,ppid,is,
     &                   lallindomain,zmmin,zmmax,lmomentum)
      use ParticleGroupmodule
      use Subtimerstop
      use Constant
      use Particles,Only: spid,ssn
      use Beam_acc
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,nn,nnpid
      real(kind=8):: x(nn),y(nn),z(nn),vx(nn),vy(nn),vz(nn),gi(nn)
      real(kind=8):: ppid(nn,nnpid)
      logical(ISZ):: lallindomain,lmomentum
      real(kind=8):: zmmin,zmmax

c Adds new particles to the simulation
c   nn: number of particles to add
c   x,y,z,vx,vy,vz,gi: coordinates and velocities
c   is: species number to add particles to
c   lallindomain: When true, all particles are assumed to be within the
c                 z extent of the domain and the scraping check is skipped
c   zmmin,zmmax: z extent of the domain
c   lmomentum: Set to false when velocites are input as velocities, and true
c              when input as massless momentum. Only used when lrelativ is
c              true.

      integer(ISZ):: ip,n1,i1,i2,i
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call setuppgroup(pgroup)
      if (pgroup%sid(is-1) == -1) pgroup%sid(is-1) = is-1

c     --- Make sure that the particles are within the z extent of the domain,
c     --- unless the input declares this to be so. Note that no check is made
c     --- transversely.
      if (.not. lallindomain) then
        n1 = 0
        do ip=1,nn
          if (zmmin <= z(ip) .and. z(ip) < zmmax) then
            n1 = n1 + 1
            if (n1 < ip) then
              x(n1) = x(ip)
              y(n1) = y(ip)
              z(n1) = z(ip)
              vx(n1) = vx(ip)
              vy(n1) = vy(ip)
              vz(n1) = vz(ip)
              gi(n1) = gi(ip)
              if (nnpid > 0) ppid(n1,:) = ppid(ip,:)
            endif
          endif
        enddo
      else
        n1 = nn
      endif

c     --- If velocities input are not momentum, and lrelativ is true,
c     --- calculate gamma inverse and convert velocities to massless
c     --- momentum.
      if (.not. lmomentum .and. lrelativ) then
        gi(1:n1) = sqrt(1.- (vx(1:n1)**2+vy(1:n1)**2+vz(1:n1)**2)/clight**2)
        vx(1:n1) = vx(1:n1)/gi(1:n1)
        vy(1:n1) = vy(1:n1)/gi(1:n1)
        vz(1:n1) = vz(1:n1)/gi(1:n1)
      endif

c     --- Make room for the particles
c     pgroup%npid = max(pgroup%npid,nnpid)
      call chckpart(pgroup,is,0,n1,.false.)

c     --- Now, copy data into arrays
      i1 = pgroup%ins(is) + pgroup%nps(is)
      i2 = i1 + n1 - 1
      pgroup%xp(i1:i2) = x(1:n1)
      pgroup%yp(i1:i2) = y(1:n1)
      pgroup%zp(i1:i2) = z(1:n1)
      pgroup%uxp(i1:i2) = vx(1:n1)
      pgroup%uyp(i1:i2) = vy(1:n1)
      pgroup%uzp(i1:i2) = vz(1:n1)
      pgroup%gaminv(i1:i2) = gi(1:n1)
      if (nnpid > 0) pgroup%pid(i1:i2,1:nnpid) = ppid(1:n1,:)
      if (spid > 0) then
c       --- Only set the ssn if it is not already set.
        do ip=i1,i2
          if (pgroup%pid(ip,spid) == 0) then
            pgroup%pid(ip,spid) = ssn
            ssn = ssn + 1
          endif
        enddo
      end if
      pgroup%nps(is) = pgroup%nps(is) + n1

!$OMP MASTER
      if (ltoptimesubs) timeaddpart = timeaddpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine clearpart(pgroup,is,fillmethod)
      use ParticleGroupmodule
      use Subtimerstop
      use InPart
      use Particles,Only: wpid
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,fillmethod

c Clears out all of the lost particles. If is < 0, the clears out all
c species, otherwise only species requested.
c The empty spaces can be fill in multiple ways, based on value of fillmethod
c   0 : don't do anything
c   1 : fill spaces with particles from the end of the list
c   2 : shift particles above down (they keep the same order)

      integer(ISZ):: is1,is2,i,ip
      integer(ISZ):: ips,ds,ipslimit
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (is < 0) then
        is1 = 1
        is2 = pgroup%ns
      else
        is1 = is
        is2 = is
      endif

      if (fillmethod == 1) then

c       --- Live particles from the end of the species block are moved
c       --- to fill in the empty spaces.
        do i=is1,is2
          ips = pgroup%ins(i) + pgroup%nps(i) - 1
          ip = pgroup%ins(i)
          pgroup%nps(i) = 0
          do while (ip <= ips)
            if (pgroup%uzp(ip) == 0.) then
              do while (pgroup%uzp(ips) == 0. .and. ips > ip)
                ips = ips - 1
              enddo
              if (ips > ip) then
                pgroup%xp(ip) = pgroup%xp(ips)
                pgroup%yp(ip) = pgroup%yp(ips)
                pgroup%zp(ip) = pgroup%zp(ips)
                pgroup%uxp(ip) = pgroup%uxp(ips)
                pgroup%uyp(ip) = pgroup%uyp(ips)
                pgroup%uzp(ip) = pgroup%uzp(ips)
                pgroup%gaminv(ip) = pgroup%gaminv(ips)
                pgroup%ex(ip) = pgroup%ex(ips)
                pgroup%ey(ip) = pgroup%ey(ips)
                pgroup%ez(ip) = pgroup%ez(ips)
                pgroup%bx(ip) = pgroup%bx(ips)
                pgroup%by(ip) = pgroup%by(ips)
                pgroup%bz(ip) = pgroup%bz(ips)
                if (pgroup%npid > 0) then
                  pgroup%pid(ip,:) = pgroup%pid(ips,:)
                  pgroup%pid(ips,:) = 0
                endif
                pgroup%uzp(ips) = 0.
                pgroup%gaminv(ips) = 0.
                if(wpid>0) pgroup%pid(ips,wpid) = 0.
                ips = ips - 1
                pgroup%nps(i) = pgroup%nps(i) + 1
              endif
            else
              pgroup%nps(i) = pgroup%nps(i) + 1
            endif
            ip = ip + 1
          enddo
        enddo

      else if (fillmethod == 2) then

c       --- Live particles are shifted downward to fill in the empty spaces
        do i=is1,is2
          ipslimit = pgroup%ins(i)+pgroup%nps(i)-1
          pgroup%nps(i) = 0
          ip = pgroup%ins(i)
          do ips=ip,ipslimit
            if (pgroup%uzp(ips) /= 0.) then
              if (ips > ip) then
                pgroup%xp(ip) = pgroup%xp(ips)
                pgroup%yp(ip) = pgroup%yp(ips)
                pgroup%zp(ip) = pgroup%zp(ips)
                pgroup%uxp(ip) = pgroup%uxp(ips)
                pgroup%uyp(ip) = pgroup%uyp(ips)
                pgroup%uzp(ip) = pgroup%uzp(ips)
                pgroup%gaminv(ip) = pgroup%gaminv(ips)
                pgroup%ex(ip) = pgroup%ex(ips)
                pgroup%ey(ip) = pgroup%ey(ips)
                pgroup%ez(ip) = pgroup%ez(ips)
                pgroup%bx(ip) = pgroup%bx(ips)
                pgroup%by(ip) = pgroup%by(ips)
                pgroup%bz(ip) = pgroup%bz(ips)
                if (pgroup%npid > 0) then
                  pgroup%pid(ip,:) = pgroup%pid(ips,:)
                  pgroup%pid(ips,:) = 0
                endif
                pgroup%uzp(ips) = 0.
                pgroup%gaminv(ips) = 0.
                if(wpid>0) pgroup%pid(ips,wpid) = 0.
              endif
              ip = ip + 1
              pgroup%nps(i) = pgroup%nps(i) + 1
            endif
          enddo
        enddo

      endif

!$OMP MASTER
      if (ltoptimesubs) timeclearpart = timeclearpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine shrinkpart(pgroup)
      use ParticleGroupmodule
      use Subtimerstop
      use InPart
      type(ParticleGroup):: pgroup

c Removes all unused space from the particle arrays. Particles are shifted
c down and then gchange is called to free the left over space.

      integer(ISZ):: is,oldnpmax
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      oldnpmax = pgroup%npmax

      call shftpart(pgroup,1,1-pgroup%ins(1))
      pgroup%ins(1) = 1
      pgroup%ipmax_s(1) = pgroup%ins(1) + pgroup%nps(1) - 1

      do is=2,pgroup%ns
        call shftpart(pgroup,is,pgroup%ins(is-1)+pgroup%nps(is-1)-pgroup%ins(is))
        pgroup%ins(is) = pgroup%ins(is-1)+pgroup%nps(is-1)
        pgroup%ipmax_s(is) = pgroup%ins(is) + pgroup%nps(is) - 1
      enddo

      pgroup%npmax = pgroup%ipmax_s(pgroup%ns)

c     --- Force npmax to change to ensure that the particles arrays are
c     --- reallocated.
      if (pgroup%npmax == oldnpmax) pgroup%npmax = pgroup%npmax + 1

      call ParticleGroupchange(pgroup)

!$OMP MASTER
      if (ltoptimesubs) timeshrinkpart = timeshrinkpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine processlostpart(pgroup,is,clearlostpart,time,zbeam)
      use ParticleGroupmodule
      use Subtimerstop
      use Constant
      use Beam_acc
      use InPart
      use InDiag
      use Particles,Only: wpid
      use LostParticles
      use Z_arrays
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,clearlostpart
      real(kind=8):: time,zbeam

c Processes lost particles, which are flagged by having gaminv set to zero.
c The particle data is optionally saved. The velocities are set to zero and
c the routine which clears out the lost particles is called.

      integer(ISZ):: nlost,i1,i2,js,isid
      integer(ISZ):: ip,iz
      real(kind=8):: usq
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (lsavelostpart) then

c       --- Make sure that if npid has changed, then npidlost is also
c       --- changed.
        if (pgroup%npid > npidlost) then
          npidlost = pgroup%npid
          call gchange("LostParticles",0)
        endif

c       --- Get the number of lost particles and make sure there is
c       --- enough space to save them in the LostParticles arrays.
        i1 = pgroup%ins(is)
        i2 = pgroup%ins(is)+pgroup%nps(is)-1
        nlost = count((pgroup%gaminv(i1:i2)==0 .and. pgroup%uzp(i1:i2)/=0.))
        call chcklostpart(is,0,nlost)

      endif

c     --- Get global species number
      isid = pgroup%sid(is-1) + 1
      if (isid == 0) return

c     --- Check lspeciesmoments, since it affects lostpars
      call checkz_arrays(ns)
      if (lspeciesmoments) then
        js = isid - 1
      else
        js = 0
      endif

c     --- Loop over the particles arrays and process the ones with gaminv
c     --- set to zero, the flag for scraped particles.
      do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
        if (pgroup%gaminv(ip) == 0.) then

c         --- Accumulate the array holding the location where particles
c         --- are lost.
          if (pgroup%uzp(ip) /= 0.) then
            iz = max(0,min(nzzarr,int((pgroup%zp(ip) - zzmin - zbeam)*dzzi)))
            lostpars(iz,js) = lostpars(iz,js) + 1
            if (nszarr > 0) lostpars(iz,nszarr) = lostpars(iz,nszarr) + 1
          endif

c         --- If requested to save lost particles, do it. Particle data is
c         --- only save if uzp /= 0, since if uzp is 0, the particle was
c         --- lost on a previous step but wasn't cleared out of the array.
          if (lsavelostpart .and. pgroup%uzp(ip) /= 0.) then
            i1 = inslost(isid) + npslost(isid)
            xplost(i1) = pgroup%xp(ip)
            yplost(i1) = pgroup%yp(ip)
            zplost(i1) = pgroup%zp(ip)
            uxplost(i1) = pgroup%uxp(ip)
            uyplost(i1) = pgroup%uyp(ip)
            uzplost(i1) = pgroup%uzp(ip)
            if (lrelativ) then
              usq = pgroup%uxp(ip)**2 + pgroup%uyp(ip)**2 + pgroup%uzp(ip)**2
              gaminvlost(i1) = 1./sqrt(1. + usq/clight**2)
            else
              gaminvlost(i1) = 1.
            endif
            tplost(i1) = time
            if (pgroup%npid > 0) then
              pidlost(i1,1:pgroup%npid) = pgroup%pid(ip,:)
            endif
            if (npidlost > pgroup%npid) pidlost(i1,pgroup%npid+1:)=0.
            npslost(isid) = npslost(isid) + 1
          endif

c         --- Zero out the velocity, the flag for a lost particle.
          pgroup%uxp(ip) = 0.0
          pgroup%uyp(ip) = 0.0
          pgroup%uzp(ip) = 0.0
          if(wpid>0) pgroup%pid(ip,wpid) = 0.

        endif
      enddo

      call clearpart(pgroup,is,clearlostpart)

!$OMP MASTER
      if (ltoptimesubs) timeprocesslostpart = timeprocesslostpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine particlesortyzwithcopy(pgroup,dy,dz,ymmin,zmmin,ny,nz)
      use ParticleGroupmodule
      type(ParticleGroup):: pgroup
      real(kind=8):: dy,dz,ymmin,zmmin
      integer(ISZ):: ny,nz

c Sorts particles, using a full size extra particle array for temporary
c space.

      integer(ISZ):: is,ip,iy,iz,pindexmin
      integer(ISZ),allocatable:: pindex(:)
      integer(ISZ):: npblock((1+ny)*(1+nz))

c     --- Treat the particles in each species separately
      do is=1,pgroup%ns
        allocate(pindex(pgroup%nps(is)))

c       --- First, find which block each particle goes into and count how
c       --- many are in each block
        pindex = 0
        do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
          iy = abs(pgroup%yp(ip) - ymmin)/dy
          iz = abs(pgroup%zp(ip) - zmmin)/dz
          pindex(ip) = iy + (ny+1)*iz
        enddo

c       --- Do the sorting
        pindexmin = minval(pindex)
        call particlesortbyindex(pgroup,pindex,pindexmin,
     &                           pgroup%ins(is),pgroup%nps(is),
     &                           (1+ny)*(1+nz),npblock)

        deallocate(pindex)

c     --- End of loop over species
      enddo

      return
      end
c=============================================================================
      subroutine particlesortbyindex(pgroup,pindex,pindexmin,ipmin,np,
     &                               nn,npblock)
      use ParticleGroupmodule
      use Beam_acc, Only: lrelativ
      type(ParticleGroup):: pgroup
      integer(ISZ):: pindexmin,ipmin,np,nn
      integer(ISZ):: npblock(nn)
      integer(ISZ):: pindex(np)

      integer(ISZ):: ip,ipb,ipid
      integer(ISZ):: ipblock(nn)
      real(kind=8),pointer:: pdata(:)
      
      npblock = 0

c     --- First, count how many are in each block
      do ip=1,np
        ipb = pindex(ip) - pindexmin + 1
        npblock(ipb) = npblock(ipb) + 1
      enddo

c     --- Calculate the starting location of each block
      ipblock(1) = 1
      do ipb=2,nn
        ipblock(ipb) = ipblock(ipb-1) + npblock(ipb-1)
      enddo

c     --- Replace pindex with the new array index.
      do ip=1,np
        ipb = pindex(ip) - pindexmin + 1
        pindex(ip) = ipblock(ipb)
        ipblock(ipb) = ipblock(ipb) + 1
      enddo

c     --- Now, the sort can be done.
      allocate(pdata(np))
      call particlecopysorteddata(np,pgroup%xp(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%yp(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%zp(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%uxp(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%uyp(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%uzp(ipmin:ipmin+np-1),pdata,pindex)
      if (lrelativ) then
        call particlecopysorteddata(np,pgroup%gaminv(ipmin:ipmin+np-1),pdata,pindex)
      endif
      call particlecopysorteddata(np,pgroup%ex(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%ey(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%ez(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%bx(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%by(ipmin:ipmin+np-1),pdata,pindex)
      call particlecopysorteddata(np,pgroup%bz(ipmin:ipmin+np-1),pdata,pindex)

      do ipid=1,pgroup%npid
        call particlecopysorteddata(np,pgroup%pid(ipmin:ipmin+np-1,ipid),pdata,
     &                              pindex)
      enddo

      deallocate(pdata)

      return
      end
c=============================================================================
      subroutine particlecopysorteddata(npmax,source,target,pindex)
      integer(ISZ):: npmax
      real(kind=8):: source(npmax),target(npmax)
      integer(ISZ):: pindex(npmax)

      integer(ISZ):: ip,ipi

      do ip=1,npmax
        ipi = pindex(ip)
        if (ipi == 0) cycle
        target(ipi) = source(ip)
      enddo
      do ip=1,npmax
        source(ip) = target(ip)
      enddo

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine alotlostpart
      use Subtimerstop
      use InPart
      use Particles,Only: npid
      use LostParticles

c Allocate space for lost particles and set up associated arrays.

      integer(ISZ):: is
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      npidlost = max(npid,npidlost)
      call gchange("LostParticles",0)

      npmaxlost_s(0) = 0
      do is=1,ns
        npmaxlost_s(is) = npmaxlost_s(is-1) + npmaxlost*sp_fract(is)
        inslost(is) = npmaxlost_s(is-1) + 1
        npslost(is) = 0
      enddo

!$OMP MASTER
      if (ltoptimesubs) timealotlostpart = timealotlostpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine chcklostpart(is,nlower,nhigher)
      use Subtimerstop
      use InPart
      use LostParticles
      use Particles,Only: npid
      integer(ISZ):: is,nlower,nhigher

c Make sure that there is enough space in the lost particle arrays for nlower
c new particles below and nhigher above the lost particles.  Returns if
c there is already enough space above and below.  If there is enough total
c space but not enough room above or below, the lost particles are shifted
c appropriately. If there is not enough space, add more to the arrays.
c Particle data is shifted appropriately.

      integer(ISZ):: nadd,i,ishift
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Make sure that npidlost is up to date.
      if (npid > npidlost) then
        npidlost = npid
        call gchange("LostParticles",0)
      endif

c     --- If there is already enough space, then return.
      if (inslost(is)-npmaxlost_s(is-1)-1 >= nlower .and.
     &    npmaxlost_s(is)-inslost(is)-npslost(is)+1 >= nhigher) return

c     --- If there is not enough space, then allocate new space.
      if (inslost(is)-npmaxlost_s(is-1)-1 +
     &    npmaxlost_s(is)-inslost(is)-npslost(is)+1 < nlower + nhigher) then

c       --- Amount of new space to add: this is calculated to be a sizable
c       --- chunk, but not too big.  By default, add 10 times the space needed.
c       --- Limit that to 10000, except, if the space needed is greater than
c       --- 10000, add the amount of space needed.
        nadd = nlower + nhigher -
     &         (inslost(is)-npmaxlost_s(is-1)-1 +
     &          npmaxlost_s(is)-inslost(is)-npslost(is)+1)
        nadd = max(lostpartchunksize, nadd)

c       --- change the particle array dimension variables and allot space
        npmaxlost = npmaxlost + nadd
        call gchange("LostParticles",0)

c       --- Loop over species species above 'is', shifting
c       --- them up to make space below for species is.
        do i=ns,is+1,-1
          npmaxlost_s(i) = npmaxlost_s(i) + nadd
          call shftlostpart(i,nadd)
        enddo

c       --- Increase size of the block for species 'is'.
        npmaxlost_s(is) = npmaxlost_s(is) + nadd

      endif

c     --- Now, there is enough room to shift the particles appropriately.
      ishift = max(0,nlower - (inslost(is)-npmaxlost_s(is-1)-1)) +
     &        min(0,(npmaxlost_s(is)-inslost(is)-npslost(is)+1) - nhigher)
      call shftlostpart(is,ishift)

!$OMP MASTER
      if (ltoptimesubs) timechcklostpart = timechcklostpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine shftlostpart(is,ishift)
      use Subtimerstop
      use InPart
      use LostParticles
      integer(ISZ):: is,ishift
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (npidlost > 0) then
        call shftpartwork(is,ishift,npmaxlost,ns,npmaxlost_s,
     &                    inslost,npslost,npidlost,
     &                    xplost,yplost,zplost,uxplost,uyplost,uzplost,
     &                    gaminvlost,
     &                    .false.,0.,0.,0.,0.,0.,0.,pidlost)
      else
        call shftpartwork(is,ishift,npmaxlost,ns,npmaxlost_s,
     &                    inslost,npslost,npidlost,
     &                    xplost,yplost,zplost,uxplost,uyplost,uzplost,
     &                    gaminvlost,
     &                    .false.,0.,0.,0.,0.,0.,0.,0.)
      endif

!$OMP MASTER
      if (ltoptimesubs) timeshftlostpart = timeshftlostpart + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine setupSubcycling(pgroup)
      use ParticleGroupmodule
      use Subcycling
      type(ParticleGroup):: pgroup

c Setups in the bookkeepping variables and arrays for particle subcycling.
c These are generic and can be used by any field solver.

      integer(ISZ):: js,i,nsndtstemp,newndts

c     --- If there are no particle species, then do nothing
      if (pgroup%ns == 0) return

c     --- Find maximum value of ndts and use it to size the ndtstorho array
      ndtsmax = max(ndtsmax,maxval(pgroup%ndts))
      nsndts = max(1,nsndts)
      call gchange("Subcycling",0)

c     --- Force the first subcycling group to be the one which advances
c     --- every time step (even if no such group exists).
      ndts(0) = 1
      ndtstorho(1) = 0

c     --- Find any new time step sizes which havn't been included and add
c     --- space to the rhopndts array.
      newndts = 0
      do js = 0, pgroup%ns-1
        if (ndtstorho(pgroup%ndts(js)) == -1) then
          ndtstorho(pgroup%ndts(js)) = -2
          newndts = newndts + 1
        endif
      enddo
      if (newndts > 0) then
        nsndtstemp = nsndts
        nsndts = nsndts + newndts
        if (ndtsmax > 1) nrhopndtscopies = 2
        call gchange("Subcycling",0)
        do js = 0, pgroup%ns-1
          if (ndtstorho(pgroup%ndts(js)) == -2) then
            ndtstorho(pgroup%ndts(js)) = nsndtstemp
            ndts(nsndtstemp) = pgroup%ndts(js)
            nsndtstemp = nsndtstemp + 1
          endif
        enddo
      endif

c     --- Make sure that nsndtsphi is set properly, and call the gchange
c     --- to allocate phipndts.
      if (ndtsaveraging == 0 .or. ndtsaveraging == 1) then
        nsndtsphi = 1
      else
        nsndtsphi = nsndts
      endif
      call gchange("Subcycling",0)

      return
      end
c=============================================================================
      integer(ISZ) function getnsndtsforsubcycling()
      use Subcycling,Only: ndtsaveraging,nsndts

      integer(ISZ):: tmpnsndts

      if (ndtsaveraging == 0 .or. ndtsaveraging == 1) then
        tmpnsndts = 1
      elseif (ndtsaveraging == 2) then
        tmpnsndts = nsndts
      elseif (ndtsaveraging == 3) then
        tmpnsndts = nsndts
      endif

      getnsndtsforsubcycling = tmpnsndts
      return
      end
c=============================================================================
c=============================================================================
      subroutine chgparticlesdts(pgroup)
      use ParticleGroupmodule
      use Particles, Only: chdtspid,spid
      type(ParticleGroup):: pgroup
      
      integer(ISZ):: ins_init(pgroup%ns),nps_init(pgroup%ns),is,ip,idest,ip0,idest0
      integer(ISZ):: ins_save(pgroup%ns),nps_save(pgroup%ns)
      real(kind=8):: pidtemp(pgroup%npid)
      real(kind=8):: temp
      integer(ISZ):: nup,ndown

c     --- first, move particles that have been flagged to the beggining or end
c     --- of the particle arrays.
      do is=1,pgroup%ns

c       --- initialize counter and save particle's start and number
        ip0 = pgroup%ins(is)
        ins_init(is) = pgroup%ins(is)
        nps_init(is) = pgroup%nps(is)

c       --- Loop over particles and pick out the ones which have been flagged
c       --- to be moved to a different ndts group.
        do while (ip0 < pgroup%ins(is) + pgroup%nps(is))

c         --- If particle is going to a faster group, switch with lowest particle
c         --- in the arrays.
          if (pgroup%pid(ip0,chdtspid) < 0.) then

            ip = ip0
            idest = pgroup%ins(is)
            pgroup%ins(is) = pgroup%ins(is) + 1
            pgroup%nps(is) = pgroup%nps(is) - 1

          elseif (pgroup%pid(ip0,chdtspid) > 0.) then

            ip = ip0
            idest = pgroup%ins(is) + pgroup%nps(is) - 1
            pgroup%nps(is) = pgroup%nps(is) - 1
c           --- Note that the particle swapped in stills need to be checked.
            ip0 = ip0 - 1

          else
            ip = -1
          endif

          if (ip > 0) then

            temp = pgroup%xp(ip); pgroup%xp(ip) = pgroup%xp(idest); pgroup%xp(idest) = temp
            temp = pgroup%yp(ip); pgroup%yp(ip) = pgroup%yp(idest); pgroup%yp(idest) = temp
            temp = pgroup%zp(ip); pgroup%zp(ip) = pgroup%zp(idest); pgroup%zp(idest) = temp

            temp = pgroup%uxp(ip); pgroup%uxp(ip) = pgroup%uxp(idest); pgroup%uxp(idest) = temp
            temp = pgroup%uyp(ip); pgroup%uyp(ip) = pgroup%uyp(idest); pgroup%uyp(idest) = temp
            temp = pgroup%uzp(ip); pgroup%uzp(ip) = pgroup%uzp(idest); pgroup%uzp(idest) = temp

            temp = pgroup%gaminv(ip); pgroup%gaminv(ip) = pgroup%gaminv(idest); pgroup%gaminv(idest) = temp

            temp = pgroup%ex(ip); pgroup%ex(ip) = pgroup%ex(idest); pgroup%ex(idest) = temp
            temp = pgroup%ey(ip); pgroup%ey(ip) = pgroup%ey(idest); pgroup%ey(idest) = temp
            temp = pgroup%ez(ip); pgroup%ez(ip) = pgroup%ez(idest); pgroup%ez(idest) = temp

            temp = pgroup%bx(ip); pgroup%bx(ip) = pgroup%bx(idest); pgroup%bx(idest) = temp
            temp = pgroup%by(ip); pgroup%by(ip) = pgroup%by(idest); pgroup%by(idest) = temp
            temp = pgroup%bz(ip); pgroup%bz(ip) = pgroup%bz(idest); pgroup%bz(idest) = temp

            pgroup%pid(ip,chdtspid) = 0.
            pidtemp = pgroup%pid(ip,:); pgroup%pid(ip,:) = pgroup%pid(idest,:); pgroup%pid(idest,:) = pidtemp

          endif

c         --- advance ip0
          ip0 = ip0 + 1

        enddo

c       --- end of loop over species
      enddo

c     --- Now, swap the particle between species.
      do is=1,pgroup%ns - 1
        ndown = pgroup%ins(is+1) - ins_init(is+1)
        nup = (ins_init(is) + nps_init(is)) - (pgroup%ins(is) + pgroup%nps(is))
        do ip0 = 0,min(ndown,nup)-1
          ip = pgroup%ins(is) + pgroup%nps(is) + ip0
          idest = pgroup%ins(is+1) - 1 - ip0

          temp = pgroup%xp(ip); pgroup%xp(ip) = pgroup%xp(idest); pgroup%xp(idest) = temp
          temp = pgroup%yp(ip); pgroup%yp(ip) = pgroup%yp(idest); pgroup%yp(idest) = temp
          temp = pgroup%zp(ip); pgroup%zp(ip) = pgroup%zp(idest); pgroup%zp(idest) = temp

          temp = pgroup%uxp(ip); pgroup%uxp(ip) = pgroup%uxp(idest); pgroup%uxp(idest) = temp
          temp = pgroup%uyp(ip); pgroup%uyp(ip) = pgroup%uyp(idest); pgroup%uyp(idest) = temp
          temp = pgroup%uzp(ip); pgroup%uzp(ip) = pgroup%uzp(idest); pgroup%uzp(idest) = temp

          temp = pgroup%gaminv(ip); pgroup%gaminv(ip) = pgroup%gaminv(idest); pgroup%gaminv(idest) = temp

          temp = pgroup%ex(ip); pgroup%ex(ip) = pgroup%ex(idest); pgroup%ex(idest) = temp
          temp = pgroup%ey(ip); pgroup%ey(ip) = pgroup%ey(idest); pgroup%ey(idest) = temp
          temp = pgroup%ez(ip); pgroup%ez(ip) = pgroup%ez(idest); pgroup%ez(idest) = temp

          temp = pgroup%bx(ip); pgroup%bx(ip) = pgroup%bx(idest); pgroup%bx(idest) = temp
          temp = pgroup%by(ip); pgroup%by(ip) = pgroup%by(idest); pgroup%by(idest) = temp
          temp = pgroup%bz(ip); pgroup%bz(ip) = pgroup%bz(idest); pgroup%bz(idest) = temp

          pidtemp = pgroup%pid(ip,:); pgroup%pid(ip,:) = pgroup%pid(idest,:); pgroup%pid(idest,:) = pidtemp

        enddo
        pgroup%nps(is) = pgroup%nps(is) + min(ndown,nup)
        pgroup%ins(is+1) = pgroup%ins(is+1) - min(ndown,nup)
        pgroup%nps(is+1) = pgroup%nps(is+1) + min(ndown,nup)

c       --- Move any particles that are left over
        if (ndown > nup) then
          idest0 = pgroup%ins(is) + pgroup%nps(is)
          idest = idest0 + (ndown - nup) - 1
          ip0 = pgroup%ins(is+1) - (ndown - nup)
          ip = pgroup%ins(is+1) - 1
          pgroup%nps(is) = pgroup%nps(is) + (ndown - nup)
        else if (nup > ndown) then
          idest0 = pgroup%ins(is+1) - (nup - ndown)
          idest = pgroup%ins(is+1) - 1
          ip0 = pgroup%ins(is) + pgroup%nps(is)
          ip = ip0 + (nup - ndown) - 1
          pgroup%ins(is+1) = pgroup%ins(is+1) - (nup - ndown)
          pgroup%nps(is+1) = pgroup%nps(is+1) + (nup - ndown)
        endif
        if (nup .ne. ndown) then
          pgroup%xp(idest0:idest) = pgroup%xp(ip0:ip)
          pgroup%yp(idest0:idest) = pgroup%yp(ip0:ip)
          pgroup%zp(idest0:idest) = pgroup%zp(ip0:ip)
  
          pgroup%uxp(idest0:idest) = pgroup%uxp(ip0:ip)
          pgroup%uyp(idest0:idest) = pgroup%uyp(ip0:ip)
          pgroup%uzp(idest0:idest) = pgroup%uzp(ip0:ip)
  
          pgroup%gaminv(idest0:idest) = pgroup%gaminv(ip0:ip)
  
          pgroup%ex(idest0:idest) = pgroup%ex(ip0:ip)
          pgroup%ey(idest0:idest) = pgroup%ey(ip0:ip)
          pgroup%ez(idest0:idest) = pgroup%ez(ip0:ip)
  
          pgroup%bx(idest0:idest) = pgroup%bx(ip0:ip)
          pgroup%by(idest0:idest) = pgroup%by(ip0:ip)
          pgroup%bz(idest0:idest) = pgroup%bz(ip0:ip)
  
          pgroup%pid(idest0:idest,:) = pgroup%pid(ip0:ip,:)
          pgroup%pid(ip0:ip,:) = 0.
        endif

        if ((pgroup%nps(is) + pgroup%nps(is+1)) > 0) then
c         --- This puts ipmax_s is a nice place between the species that allows
c         --- each to grew without shifting particles around.
c         --- Though ipmax_s will most likely be removed since it is not really
c         --- needed.
          pgroup%ipmax_s(is) = (pgroup%nps(is)*(pgroup%ins(is) + pgroup%nps(is) - 1) +
     &                          pgroup%nps(is+1)*(pgroup%ins(is+1) - 1))/
     &                         (pgroup%nps(is) + pgroup%nps(is+1))
        else
          pgroup%ipmax_s(is) = (pgroup%ins(is) + pgroup%nps(is) - 1 +
     &                          pgroup%ins(is+1) - 1)/2
        endif

      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine setupSelfB(pgroup)
      use ParticleGroupmodule
      use SelfB
      type(ParticleGroup):: pgroup

c Setups in the bookkeepping variables and arrays for particle which require
c the self B correction.
c These are generic and can be used by any field solver.

      integer(ISZ):: js1,js2,jsid

      logical(ISZ):: lnew

c     --- Find any new species with different fselfb.
      do js1=0,pgroup%ns-1
        jsid = pgroup%sid(js1)
        if (jsid == -1) cycle
        lnew = .true.
        do js2=0,nsselfb-1
          if (pgroup%fselfb(js1) == fselfb(js2)) then
            lnew = .false.
            pgroup%iselfb(js1) = js2
            iselfb(jsid) = js2
          endif
        enddo
        if (lnew) then
          nsselfb = nsselfb + 1
          call gchange("SelfB",0)
          pgroup%iselfb(js1) = nsselfb-1
          iselfb(jsid) = nsselfb-1
          fselfb(nsselfb-1) = pgroup%fselfb(js1)
        endif
      enddo

      return
      end
c=============================================================================
c Routines for handling particle groups
c=============================================================================
      subroutine checkparticlegroup(pgroup,is,nlower,nhigher)
      use Subtimerstop
      use ParticleGroupmodule
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,nlower,nhigher

c Make sure that there is enough space in the particle arrays for nlower
c new particles below and nhigher above the live particles.  Returns if
c there is already enough space above and below.  If there is enough total
c space but not enough room above or below, the particles are shifted
c appropriately.  If lfullshft is true, the particles are shifted all the
c way to the end of the species pgroup minus nhigher.  If there is not
c enough space, add more to the arrays. Particle data is shifted
c appropriately.

      integer(ISZ):: ilower,ihigher,spacebelow,spaceabove
      integer(ISZ):: naddbelow,naddabove,i,ishift
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- If there is already enough space, then return.
      if (is == 1) then
        ilower = 0
      else
        ilower = pgroup%ins(is-1) + pgroup%nps(is-1) - 1
      endif
      if (is == pgroup%ns) then
        ihigher = pgroup%npmax + 1
      else
        ihigher = pgroup%ins(is+1)
      endif
      spacebelow = pgroup%ins(is) - ilower - 1
      spaceabove = ihigher - pgroup%ins(is) - pgroup%nps(is)
      if (spacebelow >= nlower .and.
     &    spaceabove >= nhigher) return

c     --- If there is not enough space, then allocate new space.
      if (spacebelow + spaceabove < nlower + nhigher) then

c       --- Amount of new space to add: this is calculated to be a sizable
c       --- chunk, but not too big.  By default, add 10 times the space needed.
c       --- Limit that to 10000, except, if the space needed is greater than
c       --- 10000, add the amount of space needed.
        naddbelow = max(0,nlower - spacebelow)
        naddbelow = min(10*naddbelow, max(10000, naddbelow))
        naddabove = max(0,nhigher - spaceabove)
        naddabove = min(10*naddabove, max(10000, naddabove))

c       --- change the particle array dimension variables and allot space
        pgroup%npmax = pgroup%npmax + naddbelow + naddabove
        call ParticleGroupchange(pgroup)

c       --- Loop over species species above 'is', shifting
c       --- them up to make space below for species is.
        do i=pgroup%ns,is+1,-1
          call shiftparticlegroup(pgroup,i,naddbelow + naddabove)
        enddo

      endif

c     --- Now, there is enough room to shift the particles appropriately.
c     --- Note that in some cases no shift will be needed.  In those cases,
c     --- ishift comes out to be zero and shftpart returns immediately when
c     --- ishift is zero.
      if (is == 1) then
        ilower = 0
      else
        ilower = pgroup%ins(is-1) + pgroup%nps(is-1) - 1
      endif
      if (is == pgroup%ns) then
        ihigher = pgroup%npmax + 1
      else
        ihigher = pgroup%ins(is+1)
      endif
      spacebelow = pgroup%ins(is) - ilower - 1
      spaceabove = ihigher - pgroup%ins(is) - pgroup%nps(is)
      ishift = max(0,nlower - spacebelow) + min(0,spaceabove - nhigher)
      call shiftparticlegroup(pgroup,is,ishift)

!$OMP MASTER
      if (ltoptimesubs) timecheckparticlegroup = timecheckparticlegroup + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine shiftparticlegroup(pgroup,is,ishift)
      use Subtimerstop
      use ParticleGroupmodule
      type(ParticleGroup):: pgroup
      integer(ISZ):: is,ishift

c Shift particles by 'ishift'.
c Assumes that the is enough space in the particle arrays to make the
c shift without clobbering other particles or going past the ends of
c the arrays.

      integer(ISZ):: ip,ilower
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- return if shift is zero
      if (ishift == 0) return

c     --- Positive and negative shifts are treated differently since
c     --- a negative shift may result in lost data.
      if (ishift > 0) then

c       --- Loop only over live particles
        do ip=pgroup%ins(is)+pgroup%nps(is)-1, pgroup%ins(is), -1
          pgroup%xp(ip+ishift) = pgroup%xp(ip)
          pgroup%yp(ip+ishift) = pgroup%yp(ip)
          pgroup%zp(ip+ishift) = pgroup%zp(ip)
          pgroup%uxp(ip+ishift) = pgroup%uxp(ip)
          pgroup%uyp(ip+ishift) = pgroup%uyp(ip)
          pgroup%uzp(ip+ishift) = pgroup%uzp(ip)
          pgroup%gaminv(ip+ishift) = pgroup%gaminv(ip)
          pgroup%ex(ip+ishift) = pgroup%ex(ip)
          pgroup%ey(ip+ishift) = pgroup%ey(ip)
          pgroup%ez(ip+ishift) = pgroup%ez(ip)
          pgroup%bx(ip+ishift) = pgroup%bx(ip)
          pgroup%by(ip+ishift) = pgroup%by(ip)
          pgroup%bz(ip+ishift) = pgroup%bz(ip)
          if (pgroup%npid > 0) then
            pgroup%pid(ip+ishift,:) = pgroup%pid(ip,:)
c           pgroup%pid(ip,:) = 0
          endif
c         pgroup%uzp(ip) = 0.
        enddo

c       --- change lower bound of live particles
        pgroup%ins(is) = pgroup%ins(is) + ishift

      else
c       --- Loop only over live particles.
c       --- Note that particles will be clobbered if shifted onto species below.
        if (is == 1) then
          ilower = 0
        else
          ilower = pgroup%ins(is-1) + pgroup%nps(is-1) - 1
        endif
        if (pgroup%ins(is) - ilower - 1 < -ishift) then
          print*,"WARNING: shiftparticlegroup: particle data will be lost"
          print*,"         since there is  not enough space in the arrays."
        endif

        do ip=pgroup%ins(is), pgroup%ins(is)+pgroup%nps(is)-1
          pgroup%xp(ip+ishift) = pgroup%xp(ip)
          pgroup%yp(ip+ishift) = pgroup%yp(ip)
          pgroup%zp(ip+ishift) = pgroup%zp(ip)
          pgroup%uxp(ip+ishift) = pgroup%uxp(ip)
          pgroup%uyp(ip+ishift) = pgroup%uyp(ip)
          pgroup%uzp(ip+ishift) = pgroup%uzp(ip)
          pgroup%gaminv(ip+ishift) = pgroup%gaminv(ip)
          pgroup%ex(ip+ishift) = pgroup%ex(ip)
          pgroup%ey(ip+ishift) = pgroup%ey(ip)
          pgroup%ez(ip+ishift) = pgroup%ez(ip)
          pgroup%bx(ip+ishift) = pgroup%bx(ip)
          pgroup%by(ip+ishift) = pgroup%by(ip)
          pgroup%bz(ip+ishift) = pgroup%bz(ip)
          if (pgroup%npid > 0) then
            pgroup%pid(ip+ishift,:) = pgroup%pid(ip,:)
c           pgroup%pid(ip,:) = 0
          endif
c         pgroup%uzp(ip) = 0.
        enddo

c       --- change number of particles if data is lost
        if (pgroup%ins(is) - ilower - 1 < -ishift) then
          pgroup%nps(is) = pgroup%nps(is) - (-ishift - (pgroup%ins(is) - ilower - 1))
          pgroup%ins(is) = ilower + 1
        else
c         --- change lower bound only
          pgroup%ins(is) = pgroup%ins(is) + ishift
        endif

      endif

!$OMP MASTER
      if (ltoptimesubs) timeshiftparticlegroup = timeshiftparticlegroup + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine copygrouptogroup(pgroupin,nn,ii,istart,pgroupout,it)
      use ParticleGroupmodule
      use Subtimerstop
      use InPart
      integer(ISZ):: it,nn,istart
      integer(ISZ):: ii(0:nn-1)
      type(ParticleGroup):: pgroupin,pgroupout

c Copies particle data from locations given by ii to the locations starting
c at 'it', or if istart > 0, copy sequential particles from locations
c starting at 'istart' to locations starting at 'it'.

      integer(ISZ):: i,in,ij
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call setuppgroup(pgroupin)
      call setuppgroup(pgroupout)

      if (istart > 0) then
        do i=0,nn-1
          ij = istart + i
          in = it + i
          pgroupout%xp(in) = pgroupin%xp(ij)
          pgroupout%yp(in) = pgroupin%yp(ij)
          pgroupout%zp(in) = pgroupin%zp(ij)
          pgroupout%uxp(in) = pgroupin%uxp(ij)
          pgroupout%uyp(in) = pgroupin%uyp(ij)
          pgroupout%uzp(in) = pgroupin%uzp(ij)
          pgroupout%gaminv(in) = pgroupin%gaminv(ij)
          pgroupout%ex(in) = pgroupin%ex(ij)
          pgroupout%ey(in) = pgroupin%ey(ij)
          pgroupout%ez(in) = pgroupin%ez(ij)
          pgroupout%bx(in) = pgroupin%bx(ij)
          pgroupout%by(in) = pgroupin%by(ij)
          pgroupout%bz(in) = pgroupin%bz(ij)
          if (pgroupin%npid > 0) then
            pgroupout%pid(in,:) = pgroupin%pid(ij,:)
          endif
        enddo
      else
        do i=0,nn-1
          ij = ii(i)
          in = it + i
          pgroupout%xp(in) = pgroupin%xp(ij)
          pgroupout%yp(in) = pgroupin%yp(ij)
          pgroupout%zp(in) = pgroupin%zp(ij)
          pgroupout%uxp(in) = pgroupin%uxp(ij)
          pgroupout%uyp(in) = pgroupin%uyp(ij)
          pgroupout%uzp(in) = pgroupin%uzp(ij)
          pgroupout%gaminv(in) = pgroupin%gaminv(ij)
          pgroupout%ex(in) = pgroupin%ex(ij)
          pgroupout%ey(in) = pgroupin%ey(ij)
          pgroupout%ez(in) = pgroupin%ez(ij)
          pgroupout%bx(in) = pgroupin%bx(ij)
          pgroupout%by(in) = pgroupin%by(ij)
          pgroupout%bz(in) = pgroupin%bz(ij)
          if (pgroupin%npid > 0) then
            pgroupout%pid(in,:) = pgroupin%pid(ij,:)
          endif
        enddo
      endif

!$OMP MASTER
      if (ltoptimesubs) timecopyparttogroup = timecopyparttogroup + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c  VARIOUS UTILITIES
c=============================================================================
      integer(ISZ) function nextpid()
      use Particles,Only: npid

c     --- Returns the next value of pid and increments npid.
c     --- Note that this should be used and that npid should not be directly
c     --- changed.

      npid = npid + 1
      nextpid = npid

      return
      end
c=============================================================================
      subroutine copyarry (source, target, nwords)
      integer(ISZ):: nwords
      real(kind=8):: source(nwords), target(nwords)
      integer(ISZ):: i

!$OMP PARALLEL DO IF (nwords > 1000)
      do 100 i = 1, nwords
         target(i) = source(i)
  100 continue
!$OMP END PARALLEL DO

      return
      end
c=============================================================================
      subroutine sumarry (source, target, nwords)
      integer(ISZ):: nwords
      real(kind=8):: source(nwords), target(nwords)
      integer(ISZ):: i

!$OMP PARALLEL DO IF (nwords > 1000)
      do 100 i = 1, nwords
         target(i) = target(i) + source(i)
  100 continue
!$OMP END PARALLEL DO

      return
      end
c=============================================================================
      subroutine zeroarry (arry,n)
      integer(ISZ):: n
      real(kind=8):: arry(*)

c  Sets ARRY to zero
      integer(ISZ):: i
!$OMP PARALLEL DO IF (n > 1000)
      do i=1,n
        arry(i) = 0.0
      enddo
!$OMP END PARALLEL DO

      return
      end
c=============================================================================
      subroutine writarry(nn,arry,filename)
      integer(ISZ):: nn
      real(kind=8):: arry(nn)
      character(*):: filename
      integer(ISZ):: ii

      open(unit=27,file=filename,status='new')

      do ii=1,nn
        write(27,20) arry(ii)
      enddo
 20   format(1x,e10.3)

      close(27)

      return
      end
c=============================================================================
      real(kind=8) function fnice(i,e10)
      integer(ISZ):: i
      real(kind=8):: e10
c this is used to pick nice numbers; also gives next larger power of 10.
c for: i =  -4 -3 -2 -1 0  1  2  3   4   5
c  fnice = .05 .1 .2 .5 1  2  5 10  20  50
c    e10 =  .1 .1  1  1 1 10 10 10 100 100
c     i3 =   2  3  1  2 3  1  2  3   1   2
c    i10 =  -1 -1  0  0 0  1  1  1   2   2

      integer(ISZ):: i10,i3
      real(kind=8):: a(3) = (/.2,.5,1./)
      i10=i/3
      if(i > 0) i10=(i+2)/3
      i3=i+3-3*i10
      e10=10.**i10
      fnice=a(i3)*e10
      return
      end
c=============================================================================
      subroutine psumx(a,b,n)
      integer(ISZ):: n
      real(kind=8):: a(n), b(n)
c     b := partial sum of a
      integer(ISZ):: i
      b(1) = a(1)
      do i = 2, n
        b(i) = a(i) + b(i-1)
      enddo
      return
      end
c=============================================================================
      logical function inlist (list, listlen, text)
      integer(ISZ):: listlen
      character(8):: text, list(listlen)
      integer(ISZ):: i

      inlist = .false.
      do i = 1,listlen
         if (text == list(i)) inlist = .true.
      enddo

      return
      end
c=======================================================================
      real(kind=8) function wtime()
      real(kind=8):: cpu,io,sys,mem
c     real(kind=8):: walltime

c     --- Get current time
      call ostime(cpu, io, sys, mem)
      wtime = cpu
c     wtime = walltime()

      return
      end
c=======================================================================
      subroutine wtimeon()
      real(kind=8):: cpu,io,sys,mem
      common /wtimer/cpu,io,sys,mem

c     --- turns timer on
      call ostime(cpu, io, sys, mem)

      return
      end
c=======================================================================
      real(kind=8) function wtimeoff()
      real(kind=8):: cpu,io,sys,mem
      real(kind=8):: fcpu,fio,fsys,fmem
      common /wtimer/cpu,io,sys,mem

c     --- returns elapsed time in milliseconds since last wtimeon call
      call ostime(fcpu, fio, fsys, fmem)
      wtimeoff = (fcpu - cpu)*1.e3

      return
      end
c=======================================================================
      real(kind=8) function wtremain()
c Returns the time remaining in the job. Only works on the T3E. Otherwise
c returns a large number.
      real(kind=8):: py_tremain
      wtremain = py_tremain()
      return
      end
c=============================================================================
c  RANDOM NUMBER , QUIET-START, and RELATED ROUTINES
c=============================================================================
      real(kind=8) function wrandom(method,ii,idig,ifib1,ifib2)
      character(*):: method
      integer(ISZ):: ii,idig,ifib1,ifib2
      real(kind=8):: rnrev,wranf
      character(72):: errline

      if (method == "digitrev") then
        wrandom = rnrev(ii,idig)
      elseif (method == "pseudo") then
        wrandom = wranf()
      elseif (method == "fibonacc") then
        wrandom = mod((ifib1*(ii - 1.) + 0.5)/ifib2,1.0)
      else
        wrandom = 0.
        write(errline,'("ERROR: random number generator method ",a8,
     &                  "not supported")') method
        call remark (errline)
        call kaboom (0)
      endif

      return
      end
c=============================================================================
      real(kind=8) function wrandomgauss(method,ii,idig1,idig2,ifib1,ifib2,
     &                                   lefficient)
      use Constant
      character(*):: method
      integer(ISZ):: ii,idig1,idig2,ifib1,ifib2
      logical(ISZ):: lefficient

c  Gaussian random number generator: returns a value drawn from a normal
c    distribution with mean 0 and variance 1.
c  Uses Box-Muller-Marsaglia method:
c  See Knuth, vol. 2, p.104 (first ed.), p.117 (second ed.)
c  Alex Friedman, July 1989 (after an earlier version).

      real(kind=8):: wrandom
      real(kind=8):: y,v1,v2,sq,x,wranf
      integer(ISZ):: parity,i1
      integer(ISZ):: iflag = -1
c     --- The save is necessary to preserve state memory
      save y, iflag

c     --- The algorithm used calculated two random numbers at a time. If it is
c     --- not needed, there is no reason to repeat the calculation. In this
c     --- case, the extra number calculated is saved and returned on the next
c     --- call. This happens when lefficient is true - which requires that
c     --- this routine be called in the correct order. Otherwise, the work
c     --- is repeated each time and the appropriate of the two numbers is
c     --- returned.

      parity = mod(ii,2)
      if (lefficient) then

c       --- In order to have a consistent stream of random numbers, the x random
c       --- number is always taken from ii is odd, and the y when ii is even.
c       --- This is really only needed for the parallel version.
c       --- If iflag is inconsistent with the parity of ii, then adjust iflag.
c       --- If odd ii is expected but ii is even, then set calculation to use
c       --- ii-1 so that number used for v1 below is odd.

        iflag = -iflag
        i1 = ii
        if (iflag < 0 .and. parity == 0) then
          wrandomgauss = y
          return
        elseif (iflag < 0 .and. parity == 1) then
          iflag = -iflag
        elseif (iflag > 0 .and. parity == 0) then
          iflag = -iflag
          i1 = ii - 1
        endif

      else
        if (parity == 0) then
          i1 = ii - 1
        else
          i1 = ii
        endif
      endif

      v1 = wrandom(method,i1,idig1,ifib1,ifib2)
      v2 = 2.*pi*wrandom(method,i1+1,idig2,ifib1,ifib2)
      sq = sqrt(-2.*log(v1))
      x = sq*cos(v2)
      y = sq*sin(v2)
      if (parity == 1) then
        wrandomgauss = x
      else
        wrandomgauss = y
      endif

      return
      end
c=============================================================================
      real(kind=8) function r2rev(xs)
      real(kind=8):: xs
c  xs must be initialized before use (usually to 0)
c  and preserved between calls.
c  value is given by bit-reversed counter -compare sorter in cpft.
c  xs=.000,.100,.010,.110,.001,.101,.011,.111,.0001.. (binary fractions)
c  is the value of r2rev on the NEXT call. 
      real(kind=8):: xsi
      r2rev = xs
      xsi = 1.
   42  xsi = .5 * xsi
       xs = xs - xsi
      if (xs >= 0.) go to 42
      xs = xs + 2. * xsi
      return
      end
c=============================================================================
      real(kind=8) function rnrev(i,nbase)
      integer(ISZ):: i,nbase
c  Returns a fraction (0 <= rnrev <= 1) representing the
c  nbase-reversed number corresponding to i; i.e.
c  for i written in base nbase as n4n3n2n1.0, this returns
c  rnrev = 0.n1n2n3n4, where n1, n2.. are the digits in base n.
c  Disregards roundoff errors for nbase not a power of 2.
c  From N. Otani, 1985
      integer(ISZ):: j,n,jnext
      real(kind=8):: powern
      rnrev = 0.
      j = i
      n = nbase
      powern = 1.
   10 jnext = j / n
      powern = powern / n
      rnrev = rnrev + (j-jnext*n) * powern
      if (jnext == 0) return
      j = jnext
      go to 10
      end
c=============================================================================
      subroutine rnrevarray(n,x,i,nbase)
      integer(ISZ):: n,i,nbase
      real(kind=8):: x(n)
c Fills an array with uniformly distributed digit reversed numbers.
c The numbers range from 0 to 1.
      real(kind=8):: rnrev
      integer(ISZ):: j
      do j=1,n
        x(j) = rnrev(i+j-1,nbase)
      enddo
      return
      end
c=============================================================================
      subroutine sphere4(a,b,c,d,n)
      use Constant
      integer(ISZ):: n
      real(kind=8):: a(n), b(n), c(n), d(n)
c   Distribute points on surface of a 4-d unit sphere.
c   (a*a + b*b) increases, and (c*c + d*d) decreases, monotonically.
      real(kind=8):: xsrev,rr,phi,r2rev,rnrev
      integer(ISZ):: i

      xsrev = 0.                ! get same sequence if call again
      do i = 1,n
        rr = sqrt((i-.5)/n)
        phi = 2.*pi*r2rev(xsrev)
        a(i) = rr*cos(phi)
        b(i) = rr*sin(phi)
        phi = 2.*pi*rnrev(i-1,3)
        rr = sqrt(1.-rr*rr)
        c(i) = rr*cos(phi)
        d(i) = rr*sin(phi)
      enddo
      return
      end
c=============================================================================
      subroutine sphere4f(a,b,c,d,g1,g2,g3)
      use Constant
      integer(ISZ):: g1,g2,g3
      real(kind=8):: a(g1), b(g1), c(g1), d(g1)
c   Distribute points on surface of a 4-d unit sphere using fibonacci numbers.
c   (a*a + b*b) increases, and (c*c + d*d) decreases, monotonically.
      integer(ISZ):: i
      real(kind=8):: rr,phi

      do i = 1,g1
         rr = sqrt((i-.5)/g1)
         phi = 2.*pi*mod((g2*(i-1)+0.5)/g1, 1.0)
         a(i) = rr*cos(phi)
         b(i) = rr*sin(phi)
         phi = 2.*pi*mod((g3*(i-1)+0.5)/g1, 1.0)
         rr = sqrt(1.-rr*rr)
         c(i) = rr*cos(phi)
         d(i) = rr*sin(phi)
      enddo
      return
      end
c=============================================================================
      real(kind=8) function rnorm ()
      use Constant

c  Gaussian random number generator: returns a value drawn from a normal
c    distribution with mean 0 and variance 1.
c  Uses Box-Muller-Marsaglia method:
c  See Knuth, vol. 2, p.104 (first ed.), p.117 (second ed.)
c  Alex Friedman, July 1989 (after an earlier version).

c     --- This common block necessary to preserve state memory:
      real(kind=8):: y,v1,v2,sq,x,wranf
c     real(kind=8):: s
      integer(ISZ):: iflag = -1
      save y, iflag

      iflag = -iflag
      if (iflag < 0) then
         rnorm = y
         return
      endif

c 100 v1 = 2. * wranf() - 1.
c     v2 = 2. * wranf() - 1.
c     s = v1**2 + v2**2 + SMALLPOS
c     if (s > 1.) go to 100
c     sq = sqrt(-2. * log(s) / s)
c     x = v1 * sq
c     y = v2 * sq
c     rnorm = x

      v1 = wranf()
      v2 = 2.*pi*wranf()
      sq = sqrt(-2.*log(v1))
      x = sq*cos(v2)
      y = sq*sin(v2)
      rnorm = x

      return
      end
c=============================================================================
      subroutine rnormdig (i1,n,nbase1,nbase2,dx,x)
      use Constant
      integer(ISZ):: i1,n,nbase1,nbase2
      real(kind=8):: dx
      real(kind=8):: x(n)

c  Gaussian random number generator: returns a value drawn from a normal
c    distribution with mean 0 and variance 1.
c    chosen using digit reversed type random numbers
c  Uses Box-Muller-Marsaglia method:
c  See Knuth, vol. 2, p.104 (first ed.), p.117 (second ed.)
c  Alex Friedman, July 1989 (after an earlier version).
c  David P. Grote, Febuary 1990 (added digit reversed stuff)
      integer(ISZ):: i,iparity
      real(kind=8):: s,phi,sq,rnrev

c     --- Get first number separately if parity of i1 is even. Since i1
c     --- always starts with 1, the first number is always odd to keep
c     --- consistency. This primarily effects the parallel version.
c     --- The parity of i1+1 is used for convenience.
      iparity = mod(i1+1,2)
      if (iparity == 1) then
        s = rnrev(i1-1,nbase1) + dx
        phi = 2.*pi*rnrev(i1,nbase2)
        sq = sqrt(-2.*log(s))
        x(1) = sq*sin(phi)
      endif

      do i=i1+iparity,i1+n-2,2
        s = rnrev(i,nbase1) + dx
        phi = 2.*pi*rnrev(i+1,nbase2)
        sq = sqrt(-2.*log(s))
        x(i-i1+1) = sq*cos(phi)
        x(i-i1+2) = sq*sin(phi)
      enddo

c     --- Get last number if (n-iparity) is odd.
      if (mod(n-iparity,2) == 1) then
        s = rnrev(i1+n-1,nbase1) + dx
        phi = 2.*pi*rnrev(i1+n,nbase2)
        sq = sqrt(-2.*log(s))
        x(n) = sq*cos(phi)
      endif

      return
      end
c=============================================================================
      real(kind=8) function rm()

c  Crude approximation to a Gaussian with mean 0, standard deviation 1
c  Advantage is, it cuts off smoothly at 3 root 2. 
c  A similar one (used in some older AF codes) adds 12 variates, subtracts 6.

      real(kind=8):: wranf
      real(kind=8):: root2 = 1.41421356237
      rm = root2*(wranf()+wranf()+wranf()+wranf()+wranf()+wranf()-3.)
      return
      end
c=============================================================================
      subroutine rma(a,n)
      integer(ISZ):: n
      real(kind=8):: a(n)
c   From interpreter, rma(&a,n) returns n Gaussian random numbers.
      integer(ISZ):: i
      real(kind=8):: rm

      do i = 1,n
        a(i) = rm()
      enddo
      return
      end
c=============================================================================
c=============================================================================
      subroutine load2d(np,x,y,nx,ny,n,dx,dy)
      integer(ISZ):: np,nx,ny
      real(kind=8):: x(np),y(np),n(0:nx,0:ny)
      real(kind=8):: dx,dy

c Load particles onto a two dimensional distribution.  This load is only
c approximate in that when the resulting particles are loaded back onto
c the same 2-D grid from which the distribution is taken, the distribution
c is not exactly reproduced. (The exact version is still under development.)

      real(kind=8):: dyi,ni,nip1,ninti,sumny,sumnx,wy
      real(kind=8):: rnrev
      integer     :: i,ip,ix,iy
      real(kind=8):: n1x(0:nx),n1y(0:ny)
      real(kind=8):: nintx(0:nx),ninty(0:ny)

      dyi = 1./dy

c     --- Initialize xp and yp
      do i=1,np
c       x(i) = .5/np + wranf()*(np-1.)/np
        x(i) = rnrev(i,2)
        y(i) = .5/np + (i-1.)/np
      enddo

c     --- n1y = sum(n,1)
      do iy=0,ny
        n1y(iy) = sum(n(:,iy))
      enddo
      sumny = sum(n1y) - 0.5*n1y(0) - 0.5*n1y(ny)
      n1y = n1y/sumny

c     --- ninty(0:ny) = psum(n1y) - 0.5*n1y(0) - 0.5*n1y
      ninty(0) = 0.
      do iy=1,ny
        ninty(iy) = ninty(iy-1) + 0.5*(n1y(iy-1) + n1y(iy))
      enddo

c     --- Now distribute the particles in y
c     --- When the distribution has a nearly uniform section, a modified
c     --- expression is required.
      do ip=1,np
        i=0
        do while (y(ip) > ninty(i+1))
          i=i+1
        end do
        if (abs(n1y(i+1)-n1y(i))/(n1y(i+1)+n1y(i)) > 1.e-4) then
          y(ip)=((sqrt(n1y(i)**2 - 2.*(n1y(i+1)-n1y(i))*(ninty(i) - y(ip))) -
     &           n1y(i))/(n1y(i+1)-n1y(i)) + i)*dy
        else
          y(ip) = ((y(ip) - ninty(i))/n1y(i)  + i)*dy
        endif
      enddo

      do ip=1,np
        iy = y(ip)*dyi
        wy = y(ip)*dyi - iy

c       --- n1x = (n(,iy)*(1.-wy)+n(,iy+1)*wy)/sum(n(,iy)*(1.-wy)+n(,iy+1)*wy)
        n1x = n(:,iy)*(1. - wy) + n(:,iy+1)*wy
        sumnx = sum(n1x) - 0.5*n1x(0) - 0.5*n1x(nx)
        n1x = n1x/sumnx

c       --- nintx = psum(n1x) - 0.5*n1x
        nintx(0) = 0.
        do ix=1,nx
          nintx(ix) = nintx(ix-1) + 0.5*(n1x(ix-1) + n1x(ix))
        enddo

        i=0
        do while (x(ip) > nintx(i+1))
          i=i+1
        end do
        ni = n1x(i)
        nip1 = n1x(i+1)
        ninti = nintx(i)
        if (abs(nip1-ni)/(nip1+ni) > 1.e-4) then
          x(ip) = ((sqrt(max(0.,ni**2 - 2.*(nip1-ni)*(ninti - x(ip))))-ni)/
     &            (nip1-ni) + i)*dx
        else
          x(ip) = ((x(ip) - ninti)/ni  + i)*dx
        endif
      enddo

      return
      end
c=============================================================================
      integer(ISZ) function oneiftrue(l)
      logical(ISZ), intent(IN) :: l

        if(l) then
          oneiftrue = 1
        else
          oneiftrue = 0
        end if

        return
      end function oneiftrue
c=============================================================================
      subroutine reset_temperature()
      use Temperatures
      use Timers
      real(kind=8):: timetemp,wtime
      timetemp = wtime()
        pnumt = 0.
        pnumtw = 0.
        vxbart = 0.
        vybart = 0.
        vzbart = 0.
        vxsqbart = 0.
        vysqbart = 0.
        vzsqbart = 0.
        tempxz = 0.
        tempyz = 0.
        tempzz = 0.
        if(l_temp_rmcorrelations) then
          xbart = 0.
          ybart = 0.
          zbart = 0.
          xsqbart = 0.
          ysqbart = 0.
          zsqbart = 0.
          xvxbart = 0.
          yvybart = 0.
          zvzbart = 0.
        endif
      temperaturestime = temperaturestime + (wtime() - timetemp)
      return
      end
c=============================================================================
      subroutine accumulate_temperature(np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,dt,
     &                          uxpo,uypo,uzpo,is,wp,lw,lrtheta,l2symtry,l4symtry)
c     Compute temperature in Z-slices for species 'is' on a 3-D grid:
c       - the slices can have any position and thickness but cannot overlap,
c       - the min and max of each slice in x, y and z are given respectively in the arrays
c         tslicexmin, tslicexmax, tsliceymin, tsliceymax, tslicezmin and tslicezmax
c         (the reason for having slices with different thickness and dimensions is
c          to allow the temperature measurement to the shape of a distribution, like
c          for example a beam extending over several quadrupoles and accelerating gaps),
c       - the x, y and z temperatures are given in the arrays tempx, tempy and tempz,
c         while averages in each slice are given in the arrays tempxz, tempyz, tempzz.
c       - the calculation is done in three parts:
c         o reset_temperature: zero out all moments
c         o accumulate_temperature: accumulate moments from particles
c         o finalize_temperature: compute final quantities
c       - l_temp_collapseinz=.true.: collapse slices in z, i.e. align Z-locations
c                                    of particles using current velocity (uxp,uyp,yzp) and
c                                    velocity from previous time step (uxpo,uypo,uzpo)
c       - lrtheta=.true.: radial and azimuthal are computed in place of X and Y,
c       - the default temperature unit is in electron-volt. To select the units,
c         set the variable t_unit to the default integers evolt, joule or kelvin.

      use Beam_acc
      use InDiag
      use Picglb
      use ExtPart
      use Temperatures
      use Timers

      integer(ISZ):: np,is
      real(kind=8):: w,dt
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np), wp(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np) ! uxp, uyp anx uzp of previous time step
      logical(ISZ):: lw,lrtheta,l2symtry,l4symtry

      integer(ISZ):: ip,its,izl,ixt,iyt,izt
      real(kind=8):: dti,wt,ddx,ddy,oddx,oddy,xt,yt,wpp,z_local
      real(kind=8):: oneondt,clighti,vzi,zc
      real(kind=8):: xpt,ypt,zpt,vxpt,vypt,vzpt,cost,sint,rpt
      real(kind=8):: timetemp,wtime

      timetemp = wtime()

      if (np==0) return

      wpp = 1.

      oneondt = 1./dvnz(dt)

c       --- loop over particles
      do ip=1,np
        z_local = zp(ip)-zbeam
        izl  = 1+int((z_local - tloc_zmin)*tloc_dzi)

c       --- cycle if particle out of zone of calculation
        if(izl<1 .or. izl>nztlocator) cycle

c       --- loop over temperature slices
        do its = 1, ntl(izl)
          izt = tslice_locator(izl,its)

c         --- cycle if particle not in slice
          if(xp(ip)<=tslicexmin(izt) .or. xp(ip)>tslicexmax(izt) .or.
     &       yp(ip)<=tsliceymin(izt) .or. yp(ip)>tsliceymax(izt) .or.
     &       z_local<=tslicezmin(izt) .or. z_local>tslicezmax(izt)) cycle

          if(l_temp_collapseinz) then
c           --- collapse slice in z, i.e. align Z-locations of particles using
c           --- current velocity and from previous time step
            zc = 0.5*(tslicezmin(izt)+tslicezmax(izt))
            vzi = 1./(uzp(ip)*gaminv(ip)+SMALLPOS)
            dti  = (zbeam+zc-zp(ip))*vzi
            xpt  = xp(ip) + uxp(ip)*dti*gaminv(ip)
            ypt  = yp(ip) + uyp(ip)*dti*gaminv(ip)
            zpt  = z_local + uzp(ip)*dti*gaminv(ip)
            vxpt = (uxp(ip)*(1. + dti*oneondt) - uxpo(ip)*dti*oneondt) * gaminv(ip)
            vypt = (uyp(ip)*(1. + dti*oneondt) - uypo(ip)*dti*oneondt) * gaminv(ip)
            vzpt = (uzp(ip)*(1. + dti*oneondt) - uzpo(ip)*dti*oneondt) * gaminv(ip)
          else
            xpt  = xp(ip)
            ypt  = yp(ip)
            zpt  = z_local
            vxpt = uxp(ip) * gaminv(ip)
            vypt = uyp(ip) * gaminv(ip)
            vzpt = uzp(ip) * gaminv(ip)
          end if
          if(lrtheta) then
            rpt = sqrt(xpt**2+ypt**2)
            if(rpt>SMALLPOS) then
              cost = xpt/rpt
              sint = ypt/rpt
              xpt = vxpt
              vxpt = vxpt*cost + vypt*sint
              vypt = -xpt*sint + vypt*cost
            else
              vxpt = sqrt(vxpt**2+vypt**2)
              vypt = 0.
            end if
            xpt = rpt
            ypt = 0.
          else
            if(l2symtry) ypt = abs(ypt)
            if(l4symtry) then
              xpt=abs(xpt)
              ypt=abs(ypt)
            end if
          end if

c         --- compute coefficients for transverse linear deposition
          xt = (xpt-tslicexmin(izt))*dxti(izt)
          yt = (ypt-tsliceymin(izt))*dyti(izt)
          ixt = max(0,min(nxtslices-1,int(xt)))
          iyt = max(0,min(nytslices-1,int(yt)))
          ddx = xt-ixt
          ddy = yt-iyt
          oddx = 1.-ddx
          oddy = 1.-ddy

          if(lw) wpp = wp(ip)

c         --- perform transverse linear deposition
c         --- node x, y
          wt = oddx*oddy*wpp
          pnumt  (ixt,iyt,izt) = pnumt  (ixt,iyt,izt) + oddx*oddy
          pnumtw (ixt,iyt,izt) = pnumtw (ixt,iyt,izt) + wt

          vxbart  (ixt,iyt,izt) = vxbart  (ixt,iyt,izt) + wt * vxpt
          vybart  (ixt,iyt,izt) = vybart  (ixt,iyt,izt) + wt * vypt
          vzbart  (ixt,iyt,izt) = vzbart  (ixt,iyt,izt) + wt * vzpt

          vxsqbart(ixt,iyt,izt) = vxsqbart(ixt,iyt,izt) + wt * vxpt**2
          vysqbart(ixt,iyt,izt) = vysqbart(ixt,iyt,izt) + wt * vypt**2
          vzsqbart(ixt,iyt,izt) = vzsqbart(ixt,iyt,izt) + wt * vzpt**2

c         --- node x+1, y
          if(ixt+1 <= nxtslices) then
            wt = ddx*oddy*wpp
            pnumt  (ixt+1,iyt,izt) = pnumt  (ixt+1,iyt,izt) + ddx*oddy
            pnumtw (ixt+1,iyt,izt) = pnumtw (ixt+1,iyt,izt) + wt

            vxbart  (ixt+1,iyt,izt) = vxbart  (ixt+1,iyt,izt) + wt * vxpt
            vybart  (ixt+1,iyt,izt) = vybart  (ixt+1,iyt,izt) + wt * vypt
            vzbart  (ixt+1,iyt,izt) = vzbart  (ixt+1,iyt,izt) + wt * vzpt

            vxsqbart(ixt+1,iyt,izt) = vxsqbart(ixt+1,iyt,izt) + wt * vxpt**2
            vysqbart(ixt+1,iyt,izt) = vysqbart(ixt+1,iyt,izt) + wt * vypt**2
            vzsqbart(ixt+1,iyt,izt) = vzsqbart(ixt+1,iyt,izt) + wt * vzpt**2
          end if

c         --- node x, y+1
          if(iyt+1 <= nytslices) then
            wt = oddx*ddy*wpp
            pnumt  (ixt,iyt+1,izt) = pnumt  (ixt,iyt+1,izt) + oddx*ddy
            pnumtw (ixt,iyt+1,izt) = pnumtw (ixt,iyt+1,izt) + wt

            vxbart  (ixt,iyt+1,izt) = vxbart  (ixt,iyt+1,izt) + wt * vxpt
            vybart  (ixt,iyt+1,izt) = vybart  (ixt,iyt+1,izt) + wt * vypt
            vzbart  (ixt,iyt+1,izt) = vzbart  (ixt,iyt+1,izt) + wt * vzpt

            vxsqbart(ixt,iyt+1,izt) = vxsqbart(ixt,iyt+1,izt) + wt * vxpt**2
            vysqbart(ixt,iyt+1,izt) = vysqbart(ixt,iyt+1,izt) + wt * vypt**2
            vzsqbart(ixt,iyt+1,izt) = vzsqbart(ixt,iyt+1,izt) + wt * vzpt**2
          end if

c         --- node x+1, y+1
          if(ixt+1 <= nxtslices .and. iyt+1 <= nytslices) then
            wt = ddx*ddy*wpp
            pnumt  (ixt+1,iyt+1,izt) = pnumt  (ixt+1,iyt+1,izt) + ddx*ddy
            pnumtw (ixt+1,iyt+1,izt) = pnumtw (ixt+1,iyt+1,izt) + wt

            vxbart  (ixt+1,iyt+1,izt) = vxbart  (ixt+1,iyt+1,izt) + wt * vxpt
            vybart  (ixt+1,iyt+1,izt) = vybart  (ixt+1,iyt+1,izt) + wt * vypt
            vzbart  (ixt+1,iyt+1,izt) = vzbart  (ixt+1,iyt+1,izt) + wt * vzpt

            vxsqbart(ixt+1,iyt+1,izt) = vxsqbart(ixt+1,iyt+1,izt) + wt * vxpt**2
            vysqbart(ixt+1,iyt+1,izt) = vysqbart(ixt+1,iyt+1,izt) + wt * vypt**2
            vzsqbart(ixt+1,iyt+1,izt) = vzsqbart(ixt+1,iyt+1,izt) + wt * vzpt**2
          end if

          if(l_temp_rmcorrelations) then
c           --- node x, y
            wt = oddx*oddy*wpp
            xbart   (ixt,iyt,izt) = xbart   (ixt,iyt,izt) + wt * xpt
            ybart   (ixt,iyt,izt) = ybart   (ixt,iyt,izt) + wt * ypt
            zbart   (ixt,iyt,izt) = zbart   (ixt,iyt,izt) + wt * zpt

            xsqbart (ixt,iyt,izt) = xsqbart (ixt,iyt,izt) + wt * xpt**2
            ysqbart (ixt,iyt,izt) = ysqbart (ixt,iyt,izt) + wt * ypt**2
            zsqbart (ixt,iyt,izt) = zsqbart (ixt,iyt,izt) + wt * zpt**2

            xvxbart (ixt,iyt,izt) = xvxbart (ixt,iyt,izt) + wt * xpt * vxpt
            yvybart (ixt,iyt,izt) = yvybart (ixt,iyt,izt) + wt * ypt * vypt
            zvzbart (ixt,iyt,izt) = zvzbart (ixt,iyt,izt) + wt * zpt * vzpt

c           --- node x+1, y
            if(ixt+1 <= nxtslicesc) then
              wt = ddx*oddy*wpp
              xbart   (ixt+1,iyt,izt) = xbart   (ixt+1,iyt,izt) + wt * xpt
              ybart   (ixt+1,iyt,izt) = ybart   (ixt+1,iyt,izt) + wt * ypt
              zbart   (ixt+1,iyt,izt) = zbart   (ixt+1,iyt,izt) + wt * zpt

              xsqbart (ixt+1,iyt,izt) = xsqbart (ixt+1,iyt,izt) + wt * xpt**2
              ysqbart (ixt+1,iyt,izt) = ysqbart (ixt+1,iyt,izt) + wt * ypt**2
              zsqbart (ixt+1,iyt,izt) = zsqbart (ixt+1,iyt,izt) + wt * zpt**2

              xvxbart (ixt+1,iyt,izt) = xvxbart (ixt+1,iyt,izt) + wt * xpt * vxpt
              yvybart (ixt+1,iyt,izt) = yvybart (ixt+1,iyt,izt) + wt * ypt * vypt
              zvzbart (ixt+1,iyt,izt) = zvzbart (ixt+1,iyt,izt) + wt * zpt * vzpt
            end if

c           --- node x, y+1
            if(iyt+1 <= nytslicesc) then
              wt = oddx*ddy*wpp
              xbart   (ixt,iyt+1,izt) = xbart   (ixt,iyt+1,izt) + wt * xpt
              ybart   (ixt,iyt+1,izt) = ybart   (ixt,iyt+1,izt) + wt * ypt
              zbart   (ixt,iyt+1,izt) = zbart   (ixt,iyt+1,izt) + wt * zpt

              xsqbart (ixt,iyt+1,izt) = xsqbart (ixt,iyt+1,izt) + wt * xpt**2
              ysqbart (ixt,iyt+1,izt) = ysqbart (ixt,iyt+1,izt) + wt * ypt**2
              zsqbart (ixt,iyt+1,izt) = zsqbart (ixt,iyt+1,izt) + wt * zpt**2

              xvxbart (ixt,iyt+1,izt) = xvxbart (ixt,iyt+1,izt) + wt * xpt * vxpt
              yvybart (ixt,iyt+1,izt) = yvybart (ixt,iyt+1,izt) + wt * ypt * vypt
              zvzbart (ixt,iyt+1,izt) = zvzbart (ixt,iyt+1,izt) + wt * zpt * vzpt
            end if

c           --- node x+1, y+1
            if(ixt+1 <= nxtslicesc .and. iyt+1 <= nytslicesc) then
              wt = ddx*ddy*wpp
              xbart   (ixt+1,iyt+1,izt) = xbart   (ixt+1,iyt+1,izt) + wt * xpt
              ybart   (ixt+1,iyt+1,izt) = ybart   (ixt+1,iyt+1,izt) + wt * ypt
              zbart   (ixt+1,iyt+1,izt) = zbart   (ixt+1,iyt+1,izt) + wt * zpt

              xsqbart (ixt+1,iyt+1,izt) = xsqbart (ixt+1,iyt+1,izt) + wt * xpt**2
              ysqbart (ixt+1,iyt+1,izt) = ysqbart (ixt+1,iyt+1,izt) + wt * ypt**2
              zsqbart (ixt+1,iyt+1,izt) = zsqbart (ixt+1,iyt+1,izt) + wt * zpt**2

              xvxbart (ixt+1,iyt+1,izt) = xvxbart (ixt+1,iyt+1,izt) + wt * xpt * vxpt
              yvybart (ixt+1,iyt+1,izt) = yvybart (ixt+1,iyt+1,izt) + wt * ypt * vypt
              zvzbart (ixt+1,iyt+1,izt) = zvzbart (ixt+1,iyt+1,izt) + wt * zpt * vzpt

            endif
          endif
        enddo
      enddo

      temperaturestime = temperaturestime + (wtime() - timetemp)
      return
      end
c===========================================================================
      subroutine finalize_temperature(is,m)
      use Constant
      use Temperatures
      use Timers
      integer(ISZ):: is, ixt, iyt, izt
      real(kind=8):: m,timetemp,wtime
      real(kind=8):: pnumi,tfact,tottmp
      real(kind=8):: delxsq,delvxsq,delxvx
      real(kind=8):: delysq,delvysq,delyvy
      real(kind=8):: delzsq,delvzsq,delzvz
      timetemp = wtime()

#ifdef MPIPARALLEL
c     --- For slave, call routine which sums moments over processors.
      call parallel_sum_temperature
#endif

c     --- set multiplying factor for proper units
      if(t_units==evolt)  tfact = 0.5*m/echarge
      if(t_units==joule)  tfact = 0.5*m
      if(t_units==kelvin) tfact = 0.5*m/boltzmann

c     --- Complete the calculation of temperatures: divide by particle number
      do izt = 1, nztslices
        do iyt = 0, nytslices
          do ixt = 0, nxtslices

            pnumi = 1./(pnumtw(ixt,iyt,izt)+SMALLPOS)

c           --- Compute averages
            vxbart  (ixt,iyt,izt) = vxbart  (ixt,iyt,izt) * pnumi
            vybart  (ixt,iyt,izt) = vybart  (ixt,iyt,izt) * pnumi
            vzbart  (ixt,iyt,izt) = vzbart  (ixt,iyt,izt) * pnumi

            vxsqbart(ixt,iyt,izt) = vxsqbart(ixt,iyt,izt) * pnumi
            vysqbart(ixt,iyt,izt) = vysqbart(ixt,iyt,izt) * pnumi
            vzsqbart(ixt,iyt,izt) = vzsqbart(ixt,iyt,izt) * pnumi

c           --- Compute second order moments with averages subtracted
            delvxsq = vxsqbart(ixt,iyt,izt) - vxbart(ixt,iyt,izt)**2
            delvysq = vysqbart(ixt,iyt,izt) - vybart(ixt,iyt,izt)**2
            delvzsq = vzsqbart(ixt,iyt,izt) - vzbart(ixt,iyt,izt)**2

c           --- Compute linear correlations
            if(l_temp_rmcorrelations) then
              xbart   (ixt,iyt,izt) = xbart   (ixt,iyt,izt) * pnumi
              ybart   (ixt,iyt,izt) = ybart   (ixt,iyt,izt) * pnumi
              zbart   (ixt,iyt,izt) = zbart   (ixt,iyt,izt) * pnumi

              xsqbart (ixt,iyt,izt) = xsqbart (ixt,iyt,izt) * pnumi
              ysqbart (ixt,iyt,izt) = ysqbart (ixt,iyt,izt) * pnumi
              zsqbart (ixt,iyt,izt) = zsqbart (ixt,iyt,izt) * pnumi

              xvxbart (ixt,iyt,izt) = xvxbart (ixt,iyt,izt) * pnumi
              yvybart (ixt,iyt,izt) = yvybart (ixt,iyt,izt) * pnumi
              zvzbart (ixt,iyt,izt) = zvzbart (ixt,iyt,izt) * pnumi

              delxsq  = xsqbart (ixt,iyt,izt) - xbart (ixt,iyt,izt)**2
              delxvx  = xvxbart (ixt,iyt,izt) - xbart (ixt,iyt,izt)*vxbart(ixt,iyt,izt)

              delysq  = ysqbart (ixt,iyt,izt) - ybart (ixt,iyt,izt)**2
              delyvy  = yvybart (ixt,iyt,izt) - ybart (ixt,iyt,izt)*vybart(ixt,iyt,izt)

              delzsq  = zsqbart (ixt,iyt,izt) - zbart (ixt,iyt,izt)**2
              delzvz  = zvzbart (ixt,iyt,izt) - zbart (ixt,iyt,izt)*vzbart(ixt,iyt,izt)

              if(abs(delxsq)>SMALLPOS) delvxsq = max(SMALLPOS,delvxsq - delxvx**2/delxsq)
              if(abs(delysq)>SMALLPOS) delvysq = max(SMALLPOS,delvysq - delyvy**2/delysq)
              if(abs(delzsq)>SMALLPOS) delvzsq = max(SMALLPOS,delvzsq - delzvz**2/delzsq)
            end if

c           --- Compute temperatures
            tempx(ixt,iyt,izt,is) = tfact*delvxsq
            tempy(ixt,iyt,izt,is) = tfact*delvysq
            tempz(ixt,iyt,izt,is) = tfact*delvzsq

          enddo
        enddo
        tottmp = sum(pnumtw(:,:,izt))
        if(tottmp>SMALLPOS) then
          tempxz(izt,is) = sum(tempx(:,:,izt,is)*pnumtw(:,:,izt))/tottmp
          tempyz(izt,is) = sum(tempy(:,:,izt,is)*pnumtw(:,:,izt))/tottmp
          tempzz(izt,is) = sum(tempz(:,:,izt,is)*pnumtw(:,:,izt))/tottmp
        end if
      enddo
      temperaturestime = temperaturestime + (wtime() - timetemp)
      return
      end
c=============================================================================
      subroutine gett(is,lrtheta,l2symtry,l4symtry)
c     Compute temperature in Z-slices for species 'is' on a 3-D grid:
c       - the slices can have any position and thickness but cannot overlap,
c       - the min and max of each slice in x, y and z are given respectively in the arrays
c         tslicexmin, tslicexmax, tsliceymin, tsliceymax, tslicezmin and tslicezmax
c         (the reason for having slices with different thickness and dimensions is
c          to allow the temperature measurement to the shape of a distribution, like
c          for example a beam extending over several quadrupoles and accelerating gaps),
c       - the x, y and z temperatures are given in the arrays tempx, tempy and tempz,
c         while averages in each slice are given in the arrays tempxz, tempyz, tempzz.
c       - the calculation is done in three parts:
c         o reset_temperature: zero out all moments
c         o accumulate_temperature: accumulate moments from particles
c         o finalize_temperature: compute final quantities
c       - l_temp_collapseinz=.true.: collapse slices in z, i.e. align Z-locations
c                                    of particles using current velocity (uxp,uyp,yzp) and
c                                    velocity from previous time step (uxpo,uypo,uzpo)
c       - lrtheta=.true.: radial and azimuthal are computed in place of X and Y.
c       - the default temperature unit is in electron-volt. To select the units,
c         set the variable t_unit to the default integers evolt, joule or kelvin.
        use InGen
        use Particles,Only: pgroup,wpid
        integer(ISZ):: is, ipmin, itask
        logical(ISZ):: lrtheta,l2symtry,l4symtry

        ipmin = pgroup%ins(is)
        call reset_temperature()
        if(wpid>0) then
          call accumulate_temperature(pgroup%nps(is),
     &                pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     &                pgroup%gaminv(ipmin),pgroup%sw(is),dt,
     &                pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),is,
     &                pgroup%pid(ipmin,wpid),.true.,lrtheta,l2symtry,l4symtry)
        else
          call accumulate_temperature(pgroup%nps(is),
     &                pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     &                pgroup%gaminv(ipmin),pgroup%sw(is),dt,
     &                pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),is,
     &                pgroup%pid(ipmin,1),.false.,lrtheta,l2symtry,l4symtry)
        end if
        call finalize_temperature(is,pgroup%sm(is))
        return
      end
c=============================================================================
      subroutine setregulartgrid(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,dz,nzloc,lcollapse,lcorrel)
c       Setup regular grid for temperature calculation. The temperature is calculated in nz slices
c       evenly spaced of thickness dz. The spacing between each slice is given by (zmax-zmin)/(nz-1).
c       In each slice, the temperature will be computed on a nx+1*ny+1 grid of size (xmin,xmax,ymin,ymax).
        use Temperatures
        use InPart

        integer(ISZ) :: nx, ny, nz, nzloc
        logical(ISZ) :: lcollapse, lcorrel
        real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax, dz

        integer(ISZ) :: i

        nstemp = ns

        nxtslices = nx
        nytslices = ny
        nztslices = nz
        if(lcorrel) then
          l_temp_rmcorrelations = .true.
          nxtslicesc = nx
          nytslicesc = ny
          nztslicesc = nz
        else
          l_temp_rmcorrelations = .false.
        end if
        nztlocator = nzloc
        call gchange("Temperatures",0)
        tslicexmin = xmin
        tslicexmax = xmax
        tsliceymin = ymin
        tsliceymax = ymax
        dxti = nx/(xmax-xmin)
        dyti = ny/(ymax-ymin)
        if(nztslices>1) then
          do i = 1, nztslices
            tslicezmin(i) = zmin+(i-1)*(zmax-zmin)/(nztslices-1)-0.5*dz
            tslicezmax(i) = tslicezmin(i) + dz
          end do
        else
          tslicezmin = zmin-0.5*dz
          tslicezmax = zmin+0.5*dz
        end if
        tloc_dzi = nzloc/(tslicezmax(nztslices)-tslicezmin(1))
        tloc_zmin = tslicezmin(1)
        tloc_zmax = tslicezmax(nztslices)
        call set_tslice_locator()

        return
      end
c=============================================================================
      subroutine set_tslice_locator()
        use Temperatures
        integer(ISZ) :: i, ii, izmin, izmax

        integer(ISZ), allocatable :: ntsloc(:,:)
        integer(ISZ):: allocerror

        allocate(ntsloc(nztlocator,nztslices),stat=allocerror)
        if (allocerror /= 0) then
          print*,"set_tslice_locator: allocation error ",allocerror,
     &           ": could not allocate ntsloc to shape ",nztlocator,nztslices
          stop
        endif

        ntsloc = 0
        ntl = 0
        do i = 1, nztslices
          izmin = min(nztlocator,max(1,1+int((tslicezmin(i)-tloc_zmin)*tloc_dzi)))
          izmax = min(nztlocator,max(1,1+int((tslicezmax(i)-tloc_zmin)*tloc_dzi)))
          do ii = izmin,izmax
            ntl(ii) = ntl(ii)+1
            ntsloc(ii,ntl(ii)) = i
          end do
        end do
        ntlmax = maxval(ntl)
        call gchange("Temperatures",0)
        tslice_locator = ntsloc(:,:ntlmax)
c there is a problem here with the intel compiler
        deallocate(ntsloc)

        return
      end
c=============================================================================
      subroutine impact_ion(is1,is2,nbp,wp,shiftx,shifty,shiftz,deltax,deltay,deltaz,condid)
      use Constant
      use Particles,Only: pgroup,wpid,spid,ssn
      use InPart
      use LostParticles
      INTEGER, INTENT(IN) :: is1, is2, condid
      REAL(8), INTENT(IN) :: nbp,wp,shiftx,shifty,shiftz,deltax,deltay,deltaz

c Create particles of specie is2 created by impact of particles of species is1
c - nbp: number of new particles generated for each lost particle (can be fractional)
c - wp: energy in volts of new particles
c - shiftz,y,z: the location where the new particles are created relative
c               to where the original particles are lost
c - deltax,y,z: the length over which the new particles are distributed
c               relative to the location where the original particles are lost.
c               The distribution is uniform between -delta/2 and +delta/2.
c - condid: id of conductor - only particles that are lost on the specified
c           conductor generate new particles.

      INTEGER :: ip,i1,i2,ir,nbpt
      REAL(8) :: gamma, uz_tmp
      REAL(8), EXTERNAL :: wranf

      IF(nbp==0) return

      gamma = 1.+wp*echarge/(pgroup%sm(is2)*clight**2)
      IF(gamma-1.<1.e-6) then
        uz_tmp = gamma*SQRT(1.-1./(gamma**2))*clight
      else
        uz_tmp = SQRT(2.*wp*echarge/pgroup%sm(is2))
      END if

      do ip = inslost(is1), inslost(is1)+npslost(is1)-1
        if (condid > 0 .and. pidlost(ip,npidlost) .ne. condid) cycle
        nbpt = int(nbp+wranf())
        i1 = pgroup%ins(is2) + pgroup%nps(is2)
        i2 = i1 + nbpt - 1
        call chckpart(pgroup,is2,0,nbpt,.false.)
        do ir = 1, nbpt
          pgroup%xp(i1+ir-1) = xplost(ip)+shiftx+(wranf()-0.5)*deltax
          pgroup%yp(i1+ir-1) = yplost(ip)+shifty+(wranf()-0.5)*deltay
          pgroup%zp(i1+ir-1) = zplost(ip)+shiftz+(wranf()-0.5)*deltaz
          if(spid>0) then
            pgroup%pid(i1+ir-1,spid) = ssn
            ssn=ssn+1
          end if
        end do
        pgroup%uxp(i1:i2) = (wranf()-0.5)*uz_tmp
        pgroup%uyp(i1:i2) = (wranf()-0.5)*uz_tmp
c       uzp(i1:i2) = (wranf()-0.5)*uz_tmp !old
        pgroup%uzp(i1:i2) = -abs((wranf()-0.5)*uz_tmp)
        if(wpid>0) pgroup%pid(i1:i2,wpid) = pidlost(ip,wpid)
        pgroup%gaminv(i1:i2) = 1./sqrt(1. +
     &                                 (pgroup%uxp(i1:i2)**2 +
     &                                  pgroup%uyp(i1:i2)**2 +
     &                                  pgroup%uzp(i1:i2)**2)/clight**2)
        pgroup%nps(is2) = pgroup%nps(is2) + nbpt
      end do

      return
      end subroutine impact_ion
c=============================================================================
