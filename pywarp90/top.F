#include "top.h"
c=============================================================================
c@(#) File TOP.M, version $Revision: 3.1 $, $Date: 2001/04/12 23:19:26 $
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
      use Lattice

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
      use Topversion
c  Echoes code version,etc. to output files when they're created
      integer     :: iout
      call printpkgversion(iout,"Main package TOP",verstop)
      return
      end
c=============================================================================
      subroutine cigar(np,zunifrm,zpunifrm,z,zp,perpscal,straight,
     &                 szp,phi,zpm,zpms)
      use Constant
      integer     :: np
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
      real(kind=8):: pion2,oneonpi,twoonpi,zs,zend,zendr,zmap,zl,zendi,a1,a2
      real(kind=8):: third,sixth,root3,q
      integer     :: i

      pion2 = pi*0.5
      oneonpi = 1./pi
      twoonpi = 1./pion2
      zs = straight
      zend = (1.-zs)*0.5
      zendr = 1.-zend
      zmap = zs+(4./3.)*zend
      zl = zend/3.
      zendi = 1./dvnz(zend)
      a1 = -3.
      a2 = 0.

c initialize stuff for cubic solver
      third = 1./3.
      sixth = 1./6.
      root3 = sqrt(3.)
      q = a1*third-a2**2/9.0

      do i=1,np
c
c     map the center of the beam into the uniform region
c
          z(i) = zl + zunifrm(i)*zmap
          zpm(i) = 0.
c
c     solve a cubic to map the two parabolic end caps
c
          if(z(i).lt.zend) then
c  Solve for the third cubic root
c            t = q**3+r**2<0 always, so three real(kind=8):: roots
             phi(i) = atan2(sqrt(1.-(1.5*(z(i)-zend)*zendi)**2),
     &                     -1.5*(z(i)-zend)*zendi)*third
             zpm(i) = -cos(phi(i))+root3*sin(phi(i))
             z(i) = (1.+zpm(i))*zend
          elseif(z(i).gt.zendr) then
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
      use Particles

c Sets up species related quantities which were not set by the user, based
c on the rest of the user's input.    


      integer     :: is

c     First, set parameters for whole beam (for envelope solver).  If
c     they are unset, set to parameters of species 1.

c     --- Beam size in X
      if (a0 == 0.) then
        a0 = a0_s(1)
      endif
c     --- Beam divergence in X
      if (ap0 == 0.) then
        ap0 = ap0_s(1)
      endif
c     --- Beam size in Y
      if (b0 == 0.) then
        b0 = b0_s(1)
      endif
c     --- Beam divergence in Y
      if (bp0 == 0.) then
        bp0 = bp0_s(1)
      endif
c     --- Beam centroid in X
      if (x0 == 0.) then
        x0 = x0_s(1)
      endif
c     --- Beam centroid angle in X
      if (xp0 == 0.) then
        xp0 = xp0_s(1)
      endif
c     --- Beam centroid in Y
      if (y0 == 0.) then
        y0 = y0_s(1)
      endif
c     --- Beam centroid angle in Y
      if (yp0 == 0.) then
        yp0 = yp0_s(1)
      endif
c     --- Atomic number
      if (aion == 0.) then
        aion = aion_s(1)
      endif
c     --- Particle energy
      if (ekin == 0.) then
        ekin = ekin_s(1)
      endif
c     --- Emittance
      if (emit == 0.) then
        emit = emit_s(1)
      endif
      if (emitx == 0.) then
        emitx = emitx_s(1)
      endif
      if (emity == 0.) then
        emity = emity_s(1)
      endif
c     --- Normalized emittance
      if (emitn == 0.) then
        emitn = emitn_s(1)
      endif
      if (emitnx == 0.) then
        emitnx = emitnx_s(1)
      endif
      if (emitny == 0.) then
        emitny = emitny_s(1)
      endif
c     --- Current
      if (ibeam == 0.) then
        ibeam = ibeam_s(1)
      endif
c     --- Charge state
      if (zion == 0.) then
        zion = zion_s(1)
      endif
c     --- Axial velocity
      if (vbeam == 0.) then
        vbeam = vbeam_s(1)
      endif
c     --- Axial velocity tilt
      if (vtilt == 0.) then
        vtilt = vtilt_s(1)
      endif
c     --- Transverse thermal spread
      if (vthperp == 0.) then
        vthperp = vthperp_s(1)
      endif
c     --- Axial thermal spread
      if (vthz == 0.) then
        vthz = vthz_s(1)
      endif
c     --- Minimum initial z of beam
      if (zimin == 0.) then
        zimin = zimin_s(1)
      endif
c     --- Maximum initial z of beam
      if (zimax == 0.) then
        zimax = zimax_s(1)
      endif
c     --- Fractional length of uniform part of beam
      if (straight == 0.) then
        straight = straight_s(1)
      endif

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
        if (a0_s(is) == 0. .or. is == 1) then
          a0_s(is) = a0
        endif
c       --- Beam divergence in X
        if (ap0_s(is) == 0. .or. is == 1) then
          ap0_s(is) = ap0
        endif
c       --- Beam size in Y
        if (b0_s(is) == 0. .or. is == 1) then
          b0_s(is) = b0
        endif
c       --- Beam divergence in Y
        if (bp0_s(is) == 0. .or. is == 1) then
          bp0_s(is) = bp0
        endif
c       --- Beam centroid in X
        if (x0_s(is) == 0. .or. is == 1) then
          x0_s(is) = x0
        endif
c       --- Beam centroid angle in X
        if (xp0_s(is) == 0. .or. is == 1) then
          xp0_s(is) = xp0
        endif
c       --- Beam centroid in Y
        if (y0_s(is) == 0. .or. is == 1) then
          y0_s(is) = y0
        endif
c       --- Beam centroid angle in Y
        if (yp0_s(is) == 0. .or. is == 1) then
          yp0_s(is) = yp0
        endif
c       --- Atomic number
        if (aion_s(is) == 0 .or. is == 1) then
          aion_s(is) = aion
        endif
c       --- Particle energy
        if (ekin_s(is) == 0 .or. is == 1) then
          ekin_s(is) = ekin
        endif
c       --- Emittance
        if (emit_s(is) == 0 .or. is == 1) then
          emit_s(is) = emit
        endif
        if (emitx_s(is) == 0 .or. is == 1) then
          emitx_s(is) = emitx
        endif
        if (emity_s(is) == 0 .or. is == 1) then
          emity_s(is) = emity
        endif
c       --- Normalized emittance
        if (emitn_s(is) == 0 .or. is == 1) then
          emitn_s(is) = emitn
        endif
        if (emitnx_s(is) == 0 .or. is == 1) then
          emitnx_s(is) = emitnx
        endif
        if (emitny_s(is) == 0 .or. is == 1) then
          emitny_s(is) = emitny
        endif
c       --- Current
        if (ibeam_s(is) == 0 .or. is == 1) then
          ibeam_s(is) = ibeam
        endif
c       --- Charge state
        if (zion_s(is) == 0 .or. is == 1) then
          zion_s(is) = zion
        endif
c       --- Axial velocity
        if (vbeam_s(is) == 0 .or. is == 1) then
          vbeam_s(is) = vbeam
        endif
c       --- Axial velocity tilt
        if (vtilt_s(is) == 0 .or. is == 1) then
          vtilt_s(is) = vtilt
        endif
c       --- Transverse thermal spread
        if (vthperp_s(is) == 0 .or. is == 1) then
          vthperp_s(is) = vthperp
        endif
c       --- Axial thermal spread
        if (vthz_s(is) == 0 .or. is == 1) then
          vthz_s(is) = vthz
        endif
c       --- Minimum initial z of beam
        if (zimin_s(is) == 0. .or. is == 1) then
          zimin_s(is) = zimin
        endif
c       --- Maximum initial z of beam
        if (zimax_s(is) == 0. .or. is == 1) then
          zimax_s(is) = zimax
        endif
c       --- Fractional length of uniform part of beam
        if (straight_s(is) == 0. .or. is == 1) then
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
      use Particles

c  Sets globally derived quantities that need computation.


      integer     :: is
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
            if (vbeam_s(is) .eq. 0.) then

c             --- Beam energy in units of mc**2
              ke = jperev*ekin_s(is)/dvnz(aion_s(is)*amu*clight**2)
              gammabar = 1. + ke

c             --- The expression for vbeam was rewritten to avoid the
c             --- problem of taking the difference of nearly equal numbers
c             --- (i.e.  (1-1/gammabar) where gammabar ~ 1).
c             --- vbeam_s(is) = clight * sqrt(1.-1./gammbar**2)
              vbeam_s(is) = clight * sqrt((2*ke+ke**2)/gammabar**2)

            elseif(ekin_s(is) .eq. 0.) then

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
        if (vbeam .eq. 0.) then
          ke = jperev * ekin / dvnz(aion * amu * clight**2)
          gammabar = 1. + ke
          vbeam = clight * sqrt((2*ke+ke**2)/gammabar**2)
        elseif(ekin .eq. 0.) then
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
            if (vbeam_s(is) .eq. 0.) then
              vbeam_s(is) = sqrt(2.*ekin_s(is)*jperev/dvnz(aion_s(is))/amu)
            elseif(ekin_s(is) .eq. 0.) then
              ekin_s(is) = 0.5*(aion_s(is)*amu*vbeam_s(is)**2)/jperev
            endif
          endif
        enddo

c       --- for the beam in general
        if (vbeam .eq. 0.) then
          vbeam = sqrt(2.*ekin*jperev/dvnz(aion)/amu)
        elseif(ekin .eq. 0.) then
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
      if (rwall.gt.0. .and. gfactor.eq.0.) then
         gfactor = log(rwall**2/dvnz(a0*b0))
      endif

      return
      end
c============================================================================ 
      subroutine resetlat
      use Lattice

c  Resizes the lattice arrays to the final lengths needed by the actual data.
c  Sets dipole field, lattice element 0, flags for bends, etc.


      integer     :: idrft,ib,id,iq,is,ih,ii,ia,im,ip
      integer     :: nndrft,nnbend,nndipo,nnquad,nnsext 
      integer     :: nnhele,nnhmlt,nnaccl,nnemlt,nnmmlt,nnmmlt2,nnbgrd,nnbgrd2,nnpgrd

c     Find size of lattice element arrays for resizing 
c     --- initialize counting variables  
      nndrft = 0
      nnbend = 0
      nndipo = 0
      nnquad = 0
      nnsext = 0
      nnhele = 0
      nnhmlt = 0
      nnaccl = 0
      nnemlt = 0
      nnmmlt = 0
      nnmmlt2 = 0 
      nnbgrd = 0
      nnbgrd2 = 0 
      nnpgrd = 0
c     --- count drfts 
      do idrft = 0, ndrft
         if (drftzs(idrft) .ne. drftze(idrft)) nndrft = idrft
      enddo
c     --- count bends 
      do ib = 0, nbend
         if (bendzs(ib) .ne. bendze(ib)) nnbend = ib
      enddo
c     --- count dipoles 
      do id = 0, ndipo
         if (dipozs(id) .ne. dipoze(id)) nndipo = id
      enddo
c     --- count quads 
      do iq = 0, nquad
         if (quadzs(iq) .ne. quadze(iq)) nnquad = iq
      enddo
c     --- count sextupoles 
      do is = 0, nsext
         if (sextzs(is) .ne. sextze(is)) nnsext = is
      enddo
c     --- count hard-edge multipole elements and number of multipole 
c     --- components in each element
      do ih = 0, nhele
         if (helezs(ih) .ne. heleze(ih)) nnhele = ih
         do ii = 1, nhmlt 
            if (heleae(ii,ih) .ne. 0. .or. heleep(ii,ih) .ne. 0.) helene(ih)=ii
            if (heleam(ii,ih) .ne. 0. .or. helemp(ii,ih) .ne. 0.) helenm(ih)=ii
         enddo
         if (nnhmlt .lt. helene(ih)) nnhmlt = helene(ih)
         if (nnhmlt .lt. helenm(ih)) nnhmlt = helenm(ih)           
      enddo
c     --- count acceleration gaps 
      do ia = 0, naccl
         if (acclzs(ia) .ne. acclze(ia)) nnaccl = ia
      enddo
c     --- count multipole data sets 
      do im = 0, nemlt
         if (emltid(im) .gt. 0) nnemlt = im
      enddo
      do im = 0, nmmlt
         if (mmltid(im) .gt. 0) nnmmlt = im
      enddo
      do im = 0, nmmlt2
         if (mmlt2id(im) .gt. 0) nnmmlt2 = im
      enddo
c     --- count bgridded field data sets 
      do ib = 0, nbgrd
         if (bgrdid(ib) .gt. 0) nnbgrd = ib
      enddo
      do ib = 0, nbgrd2
         if (bgrd2id(ib) .gt. 0) nnbgrd2 = ib
      enddo
      do ip = 0, npgrd
         if (pgrdid(ip) .gt. 0) nnpgrd = ip
      enddo
c     --- set lengths equal to counted size 
      ndrft = nndrft
      nbend = nnbend
      ndipo = nndipo
      nquad = nnquad
      nsext = nnsext
      nhele = nnhele
      nhmlt = nnhmlt  
      naccl = nnaccl
      nemlt = nnemlt
      nmmlt = nnmmlt
      nmmlt2 = nnmmlt2 
      nbgrd = nnbgrd
      nbgrd2 = nnbgrd2
      npgrd = nnpgrd

c     --- resize lattice arrays (use corrected ndipo size if autoset dipoles, 
c     --- bends imply dipoles for autoset dipoles, see following code)
      if (diposet) ndipo = max(ndipo, nbend)    
      call gchange("Lattice", 0)

c     --- the rest of the lattice reset is carried out in a work routine 
c     --- to avoid possibilites of a bug associated with resizing 
c     --- multidimensional dynamics arrays (used for hard-edge multipole 
c     --- elements) and then using them in the same routine  

      call resetlat_work

      return 
      end 
c============================================================================ 
      subroutine resetlat_work
      use InGen
      use Lattice
      use Beam_acc
      use Constant
      use Mult_data
      use BGRDdata
      use PGRDdata

c  work routine for the resizing of lattice arrays by routine resetlat  

      integer     :: i,idrft,ib,id,iq,is,ih,ii,ia,ie,im,iz,l
      real(kind=8):: rnorm,emax,mmax,ephmax,mphmax,dzio2

c     --- if autoset dipoles are specified, then derive dipole locations 
c     --- consistently from bend locations and radii 
      if (diposet) then 
        do id = 0, ndipo
           if (id .le. nbend) then
              if (bendrc(id) .ne. 0.) then
c                --- auto-set dipole field to match bend radius of curvature
                 if (dipoby(id) .eq. 0. .and. dipoex(id) .eq. 0.) then 
                    dipoby(id) = aion * amu / (zion * echarge)
     &                            * gammabar * vbeam / dvnz(bendrc(id))
                 endif
c                --- auto-set dipole starts and ends to match those of bends
                 if (dipozs(id) .eq. dipoze(id)) then
                    dipozs(id) = bendzs(id)
                    dipoze(id) = bendze(id)
                 endif
c                --- auto-set dipole entrance & exit angles for "box dipoles"
                 if (dipotype .eq. "box") then
                  dipota(id) = tan(.5*(dipoze(id)-dipozs(id))/dvnz(bendrc(id)))
                  dipotb(id) = tan(.5*(dipoze(id)-dipozs(id))/dvnz(bendrc(id)))
                 endif
              endif
           endif
        enddo
      endif 

c     --- The following is done only when zlatperi > 0. Otherwise, the lattice
c     --- is assumed not to repeat so the first element of each type must
c     --- already be defined appropriately.
c     --- Set lattice element zero to be a periodic repeat of the last element.
c     --- (Do this only if element hasn't been set already).
c     --- NOTE: all lattice elements must have their starts in [0,zlatperi).
c     --- The end of the final element may be beyond zlatperi, but should not
c     --- be so far beyond that it overlaps the periodic repeat of element 1.
      if (zlatperi > 0.) then
        if (drftzs(0) .eq. drftze(0)) then
           drftzs(0) = drftzs(ndrft) - zlatperi
           drftze(0) = drftze(ndrft) - zlatperi
           drftap(0) = drftap(ndrft)
           drftox(0) = drftox(ndrft)
           drftoy(0) = drftoy(ndrft)
        endif
        if (bendzs(0) .eq. bendze(0)) then
           bendzs(0) = bendzs(nbend) - zlatperi
           bendze(0) = bendze(nbend) - zlatperi
           bendrc(0) = bendrc(nbend)
        endif
        if (dipozs(0) .eq. dipoze(0)) then
           dipozs(0) = dipozs(ndipo) - zlatperi
           dipoze(0) = dipoze(ndipo) - zlatperi
           dipoby(0) = dipoby(ndipo)
           dipobx(0) = dipobx(ndipo)
           dipoex(0) = dipoex(ndipo)
           dipoey(0) = dipoey(ndipo)
           dipox1(0) = dipox1(ndipo)
           dipox2(0) = dipox2(ndipo)
           dipov1(0) = dipov1(ndipo)
           dipov2(0) = dipov2(ndipo)
           dipol1(0) = dipol1(ndipo)
           dipol2(0) = dipol2(ndipo)
           dipow1(0) = dipow1(ndipo)
           dipow2(0) = dipow2(ndipo)
        endif
        if (quadzs(0) .eq. quadze(0)) then
           quadzs(0) = quadzs(nquad) - zlatperi
           quadze(0) = quadze(nquad) - zlatperi
           quaddb(0) = quaddb(nquad)
           quadde(0) = quadde(nquad)
           quadts(0) = quadts(nquad)
           quaddt(0) = quaddt(nquad)
           do l = 0,ntquad
             quadet(l,0) = quadet(l,nquad)
             quadbt(l,0) = quadbt(l,nquad)
           enddo
           quadvx(0) = quadvx(nquad)
           quadvy(0) = quadvy(nquad)
           quadap(0) = quadap(nquad)
           quadrr(0) = quadrr(nquad)
           quadrl(0) = quadrl(nquad)
           quadgl(0) = quadgl(nquad)
           quadgp(0) = quadgp(nquad)
           quadpw(0) = quadpw(nquad)
           quadpa(0) = quadpa(nquad)
           quadpr(0) = quadpr(nquad)
           quadsl(0) = quadsl(nquad)
           qdelglx(0) = qdelglx(nquad)
           qdelgly(0) = qdelgly(nquad)
           qdelaxp(0) = qdelaxp(nquad)
           qdelaxm(0) = qdelaxm(nquad)
           qdelayp(0) = qdelayp(nquad)
           qdelaym(0) = qdelaym(nquad)
           qdelrxp(0) = qdelrxp(nquad)
           qdelrxm(0) = qdelrxm(nquad)
           qdelryp(0) = qdelryp(nquad)
           qdelrym(0) = qdelrym(nquad)
           qdelvxp(0) = qdelvxp(nquad)
           qdelvxm(0) = qdelvxm(nquad)
           qdelvyp(0) = qdelvyp(nquad)
           qdelvym(0) = qdelvym(nquad)
           qdeloxp(0) = qdeloxp(nquad)
           qdeloxm(0) = qdeloxm(nquad)
           qdeloyp(0) = qdeloyp(nquad)
           qdeloym(0) = qdeloym(nquad)
           qdelpwl(0) = qdelpwl(nquad)
           qdelpwr(0) = qdelpwr(nquad)
           qdelpal(0) = qdelpal(nquad)
           qdelpar(0) = qdelpar(nquad)
           qdelprl(0) = qdelprl(nquad)
           qdelprr(0) = qdelprr(nquad)
        endif
        if (sextzs(0) .eq. sextze(0)) then
           sextzs(0) = sextzs(nsext) - zlatperi
           sextze(0) = sextze(nsext) - zlatperi
           sextdb(0) = sextdb(nsext)
           sextde(0) = sextde(nsext)
        endif
        if (helezs(0) .eq. heleze(0)) then
           helezs(0) = helezs(nhele) - zlatperi
           heleze(0) = heleze(nhele) - zlatperi
           do ii = 1, nhmlt 
              heleae(ii,0) = heleae(ii,nhele)
              heleep(ii,0) = heleep(ii,nhele)
              heleam(ii,0) = heleam(ii,nhele)
              helemp(ii,0) = helemp(ii,nhele)
              hele_n(ii,0) = hele_n(ii,nhele)
              hele_v(ii,0) = hele_v(ii,nhele)
              helepe(ii,0) = helepe(ii,nhele)
              helepm(ii,0) = helepm(ii,nhele)
           enddo 
           helene(0) = helene(nhele)  
           helenm(0) = helenm(nhele)  
           heleox(0) = heleox(nhele)   
           heleoy(0) = heleoy(nhele) 
           helerr(0) = helerr(nhele)
           helerl(0) = helerl(nhele)
           helegl(0) = helegl(nhele)
           helegp(0) = helegp(nhele)
           helepw(0) = helepw(nhele)
           helepa(0) = helepa(nhele)
        endif
        if (acclzs(0) .eq. acclze(0)) then
           acclzs(0) = acclzs(naccl) - zlatperi
           acclze(0) = acclze(naccl) - zlatperi
           acclez(0) = acclez(naccl)
           acclxw(0) = acclxw(naccl)
           acclsw(0) = acclsw(naccl)
           acclts(0) = acclts(naccl)
           accldt(0) = accldt(naccl)
           do l = 0,ntaccl
             acclet(l,0) = acclet(l,naccl)
           enddo
        endif
c       --- Note that these elemente are different since only the centers are
c       --- specified.
        if (emltid(0) .eq. 0) then
           emltzs(0) = emltzs(nemlt) - zlatperi
           emltze(0) = emltze(nemlt) - zlatperi
           emltph(0) = emltph(nemlt)
           emltsf(0) = emltsf(nemlt)
           emltsc(0) = emltsc(nemlt)
           emltid(0) = emltid(nemlt)
           emltap(0) = emltap(nemlt) 
           emltrr(0) = emltrr(nemlt)
           emltrl(0) = emltrl(nemlt)
           emltgl(0) = emltgl(nemlt)
           emltgp(0) = emltgp(nemlt)
           emltpw(0) = emltpw(nemlt)
           emltpa(0) = emltpa(nemlt)
        endif
        if (mmltid(0) .eq. 0) then
           mmltzs(0) = mmltzs(nmmlt) - zlatperi
           mmltze(0) = mmltze(nmmlt) - zlatperi
           mmltph(0) = mmltph(nmmlt)
           mmltsf(0) = mmltsf(nmmlt)
           mmltsc(0) = mmltsc(nmmlt)
           mmltid(0) = mmltid(nmmlt)
           mmltap(0) = mmltap(nmmlt) 
        endif
        if (mmlt2id(0) .eq. 0) then
           mmlt2zs(0) = mmlt2zs(nmmlt2) - zlatperi
           mmlt2ze(0) = mmlt2ze(nmmlt2) - zlatperi
           mmlt2ph(0) = mmlt2ph(nmmlt2)
           mmlt2sf(0) = mmlt2sf(nmmlt2)
           mmlt2sc(0) = mmlt2sc(nmmlt2)
           mmlt2id(0) = mmlt2id(nmmlt2)
           mmlt2ap(0) = mmlt2ap(nmmlt2) 
        endif
        if (bgrdid(0) .eq. 0) then
           bgrdzs(0) = bgrdzs(nbgrd) - zlatperi
           bgrdze(0) = bgrdze(nbgrd) - zlatperi
           bgrdxs(0) = bgrdxs(nbgrd)
           bgrdys(0) = bgrdys(nbgrd)
           bgrdid(0) = bgrdid(nbgrd)
           bgrdsf(0) = bgrdsf(nbgrd)
           bgrdsc(0) = bgrdsc(nbgrd)
           bgrdsy(0) = bgrdsy(nbgrd)
           bgrdox(0) = bgrdox(nbgrd)
           bgrdoy(0) = bgrdoy(nbgrd)
           bgrdph(0) = bgrdph(nbgrd) 
           bgrdap(0) = bgrdap(nbgrd)
        endif
        if (bgrd2id(0) .eq. 0) then
           bgrd2zs(0) = bgrd2zs(nbgrd2) - zlatperi
           bgrd2ze(0) = bgrd2ze(nbgrd2) - zlatperi
           bgrd2xs(0) = bgrd2xs(nbgrd2)
           bgrd2ys(0) = bgrd2ys(nbgrd2)
           bgrd2id(0) = bgrd2id(nbgrd2)
           bgrd2sf(0) = bgrd2sf(nbgrd2)
           bgrd2sc(0) = bgrd2sc(nbgrd2)
           bgrd2sy(0) = bgrd2sy(nbgrd2)
           bgrd2ox(0) = bgrd2ox(nbgrd2)
           bgrd2oy(0) = bgrd2oy(nbgrd2)
           bgrd2ph(0) = bgrd2ph(nbgrd2) 
           bgrd2ap(0) = bgrd2ap(nbgrd2)
        endif
        if (pgrdid(0) .eq. 0) then
           pgrdzs(0) = pgrdzs(npgrd) - zlatperi
           pgrdze(0) = pgrdze(npgrd) - zlatperi
           pgrdxs(0) = pgrdxs(npgrd)
           pgrdys(0) = pgrdys(npgrd)
           pgrdid(0) = pgrdid(npgrd)
           pgrdsf(0) = pgrdsf(npgrd)
           pgrdsc(0) = pgrdsc(npgrd)
           pgrdox(0) = pgrdox(npgrd)
           pgrdoy(0) = pgrdoy(npgrd)
           pgrdph(0) = pgrdph(npgrd) 
           pgrdap(0) = pgrdap(npgrd)
           pgrdrr(0) = pgrdrr(npgrd)
           pgrdrl(0) = pgrdrl(npgrd)
           pgrdgl(0) = pgrdgl(npgrd)
           pgrdgp(0) = pgrdgp(npgrd)
           pgrdpw(0) = pgrdpw(npgrd)
           pgrdpa(0) = pgrdpa(npgrd)
        endif
      endif

c     --- Set flags for existence of bends, dipos, quads, and mults.
      drfts = .false.
      do idrft = 0, ndrft
         if (drftzs(idrft) .ne. drftze(idrft)) drfts = .true.
      enddo
      bends = .false.
      do ib = 0, nbend
         if (bendrc(ib) .ne. 0.) bends = .true.
      enddo
      dipos = .false.
      do id = 0, ndipo
         if (dipoby(id).ne.0. .or. dipoex(id).ne.0.) dipos = .true.
         if (dipobx(id).ne.0. .or. dipoey(id).ne.0.) dipos = .true.
      enddo
      quads = .false.
      if (ntquad .gt. 0) quads = .true.
      do iq = 0, nquad
         if (quaddb(iq).ne.0. .or. quadde(iq).ne.0.) quads = .true.
      enddo
      sexts = .false.
      do is = 0, nsext
         if (sextdb(is).ne.0. .or. sextde(is).ne.0.) sexts = .true.
      enddo
      heles = .false.
      do ih = 0, nhele
         do ii = 1,nhmlt 
            if (heleae(ii,ih).ne.0. .or. heleam(ii,ih).ne.0. .or.
     &          heleep(ii,ih).ne.0. .or. helemp(ii,ih).ne.0.) heles = .true.
         enddo 
      enddo
      accls = .false.
      if (ntaccl .gt. 0) accls = .true.
      do ia = 0, naccl
         if (acclez(ia).ne.0) accls = .true.
      enddo
      emlts = .false.
      do im = 0, nemlt
         if (emltid(im).ne.0) emlts = .true.
      enddo
      mmlts = .false.
      do im = 0, nmmlt
         if (mmltid(im).ne.0) mmlts = .true.
      enddo
      mmlt2s = .false.
      do im = 0, nmmlt2
         if (mmlt2id(im).ne.0) mmlt2s = .true.
      enddo
      bgrds = .false.
      do im = 0, nbgrd
         if (bgrdid(im).ne.0) bgrds = .true.
      enddo
      bgrd2s = .false.
      do im = 0, nbgrd2
         if (bgrd2id(im).ne.0) bgrd2s = .true.
      enddo
      pgrds = .false.
      do im = 0, npgrd
         if (pgrdid(im).ne.0) pgrds = .true.
      enddo

c     --- If qerrxrms or qerryrms not zero, then set qoff arrays using cutoff
c     --- Gaussian distribution of offsets
      if (qerrxrms .ne. 0. .or. qerryrms .ne. 0.) then
        do l = 0,iqerr-1
          qoffx(l) = 0.
          qoffy(l) = 0.
        enddo
        do l = iqerr, nqerr
          qoffx(l) = qerrxrms*rnorm()
          qoffy(l) = qerryrms*rnorm()
        enddo
      endif

c     --- Calculate 1/bgrddx etc. for efficiency (minimizes divides)
c     --- Note that this coding is not executed if bgrdns=0.
      do i=1,bgrdns
        if (bgrddxi(i) .eq. 0.) bgrddxi(i) = 1./bgrddx(i)
        if (bgrddyi(i) .eq. 0.) bgrddyi(i) = 1./bgrddy(i)
        if (bgrddzi(i) .eq. 0.) bgrddzi(i) = 1./bgrddz(i)
      enddo

c     --- Calculate sines and cosines of bgrdph for efficiency
      do i=0,nbgrd
        bgrdsp(i) = sin(bgrdph(i)) 
        bgrdcp(i) = cos(bgrdph(i)) 
      enddo 

c     --- Calculate sines and cosines of bgrd2ph for efficiency
      do i=0,nbgrd2
        bgrd2sp(i) = sin(bgrd2ph(i)) 
        bgrd2cp(i) = cos(bgrd2ph(i)) 
      enddo 

c     --- Calculate 1/pgrddx etc. for efficiency (minimizes divides)
c     --- Note that this coding is not executed if pgrdns=0.
      do i=1,pgrdns
        if (pgrddxi(i) .eq. 0.) pgrddxi(i) = 1./pgrddx(i)
        if (pgrddyi(i) .eq. 0.) pgrddyi(i) = 1./pgrddy(i)
        if (pgrddzi(i) .eq. 0.) pgrddzi(i) = 1./pgrddz(i)
      enddo

c     --- Calculate sines and cosines of pgrdph for efficiency
      do i=0,npgrd
        pgrdsp(i) = sin(pgrdph(i)) 
        pgrdcp(i) = cos(pgrdph(i)) 
      enddo 

c     --- Precalculate the axial derivatives of the multipole moments
c     --- if they were not supplied by the user.

c     --- Electric moments first.
      do i=1,nemltsets
        dzio2 = .5/dzemlt(i)
        do ie=1,nesmult
c         --- Axial derivative of the moment strengths
c         --- Find the maximum value.
          emax = 0.
          do iz=0,nzemltmax
            if (abs(esemltp(iz,ie,i)) .gt. emax) emax = abs(esemltp(iz,ie,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (emax .eq. 0.) then
            esemltp(0,ie,i) = 0.
            do iz=1,nzemltmax-1
              esemltp(iz,ie,i) = (esemlt(iz+1,ie,i) - esemlt(iz-1,ie,i))*dzio2
            enddo
            esemltp(nzemltmax,ie,i) = 0.
          endif
c         --- Axial derivative of the phase angle
c         --- Find the maximum value.
          ephmax = 0.
          do iz=0,nzemltmax
            if (abs(esemltphp(iz,ie,i)) .gt. ephmax)
     &        ephmax = abs(esemltphp(iz,ie,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (ephmax .eq. 0.) then
            esemltphp(0,ie,i) = 0.
            do iz=1,nzemltmax-1
              esemltphp(iz,ie,i)=(esemltph(iz+1,ie,i)-esemltph(iz-1,ie,i))*dzio2
            enddo
            esemltphp(nzemltmax,ie,i) = 0.
          endif
        enddo
      enddo
c     --- Magnetic moments second.
      do i=1,nmmltsets
        dzio2 = .5/dzmmlt(i)
        do im=1,nmsmult
c         --- Axial derivative of the moment strengths
c         --- Find the maximum value.
          mmax = 0.
          do iz=0,nzmmltmax
            if (abs(msmmltp(iz,im,i)) .gt. mmax) mmax = abs(msmmltp(iz,im,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (mmax .eq. 0.) then
            msmmltp(0,im,i) = 0.
            do iz=1,nzmmltmax-1
              msmmltp(iz,im,i) = (msmmlt(iz+1,im,i) - msmmlt(iz-1,im,i))*dzio2
            enddo
            msmmltp(nzmmltmax,im,i) = 0.
          endif
c         --- Axial derivative of the phase angle
c         --- Find the maximum value.
          mphmax = 0.
          do iz=0,nzmmltmax
            if (abs(msmmltphp(iz,im,i)) .gt. mphmax)
     &        mphmax = abs(msmmltphp(iz,im,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (mphmax .eq. 0.) then
            msmmltphp(0,im,i) = 0.
            do iz=1,nzmmltmax-1
              msmmltphp(iz,im,i)=(msmmltph(iz+1,im,i)-msmmltph(iz-1,im,i))*dzio2
            enddo
            msmmltphp(nzmmltmax,im,i) = 0.
          endif
        enddo
      enddo

      return
      end
c============================================================================ 
      subroutine setlatt
      use InGen
      use Picglb
c  Sets the lattice data for the current beam location.
c  These ride with the beam, just as the self-field does, and
c  should be loaded at the same time.
      call setlattzt(zbeam,time,fstype)
      return
      end
c============================================================================ 
      subroutine setlattzt(zbeam,time,fstype)
      use Lattice
      use LatticeInternal
      real(kind=8):: zbeam,time
      integer     :: fstype

c  Sets the lattice data at zbeam and time into 1d arrays.


c  If the transition from one element to the next falls to left of cell 
c  center, associate the "right" element with the cell; here, idipo(j)=i+1.

c                                                            (dipole elements
c                                                             begin and end
c            |<- cell j-1 ->|<-- cell j -->|<- cell j+1 ->|   at midpoints of 
c                                                             dipo-drifts)
c        |                      |                                |
c        |<----- element i ---->|<-------- element i+1 --------->|
c        |                      |                                |
c |             |xxxxxxx|               |xxxxxxxxxxxxxxx|                 |
c | dipo-drift  |dipo i |  dipo-drift   |   dipo i+1    |   dipo-drift    |
c |             |xxxxxxx|               |xxxxxxxxxxxxxxx|                 |
c              /         \             /                 \              z --->
c          dipozs(i)  dipoze(i)    dipozs(i+1)         dipoze(i+1)     
c       cdipozs(j-1)  cdipoze(j-1) cdipozs(j)          cdipoze(j)
c                                  cdipozs(j+1)        cdipoze(j+1)

      integer     :: iz,iquadert,iquaderu,iquaderr,ii,imultert,imulteru,imulterr
      integer     :: id
      real(kind=8):: zcellc,offset,wi

c     --- Save the z location and time at which the internal lattice frame
c     --- is set.
      zlframe = zbeam
      zltime = time

c     --- Assume at the start that there are not any of the following elements
c     --- in the mesh.
      linbend = .false.
      linemlt = .false.
      linmmlt = .false.
      linmmlt2 = .false. 
      linbgrd = .false.
      linbgrd2 = .false. 
      linpgrd = .false.

      do iz = 0, nzl
c        --- compute z of center of cell (may be off lattice)
c        --- If nzl is zero, then find the elements at zbeam. Otherwise,
c        --- find the elements at each of the grid points.
         if (nzl > 0) then
           zcellc = iz*dzl + zlmin + zbeam + 0.5*dzl
         else
           zcellc = zbeam
         endif

c        --- find index of nearest drft element
         if (drfts) then
           call getelemid(zcellc,offset,ndrft,drftzs,drftze,cdrftid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cdrftzs(iz) = drftzs(cdrftid(iz)) + offset
           cdrftze(iz) = drftze(cdrftid(iz)) + offset
         endif

c        --- find index of nearest bend element
         if (bends) then
           call getelemid(zcellc,offset,nbend,bendzs,bendze,cbendid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cbendzs(iz) = bendzs(cbendid(iz)) + offset
           cbendze(iz) = bendze(cbendid(iz)) + offset
           cbendrc(iz) = bendrc(cbendid(iz))
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cbendzs(iz) <= zcellc+dzl .and. zcellc-dzl <= cbendze(iz))
     &         linbend = .true.
         endif

c        --- find index of nearest dipole element
c        --- The cdipo arrays calculated here are needed for residence
c        --- corrections (when dipos is true), and internal conductors (when
c        --- fstype is 3).
         if (dipos .or. fstype .eq. 3) then
           call getelemid(zcellc,offset,ndipo,dipozs,dipoze,cdipoid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cdipozs(iz) = dipozs(cdipoid(iz)) + offset
           cdipoze(iz) = dipoze(cdipoid(iz)) + offset
           cdipoby(iz) = dipoby(cdipoid(iz))
           cdipobx(iz) = dipobx(cdipoid(iz))
           cdipota(iz) = dipota(cdipoid(iz))
           cdipotb(iz) = dipotb(cdipoid(iz))
           cdipoex(iz) = dipoex(cdipoid(iz))
           cdipoey(iz) = dipoey(cdipoid(iz))
         endif

c        --- find index of nearest quad element
c        --- The cquad arrays calculated here are needed for residence
c        --- corrections (when quads is true) and the internal conductors (when
c        --- fstype is either 2 or 3).
         if (quads .or. fstype .eq. 2 .or. fstype .eq. 3) then
           call getelemid(zcellc,offset,nquad,quadzs,quadze,cquadid(iz))
           id = cquadid(iz)
c          --- load lattice info into comoving 1d array for easy accss by ptcls
           cquadzs(iz) = quadzs(id) + offset
           cquadze(iz) = quadze(id) + offset
           cquaddb(iz) = quaddb(id)
           cquadde(iz) = quadde(id)
           cquadvx(iz) = quadvx(id)
           cquadvy(iz) = quadvy(id)
           if (quadts(id)<=time .and. time<quadts(id)+ntquad*quaddt(id)) then
             ii = int((time - quadts(id))/quaddt(id))
             wi = (time - quadts(id))/quaddt(id) - ii
             cquadde(iz) = cquadde(iz) + quadet(ii  ,id)*(1. - wi) +
     &                                   quadet(ii+1,id)*      wi 
             cquaddb(iz) = cquaddb(iz) + quadbt(ii  ,id)*(1. - wi) +
     &                                   quadbt(ii+1,id)*      wi 
           endif
c          --- set index for offset quads.  NOTE that this uses "nqerr"
c          --- and not "nquad".  This gives us a choice of either:
c          --- periodic errors (recirculator or ring) via nqerr=nquad, or
c          --- aperiodic errors (straight lattice w/ repeat).
c          --- THIS NEEDS CHECKING, especially when iqerr is nonzero.
           if (zlatperi > 0.) then
             iquadert = id + (int(10.+(zcellc-zlatstrt)/zlatperi)-10)*nquad
           else
             iquadert = id
           endif
           iquaderu = max(iquadert, 0)
           iquaderr = mod(iquaderu, nqerr) 
           cqoffx(iz) = qoffx(iquaderr)
           cqoffy(iz) = qoffy(iquaderr)
         endif

c        --- find index of nearest sext element
c        --- The csext arrays calculated here are needed for residence
c        --- corrections (when sexts is true).
         if (sexts) then
           call getelemid(zcellc,offset,nsext,sextzs,sextze,csextid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls
           csextzs(iz) = sextzs(csextid(iz)) + offset
           csextze(iz) = sextze(csextid(iz)) + offset
           csextdb(iz) = sextdb(csextid(iz))
           csextde(iz) = sextde(csextid(iz))
         endif

c        --- find index of nearest hard-edge multipole element
         if (heles) then
           call getelemid(zcellc,offset,nhele,helezs,heleze,cheleid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls
           chelezs(iz) = helezs(cheleid(iz)) + offset
           cheleze(iz) = heleze(cheleid(iz)) + offset
         endif

c        --- find index of nearest acceleration element
         if (accls) then
           call getelemid(zcellc,offset,naccl,acclzs,acclze,cacclid(iz))
           id = cacclid(iz)
           if (acclze(id) + offset .gt. acclzstt) then
             cacclzs(iz) = acclzs(id) + offset
             cacclze(iz) = acclze(id) + offset
             cacclez(iz) = acclez(id)
             cacclxw(iz) = acclxw(id) 
             cacclsw(iz) = acclsw(id) 
             if (acclts(id) <= time .and.
     &                         time < acclts(id)+ntaccl*accldt(id)) then
               ii = int((time - acclts(id))/accldt(id))
               wi = (time - acclts(id))/accldt(id) - ii
               cacclez(iz) = cacclez(iz) + acclet(ii  ,id)*(1. - wi) +
     &                                     acclet(ii+1,id)*      wi 
             endif
           else
             cacclzs(iz) = 0.
             cacclze(iz) = 0.
             cacclez(iz) = 0.
             cacclxw(iz) = 0.
             cacclsw(iz) = 0.  
           endif
         endif

c        --- find index of nearest emlt element
c        --- The cemlt arrays calculated here are needed for the elements
c        --- described in terms of their multipole components.
         if (emlts) then
           call getelemid(zcellc,offset,nemlt,emltzs,emltze,cemltid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cemltzs(iz) = emltzs(cemltid(iz)) + offset
           cemltze(iz) = emltze(cemltid(iz)) + offset
           cemltph(iz) = emltph(cemltid(iz))
           cemltsf(iz) = emltsf(cemltid(iz))
           cemltsc(iz) = emltsc(cemltid(iz))
           cemltim(iz) = emltid(cemltid(iz))
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cemltzs(iz) <= zcellc+dzl .and. zcellc-dzl <= cemltze(iz))
     &         linemlt = .true.
c          --- set index for offset emlts.  NOTE that this uses "neerr"
c          --- and not "nemlt".  This gives us a choice of either:
c          --- periodic errors (recirculator or ring) via neerr=nemlt, or
c          --- aperiodic errors (straight lattice w/ repeat).
c          --- THIS NEEDS CHECKING
           if (zlatperi > 0.) then
             imultert=cemltid(iz)+(int(10.+(zcellc-zlatstrt)/zlatperi)-10)*nemlt
           else
             imultert = cemltid(iz)
           endif
           imulteru = max(imultert, 0)
           imulterr = mod(imulteru, neerr)
           cemltox(iz) = emltox(imulterr)
           cemltoy(iz) = emltoy(imulterr)

         endif

c        --- find index of nearest mmlt element
c        --- The cmmlt arrays calculated here are needed for the elements
c        --- described in terms of their multipole components.
         if (mmlts) then
           call getelemid(zcellc,offset,nmmlt,mmltzs,mmltze,cmmltid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cmmltzs(iz) = mmltzs(cmmltid(iz)) + offset
           cmmltze(iz) = mmltze(cmmltid(iz)) + offset
           cmmltph(iz) = mmltph(cmmltid(iz))
           cmmltsf(iz) = mmltsf(cmmltid(iz))
           cmmltsc(iz) = mmltsc(cmmltid(iz))
           cmmltim(iz) = mmltid(cmmltid(iz))
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cmmltzs(iz) <= zcellc+dzl .and. zcellc-dzl <= cmmltze(iz))
     &         linmmlt = .true.
c          --- set index for offset mmlts.  NOTE that this uses "nmerr"
c          --- and not "nmmlt".  This gives us a choice of either:
c          --- periodic errors (recirculator or ring) via nmerr=nmmlt, or
c          --- aperiodic errors (straight lattice w/ repeat).
c          --- THIS NEEDS CHECKING
           if (zlatperi > 0.) then
             imultert=cmmltid(iz)+(int(10.+(zcellc-zlatstrt)/zlatperi)-10)*nmmlt
           else
             imultert = cmmltid(iz)
           endif
           imulteru = max(imultert, 0)
           imulterr = mod(imulteru, nmerr)
           cmmltox(iz) = mmltox(imulterr)
           cmmltoy(iz) = mmltoy(imulterr)
         endif

c        --- find index of nearest mmlt2 element
c        --- The cmmlt2 arrays calculated here are needed for the elements
c        --- described in terms of their multipole components.
         if (mmlt2s) then
           call getelemid(zcellc,offset,nmmlt2,mmlt2zs,mmlt2ze,cmmlt2id(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cmmlt2zs(iz) = mmlt2zs(cmmlt2id(iz)) + offset
           cmmlt2ze(iz) = mmlt2ze(cmmlt2id(iz)) + offset
           cmmlt2ph(iz) = mmlt2ph(cmmlt2id(iz))
           cmmlt2sf(iz) = mmlt2sf(cmmlt2id(iz))
           cmmlt2sc(iz) = mmlt2sc(cmmlt2id(iz))
           cmmlt2im(iz) = mmlt2id(cmmlt2id(iz))
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cmmlt2zs(iz) <= zcellc+dzl .and. zcellc-dzl <= cmmlt2ze(iz))
     &         linmmlt2 = .true.
c          --- set index for offset mmlt2s.  NOTE that this uses "nmerr"
c          --- and not "nmmlt2".  This gives us a choice of either:
c          --- periodic errors (recirculator or ring) via nmerr=nmmlt2, or
c          --- aperiodic errors (straight lattice w/ repeat).
c          --- THIS NEEDS CHECKING
           if (zlatperi > 0.) then
             imultert=cmmlt2id(iz)+(int(10.+(zcellc-zlatstrt)/zlatperi)-10)*nmmlt2
           else
             imultert = cmmlt2id(iz)
           endif
           imulteru = max(imultert, 0)
           imulterr = mod(imulteru, nmerr)
           cmmlt2ox(iz) = mmlt2ox(imulterr)
           cmmlt2oy(iz) = mmlt2oy(imulterr)
         endif

c        --- find index of nearest bgrd element
c        --- The cbgrdid array calculated here is needed for the magnetic
c        --- elements described by B field data on a 3-D grid.
         if (bgrds) then
           call getelemid(zcellc,offset,nbgrd,bgrdzs,bgrdze,cbgrdid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cbgrdzs(iz) = bgrdzs(cbgrdid(iz)) + offset
           cbgrdze(iz) = bgrdze(cbgrdid(iz)) + offset
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cbgrdzs(iz) <= zcellc+dzl .and. zcellc-dzl <= cbgrdze(iz))
     &         linbgrd = .true.
         endif

c        --- find index of nearest bgrd2 element
c        --- The cbgrd2id array calculated here is needed for the magnetic
c        --- elements described by B field data on a 3-D grid.
         if (bgrd2s) then
           call getelemid(zcellc,offset,nbgrd2,bgrd2zs,bgrd2ze,cbgrd2id(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cbgrd2zs(iz) = bgrd2zs(cbgrd2id(iz)) + offset
           cbgrd2ze(iz) = bgrd2ze(cbgrd2id(iz)) + offset
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cbgrd2zs(iz) <= zcellc+dzl .and. zcellc-dzl <= cbgrd2ze(iz))
     &         linbgrd2 = .true.
         endif

c        --- find index of nearest pgrd element
c        --- The cpgrdid array calculated here is needed for the electrostatic
c        --- elements described by potential data on a 3-D grid.
         if (pgrds) then
           call getelemid(zcellc,offset,npgrd,pgrdzs,pgrdze,cpgrdid(iz))
c          --- load lattice info into comoving 1d array for easy accss by ptcls.
           cpgrdzs(iz) = pgrdzs(cpgrdid(iz)) + offset
           cpgrdze(iz) = pgrdze(cpgrdid(iz)) + offset
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cpgrdzs(iz) <= zcellc+dzl .and. zcellc-dzl <= cpgrdze(iz))
     &         linpgrd = .true.
         endif

      enddo

c     --- All parallel slave processors must know if any processors
c     --- are in a bend.  There is probably a better way of doing this.
#ifdef PARALLEL
      if (bends) then
        wi = 0.
        if (linbend) wi = 1.
        call parallelmax(wi,1)
        if (wi .gt. 0.5) linbend = .true.
      endif
#endif

      return
      end
c=============================================================================
      subroutine getelemid(z,offset,nelem,elemzs,elemze,id)
      use Lattice
      integer     :: nelem,id
      real(kind=8):: z,offset,elemzs(0:nelem),elemze(0:nelem)

c Given information about the starts and stops of lattice elements, this
c routine finds the number of the element closest to the specified z value.
c It assumes that some previous value of the element number is passed in
c through id. This is used as a starting point for the search.

c  If the transition from one element to the next falls to left of z,
c  associate the "right" element with z; here, ielem(j)=i+1.

c                                                            (elements
c                                                             begin and end
c            |<- cell j-1 ->|<-- cell j -->|<- cell j+1 ->|   at midpoints of 
c                                                             elem-drifts)
c        |                      |                                |
c        |<----- element i ---->|<-------- element i+1 --------->|
c        |                      |                                |
c |             |xxxxxxx|               |xxxxxxxxxxxxxxx|                 |
c | elem-drift  |elem i |  elem-drift   |   elem i+1    |   elem-drift    |
c |             |xxxxxxx|               |xxxxxxxxxxxxxxx|                 |
c              /         \             /                 \              z --->
c          elemzs(i)  elemze(i)    elemzs(i+1)         elemze(i+1)     
c       celemzs(j-1)  celemze(j-1) celemzs(j)          celemze(j)
c                                  celemzs(j+1)        celemze(j+1)

c The Lattice is only needed for zlatstrt and zlatperi.

c --- Offset of the lattice, for periodic repeat.
      if (zlatperi > 0.) then
        offset = (int(10.+(z-zlatstrt)/zlatperi)-10)*zlatperi+zlatstrt
      else
        offset = zlatstrt
      endif

c --- Check if there are more than one element. If not, return id=0.
      if (nelem == 0) then
        id = 0
        return
      endif

c --- Check if z is somehow greater than the current element (id) location.
c --- If so, reset it to 0 and start from the beginning. This happens if
c --- during a periodic repeat of the lattice or if setlatt was called out
c --- of order. Also, make sure that id is within the proper bounds.
      if (id > 0) then
        if (id > nelem) then
          id = 0
        else if (z < 0.5*(elemze(id-1) + elemzs(id)) + offset) then
          id = 0
        endif
      else if (id < 0) then
        id = 0
      endif

c --- Scan sequentially through the list until the element is found.
c --- Note that there is a potential memory access error here when
c --- id=nelem, but the loop will exit anyway.
      do while (id < nelem .and. z > 0.5*(elemze(id) + elemzs(id+1)) + offset)
        id = id + 1
      end do

      return
      end
c=============================================================================
      subroutine applyaccl(np,xp,zp,uzp,gaminv,dtl,dtr,dt,qoverm,ez)
      use Lattice
      use LatticeInternal
      integer     :: np
      real(kind=8):: xp(np),zp(np),uzp(np),gaminv(np)
      real(kind=8):: dtl,dtr,dt,qoverm
      real(kind=8):: ez(np)

c Apply the accelerating gap element accl.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c   qoverm   charge over mass
c Output:
c   ez       axial electric field


      real(kind=8):: dti,vz,gapez,frac,z1,z2,vn,zl,zr,cacclz
      integer     :: ip,j

      if (.not. accls) return

      dti = 1./(dtr-dtl)

c     --- Finite length gaps
      if (.not. lacclzl) then

c       --- For velocity correction, first calculate the velocity at
c       --- time level n.
        do ip=1,np
          vz = uzp(ip)*gaminv(ip)
          if (vz .eq. 0.) vz = LARGEPOS
c         --- find z-cell in which particle lies
          j = int(max(0., (zp(ip) - zlmin - zlframe) * dzli + 0.5))
          gapez = cacclez(j) + cacclxw(j)*xp(ip)
          z1 = zp(ip) + dtl*vz - 0.5*gapez*qoverm*dtl**2
          z2 = zp(ip) + dtl*vz
          if (zp(ip) .le. cacclzs(j)) then
            vn = vz
          elseif (z1 .lt. cacclzs(j)) then
            vn = sqrt(vz**2 + 2.*gapez*qoverm*(zp(ip) - cacclzs(j)))
          elseif (zp(ip) .le. cacclze(j)) then
            vn = vz - gapez*qoverm*dtl
          elseif (z2 .lt. cacclze(j)) then
            vn = 0.5*(vz + 0.5*gapez*qoverm*dt + sqrt((vz +
     &           0.5*gapez*qoverm*dt)**2 -
     &           4.*gapez*qoverm*(zp(ip) - cacclze(j))))
          else
            vn = vz
          endif
c         --- Calculate the fraction of time in the gap.  Cases inside
c         --- and outside are included implicitly in the max and min calls.
c         --- Note that the max's inside the sqrt are for idiot proofing.
          if (zp(ip) .le. cacclzs(j)) then
            frac = max((dtr + (zp(ip) - cacclzs(j))/vn)*dti , 0.)
          elseif (zp(ip) .le. (cacclzs(j) + cacclze(j))*0.5) then
            frac = min((dtr + 2.*(zp(ip) - cacclzs(j))/
     &             (sqrt(max(0.,vn**2 - 2.*gapez*qoverm*
     &             (zp(ip) - cacclzs(j)))) + vn))*dti, 1.)
          elseif (zp(ip) .le. cacclze(j)) then
            frac = min((-dtl + 2.*(cacclze(j) - zp(ip))/
     &             (sqrt(max(0.,vn**2 + 2.*gapez*qoverm*
     &             (cacclze(j) - zp(ip)))) + vn))*dti, 1.)
          else
            frac = max((-dtl - (zp(ip) - cacclze(j))/vn)*dti , 0.)
          endif
c         --- add acceleration field to Ez field
          ez(ip) = ez(ip) + gapez*frac
        enddo

      else

c       --- Zero length gaps
        do ip = 1, np
          vz = uzp(ip)*gaminv(ip)
          j = max(0., (zp(ip) - zlmin - zlframe) * dzli + 0.5)
c         --- "left" end of velocity advance step
          zl = zp(ip) + vz*dtl
c         --- "right" end of velocity advance step
          zr = zp(ip) + vz*dtr
c         --- Gap center
          cacclz = 0.5*(cacclzs(j) + cacclze(j))
c         --- Add acceleration to velocity.
c         --- Calculate the change in velocity and convert that into
c         --- and Ez - the expression reduces to q*V/m for small V.
          if (zl .le. cacclz .and. cacclz .lt. zr) then
            gapez = cacclez(j) + cacclxw(j)*xp(ip)
            ez(ip) = ez(ip) + (sqrt(vz**2 +
     &          2.*qoverm*gapez*(cacclze(j)-cacclzs(j))) - vz)/qoverm*dti
          endif
        enddo

      endif

      return
      end
c=============================================================================
      subroutine getzmmnt(np,xp,yp,zp,uxp,uyp,uzp,gaminv,q,m,w,dt,itask,nplive,
     &                    uxpo,uypo,uzpo,is,ns)
      use Constant
      use Beam_acc
      use InDiag
      use Z_Moments
      use Win_Moments
      use Moments
      use Picglb
      use ExtPart
      use Timers
      integer     :: np,itask,nplive,is,ns
      real(kind=8):: q,m,w,dt
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np)

c  Sets moments for species 1 in as a function of z.
c  Interpolation is done to grid centers
c  Note that dead particles are to the left of the z grid.
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


      real(kind=8):: pp(0:16,NPARPGRP)
      integer     :: iz(NPARPGRP)
      integer     :: im,izm,iw,ip,icell,i
      real(kind=8):: dti,dtip1,wz1,wz0,zwin,tg,th
      real(kind=8):: oneondt,clighti,vzi,pnumi,vbeami
      real(kind=8):: delxsq,delxxp,delxpsq,delysq,delyyp,delypsq
      real(kind=8):: delxy,delxyp,delyxp,delxpyp
      real(kind=8):: delzsq,delvzsq,delzvz
      real(kind=8):: delvxsq,delvysq
      real(kind=8):: gamma
      real(kind=8):: timetemp,wtime
#if SYSTEM == J90
      real(kind=8):: pmmnts(NUMZMMNT,NPARPGRP)
#endif

      timetemp = wtime()

      if (ifzmmnt .eq. 0) return

      oneondt = 1./dvnz(dt)

      if (itask .eq. 1 .and. .not. laccumulate_zmoments) then
         xmaxp = -LARGEPOS
         xminp = LARGEPOS
         ymaxp = -LARGEPOS
         yminp = LARGEPOS
         zmaxp = -LARGEPOS
         zminp = LARGEPOS
         vxmaxp = -LARGEPOS
         vxminp = LARGEPOS
         vymaxp = -LARGEPOS
         vyminp = LARGEPOS
         vzmaxp = -LARGEPOS
         vzminp = LARGEPOS
         do im=1,NUMZMMNT
           zmmnts0(im) = 0.
         enddo
         do im=1,NUMZMMNT
           do izm=0,nzmmnt
             zmmnts(izm,im) = 0.
           enddo
         enddo
         do im=1,ns
           do iw=1,nepwin
             nep(iw,im) = 0
           enddo
         enddo
      endif

      if (itask .eq. 2) then
c       --- Set maximum and minimum particle coordinate variables
c       --- set multiplier so dead particles don't contribute
        do ip = 1, np
          if (uzp(ip) .ne. 0.) then
            xmaxp  = max(xmaxp, xp(ip))
            xminp  = min(xminp, xp(ip))
            ymaxp  = max(ymaxp, yp(ip))
            yminp  = min(yminp, yp(ip))
            zmaxp  = max(zmaxp, zp(ip))
            zminp  = min(zminp, zp(ip))
            vxmaxp = max(vxmaxp, gaminv(ip)*uxp(ip))
            vxminp = min(vxminp, gaminv(ip)*uxp(ip))
            vymaxp = max(vymaxp, gaminv(ip)*uyp(ip))
            vyminp = min(vyminp, gaminv(ip)*uyp(ip))
            vzmaxp = max(vzmaxp, gaminv(ip)*uzp(ip))
            vzminp = min(vzminp, gaminv(ip)*uzp(ip))

            pp(0,ip) = 1.
            if (ifzmmnt .eq. 1) then
              if ((zp(ip) - zbeam) .lt. zmmntmin .or.
     &            (zp(ip) - zbeam) .gt. zmmntmax) pp(0,ip) = 0.
            else
              if ((zp(ip) - zbeam) .lt. zmmntmin .or.
     &            (zp(ip) - zbeam) .ge. zmmntmax) pp(0,ip) = 0.
            endif

          else
            pp(0,ip) = 0.
          endif
        enddo

c       --- set multiplier so dead particles don't contribute
c       ---   not added to total moments
c       ---   weights set to zero for z moments
c       if (ifzmmnt .eq. 1) then
c         do ip=1,np
c           pp(0,ip) = 1.
c           if (uzp(ip) .eq. 0. .or.
c    &          (zp(ip) - zbeam) .lt. zmmntmin .or.
c    &          (zp(ip) - zbeam) .gt. zmmntmax) pp(0,ip) = 0.
c         enddo
c       else
c         do ip=1,np
c           pp(0,ip) = 1.
c           if (uzp(ip) .eq. 0. .or.
c    &          (zp(ip) - zbeam) .lt. zmmntmin .or.
c    &          (zp(ip) - zbeam) .ge. zmmntmax) pp(0,ip) = 0.
c         enddo
c       endif

c       --- Accumate global moments first.

c-------------------------------------------------------------------------
#if SYSTEM == J90
c       --- This method is faster for vectorized machines, primarily the J90.
c       --- put particle moments into pmmnts array for total moments
        do ip=1,np
          vzi = 1./(uzp(ip)+SMALLPOS)
          if (pp(0,ip) > 1.) then
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
c         do im=1,NUMZMMNT
c           zmmnts0(im) = zmmnts0(im) + pmmnts(ip,im)
c         enddo
          zmmnts0(01) = zmmnts0(01) + pmmnts(01,ip)
          zmmnts0(02) = zmmnts0(02) + pmmnts(02,ip)
          zmmnts0(03) = zmmnts0(03) + pmmnts(03,ip)
          zmmnts0(04) = zmmnts0(04) + pmmnts(04,ip)
          zmmnts0(05) = zmmnts0(05) + pmmnts(05,ip)
          zmmnts0(06) = zmmnts0(06) + pmmnts(06,ip)
          zmmnts0(07) = zmmnts0(07) + pmmnts(07,ip)
          zmmnts0(08) = zmmnts0(08) + pmmnts(08,ip)
          zmmnts0(09) = zmmnts0(09) + pmmnts(09,ip)
          zmmnts0(10) = zmmnts0(10) + pmmnts(10,ip)
          zmmnts0(11) = zmmnts0(11) + pmmnts(11,ip)
          zmmnts0(12) = zmmnts0(12) + pmmnts(12,ip)
          zmmnts0(13) = zmmnts0(13) + pmmnts(13,ip)
          zmmnts0(14) = zmmnts0(14) + pmmnts(14,ip)
          zmmnts0(15) = zmmnts0(15) + pmmnts(15,ip)
          zmmnts0(16) = zmmnts0(16) + pmmnts(16,ip)
          zmmnts0(17) = zmmnts0(17) + pmmnts(17,ip)
          zmmnts0(18) = zmmnts0(18) + pmmnts(18,ip)
          zmmnts0(19) = zmmnts0(19) + pmmnts(19,ip)
          zmmnts0(20) = zmmnts0(20) + pmmnts(20,ip)
          zmmnts0(21) = zmmnts0(21) + pmmnts(21,ip)
          zmmnts0(22) = zmmnts0(22) + pmmnts(22,ip)
          zmmnts0(23) = zmmnts0(23) + pmmnts(23,ip) 
          zmmnts0(24) = zmmnts0(24) + pmmnts(24,ip) 
          zmmnts0(25) = zmmnts0(25) + pmmnts(25,ip) 
          zmmnts0(26) = zmmnts0(26) + pmmnts(26,ip) 
          zmmnts0(27) = zmmnts0(27) + pmmnts(27,ip) 
          zmmnts0(28) = zmmnts0(28) + pmmnts(28,ip) 
        enddo

#else

c-------------------------------------------------------------------------
c       --- This method is faster for RISC and pentium machines.
c       --- sum particle moments into total moments
        do ip=1,np
          vzi = 1./(uzp(ip)+SMALLPOS)
          if (pp(0,ip) > 0) then
            zmmnts0(01) = zmmnts0(01) + 1.
            zmmnts0(02) = zmmnts0(02) + xp(ip)
            zmmnts0(03) = zmmnts0(03) + yp(ip)
            zmmnts0(04) = zmmnts0(04) + zp(ip)
            zmmnts0(07) = zmmnts0(07) + uxp(ip)*gaminv(ip)
            zmmnts0(08) = zmmnts0(08) + uyp(ip)*gaminv(ip)
            zmmnts0(09) = zmmnts0(09) + uzp(ip)*gaminv(ip)
            zmmnts0(10) = zmmnts0(10) + xp(ip)*yp(ip)
            zmmnts0(12) = zmmnts0(12) + xp(ip)**2
            zmmnts0(13) = zmmnts0(13) + yp(ip)**2
            zmmnts0(14) = zmmnts0(14) + zp(ip)**2
            zmmnts0(17) = zmmnts0(17) + (gaminv(ip)*uxp(ip))**2
            zmmnts0(18) = zmmnts0(18) + (gaminv(ip)*uyp(ip))**2
            zmmnts0(19) = zmmnts0(19) + (gaminv(ip)*uzp(ip))**2
            zmmnts0(22) = zmmnts0(22) + zp(ip)*gaminv(ip)*uzp(ip)
            zmmnts0(25) = zmmnts0(25) + xp(ip)*uzp(ip)*gaminv(ip)
            zmmnts0(26) = zmmnts0(26) + yp(ip)*uzp(ip)*gaminv(ip)
            zmmnts0(27) = zmmnts0(27) + uxp(ip)*uzp(ip)*gaminv(ip)**2
            zmmnts0(28) = zmmnts0(28) + uyp(ip)*uzp(ip)*gaminv(ip)**2
            zmmnts0(05) = zmmnts0(05) + uxp(ip)*vzi
            zmmnts0(06) = zmmnts0(06) + uyp(ip)*vzi
            zmmnts0(11) = zmmnts0(11) + uxp(ip)*uyp(ip)*vzi**2
            zmmnts0(15) = zmmnts0(15) + (uxp(ip)*vzi)**2
            zmmnts0(16) = zmmnts0(16) + (uyp(ip)*vzi)**2
            zmmnts0(20) = zmmnts0(20) + xp(ip)*uxp(ip)*vzi
            zmmnts0(21) = zmmnts0(21) + yp(ip)*uyp(ip)*vzi
            zmmnts0(23) = zmmnts0(23) + xp(ip)*uyp(ip)*vzi
            zmmnts0(24) = zmmnts0(24) + yp(ip)*uxp(ip)*vzi
          endif
        enddo

c-------------------------------------------------------------------------
#endif


c       --- Now, accumulate moments onto Z grid.
        if (ifzmmnt .eq. 2) then

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
            iz(ip)  = (zp(ip) - zbeam - zmmntmin)*dzmi*pp(0,ip)
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
            pp(14,ip) = (1. - (zp(ip)-zbeam-zmmntmin)*dzmi + iz(ip))*pp(0,ip)
            pp(15,ip) = ((zp(ip) - zbeam - zmmntmin)*dzmi - iz(ip))*pp(0,ip)
          enddo

c         --- Add particles to ExtPart arrays (skips if nepwin = 0).
c         --- This allows access to the particle data when it is extrapolated
c         --- to the moments grid cell centers.
          do iw=1,nepwin
            do ip=1,np
              if (iz(ip) .eq. izepwin(iw) .and. nep(iw,is) .lt. nepmax .and.
     &            pp(0,ip) > 0.) then
                nep(iw,is) = nep(iw,is) + 1
                xep(nep(iw,is),iw,is) = pp(1,ip)
                yep(nep(iw,is),iw,is) = pp(5,ip)
                uxep(nep(iw,is),iw,is) = pp(3,ip)
                uyep(nep(iw,is),iw,is) = pp(7,ip)
                uzep(nep(iw,is),iw,is) = pp(11,ip)
              else if (iz(ip)+1 .eq. izepwin(iw) .and.
     &                 nep(iw,is) .lt. nepmax .and.
     &                 pp(0,ip) > 0.) then
                nep(iw,is) = nep(iw,is) + 1
                xep(nep(iw,is),iw,is) = pp(2,ip)
                yep(nep(iw,is),iw,is) = pp(6,ip)
                uxep(nep(iw,is),iw,is) = pp(4,ip)
                uyep(nep(iw,is),iw,is) = pp(8,ip)
                uzep(nep(iw,is),iw,is) = pp(12,ip)
              endif
            enddo
          enddo

c-------------------------------------------------------------------------
#if SYSTEM == J90
c       --- This method is faster for vectorized machines, primarily the J90.
c         --- put particle moments into pmmnts array for zmmnts(iz,)
          do ip=1,np
            vzi = 1./(pp(11,ip)+SMALLPOS)
            wz0 = pp(14,ip)
            if (ifzmmnt .eq. 1) wz0 = pp(0,ip)
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
              zmmnts(iz(ip),im) = zmmnts(iz(ip),im) + pmmnts(im,ip)
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
              zmmnts(iz(ip)+1,im) = zmmnts(iz(ip)+1,im) + pmmnts(im,ip)
            enddo
          enddo

#else

c-------------------------------------------------------------------------
c       --- This method is faster for RISC and pentium machines.

          do ip=1,np
c           --- sum particle moments into zmmnts(iz,)
            vzi = 1./(pp(11,ip)+SMALLPOS)
            wz0 = pp(14,ip)
            if (ifzmmnt .eq. 1) wz0 = pp(0,ip)
            i = iz(ip)
            zmmnts(i, 1) = zmmnts(i, 1) + wz0
            zmmnts(i, 2) = zmmnts(i, 2) + pp(1,ip)*wz0
            zmmnts(i, 3) = zmmnts(i, 3) + pp(5,ip)*wz0
            zmmnts(i, 4) = zmmnts(i, 4) + pp(9,ip)*wz0
            zmmnts(i, 7) = zmmnts(i, 7) + pp(3,ip)*pp(13,ip)*wz0
            zmmnts(i, 8) = zmmnts(i, 8) + pp(7,ip)*pp(13,ip)*wz0
            zmmnts(i, 9) = zmmnts(i, 9) + pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,10) = zmmnts(i,10) + pp(1,ip)*pp(5,ip)*wz0
            zmmnts(i,12) = zmmnts(i,12) + pp(1,ip)**2*wz0
            zmmnts(i,13) = zmmnts(i,13) + pp(5,ip)**2*wz0
            zmmnts(i,14) = zmmnts(i,14) + pp(9,ip)**2*wz0
            zmmnts(i,17) = zmmnts(i,17) + (pp(13,ip)*pp(3,ip))**2*wz0
            zmmnts(i,18) = zmmnts(i,18) + (pp(13,ip)*pp(7,ip))**2*wz0
            zmmnts(i,19) = zmmnts(i,19) + (pp(13,ip)*pp(11,ip))**2*wz0
            zmmnts(i,22) = zmmnts(i,22) + pp(9,ip)*pp(13,ip)*pp(11,ip)*wz0
            zmmnts(i,25) = zmmnts(i,25) + pp(1,ip)*pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,26) = zmmnts(i,26) + pp(5,ip)*pp(11,ip)*pp(13,ip)*wz0
            zmmnts(i,27) = zmmnts(i,27) + pp(3,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            zmmnts(i,28) = zmmnts(i,28) + pp(7,ip)*pp(11,ip)*pp(13,ip)**2*wz0
            zmmnts(i, 5) = zmmnts(i, 5) + pp(3,ip)*vzi*wz0
            zmmnts(i, 6) = zmmnts(i, 6) + pp(7,ip)*vzi*wz0
            zmmnts(i,11) = zmmnts(i,11) + pp(3,ip)*pp(7,ip)*vzi**2*wz0
            zmmnts(i,15) = zmmnts(i,15) + (pp(3,ip)*vzi)**2*wz0
            zmmnts(i,16) = zmmnts(i,16) + (pp(7,ip)*vzi)**2*wz0
            zmmnts(i,20) = zmmnts(i,20) + pp(1,ip)*pp(3,ip)*vzi*wz0
            zmmnts(i,21) = zmmnts(i,21) + pp(5,ip)*pp(7,ip)*vzi*wz0
            zmmnts(i,23) = zmmnts(i,23) + pp(1,ip)*pp(7,ip)*vzi*wz0
            zmmnts(i,24) = zmmnts(i,24) + pp(5,ip)*pp(3,ip)*vzi*wz0

c           --- sum particle moments into zmmnts(iz+1,)
            vzi = 1./(pp(12,ip)+SMALLPOS)
            wz1 = pp(15,ip)
            i = iz(ip) + 1
            zmmnts(i, 1) = zmmnts(i, 1) + wz1
            zmmnts(i, 2) = zmmnts(i, 2) + pp(2,ip)*wz1
            zmmnts(i, 3) = zmmnts(i, 3) + pp(6,ip)*wz1
            zmmnts(i, 4) = zmmnts(i, 4) + pp(10,ip)*wz1
            zmmnts(i, 7) = zmmnts(i, 7) + pp(4,ip)*pp(13,ip)*wz1
            zmmnts(i, 8) = zmmnts(i, 8) + pp(8,ip)*pp(13,ip)*wz1
            zmmnts(i, 9) = zmmnts(i, 9) + pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,10) = zmmnts(i,10) + pp(2,ip)*pp(6,ip)*wz1
            zmmnts(i,12) = zmmnts(i,12) + pp(2,ip)**2*wz1
            zmmnts(i,13) = zmmnts(i,13) + pp(6,ip)**2*wz1
            zmmnts(i,14) = zmmnts(i,14) + pp(10,ip)**2*wz1
            zmmnts(i,17) = zmmnts(i,17) + (pp(13,ip)*pp(4,ip))**2*wz1
            zmmnts(i,18) = zmmnts(i,18) + (pp(13,ip)*pp(8,ip))**2*wz1
            zmmnts(i,19) = zmmnts(i,19) + (pp(13,ip)*pp(12,ip))**2*wz1
            zmmnts(i,22) = zmmnts(i,22) + pp(10,ip)*pp(13,ip)*pp(12,ip)*wz1
            zmmnts(i,25) = zmmnts(i,25) + pp(2,ip)*pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,26) = zmmnts(i,26) + pp(6,ip)*pp(12,ip)*pp(13,ip)*wz1
            zmmnts(i,27) = zmmnts(i,27) + pp(4,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            zmmnts(i,28) = zmmnts(i,28) + pp(8,ip)*pp(12,ip)*pp(13,ip)**2*wz1
            zmmnts(i, 5) = zmmnts(i, 5) + pp(4,ip)*vzi*wz1
            zmmnts(i, 6) = zmmnts(i, 6) + pp(8,ip)*vzi*wz1
            zmmnts(i,11) = zmmnts(i,11) + pp(4,ip)*pp(8,ip)*vzi**2*wz1
            zmmnts(i,15) = zmmnts(i,15) + (pp(4,ip)*vzi)**2*wz1
            zmmnts(i,16) = zmmnts(i,16) + (pp(8,ip)*vzi)**2*wz1
            zmmnts(i,20) = zmmnts(i,20) + pp(2,ip)*pp(4,ip)*vzi*wz1
            zmmnts(i,21) = zmmnts(i,21) + pp(6,ip)*pp(8,ip)*vzi*wz1
            zmmnts(i,23) = zmmnts(i,23) + pp(2,ip)*pp(8,ip)*vzi*wz1
            zmmnts(i,24) = zmmnts(i,24) + pp(6,ip)*pp(4,ip)*vzi*wz1

          enddo

c-------------------------------------------------------------------------
#endif

c         --- End of Z moments accumulation
        endif

      endif

      if (itask .eq. 3 .and. .not. laccumulate_zmoments) then
        vbeami = 1./dvnz(vbeam)
        clighti = 1./clight
        gamma = 1./sqrt(1. - (vbeam*clighti)**2)

c       --- For slave, call routine which sums moments over processors.
#ifdef PARALLEL
        call sum_mmnts
#endif

c       --- Complete the calculation of moments: divide by particle number

c       --- Global moments first
        pnumi = 1./(zmmnts0( 1)+SMALLPOS)
        pnum(0)    = zmmnts0( 1)
        xbar(0)    = zmmnts0( 2)*pnumi
        ybar(0)    = zmmnts0( 3)*pnumi
        zbar(0)    = zmmnts0( 4)*pnumi
        xpbar(0)   = zmmnts0( 5)*pnumi
        ypbar(0)   = zmmnts0( 6)*pnumi
        vxbar(0)   = zmmnts0( 7)*pnumi
        vybar(0)   = zmmnts0( 8)*pnumi
        vzbar(0)   = zmmnts0( 9)*pnumi
        xybar(0)   = zmmnts0(10)*pnumi
        xpypbar(0) = zmmnts0(11)*pnumi
        xsqbar(0)  = zmmnts0(12)*pnumi
        ysqbar(0)  = zmmnts0(13)*pnumi
        zsqbar(0)  = zmmnts0(14)*pnumi
        xpsqbar(0) = zmmnts0(15)*pnumi
        ypsqbar(0) = zmmnts0(16)*pnumi
        vxsqbar(0) = zmmnts0(17)*pnumi
        vysqbar(0) = zmmnts0(18)*pnumi
        vzsqbar(0) = zmmnts0(19)*pnumi
        xxpbar(0)  = zmmnts0(20)*pnumi
        yypbar(0)  = zmmnts0(21)*pnumi
        zvzbar(0)  = zmmnts0(22)*pnumi
        xypbar(0)  = zmmnts0(23)*pnumi 
        yxpbar(0)  = zmmnts0(24)*pnumi 
        xvzbar(0)  = zmmnts0(25)*pnumi
        yvzbar(0)  = zmmnts0(26)*pnumi
        vxvzbar(0) = zmmnts0(27)*pnumi
        vyvzbar(0) = zmmnts0(28)*pnumi

c       --- Compute second order moments with averages subtracted
        delxsq = xsqbar(0) - xbar(0)*xbar(0)
        delxxp = xxpbar(0) - xbar(0)*xpbar(0)
        delxpsq = xpsqbar(0) - xpbar(0)*xpbar(0)
        delysq = ysqbar(0) - ybar(0)*ybar(0)
        delyyp = yypbar(0) - ybar(0)*ypbar(0)
        delypsq = ypsqbar(0) - ypbar(0)*ypbar(0)
        delxy = xybar(0) - xbar(0)*ybar(0)
        delxyp = xypbar(0) - xbar(0)*ypbar(0)
        delyxp = yxpbar(0) - ybar(0)*xpbar(0)
        delxpyp = xpypbar(0) - xpbar(0)*ypbar(0)

        delvxsq = vxsqbar(0) - vxbar(0)**2
        delvysq = vysqbar(0) - vybar(0)**2

        delzsq = zsqbar(0) - zbar(0)**2
        delvzsq = vzsqbar(0) - vzbar(0)**2
        delzvz = zvzbar(0) - zbar(0)*vzbar(0)

c       --- Compute overall true RMS positions
        xrms(0) = sqrt(max(SMALLPOS,delxsq))
        yrms(0) = sqrt(max(SMALLPOS,delysq))
        zrms(0) = sqrt(max(SMALLPOS,delzsq))

c       --- Compute overall true RMS velocities
        xprms(0) = sqrt(max(SMALLPOS,delxpsq))
        yprms(0) = sqrt(max(SMALLPOS,delypsq))
        vxrms(0) = sqrt(max(SMALLPOS,delvxsq))
        vyrms(0) = sqrt(max(SMALLPOS,delvysq))
        vzrms(0) = sqrt(max(SMALLPOS,delvzsq))

c       --- Compute overall emittances
        epsx(0) = 4.*sqrt(max(SMALLPOS,delxsq*delxpsq - delxxp**2))
        epsy(0) = 4.*sqrt(max(SMALLPOS,delysq*delypsq - delyyp**2))
        epsz(0) = 4.*sqrt(max(SMALLPOS,delzsq*delvzsq - delzvz**2))
        epsnx(0) = epsx(0)*gamma*vzbar(0)*clighti*1.e6
        epsny(0) = epsy(0)*gamma*vzbar(0)*clighti*1.e6
        epsnz(0) = epsz(0)*gamma*vzbar(0)*clighti*1.e6

c       --- Compute generalized emittances
        tg = 0.5*(epsx(0)**2 + epsy(0)**2) + 16*(delxy*delxpyp - delxyp*delyxp)
        epsg(0) = sqrt(max(SMALLPOS,tg))
        th = epsx(0)**2*epsy(0)**2 +
     &       256*((delxy*delxpyp)**2 + (delxyp*delyxp)**2 -
     &       delxsq*delysq*(delxpyp)**2 - delxsq*delypsq*(delyxp)**2 -
     &       delxpsq*delysq*(delxyp)**2 - delxpsq*delypsq*(delxy)**2 -
     &       2*delxy*delxyp*delyxp*delxpyp + 2*delxxp*delypsq*delxy*delyxp -
     &       2*delxxp*delyyp*delxy*delxpyp - 2*delxxp*delyyp*delxyp*delyxp +
     &       2*delxpsq*delyyp*delxy*delxyp + 2*delxsq*delyyp*delyxp*delxpyp +
     &       2*delxxp*delysq*delxpyp*delxyp)
        epsh(0) = sqrt(sqrt(max(SMALLPOS,th)))
        epsng(0) = epsg(0)*gamma*vzbar(0)*clighti*1.e6
        epsnh(0) = epsh(0)*gamma*vzbar(0)*clighti*1.e6

c       ---  Number of live particles (integer)
        nplive = pnum(0) + 0.5

c       ---  Axial z momentum in vbeam frame (classical)
        pz = m*w*pnum(0) * (vzbar(0)-vbeam)

c       ---  Total Z kinetic energy in lab frame minus beam energy (classical)
c       ---  using <vz^2>-vbeam^2;
        ekzmbe = 0.5*m*w*pnum(0) * (vzsqbar(0) - vbeam**2)

c       ---  Z kinetic energy in beam frame (classical), using <(vz-vbeam)^2>;
        ekzbeam = 0.5*m*w*pnum(0)*(vzsqbar(0) - 2.*vzbar(0)*vbeam + vbeam**2)

c       ---  Perp kinetic energy 
        ekperp =  0.5*m*w*pnum(0) * (vxsqbar(0) + vysqbar(0))

c       --- Total kinetic energy in beam frame
        ek = ekzbeam + ekperp

c       --- Now scale the full Z moments by particle number
        if (ifzmmnt .eq. 2) then
          do izm = 0, nzmmnt
            pnumi = 1./(zmmnts(izm, 1)+SMALLPOS)
            pnumz(izm)    = zmmnts(izm, 1)
            xbarz(izm)    = zmmnts(izm, 2)*pnumi
            ybarz(izm)    = zmmnts(izm, 3)*pnumi
            zbarz(izm)    = zmmnts(izm, 4)*pnumi
            xpbarz(izm)   = zmmnts(izm, 5)*pnumi
            ypbarz(izm)   = zmmnts(izm, 6)*pnumi
            vxbarz(izm)   = zmmnts(izm, 7)*pnumi
            vybarz(izm)   = zmmnts(izm, 8)*pnumi
            vzbarz(izm)   = zmmnts(izm, 9)*pnumi
            xybarz(izm)   = zmmnts(izm,10)*pnumi
            xpypbarz(izm) = zmmnts(izm,11)*pnumi
            xsqbarz(izm)  = zmmnts(izm,12)*pnumi
            ysqbarz(izm)  = zmmnts(izm,13)*pnumi
            zsqbarz(izm)  = zmmnts(izm,14)*pnumi
            xpsqbarz(izm) = zmmnts(izm,15)*pnumi
            ypsqbarz(izm) = zmmnts(izm,16)*pnumi
            vxsqbarz(izm) = zmmnts(izm,17)*pnumi
            vysqbarz(izm) = zmmnts(izm,18)*pnumi
            vzsqbarz(izm) = zmmnts(izm,19)*pnumi
            xxpbarz(izm)  = zmmnts(izm,20)*pnumi
            yypbarz(izm)  = zmmnts(izm,21)*pnumi
            zvzbarz(izm)  = zmmnts(izm,22)*pnumi
            xypbarz(izm)  = zmmnts(izm,23)*pnumi 
            yxpbarz(izm)  = zmmnts(izm,24)*pnumi 
            xvzbarz(izm)  = zmmnts(izm,25)*pnumi
            yvzbarz(izm)  = zmmnts(izm,26)*pnumi
            vxvzbarz(izm) = zmmnts(izm,27)*pnumi
            vyvzbarz(izm) = zmmnts(izm,28)*pnumi

c           --- Compute second order moments with averages subtracted
            delxsq = xsqbarz(izm) - xbarz(izm)*xbarz(izm)
            delxxp = xxpbarz(izm) - xbarz(izm)*xpbarz(izm)
            delxpsq = xpsqbarz(izm) - xpbarz(izm)*xpbarz(izm)
            delysq = ysqbarz(izm) - ybarz(izm)*ybarz(izm)
            delyyp = yypbarz(izm) - ybarz(izm)*ypbarz(izm)
            delypsq = ypsqbarz(izm) - ypbarz(izm)*ypbarz(izm)
            delxy = xybarz(izm) - xbarz(izm)*ybarz(izm)
            delxyp = xypbarz(izm) - xbarz(izm)*ypbarz(izm)
            delyxp = yxpbarz(izm) - ybarz(izm)*xpbarz(izm)
            delxpyp = xpypbarz(izm) - xpbarz(izm)*ypbarz(izm)

            delvxsq = vxsqbarz(izm) - vxbarz(izm)**2
            delvysq = vysqbarz(izm) - vybarz(izm)**2

            delzsq = zsqbarz(izm) - zbarz(izm)**2
            delvzsq = vzsqbarz(izm) - vzbarz(izm)**2
            delzvz = zvzbarz(izm) - zbarz(izm)*vzbarz(izm)

c           --- Compute overall true RMS positions
            xrmsz(izm) = sqrt(max(SMALLPOS,delxsq))
            yrmsz(izm) = sqrt(max(SMALLPOS,delysq))
            zrmsz(izm) = sqrt(max(SMALLPOS,delzsq))

c           --- Compute overall true RMS velocities
            xprmsz(izm) = sqrt(max(SMALLPOS,delxpsq))
            yprmsz(izm) = sqrt(max(SMALLPOS,delypsq))
            vxrmsz(izm) = sqrt(max(SMALLPOS,delvxsq))
            vyrmsz(izm) = sqrt(max(SMALLPOS,delvysq))
            vzrmsz(izm) = sqrt(max(SMALLPOS,delvzsq))

c           --- Compute overall emittances
            epsxz(izm) = 4.*sqrt(max(SMALLPOS,delxsq*delxpsq - delxxp**2))
            epsyz(izm) = 4.*sqrt(max(SMALLPOS,delysq*delypsq - delyyp**2))
            epszz(izm) = 4.*sqrt(max(SMALLPOS,delzsq*delvzsq - delzvz**2))
            epsnxz(izm) = epsxz(izm)*gamma*vzbarz(izm)*clighti*1.e6
            epsnyz(izm) = epsyz(izm)*gamma*vzbarz(izm)*clighti*1.e6
            epsnzz(izm) = epszz(izm)*gamma*vzbarz(izm)*clighti*1.e6

c           --- Compute generalized emittances
            tg = 0.5*(epsxz(izm)**2 + epsyz(izm)**2) +
     &           16*(delxy*delxpyp - delxyp*delyxp)
            epsgz(izm) = sqrt(max(SMALLPOS,tg))
            th = epsxz(izm)**2*epsyz(izm)**2 +
     &          256*((delxy*delxpyp)**2 + (delxyp*delyxp)**2 -
     &          delxsq*delysq*(delxpyp)**2 - delxsq*delypsq*(delyxp)**2 -
     &          delxpsq*delysq*(delxyp)**2 - delxpsq*delypsq*(delxy)**2 -
     &          2*delxy*delxyp*delyxp*delxpyp + 2*delxxp*delypsq*delxy*delyxp -
     &          2*delxxp*delyyp*delxy*delxpyp - 2*delxxp*delyyp*delxyp*delyxp +
     &          2*delxpsq*delyyp*delxy*delxyp + 2*delxsq*delyyp*delyxp*delxpyp +
     &          2*delxxp*delysq*delxpyp*delxyp)
            epshz(izm) = sqrt(sqrt(max(SMALLPOS,th)))
            epsngz(izm) = epsgz(izm)*gamma*vzbarz(izm)*clighti*1.e6
            epsnhz(izm) = epshz(izm)*gamma*vzbarz(izm)*clighti*1.e6

          enddo

c         --- Save window moments.  Linearly interpolate the Z moments data
c         --- to the Z window center.
c         --- If the zwindow is outside the zmoments range, don't set anything.
          do iw = 1,nzwind
            zwin = 0.5*(zwindows(1,iw) + zwindows(2,iw))
            if (zmmntmin .le. zwin .and. zwin .lt. zmmntmax) then
              icell = int((zwin - zmmntmin)*dzmi)
              wz0 = (zwin - zmmntmin)*dzmi - icell
              wz1 = 1. - wz0
              pnum(iw)    = pnumz(icell)*wz1    + pnumz(icell+1)*wz0
              xbar(iw)    = xbarz(icell)*wz1    + xbarz(icell+1)*wz0
              ybar(iw)    = ybarz(icell)*wz1    + ybarz(icell+1)*wz0
              zbar(iw)    = zbarz(icell)*wz1    + zbarz(icell+1)*wz0
              xpbar(iw)   = xpbarz(icell)*wz1   + xpbarz(icell+1)*wz0
              ypbar(iw)   = ypbarz(icell)*wz1   + ypbarz(icell+1)*wz0
              vxbar(iw)   = vxbarz(icell)*wz1   + vxbarz(icell+1)*wz0
              vybar(iw)   = vybarz(icell)*wz1   + vybarz(icell+1)*wz0
              vzbar(iw)   = vzbarz(icell)*wz1   + vzbarz(icell+1)*wz0
              xybar(iw)   = xybarz(icell)*wz1   + xybarz(icell+1)*wz0
              xypbar(iw)  = xypbarz(icell)*wz1  + xypbarz(icell+1)*wz0
              yxpbar(iw)  = yxpbarz(icell)*wz1  + yxpbarz(icell+1)*wz0
              xpypbar(iw) = xpypbarz(icell)*wz1 + xpypbarz(icell+1)*wz0
              xsqbar(iw)  = xsqbarz(icell)*wz1  + xsqbarz(icell+1)*wz0
              ysqbar(iw)  = ysqbarz(icell)*wz1  + ysqbarz(icell+1)*wz0
              zsqbar(iw)  = zsqbarz(icell)*wz1  + zsqbarz(icell+1)*wz0
              xpsqbar(iw) = xpsqbarz(icell)*wz1 + xpsqbarz(icell+1)*wz0
              ypsqbar(iw) = ypsqbarz(icell)*wz1 + ypsqbarz(icell+1)*wz0
              vxsqbar(iw) = vxsqbarz(icell)*wz1 + vxsqbarz(icell+1)*wz0
              vysqbar(iw) = vysqbarz(icell)*wz1 + vysqbarz(icell+1)*wz0
              vzsqbar(iw) = vzsqbarz(icell)*wz1 + vzsqbarz(icell+1)*wz0
              xxpbar(iw)  = xxpbarz(icell)*wz1  + xxpbarz(icell+1)*wz0
              yypbar(iw)  = yypbarz(icell)*wz1  + yypbarz(icell+1)*wz0
              xvzbar(iw)  = xvzbarz(icell)*wz1  + xvzbarz(icell+1)*wz0
              yvzbar(iw)  = yvzbarz(icell)*wz1  + yvzbarz(icell+1)*wz0
              zvzbar(iw)  = zvzbarz(icell)*wz1  + zvzbarz(icell+1)*wz0
              vxvzbar(iw) = vxvzbarz(icell)*wz1 + vxvzbarz(icell+1)*wz0
              vyvzbar(iw) = vyvzbarz(icell)*wz1 + vyvzbarz(icell+1)*wz0
              xrms(iw)    = xrmsz(icell)*wz1    + xrmsz(icell+1)*wz0
              yrms(iw)    = yrmsz(icell)*wz1    + yrmsz(icell+1)*wz0
              zrms(iw)    = zrmsz(icell)*wz1    + zrmsz(icell+1)*wz0
              xprms(iw)   = xprmsz(icell)*wz1   + xprmsz(icell+1)*wz0
              yprms(iw)   = yprmsz(icell)*wz1   + yprmsz(icell+1)*wz0
              epsx(iw)    = epsxz(icell)*wz1    + epsxz(icell+1)*wz0
              epsy(iw)    = epsyz(icell)*wz1    + epsyz(icell+1)*wz0
              epsz(iw)    = epszz(icell)*wz1    + epszz(icell+1)*wz0
              epsnx(iw)   = epsnxz(icell)*wz1   + epsnxz(icell+1)*wz0
              epsny(iw)   = epsnyz(icell)*wz1   + epsnyz(icell+1)*wz0
              epsnz(iw)   = epsnzz(icell)*wz1   + epsnzz(icell+1)*wz0
              epsg(iw)    = epsgz(icell)*wz1    + epsgz(icell+1)*wz0
              epsh(iw)    = epshz(icell)*wz1    + epshz(icell+1)*wz0
              epsng(iw)   = epsngz(icell)*wz1   + epsngz(icell+1)*wz0
              epsnh(iw)   = epsnhz(icell)*wz1   + epsnhz(icell+1)*wz0
              vxrms(iw)   = vxrmsz(icell)*wz1   + vxrmsz(icell+1)*wz0
              vyrms(iw)   = vyrmsz(icell)*wz1   + vyrmsz(icell+1)*wz0
              vzrms(iw)   = vzrmsz(icell)*wz1   + vzrmsz(icell+1)*wz0
            endif
          enddo

        endif

      endif

      momentstime = momentstime + (wtime() - timetemp)

      return
      end
c=============================================================================
      subroutine getlabwn()
      use Beam_acc
      use InDiag
      use Picglb
      use Z_Moments
      use Lab_Moments
 
c  Sets moments for species 1 in all axial windows in the lab frame
c  Also saves appropriate moments in time history arrays.
c  Get info from z moments
c  The data is only saved every itlabwn steps, even if the moments
c  are actually calculated more or less often.

      real(kind=8):: w0,w1
      integer     :: i,ilw,iz

      if (iflabwn .eq. 0) return
      if (mod(it,itlabwn) .gt. 0) return
 
      do i = 1, nlabwn
        if (zlw(i) .gt. (zbeam+zmmntmin) .and.
     &      zlw(i) .lt. (zbeam+zmmntmax) .and.
     &      ilabwn(i) .lt. ntlabwn) then
          ilabwn(i) = ilabwn(i) + 1
          ilw = ilabwn(i)
          iz = (zlw(i) - zbeam - zmmntmin)*dzmi
          w0 = (zlw(i) - zbeam - zmmntmin)*dzmi - iz
          w1 = 1. - w0
          timelw(ilw,i) = time
          pnumlw(ilw,i)  = w1*pnumz(iz)  + w0*pnumz(iz+1)
          xbarlw(ilw,i)  = w1*xbarz(iz)  + w0*xbarz(iz+1)
          ybarlw(ilw,i)  = w1*ybarz(iz)  + w0*ybarz(iz+1)
          vzbarlw(ilw,i) = w1*vzbarz(iz) + w0*vzbarz(iz+1)
          epsxlw(ilw,i)  = w1*epsxz(iz)  + w0*epsxz(iz+1)
          epsylw(ilw,i)  = w1*epsyz(iz)  + w0*epsyz(iz+1)
          epszlw(ilw,i)  = w1*epszz(iz)  + w0*epszz(iz+1)
          vxrmslw(ilw,i) = w1*vxrmsz(iz) + w0*vxrmsz(iz+1)
          vyrmslw(ilw,i) = w1*vyrmsz(iz) + w0*vyrmsz(iz+1)
          vzrmslw(ilw,i) = w1*vzrmsz(iz) + w0*vzrmsz(iz+1)
          xrmslw(ilw,i)  = w1*xrmsz(iz)  + w0*xrmsz(iz+1)
          yrmslw(ilw,i)  = w1*yrmsz(iz)  + w0*yrmsz(iz+1)
        endif
      enddo

      return
      end
c=============================================================================
      subroutine gammaadv(np,gaminv,uxp,uyp,uzp,gamadv,lrelativ)
      use Constant
      integer     :: np
      real(kind=8):: gaminv(np),uxp(np),uyp(np),uzp(np)
      character*8 :: gamadv
      logical     :: lrelativ

c  Advance Gamma Inverse from it-3/2 to it-1/2 using method given by gamadv


      real(kind=8):: clghtisq,usq
      integer     :: ip

c  If not using relativity then don't do anything
      if (lrelativ) then

      clghtisq = 1./clight**2

c  Standard gamma advance
      if (gamadv .eq. "stndrd") then
        do ip=1,np
          usq = (uxp(ip)**2 + uyp(ip)**2 + uzp(ip)**2)*clghtisq
          gaminv(ip) = 1./sqrt(1. + usq)
        enddo
      endif

c  Fast gamma advance, version 1
      if (gamadv .eq. "fast 1") then
        do ip=1,np
          usq = (uxp(ip)**2 + uyp(ip)**2 + uzp(ip)**2)*clghtisq
          gaminv(ip) = gaminv(ip)*1.5 - (1.+usq)/gaminv(ip)*0.5
        enddo
      endif

c  Fast gamma advance, version 2, fastest, but least accurate
      if (gamadv .eq. "fast 2") then
        do ip=1,np
          gaminv(ip)  = 1. - 0.5*(uxp(ip)**2+uyp(ip)**2+uzp(ip)**2)*clghtisq
        enddo
      endif

      endif

      return
      end
c=============================================================================
      subroutine setcurr (curr,np,zp,uzp,gaminv,q,wght,zbeam,dzz,zzmin)
      integer     :: np
      real(kind=8):: q,wght,zbeam,dzz,zzmin
      real(kind=8):: curr(0:*), zp(np), uzp(np), gaminv(np)

c  Sets 1d beam current as a function of z, directly from particle data

      real(kind=8):: s0(NPARPGRP), s1(NPARPGRP)
      integer     :: kk(NPARPGRP)
      real(kind=8):: g,w1,w0,dzzi
      integer     :: ip

      g = wght * q / dzz
      dzzi = 1. / dzz
c     --- vector loop
      do ip = 1,np
         kk(ip) = (zp(ip) - zbeam - zzmin) * dzzi
         w1 =  (zp(ip) - zbeam - zzmin) * dzzi - kk(ip)
         w0 = 1. - w1
         s0(ip) = w0 * g * gaminv(ip) * uzp(ip)
         s1(ip) = w1 * g * gaminv(ip) * uzp(ip)
      enddo
c     --- scalar loop
      do ip = 1,np
         if (uzp(ip) .gt. 0) then
            curr(kk(ip)) = curr(kk(ip)) + s0(ip)
            curr(kk(ip)+1) = curr(kk(ip)+1) + s1(ip)
         endif
      enddo

      return
      end
c=============================================================================
      subroutine percurr(curr,nzzarr,periinz)
      integer     :: nzzarr
      real(kind=8):: curr(0:nzzarr)
      logical     :: periinz
c Make current periodic.

c Parallel version makes call to special routine percurr_slave
#ifdef PARALLEL
      call percurr_slave(curr,nzzarr,periinz)
#else
      if (periinz) then
        curr(0) = curr(0) + curr(nzzarr)
        curr(nzzarr) = curr(0)
      endif
#endif

      return
      end
c=============================================================================
      subroutine getvzofz
      use Z_arrays

c  Computes the mean beam z velocity from the current and line charge density

      integer     :: i

      do i = 0, nzzarr
         vzofz(i) = curr(i) / dvnz(linechg(i))
      enddo

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

      integer     :: i
      real(kind=8):: d

      if (ifgap) then

        do i = nzzarr-1, 1, -1
           d = dvnz(2.*rgap*dzz)
           d = dvnz(cgap*vbeam**2/dzz**2 + vbeam/d)
           egap(i-1) = (- cgap * vbeam**2 * (-2.*egap(i)+egap(i+1)) / dzz**2
     &                  + vbeam * egap(i+1) / dvnz(2.*rgap*dzz)
     &                  - egap(i) / dvnz(lgap)
     &                  + vbeam * (curr(i+1) - curr(i-1)) / (2.*dzz)
     &                 ) / d
        enddo

      endif

      return
      end
c=============================================================================
      subroutine gapfield (np,zp,ez,zbeam,zzmin,egap,dzz)
      use InGaps
      integer     :: np
      real(kind=8):: zbeam,zzmin,dzz
      real(kind=8):: zp(np),ez(np),egap(0:*)

c  Interpolates the spatially averaged gap electric field onto the particle
c  positions, then uses it to augment the individual particle Ez's.

      real(kind=8):: dzzi,w1,w0
      integer     :: ip,kk

      if (ifgap) then

        dzzi = 1./dzz
        do ip=1,np
           kk = (zp(ip) - zbeam - zzmin) * dzzi
           w1 = (zp(ip) - zbeam - zzmin) * dzzi - kk
           w0 = 1. - w1
           ez(ip) = ez(ip) + w0*egap(kk) + w1*egap(kk+1)
        enddo

      endif

      return
      end
c=============================================================================
      subroutine setgamma(lrelativ)
      use Constant
      use InPart
      use Particles
      logical     :: lrelativ
c  Converts v to u, sets inverse gamma factor for all particles


      real(kind=8):: clghtisq,gamma
      integer     :: is,ip

      if (lrelativ) then
        clghtisq = 1./clight**2
        do is=1,ns
          do ip = ins(is),ins(is)+nps(is)-1
             gamma = 1. / sqrt ( 1. -
     &        (uxp(ip)**2 + uyp(ip)**2 + uzp(ip)**2) * clghtisq )
             uxp(ip) = gamma * uxp(ip)
             uyp(ip) = gamma * uyp(ip)
             uzp(ip) = gamma * uzp(ip)
             gaminv(ip) = 1. / gamma
          enddo
        enddo
      else
        do is=1,ns
          do ip=ins(is),ins(is)+nps(is)-1
             gaminv(ip) = 1.0
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine zpartbnd(zmmax,zmmin,dz,zgrid)
      use InGen
      use InPart
      use Particles
      real(kind=8):: zmmax,zmmin,dz,zgrid

c  Impose boundary conditions on zp

      real(kind=8):: syslen
      integer     :: is,ip

c Parallel version makes call to special routine zpartbnd_slave
#ifdef PARALLEL
      call zpartbnd_slave(zmmax,zmmin,dz,zgrid)
#else

      if (periinz) then
c     --- periodic boundary condition
        syslen = zmmax - zmmin
        do is=1,ns
          do ip=ins(is),ins(is)+nps(is)-1
            if (zp(ip) .lt. zgrid + zmmin)
     &        zp(ip) = zp(ip) + syslen*int((zgrid + zmmax - zp(ip))/syslen)
            if (zp(ip) .gt. zgrid + zmmax)
     &        zp(ip) = zp(ip) - syslen*int((zp(ip) - zgrid - zmmin)/syslen)
          enddo
        enddo

      elseif (stickyz) then
c     --- sticky boundary condition
c     --- Particles stick at the extreme edges of the grid.
        do is=1,ns
          do ip=ins(is),ins(is)+nps(is)-1
            if ((zp(ip)-zgrid).ge.zmmax) then
              zp(ip) = zmmin - .9999*dz + zgrid
              uxp(ip) = 0.0
              uyp(ip) = 0.0
              uzp(ip) = 0.0
            elseif ((zp(ip)-zgrid).le.zmmin) then
              zp(ip) = zmmin - .9999*dz + zgrid
              uxp(ip) = 0.0
              uyp(ip) = 0.0
              uzp(ip) = 0.0
            elseif (uzp(ip).eq.0.) then
c             --- Make sure dead particles stay dead and are not ressurrected
c                 (only need if moving backward in time).
              zp(ip) = zmmin - .9999*dz + zgrid
            endif
          enddo
        enddo
      endif

#endif
c --- end of ifelse above for paralle version.

      return
      end
c=============================================================================
c The following two subroutines are made obsolete by zpartbnd above, but
c are maintained for use by the wrz and wxy packages.
c=============================================================================
      subroutine periz(np,zp,zgrid,zmmax,zmmin)
      integer     :: np
      real(kind=8):: zgrid,zmmax,zmmin
      real(kind=8):: zp(np)

c  Imposes periodic boundary conditions on zp

      real(kind=8):: syslen,sysleni
      integer     :: ip

      syslen = zmmax - zmmin
      sysleni = 1./syslen
      do ip=1,np
        if (zp(ip) .lt. zgrid + zmmin)
     &    zp(ip) = zp(ip) + syslen*int((zgrid + zmmax - zp(ip))*sysleni)
        if (zp(ip) .gt. zgrid + zmmax)
     &    zp(ip) = zp(ip) - syslen*int((zp(ip) - zgrid - zmmin)*sysleni)
      enddo

      return
      end
c=============================================================================
      subroutine stckyz(np,zp,zmmax,zmmin,dz,uxp,uyp,uzp,zgrid)
      integer     :: np
      real(kind=8):: zmmax,zmmin,dz,zgrid
      real(kind=8):: zp(np)
      real(kind=8):: uxp(np),uyp(np),uzp(np)

c  Enforces sticky b.c.'s on the z walls.
c  Particles stick at the extreme edges of the grid.

      integer     :: ip

      do ip=1,np
        if ((zp(ip)-zgrid).ge.zmmax) then
          zp(ip) = zmmin - .9999*dz + zgrid
          uxp(ip) = 0.0
          uyp(ip) = 0.0
          uzp(ip) = 0.0
        elseif ((zp(ip)-zgrid).le.zmmin) then
          zp(ip) = zmmin - .9999*dz + zgrid
          uxp(ip) = 0.0
          uyp(ip) = 0.0
          uzp(ip) = 0.0
        elseif (uzp(ip).eq.0.) then
c         --- Make sure dead particles stay dead and are not ressurrected
c             (only need if moving backward in time).
          zp(ip) = zmmin - .9999*dz + zgrid
        endif
      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine zbeamcom(zbeam)
      use InPart
      use Particles
      real(kind=8):: zbeam

c Calculate the center of mass of the beam along the z axis from the particles


      integer     :: is,ip
      real(kind=8):: totz,totmass

      totz = 0.
      totmass = 0.
      do is=1,ns
        do ip=ins(is),ins(is)+nps(is)-1
          if (uzp(ip) .ne. 0.) then
            totz = totz + zp(ip)*sm(is)
            totmass = totmass + sm(is)
          endif
        enddo
      enddo

      zbeam = totz/dvnz(totmass)

      return
      end
c=============================================================================
      subroutine alotpart
      use InPart
      use Particles

c Allocate space for particles and set up associated arrays.


      integer     :: is

c     --- Set one of npmax and np_s from the other. If they are both set,
c     --- then do nothing. That is not an error since the two variables
c     --- do not need to be consistent. For example, the parallel 3-D code
c     --- sets npmax to the local size of the arrays and np_s to the total
c     --- number of particles across all of the processors, leaving
c     --- np_s > npmax.
c     --- Note that here npmax is not set in the parallel code since it is
c     --- assumed that the parallel routines will set npmax properly to
c     --- distribute the particles among the processors.
      do is=1,ns
        if (np_s(is) .eq. 0) np_s(is) = npmax*sp_fract(is)
      enddo
#ifndef PARALLEL
      npmax = 0
      do is=1,ns
        npmax = npmax + np_s(is)
      enddo
      np = npmax
#endif
      npmaxb = npmax

c     --- Allocate space. A gchange is done primarily so that any data
c     --- in npmax_s is preserved.
      call gchange("Particles",0)

c     --- setup npmax_s array
      npmax_s(0) = 0
      do is=1,ns
        if (npmax_s(is) .eq. 0) npmax_s(is) = npmax_s(is-1) + np_s(is)
        ins(is) = npmax_s(is-1) + 1
        nps(is) = 0
      enddo

      return
      end
c=============================================================================
      subroutine chckpart(is,nlower,nhigher,lfullshft)
      use InPart
      use Particles
      integer     :: is,nlower,nhigher
      logical     :: lfullshft

c Make sure that there is enough space in the particle arrays for nlower
c new particles below and nhigher above the live particles.  Returns if
c there is already enough space above and below.  If there is enough total
c space but not enough room above or below, the particles are shifted
c appropriately.  If lfullshft is true, the particles are shifted all the
c way to the end of the species block minus nhigher.  If there is not
c enough space, add more to the arrays. Particle data is shifted
c appropriately.


      integer     :: nadd,i,ishft

c     --- If there is already enough space, then return.
      if (ins(is)-npmax_s(is-1)-1 .ge. nlower .and.
     &    npmax_s(is)-ins(is)-nps(is)+1 .ge. nhigher) return

c     --- If there is not enough space, then allocate new space.
      if (ins(is)-npmax_s(is-1)-1 + npmax_s(is)-ins(is)-nps(is)+1 .lt.
     &    nlower + nhigher) then

c       --- Amount of new space to add: this is calculated to be a sizable
c       --- chunk, but not too big.  By default, add 10 times the space needed.
c       --- Limit that to 10000, except, if the space needed is greater than
c       --- 10000, add the amount of space needed.
        nadd = nlower + nhigher -
     &           (ins(is)-npmax_s(is-1)-1 + npmax_s(is)-ins(is)-nps(is)+1)
        nadd = min(10*nadd, max(10000, nadd))

c       --- change the particle array dimension variables and allot space
        npmax = npmax + nadd
        npmaxb = npmax
        call gchange("Particles",0)

c       --- Loop over species species above 'is', shifting
c       --- them up to make space below for species is.
        do i=ns,is+1,-1
          npmax_s(i) = npmax_s(i) + nadd
          call shftpart(i,nadd)
        enddo

c       --- Increase size of the block for species 'is'.
        npmax_s(is) = npmax_s(is) + nadd

      endif

c     --- Now, there is enough room to shift the particles appropriately.
c     --- Even though there may be enough space needed for the new particles,
c     --- the shift should still be done for the case when lfullshft is true.
c     --- Note that is some cases no shift will be needed.  In those cases,
c     --- ishft comes out to be zero and shftpart returns immediately when
c     --- ishft is zero.
      if (lfullshft) then
        ishft = npmax_s(is) - ins(is) - nps(is) + 1 - nhigher
      else
        ishft = max(0,nlower - (ins(is)-npmax_s(is-1)-1)) +
     &          min(0,(npmax_s(is)-ins(is)-nps(is)+1) - nhigher)
      endif
      call shftpart(is,ishft)

      return
      end
c=============================================================================
      subroutine shftpart(is,ishft)
      use InPart
      use Particles
      integer     :: is,ishft

c Shift particles by 'ishft'.  All of the data from the low end of the
c species block to the top of the live particles are moved.  The area
c below the live particles may include newly injected particles which
c have not been included into the simulation yet.
c Assumes that the is enought space in the particle arrays to make the
c shift without clobbering other particles or going past the ends of
c the arrays.


      integer     :: ip

c     --- return if shift is zero
      if (ishft .eq. 0) return

c     --- Positive and negative shifts are treated differently since
c     --- a negative shift may result in lost data.
      if (ishft .gt. 0) then

c       --- Loop from top of live particles to bottom of species block.
        do ip=ins(is)+nps(is)-1, npmax_s(is-1)+1, -1
          xp(ip+ishft) = xp(ip)
          yp(ip+ishft) = yp(ip)
          zp(ip+ishft) = zp(ip)
          uxp(ip+ishft) = uxp(ip)
          uyp(ip+ishft) = uyp(ip)
          uzp(ip+ishft) = uzp(ip)
          gaminv(ip+ishft) = gaminv(ip)
        enddo

c       --- change lower bound of live particles
        ins(is) = ins(is) + ishft

      else
c       --- Loop from -ishft above bottom of species block to top of live
c       --- particles.  Note that particles will be clobbered if
c       --- (ins(is) - npmax_s(is-1)) < -ishft
        if (ins(is) - npmax_s(is-1) .lt. -ishft) then
          print*,"WARNING: shftpart: particle data will be lost since there is"
          print*,"                   not enough space in the arrays."
        endif

        do ip=npmax_s(is-1)+1-ishft, ins(is)+nps(is)-1
          xp(ip+ishft) = xp(ip)
          yp(ip+ishft) = yp(ip)
          zp(ip+ishft) = zp(ip)
          uxp(ip+ishft) = uxp(ip)
          uyp(ip+ishft) = uyp(ip)
          uzp(ip+ishft) = uzp(ip)
          gaminv(ip+ishft) = gaminv(ip)
        enddo

c       --- change number of particles if data is lost
        if (ins(is) - npmax_s(is-1) .lt. -ishft) then
          nps(is) = nps(is) - (-ishft - (ins(is) - npmax_s(is-1)))
          ins(is) = npmax_s(is-1) + 1
        else
c         --- change lower bound only
          ins(is) = ins(is) + ishft
        endif

      endif

      return
      end
c=============================================================================
      subroutine copypart(it,nn,ii,istart)
      use InPart
      use Particles
      integer     :: it,nn,istart
      integer     :: ii(0:nn-1)
c Copies particle data from locations given by ii to the locations starting
c at 'it', or if istart > 0, copy sequential particles from locations
c starting at 'it' to locations starting at 'istart'.
      integer     :: i,in,ij
      if (istart > 0) then
        do i=0,nn-1
          ij = istart + i
          in = it + i
          xp(in) = xp(ij)
          yp(in) = yp(ij)
          zp(in) = zp(ij)
          uxp(in) = uxp(ij)
          uyp(in) = uyp(ij)
          uzp(in) = uzp(ij)
          gaminv(in) = gaminv(ij)
        enddo
      else
        do i=0,nn-1
          ij = ii(i)
          in = it + i
          xp(in) = xp(ij)
          yp(in) = yp(ij)
          zp(in) = zp(ij)
          uxp(in) = uxp(ij)
          uyp(in) = uyp(ij)
          uzp(in) = uzp(ij)
          gaminv(in) = gaminv(ij)
        enddo
      endif
      return
      end
c=============================================================================
c=============================================================================
c  VARIOUS UTILITIES
c=============================================================================
      subroutine copyarry (source, target, nwords)
      integer     :: nwords
      real(kind=8):: source(nwords), target(nwords)
      integer     :: i

      do 100 i = 1, nwords
         target(i) = source(i)
  100 continue

      return
      end
c=============================================================================
      subroutine sumarry (source, target, nwords)
      integer     :: nwords
      real(kind=8):: source(nwords), target(nwords)
      integer     :: i

      do 100 i = 1, nwords
         target(i) = target(i) + source(i)
  100 continue

      return
      end
c=============================================================================
      subroutine zeroarry (arry,n)
      integer     :: n
      real(kind=8):: arry(*)

c  Sets ARRY to zero
      integer     :: i
      do i=1,n
        arry(i) = 0.0
      enddo

      return
      end
c=============================================================================
      subroutine writarry(nn,arry,filename)
      integer     :: nn
      real(kind=8):: arry(nn)
      character*(*):: filename
      integer     :: ii

      open(unit=27,file=filename,status='new')

      do ii=1,nn
        write(27,20) arry(ii)
      enddo
 20   format(x,e10.3)

      close(27)

      return
      end
c=============================================================================
      real(kind=8) function fnice(i,e10)
      integer     :: i
      real(kind=8):: e10
c this is used to pick nice numbers; also gives next larger power of 10.
c for: i =  -4 -3 -2 -1 0  1  2  3   4   5
c  fnice = .05 .1 .2 .5 1  2  5 10  20  50
c    e10 =  .1 .1  1  1 1 10 10 10 100 100
c     i3 =   2  3  1  2 3  1  2  3   1   2
c    i10 =  -1 -1  0  0 0  1  1  1   2   2

      integer     :: i10,i3
      real(kind=8):: a(3) = (/.2,.5,1./)
      i10=i/3
      if(i.gt.0) i10=(i+2)/3
      i3=i+3-3*i10
      e10=10.**i10
      fnice=a(i3)*e10
      return
      end
c=============================================================================
      subroutine psumx(a,b,n)
      integer     :: n
      real(kind=8):: a(n), b(n)
c     b := partial sum of a
      integer     :: i
      b(1) = a(1)
      do i = 2, n
        b(i) = a(i) + b(i-1)
      enddo
      return
      end
c=============================================================================
      logical function inlist (list, listlen, text)
      integer     :: listlen
      character*8 :: text, list(listlen)
      integer     :: i

      inlist = .false.
      do i = 1,listlen
         if (text .eq. list(i)) inlist = .true.
      enddo

      return
      end
c=======================================================================
      real(kind=8) function wtime()
      real(kind=4):: cpu,io,sys,mem

c     --- Get current time
      call ostime(cpu, io, sys, mem)
      wtime = cpu

      return
      end
c=======================================================================
      subroutine wtimeon()
      real(kind=4):: cpu,io,sys,mem
      common /wtimer/cpu,io,sys,mem

c     --- turns timer on
      call ostime(cpu, io, sys, mem)

      return
      end
c=======================================================================
      real(kind=8) function wtimeoff()
      real(kind=4):: cpu,io,sys,mem
      real(kind=4):: fcpu,fio,fsys,fmem
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
      if (xs .ge. 0.) go to 42
      xs = xs + 2. * xsi
      return
      end
c=============================================================================
      real(kind=8) function rnrev(i,nbase)
      integer     :: i,nbase
c  Returns a fraction (0 <= rnrev <= 1) representing the
c  nbase-reversed number corresponding to i; i.e.
c  for i written in base nbase as n4n3n2n1.0, this returns
c  rnrev = 0.n1n2n3n4, where n1, n2.. are the digits in base n.
c  Disregards roundoff errors for nbase not a power of 2.
c  From N. Otani, 1985
      integer     :: j,n,jnext
      real(kind=8):: powern
      rnrev = 0.
      j = i
      n = nbase
      powern = 1.
   10 jnext = j / n
      powern = powern / n
      rnrev = rnrev + (j-jnext*n) * powern
      if (jnext .eq. 0) return
      j = jnext
      go to 10
      end
c=============================================================================
      subroutine rnrevarray(n,x,i,nbase)
      integer     :: n,i,nbase
      real(kind=8):: x(n)
c Fills an array with uniformly distributed digit reversed numbers.
c The numbers range from 0 to 1.
      real(kind=8):: rnrev
      integer     :: j
      do j=1,n
        x(j) = rnrev(i+j-1,nbase)
      enddo
      return
      end
c=============================================================================
      subroutine sphere4(a,b,c,d,n)
      use Constant
      integer     :: n
      real(kind=8):: a(n), b(n), c(n), d(n)
c   Distribute points on surface of a 4-d unit sphere.
c   (a*a + b*b) increases, and (c*c + d*d) decreases, monotonically.
      real(kind=8):: xsrev,rr,phi,r2rev,rnrev
      integer     :: i

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
      integer     :: g1,g2,g3
      real(kind=8):: a(g1), b(g1), c(g1), d(g1)
c   Distribute points on surface of a 4-d unit sphere using fibonacci numbers.
c   (a*a + b*b) increases, and (c*c + d*d) decreases, monotonically.
      integer     :: i
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
      integer     :: iflag = -1
      save y, iflag

      iflag = -iflag
      if (iflag .lt. 0) then
         rnorm = y
         return
      endif

c 100 v1 = 2. * wranf() - 1.
c     v2 = 2. * wranf() - 1.
c     s = v1**2 + v2**2 + SMALLPOS
c     if (s .gt. 1.) go to 100
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
      integer     :: i1,n,nbase1,nbase2
      real(kind=8):: dx
      real(kind=8):: x(n)

c  Gaussian random number generator: returns a value drawn from a normal
c    distribution with mean 0 and variance 1.
c    chosen using digit reversed type random numbers
c  Uses Box-Muller-Marsaglia method:
c  See Knuth, vol. 2, p.104 (first ed.), p.117 (second ed.)
c  Alex Friedman, July 1989 (after an earlier version).
c  David P. Grote, Febuary 1990 (added digit reversed stuff)
      integer     :: i
      real(kind=8):: s,phi,sq,rnrev

      do i=i1,i1+n-2,2
        s = rnrev(i,nbase1) + dx
        phi = 2.*pi*rnrev(i+1,nbase2)
        sq = sqrt(-2.*log(s))
        x(i-i1+1) = sq*cos(phi)
        x(i-i1+2) = sq*sin(phi)
      enddo

c     --- Get last number if n is odd.
      if (mod(n,2) .eq. 1) then
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
      integer     :: n
      real(kind=8):: a(n)
c   From interpreter, rma(&a,n) returns n Gaussian random numbers.
      integer     :: i
      real(kind=8):: rm

      do i = 1,n
        a(i) = rm()
      enddo
      return
      end
c=============================================================================
c=============================================================================
      subroutine load2d(np,x,y,nx,ny,n,dx,dy)
      integer     :: np,nx,ny
      real(kind=8):: x(np),y(np),n(0:nx,0:ny)
      real(kind=8):: dx,dy

c Load particles onto a two dimensional distribution.  This load is only
c approximate in that when the resulting particles are loaded back onto
c the same 2-D grid from which the distribution is taken, the distribution
c is not exactly reproduced. (The exact version is still under development.)

      real(kind=8):: dxi,dyi,ni,nip1,ninti,sumn,sumnx,wy
      real(kind=8):: rnrev
      integer     :: i,ip,ix,iy

c Note that this relies on automatic allocation of local arrays.
c ifelse(COMPILER,PGF77,[
#ifdef PGF77
      real(kind=8):: n1x(0:1024),n1y(0:1024)
      real(kind=8):: nintx(0:1024),ninty(0:1024)
#else
      real(kind=8):: n1x(0:nx),n1y(0:ny)
      real(kind=8):: nintx(0:nx),ninty(0:ny)
#endif

      dxi = 1./dx
      dyi = 1./dy

c     --- Initialize xp and yp
      do i=1,np
c       x(i) = .5/np + wranf()*(np-1.)/np
        x(i) = rnrev(i,2)
        y(i) = .5/np + (i-1.)/np
      enddo

c     --- n = n/sum(n)
      sumn = 0.
      do ix=0,nx
        do iy=0,ny
          sumn = sumn + n(ix,iy)
        enddo
      enddo
      do ix=0,nx
        do iy=0,ny
          n(ix,iy) = n(ix,iy)/sumn
        enddo
      enddo

c     --- n1y = sum(n,1)
      do iy=0,ny
        n1y(iy) = 0.
        do ix=0,nx
          n1y(iy) = n1y(iy) + n(ix,iy)
        enddo
      enddo

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
        sumnx = 0.
        do ix=0,nx
          sumnx = sumnx + n(ix,iy)*(1.-wy)+n(ix,iy+1)*wy
        enddo
        do ix=0,nx
          n1x(ix) = (n(ix,iy)*(1.-wy)+n(ix,iy+1)*wy)/sumnx
        enddo

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

c     --- Return n to it original value
      do ix=0,nx
        do iy=0,ny
          n(ix,iy) = n(ix,iy)*sumn
        enddo
      enddo

      return
      end
c=============================================================================
