#include "top.h"
c=============================================================================
c@(#) File TOP.F, version $Revision: 3.118 $, $Date: 2004/02/06 21:37:03 $
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
      use Particles

c Sets up species related quantities which were not set by the user, based
c on the rest of the user's input.    


      integer(ISZ):: is

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
c============================================================================ 
      subroutine resetlat
      use Lattice

c  Resizes the lattice arrays to the final lengths needed by the actual data.
c  Sets dipole field, lattice element 0, flags for bends, etc.

      integer(ISZ):: idrft,ib,id,iq,is,ih,ii,ia,im,ip
      integer(ISZ):: nndrft,nnbend,nndipo,nnquad,nnsext 
      integer(ISZ):: nnhele,nnhmlt,nnaccl,nnemlt,nnmmlt
      integer(ISZ):: nnbgrd,nnpgrd

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
      nnbgrd = 0
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
            if (nnhmlt < helene(ih)) nnhmlt = helene(ih)
            if (nnhmlt < helenm(ih)) nnhmlt = helenm(ih)           
         enddo
      enddo
c     --- count acceleration gaps 
      do ia = 0, naccl
         if (acclzs(ia) .ne. acclze(ia)) nnaccl = ia
      enddo
c     --- count multipole data sets 
      do im = 0, nemlt
         if (emltid(im) > 0) nnemlt = im
      enddo
      do im = 0, nmmlt
         if (mmltid(im) > 0) nnmmlt = im
      enddo
c     --- count bgridded field data sets 
      do ib = 0, nbgrd
         if (bgrdid(ib) > 0) nnbgrd = ib
      enddo
      do ip = 0, npgrd
         if (pgrdid(ip) > 0) nnpgrd = ip
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
      nbgrd = nnbgrd
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

c     --- Get the overlap level of the lattice elements
      call getoverlaps()

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

      integer(ISZ):: i,idrft,ib,id,iq,is,ih,ii,ia,ie,im,iz,l
      real(kind=8):: rnorm,emax,mmax,ephmax,mphmax,dzio2

c     --- if autoset dipoles are specified, then derive dipole locations 
c     --- consistently from bend locations and radii 
      if (diposet) then 
        do id = 0, ndipo
           if (id <= nbend) then
              if (bendrc(id) .ne. 0.) then
c                --- auto-set dipole field to match bend radius of curvature
                 if (dipoby(id) == 0. .and. dipoex(id) == 0.) then 
                    dipoby(id) = aion * amu / (zion * echarge)
     &                            * gammabar * vbeam / dvnz(bendrc(id))
                 endif
c                --- auto-set dipole starts and ends to match those of bends
                 if (dipozs(id) == dipoze(id)) then
                    dipozs(id) = bendzs(id)
                    dipoze(id) = bendze(id)
                 endif
c                --- auto-set dipole entrance & exit angles for "box dipoles"
                 if (dipotype == "box") then
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
        if (drftzs(0) == drftze(0)) then
           drftzs(0) = drftzs(ndrft) - zlatperi
           drftze(0) = drftze(ndrft) - zlatperi
           drftap(0) = drftap(ndrft)
           drftox(0) = drftox(ndrft)
           drftoy(0) = drftoy(ndrft)
        endif
        if (bendzs(0) == bendze(0)) then
           bendzs(0) = bendzs(nbend) - zlatperi
           bendze(0) = bendze(nbend) - zlatperi
           bendrc(0) = bendrc(nbend)
        endif
        if (dipozs(0) == dipoze(0)) then
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
        if (quadzs(0) == quadze(0)) then
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
        if (sextzs(0) == sextze(0)) then
           sextzs(0) = sextzs(nsext) - zlatperi
           sextze(0) = sextze(nsext) - zlatperi
           sextdb(0) = sextdb(nsext)
           sextde(0) = sextde(nsext)
        endif
        if (helezs(0) == heleze(0)) then
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
        if (acclzs(0) == acclze(0)) then
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
        if (emltid(0) == 0) then
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
        if (mmltid(0) == 0) then
           mmltzs(0) = mmltzs(nmmlt) - zlatperi
           mmltze(0) = mmltze(nmmlt) - zlatperi
           mmltph(0) = mmltph(nmmlt)
           mmltsf(0) = mmltsf(nmmlt)
           mmltsc(0) = mmltsc(nmmlt)
           mmltid(0) = mmltid(nmmlt)
           mmltap(0) = mmltap(nmmlt) 
        endif
        if (bgrdid(0) == 0) then
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
        if (pgrdid(0) == 0) then
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
      if (ntquad > 0) quads = .true.
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
      if (ntaccl > 0) accls = .true.
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
      bgrds = .false.
      do im = 0, nbgrd
         if (bgrdid(im).ne.0) bgrds = .true.
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
        if (bgrddxi(i) == 0.) bgrddxi(i) = 1./bgrddx(i)
        if (bgrddyi(i) == 0.) bgrddyi(i) = 1./bgrddy(i)
        if (bgrddzi(i) == 0.) bgrddzi(i) = 1./bgrddz(i)
      enddo

c     --- Calculate sines and cosines of bgrdph for efficiency
      do i=0,nbgrd
        bgrdsp(i) = sin(bgrdph(i)) 
        bgrdcp(i) = cos(bgrdph(i)) 
      enddo 

c     --- Calculate 1/pgrddx etc. for efficiency (minimizes divides)
c     --- Note that this coding is not executed if pgrdns=0.
      do i=1,pgrdns
        if (pgrddxi(i) == 0.) pgrddxi(i) = 1./pgrddx(i)
        if (pgrddyi(i) == 0.) pgrddyi(i) = 1./pgrddy(i)
        if (pgrddzi(i) == 0.) pgrddzi(i) = 1./pgrddz(i)
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
            if (abs(esemltp(iz,ie,i)) > emax) emax = abs(esemltp(iz,ie,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (emax == 0.) then
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
            if (abs(esemltphp(iz,ie,i)) > ephmax)
     &        ephmax = abs(esemltphp(iz,ie,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (ephmax == 0.) then
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
            if (abs(msmmltp(iz,im,i)) > mmax) mmax = abs(msmmltp(iz,im,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (mmax == 0.) then
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
            if (abs(msmmltphp(iz,im,i)) > mphmax)
     &        mphmax = abs(msmmltphp(iz,im,i))
          enddo
c         --- If the max is zero, i.e. the array is unset, then calculate
c         --- the derivative by finite differencing the data.
          if (mphmax == 0.) then
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
      subroutine getoverlaps()
      use Lattice
      use LatticeInternal

c For each of the elements types, get the level of overlaps
c zlatbuffer is the amount of space (on meters) around each element. It is
c used so that elements right next to each other but not strictly overlapping
c are considered to be overlapping.

      call getelemoverlaps(ndrft,drftzs,drftze,drftol,zlatbuffer,ndrftol)
      call getelemoverlaps(nbend,bendzs,bendze,bendol,zlatbuffer,nbendol)
      call getelemoverlaps(ndipo,dipozs,dipoze,dipool,zlatbuffer,ndipool)
      call getelemoverlaps(nquad,quadzs,quadze,quadol,zlatbuffer,nquadol)
      call getelemoverlaps(nsext,sextzs,sextze,sextol,zlatbuffer,nsextol)
      call getelemoverlaps(nhele,helezs,heleze,heleol,zlatbuffer,nheleol)
      call getelemoverlaps(naccl,acclzs,acclze,acclol,zlatbuffer,nacclol)
      call getelemoverlaps(nemlt,emltzs,emltze,emltol,zlatbuffer,nemltol)
      call getelemoverlaps(nmmlt,mmltzs,mmltze,mmltol,zlatbuffer,nmmltol)
      call getelemoverlaps(nbgrd,bgrdzs,bgrdze,bgrdol,zlatbuffer,nbgrdol)
      call getelemoverlaps(npgrd,pgrdzs,pgrdze,pgrdol,zlatbuffer,npgrdol)

c Now, setup indicing arrays in LatticeInternal
      call gchange("LatticeInternal",0)

      call getelemoverlapindices(ndrft,drftol,ndrftol,odrftoi,odrftio,odrftii,
     &                           odrftnn)
      call getelemoverlapindices(nbend,bendol,nbendol,obendoi,obendio,obendii,
     &                           obendnn)
      call getelemoverlapindices(ndipo,dipool,ndipool,odipooi,odipoio,odipoii,
     &                           odiponn)
      call getelemoverlapindices(nquad,quadol,nquadol,oquadoi,oquadio,oquadii,
     &                           oquadnn)
      call getelemoverlapindices(nsext,sextol,nsextol,osextoi,osextio,osextii,
     &                           osextnn)
      call getelemoverlapindices(nhele,heleol,nheleol,oheleoi,oheleio,oheleii,
     &                           ohelenn)
      call getelemoverlapindices(naccl,acclol,nacclol,oaccloi,oacclio,oacclii,
     &                           oacclnn)
      call getelemoverlapindices(nemlt,emltol,nemltol,oemltoi,oemltio,oemltii,
     &                           oemltnn)
      call getelemoverlapindices(nmmlt,mmltol,nmmltol,ommltoi,ommltio,ommltii,
     &                           ommltnn)
      call getelemoverlapindices(nbgrd,bgrdol,nbgrdol,obgrdoi,obgrdio,obgrdii,
     &                           obgrdnn)
      call getelemoverlapindices(npgrd,pgrdol,npgrdol,opgrdoi,opgrdio,opgrdii,
     &                           opgrdnn)

      return
      end
c============================================================================ 
      subroutine getelemoverlaps(nelem,elemzs,elemze,elemol,zlatbuffer,nelemol)
      integer(ISZ):: nelem,nelemol
      real(kind=8):: elemzs(0:nelem),elemze(0:nelem)
      integer(ISZ):: elemol(0:nelem)
      real(kind=8):: zlatbuffer

c Loops over lattice elements and for elements of the same type that overlap,
c keep track of the level of overlaps.

      integer(ISZ):: ie,ii
      logical(ISZ):: overlaps(1+nelem)

c     --- The first element is always at the lowest overlap level.
      elemol(0) = 1

c     --- Loop over the rest of the elements...
      do ie=1,nelem
      
c       --- First assume that this element doesn't overlap any others.
        overlaps = .false.

c       --- Loop over the previous elements and check which ones overlap the
c       --- current one. For each overlap, set the corresponding overlap level
c       --- to true.
        do ii=0,ie-1
          if (elemze(ie)+zlatbuffer > elemzs(ii) .and.
     &        elemze(ii)+zlatbuffer > elemzs(ie)) then
            overlaps(elemol(ii)) = .true.
          endif
        enddo

c       --- Loop over the list of overlap levels, looking for the lowest
c       --- level which has no overlap
        ii = 1
        do while (overlaps(ii) .and. ii <= 1+nelem)
          ii = ii + 1
        enddo
        elemol(ie) = ii

      enddo

      nelemol = maxval(elemol)

      return
      end
c============================================================================ 
      subroutine getelemoverlapindices(nelem,elemol,
     &                                 nelemol,oelemoi,oelemio,oelemii,oelemnn)
      use Lattice
      use LatticeInternal
      integer(ISZ),intent(in):: nelem,nelemol
      integer(ISZ),intent(in):: elemol(0:nelem)
      integer(ISZ),intent(out):: oelemoi(0:nelem),oelemio(0:nelem)
      integer(ISZ),intent(out):: oelemii(nelemol),oelemnn(nelemol)

c Generate indicing information which is used to make finding elements
c efficient.
c    nelem    the number of elements
c    elemol   the overlap level of each element
c    nelemol  the number of overlapping levels
c The output is
c    oelemoi   a list of the element indices group by overlap level
c    oelemio   revervse list of the element indices group by overlap level
c    oelemii   the starting index in that list for each overlap level
c    oelemnn   the number of elements in each overlap level

      integer(ISZ):: io,ie,ii

c     --- Loop over the number of levels of overlap. For each level, get
c     --- the element indices that are at that level and increment
c     --- the number of elements at that level.
      ii = 0
      do io=1,nelemol
        oelemnn(io) = 0
        oelemii(io) = ii
        do ie=0,nelem
          if (elemol(ie) == io) then
            oelemoi(ii) = ie
            oelemio(ie) = ii - oelemii(io) + 1
            oelemnn(io) = oelemnn(io) + 1
            ii = ii + 1
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
      integer(ISZ):: fstype

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

      integer(ISZ):: io,iz
      integer(ISZ):: iquadert,iquaderu,iquaderr,ii,imultert,imulteru,imulterr
      integer(ISZ):: id
      real(kind=8):: zz,offset,wi

c     --- Save the z location and time at which the internal lattice frame
c     --- is set.
      zlframe = zbeam
      zltime = time

c     --- Assume at the start that there are not any elements in the mesh.
      lindrft = .false.
      linbend = .false.
      lindipo = .false.
      linquad = .false.
      linsext = .false.
      linhele = .false.
      linemlt = .false.
      linmmlt = .false.
      linaccl = .false.
      linbgrd = .false.
      linpgrd = .false.

      do iz = 0, nzl
c        --- compute z of center of cell (may be off lattice)
c        --- If nzl is zero, then find the elements at zbeam. Otherwise,
c        --- find the elements at each of the grid points.
         if (nzl > 0) then
           zz = iz*dzl + zlmin + zbeam + 0.5*dzl
         else
           zz = zbeam
         endif

c        --- find index of nearest drft element
         if (drfts) then
           do io=1,ndrftol
             call getelemid(zz,offset,ndrft,drftzs(0),drftze(0),drftol(0),
     &                      odrftnn(io),odrftoi(odrftii(io)),odrftio,
     &                      cdrftid(iz,io),io)
c            --- load lattice info into comoving 1d array
             cdrftzs(iz,io) = drftzs(cdrftid(iz,io)) + offset
             cdrftze(iz,io) = drftze(cdrftid(iz,io)) + offset
             if (cdrftzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cdrftze(iz,io))
     &           lindrft(io) = .true.
           enddo
         endif

c        --- find index of nearest bend element
c        --- Note that while some setup is done for overlapping bends, they
c        --- are actually not supported since the idea of overlapping bends
c        --- doesn't make much sense. The extra indexing is done so that the
c        --- getelemid routine can still be used.
         if (bends) then
           call getelemid(zz,offset,nbend,bendzs(0),bendze(0),bendol(0),
     &                    obendnn(1),obendoi(obendii(1)),obendio,
     &                    cbendid(iz),1)
c          --- load lattice info into comoving 1d array
           cbendzs(iz) = bendzs(cbendid(iz)) + offset
           cbendze(iz) = bendze(cbendid(iz)) + offset
           cbendrc(iz) = bendrc(cbendid(iz))
c          --- Check whether this grid point is in the element.
c          --- The grid points intentionally overlap to avoid roundoff errors.
           if (cbendzs(iz) <= zlmax+zbeam+zlatbuffer .and.
     &                        zlmin+zbeam-zlatbuffer <= cbendze(iz))
     &         linbend = .true.
         endif

c        --- find index of nearest dipole element
c        --- The cdipo arrays calculated here are needed for residence
c        --- corrections (when dipos is true), and internal conductors (when
c        --- fstype is 3).
         if (dipos .or. fstype == 3) then
           do io=1,ndipool
             call getelemid(zz,offset,ndipo,dipozs(0),dipoze(0),dipool(0),
     &                      odiponn(io),odipooi(odipoii(io)),odipoio,
     &                      cdipoid(iz,io),io)
c            --- load lattice info into comoving 1d array
             cdipozs(iz,io) = dipozs(cdipoid(iz,io)) + offset
             cdipoze(iz,io) = dipoze(cdipoid(iz,io)) + offset
             cdipoby(iz,io) = dipoby(cdipoid(iz,io))
             cdipobx(iz,io) = dipobx(cdipoid(iz,io))
             cdipota(iz,io) = dipota(cdipoid(iz,io))
             cdipotb(iz,io) = dipotb(cdipoid(iz,io))
             cdipoex(iz,io) = dipoex(cdipoid(iz,io))
             cdipoey(iz,io) = dipoey(cdipoid(iz,io))
             if (cdipozs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cdipoze(iz,io))
     &           lindipo(io) = .true.
           enddo
         endif

c        --- find index of nearest quad element
c        --- The cquad arrays calculated here are needed for residence
c        --- corrections (when quads is true) and the internal conductors (when
c        --- fstype is either 2 or 3).
         if (quads .or. fstype == 2 .or. fstype == 3) then
           do io=1,nquadol
             call getelemid(zz,offset,nquad,quadzs(0),quadze(0),quadol(0),
     &                      oquadnn(io),oquadoi(oquadii(io)),oquadio(0),
     &                      cquadid(iz,io),io)
             id = cquadid(iz,io)
c            --- load lattice info into comoving 1d array
             cquadzs(iz,io) = quadzs(id) + offset
             cquadze(iz,io) = quadze(id) + offset
             cquaddb(iz,io) = quaddb(id)
             cquadde(iz,io) = quadde(id)
             cquadvx(iz,io) = quadvx(id)
             cquadvy(iz,io) = quadvy(id)
             if (ntquad > 0 .and. 
     &           quadts(id)<=time .and. time<quadts(id)+ntquad*quaddt(id)) then
               ii = int((time - quadts(id))/quaddt(id))
               wi = (time - quadts(id))/quaddt(id) - ii
               cquadde(iz,io) = cquadde(iz,io) + quadet(ii  ,id)*(1. - wi) +
     &                                           quadet(ii+1,id)*      wi 
               cquaddb(iz,io) = cquaddb(iz,io) + quadbt(ii  ,id)*(1. - wi) +
     &                                           quadbt(ii+1,id)*      wi 
             endif
c            --- set index for offset quads.  NOTE that this uses "nqerr"
c            --- and not "nquad".  This gives us a choice of either:
c            --- periodic errors (recirculator or ring) via nqerr=nquad, or
c            --- aperiodic errors (straight lattice w/ repeat).
c            --- THIS NEEDS CHECKING, especially when iqerr is nonzero.
             if (zlatperi > 0.) then
               iquadert = id + (int(10.+(zz-zlatstrt)/zlatperi)-10)*nquad
             else
               iquadert = id
             endif
             iquaderu = max(iquadert, 0)
             if (nqerr == 0) then
               iquaderr = 0
             else
               iquaderr = mod(iquaderu, nqerr) 
             endif
             cqoffx(iz,io) = qoffx(iquaderr)
             cqoffy(iz,io) = qoffy(iquaderr)
             if (cquadzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cquadze(iz,io))
     &           linquad(io) = .true.
           enddo
         endif

c        --- find index of nearest sext element
c        --- The csext arrays calculated here are needed for residence
c        --- corrections (when sexts is true).
         if (sexts) then
           do io=1,nsextol
             call getelemid(zz,offset,nsext,sextzs(0),sextze(0),sextol(0),
     &                      osextnn(io),osextoi(osextii(io)),osextio(0),
     &                      csextid(iz,io),io)
c            --- load lattice info into comoving 1d array
             csextzs(iz,io) = sextzs(csextid(iz,io)) + offset
             csextze(iz,io) = sextze(csextid(iz,io)) + offset
             csextdb(iz,io) = sextdb(csextid(iz,io))
             csextde(iz,io) = sextde(csextid(iz,io))
             if (csextzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= csextze(iz,io))
     &           linsext(io) = .true.
           enddo
         endif

c        --- find index of nearest hard-edge multipole element
         if (heles) then
           do io=1,nheleol
             call getelemid(zz,offset,nhele,helezs(0),heleze(0),heleol(0),
     &                      ohelenn(io),oheleoi(oheleii(io)),oheleio(0),
     &                      cheleid(iz,io),io)
c            --- load lattice info into comoving 1d array
             chelezs(iz,io) = helezs(cheleid(iz,io)) + offset
             cheleze(iz,io) = heleze(cheleid(iz,io)) + offset
             if (chelezs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cheleze(iz,io))
     &           linhele(io) = .true.
           enddo
         endif

c        --- find index of nearest acceleration element
         if (accls) then
           do io=1,nacclol
             call getelemid(zz,offset,naccl,acclzs(0),acclze(0),acclol(0),
     &                      oacclnn(io),oaccloi(oacclii(io)),oacclio(0),
     &                      cacclid(iz,io),io)
             id = cacclid(iz,io)
             if (acclze(id) + offset > acclzstt) then
               cacclzs(iz,io) = acclzs(id) + offset
               cacclze(iz,io) = acclze(id) + offset
               cacclez(iz,io) = acclez(id)
               cacclxw(iz,io) = acclxw(id) 
               cacclsw(iz,io) = acclsw(id) 
               if (acclts(id) <= time .and.
     &                           time < acclts(id)+ntaccl*accldt(id)) then
                 ii = int((time - acclts(id))/accldt(id))
                 wi = (time - acclts(id))/accldt(id) - ii
                 cacclez(iz,io) = acclez(id) + acclet(ii  ,id)*(1. - wi) +
     &                                         acclet(ii+1,id)*      wi 
               endif
               if (cacclzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                               zlmin+zbeam-zlatbuffer <= cacclze(iz,io))
     &             linaccl(io) = .true.
             else
               cacclzs(iz,io) = 0.
               cacclze(iz,io) = 0.
               cacclez(iz,io) = 0.
               cacclxw(iz,io) = 0.
               cacclsw(iz,io) = 0.  
             endif
           enddo
         endif

c        --- find index of nearest emlt element
c        --- The cemlt arrays calculated here are needed for the elements
c        --- described in terms of their multipole components.
         if (emlts) then
           do io=1,nemltol
             call getelemid(zz,offset,nemlt,emltzs(0),emltze(0),emltol(0),
     &                      oemltnn(io),oemltoi(oemltii(io)),oemltio(0),
     &                      cemltid(iz,io),io)
c            --- load lattice info into comoving 1d array
             id = cemltid(iz,io)
             cemltzs(iz,io) = emltzs(id) + offset
             cemltze(iz,io) = emltze(id) + offset
             cemltph(iz,io) = emltph(id)
             cemltsf(iz,io) = emltsf(id)
             cemltsc(iz,io) = emltsc(id)
             cemltim(iz,io) = emltid(id)
c            --- Check whether this grid point is in the element.
c            --- The grid points intentionally overlap to avoid roundoff errors.
             if (cemltzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cemltze(iz,io))
     &           linemlt(io) = .true.
c            --- set index for offset emlts.  NOTE that this uses "neerr"
c            --- and not "nemlt".  This gives us a choice of either:
c            --- periodic errors (recirculator or ring) via neerr=nemlt, or
c            --- aperiodic errors (straight lattice w/ repeat).
c            --- THIS NEEDS CHECKING
             if (zlatperi > 0.) then
               imultert = id+(int(10.+(zz-zlatstrt)/zlatperi)-10)*nemlt
             else
               imultert = id
             endif
             imulteru = max(imultert, 0)
             if (neerr == 0) then
               imulterr = 0
             else
               imulterr = mod(imulteru, neerr)
             endif
             cemltox(iz,io) = emltox(imulterr)
             cemltoy(iz,io) = emltoy(imulterr)
           enddo
         endif

c        --- find index of nearest mmlt element
c        --- The cmmlt arrays calculated here are needed for the elements
c        --- described in terms of their multipole components.
         if (mmlts) then
           do io=1,nmmltol
             call getelemid(zz,offset,nmmlt,mmltzs(0),mmltze(0),mmltol(0),
     &                      ommltnn(io),ommltoi(ommltii(io)),ommltio(0),
     &                      cmmltid(iz,io),io)
c            --- load lattice info into comoving 1d array
             id = cmmltid(iz,io)
             cmmltzs(iz,io) = mmltzs(id) + offset
             cmmltze(iz,io) = mmltze(id) + offset
             cmmltph(iz,io) = mmltph(id)
             cmmltsf(iz,io) = mmltsf(id)
             cmmltsc(iz,io) = mmltsc(id)
             cmmltim(iz,io) = mmltid(id)
c            --- Check whether this grid point is in the element.
c            --- The grid points intentionally overlap to avoid roundoff errors.
             if (cmmltzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cmmltze(iz,io))
     &           linmmlt(io) = .true.
c            --- set index for offset mmlts.  NOTE that this uses "nmerr"
c            --- and not "nmmlt".  This gives us a choice of either:
c            --- periodic errors (recirculator or ring) via nmerr=nmmlt, or
c            --- aperiodic errors (straight lattice w/ repeat).
c            --- THIS NEEDS CHECKING
             if (zlatperi > 0.) then
               imultert = id+(int(10.+(zz-zlatstrt)/zlatperi)-10)*nmmlt
             else
               imultert = id
             endif
             imulteru = max(imultert, 0)
             if (nmerr == 0) then
               imulterr = 0
             else
               imulterr = mod(imulteru, nmerr)
             endif
             cmmltox(iz,io) = mmltox(imulterr)
             cmmltoy(iz,io) = mmltoy(imulterr)
           enddo
         endif

c        --- find index of nearest bgrd element
c        --- The cbgrdid array calculated here is needed for the magnetic
c        --- elements described by B field data on a 3-D grid.
         if (bgrds) then
           do io=1,nbgrdol
             call getelemid(zz,offset,nbgrd,bgrdzs(0),bgrdze(0),bgrdol(0),
     &                      obgrdnn(io),obgrdoi(obgrdii(io)),obgrdio(0),
     &                      cbgrdid(iz,io),io)
c            --- load lattice info into comoving 1d array
             cbgrdzs(iz,io) = bgrdzs(cbgrdid(iz,io)) + offset
             cbgrdze(iz,io) = bgrdze(cbgrdid(iz,io)) + offset
c            --- Check whether this grid point is in the element.
c            --- The grid points intentionally overlap to avoid roundoff errors.
             if (cbgrdzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cbgrdze(iz,io))
     &           linbgrd(io) = .true.
           enddo
         endif

c        --- find index of nearest pgrd element
c        --- The cpgrdid array calculated here is needed for the electrostatic
c        --- elements described by potential data on a 3-D grid.
         if (pgrds) then
           do io=1,npgrdol
             call getelemid(zz,offset,npgrd,pgrdzs(0),pgrdze(0),pgrdol(0),
     &                      opgrdnn(io),opgrdoi(opgrdii(io)),opgrdio(0),
     &                      cpgrdid(iz,io),io)
c            --- load lattice info into comoving 1d array
             cpgrdzs(iz,io) = pgrdzs(cpgrdid(iz,io)) + offset
             cpgrdze(iz,io) = pgrdze(cpgrdid(iz,io)) + offset
c            --- Check whether this grid point is in the element.
c            --- The grid points intentionally overlap to avoid roundoff errors.
             if (cpgrdzs(iz,io) <= zlmax+zbeam+zlatbuffer .and.
     &                             zlmin+zbeam-zlatbuffer <= cpgrdze(iz,io))
     &           linpgrd(io) = .true.
           enddo
         endif

      enddo

c     --- For each element type, check to see if any elements are in the mesh.
      lindrft(0) = any(lindrft)
      lindipo(0) = any(lindipo)
      linquad(0) = any(linquad)
      linsext(0) = any(linsext)
      linhele(0) = any(linhele)
      linemlt(0) = any(linemlt)
      linmmlt(0) = any(linmmlt)
      linaccl(0) = any(linaccl)
      linbgrd(0) = any(linbgrd)
      linpgrd(0) = any(linpgrd)

c     --- All parallel slave processors must know if any processors
c     --- are in a bend.  There is probably a better way of doing this.
#ifdef MPIPARALLEL
      if (bends) then
        wi = 0.
        if (linbend) wi = 1.
        call parallelmaxrealarray(wi,1)
        if (wi > 0.5) linbend = .true.
      endif
#endif

      return
      end
c=============================================================================
      subroutine getelemid(z,offset,nelem,elemzs,elemze,elemol,
     &                     oelemnn,oelemoi,oelemio,id,io)
      use Lattice
      integer(ISZ):: nelem,oelemnn,id,io
      real(kind=8):: z,offset,elemzs(0:nelem),elemze(0:nelem)
      integer(ISZ):: elemol(0:nelem),oelemoi(oelemnn),oelemio(0:nelem)

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

      integer(ISZ):: ii
      real(kind=8):: zc

c The Lattice is only needed for zlatstrt and zlatperi.

c --- Offset of the lattice, for periodic repeat.
      if (zlatperi > 0.) then
        offset = (int(10.+(z-zlatstrt)/zlatperi)-10)*zlatperi+zlatstrt
      else
        offset = zlatstrt
      endif

c --- Check if there is only one element. If so, return id is that element.
      if (oelemnn == 1) then
        id = oelemoi(1)
        return
      endif

c --- Make sure id actually refers to an element of the overlap level.
      if (elemol(id) .ne. io) id = oelemoi(1)

c --- Check if z is somehow greater than the current element (id) location.
c --- If so, reset it to 0 and start from the beginning. This happens if
c --- during a periodic repeat of the lattice or if setlatt was called out
c --- of order. Also, make sure that id is within the proper bounds.
      if (id > oelemoi(oelemnn)) then
        id = oelemoi(1)
      else if (id > oelemoi(1)) then
        if (id > oelemoi(oelemnn)) then
          id = oelemoi(1)
        else if (z < .5*(elemze(oelemoi(oelemio(id)-1))+elemzs(id))+offset) then
          id = oelemoi(1)
        endif
      else if (id < oelemoi(1)) then
        id = oelemoi(1)
      endif

c --- Make ii consistent with id
      ii = oelemio(id)

c --- Scan sequentially through the list until the element is found.
c --- Note that loop is written in this way so that if ii reaches oelemnn
c --- there is not an out of bounds reference (oelemoi(ii+1)). The initial
c --- statement, zc = 0., is just to insure the zc has value.
      zc = 0.
      if (ii < oelemnn) zc = 0.5*(elemze(id) + elemzs(oelemoi(ii+1))) + offset
      do while (ii < oelemnn .and. z > zc)
        ii = ii + 1
        id = oelemoi(ii)
        if (ii < oelemnn) zc = 0.5*(elemze(id) + elemzs(oelemoi(ii+1))) + offset
      end do

      return
      end
c=============================================================================
      subroutine getdipoid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the dipole at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,ndipo,dipozs(0),dipoze(0),dipool(0),
     &                 odiponn(io),odipooi(odipoii(io)),odipoio(0),id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,ndipo,dipozs(0),dipoze(0),dipool(0),
     &                 odiponn(1),odipooi(odipoii(1)),odipoio(0),id,1)
c       --- Loop over rest of overlap levels looking for a dipo with lower
c       --- z start.
        id_temp = id
        do ii=2,ndipool
          call getelemid(zz,offset_temp,ndipo,dipozs(0),dipoze(0),dipool(0),
     &                   odiponn(ii),odipooi(odipoii(ii)),odipoio(0),id_temp,ii)
          if (dipozs(id_temp)+offset_temp < dipozs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getquadid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the quadrupole at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,nquad,quadzs(0),quadze(0),quadol(0),
     &                 oquadnn(io),oquadoi(oquadii(io)),oquadio(0),id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,nquad,quadzs(0),quadze(0),quadol(0),
     &                 oquadnn(1),oquadoi(oquadii(1)),oquadio(0),id,1)
c       --- Loop over rest of overlap levels looking for a quad with lower
c       --- z start.
        id_temp = id
        do ii=2,nquadol
          call getelemid(zz,offset_temp,nquad,quadzs(0),quadze(0),quadol(0),
     &                   oquadnn(ii),oquadoi(oquadii(ii)),oquadio(0),id_temp,ii)
          if (quadzs(id_temp)+offset_temp < quadzs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getheleid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the hele at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,nhele,helezs(0),heleze(0),heleol(0),
     &                 ohelenn(io),oheleoi(oheleii(io)),oheleio(0),id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,nhele,helezs(0),heleze(0),heleol(0),
     &                 ohelenn(1),oheleoi(oheleii(1)),oheleio(0),id,1)
c       --- Loop over rest of overlap levels looking for a hele with lower
c       --- z start.
        id_temp = id
        do ii=2,nheleol
          call getelemid(zz,offset_temp,nhele,helezs(0),heleze(0),heleol(0),
     &                   ohelenn(ii),oheleoi(oheleii(ii)),oheleio(0),id_temp,ii)
          if (helezs(id_temp)+offset_temp < helezs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getemltid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the emlt at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,nemlt,emltzs(0),emltze(0),emltol(0),
     &                 oemltnn(io),oemltoi(oemltii(io)),oemltio(0),id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,nemlt,emltzs(0),emltze(0),emltol(0),
     &                 oemltnn(1),oemltoi(oemltii(1)),oemltio(0),id,1)
c       --- Loop over rest of overlap levels looking for a emlt with lower
c       --- z start.
        id_temp = id
        do ii=2,nemltol
          call getelemid(zz,offset_temp,nemlt,emltzs(0),emltze(0),emltol(0),
     &                   oemltnn(ii),oemltoi(oemltii(ii)),oemltio(0),id_temp,ii)
          if (emltzs(id_temp)+offset_temp < emltzs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getmmltid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the mmlt at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,nmmlt,mmltzs(0),mmltze(0),mmltol(0),
     &                 ommltnn(io),ommltoi(ommltii(io)),ommltio(0),id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,nmmlt,mmltzs(0),mmltze(0),mmltol(0),
     &                 ommltnn(1),ommltoi(ommltii(1)),ommltio(0),id,1)
c       --- Loop over rest of overlap levels looking for a mmlt with lower
c       --- z start.
        id_temp = id
        do ii=2,nmmltol
          call getelemid(zz,offset_temp,nmmlt,mmltzs(0),mmltze(0),mmltol(0),
     &                   ommltnn(ii),ommltoi(ommltii(ii)),ommltio(0),id_temp,ii)
          if (mmltzs(id_temp)+offset_temp < mmltzs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine getacclid(zz,offset,id,io)
      use Lattice
      use LatticeInternal
      real(kind=8):: zz,offset
      integer(ISZ):: id,io

c Gets the id and offset for the accl at the z location zz and overlap
c level io. If the overlap level is 0, finds the element with the lowest
c z start.

      integer(ISZ):: ii,id_temp
      real(kind=8):: offset_temp

      if (io > 0) then
c       --- Just call getelemid.
        call getelemid(zz,offset,naccl,acclzs(0),acclze(0),acclol(0),
     &                 oacclnn(io),oaccloi(oacclii(io)),oacclio(0),id,io)
      else
c       --- Call getelemid for overlap level 1 to get first guess.
        call getelemid(zz,offset,naccl,acclzs(0),acclze(0),acclol(0),
     &                 oacclnn(1),oaccloi(oacclii(1)),oacclio(0),id,1)
c       --- Loop over rest of overlap levels looking for a accl with lower
c       --- z start.
        id_temp = id
        do ii=2,nacclol
          call getelemid(zz,offset_temp,naccl,acclzs(0),acclze(0),acclol(0),
     &                   oacclnn(ii),oaccloi(oacclii(ii)),oacclio(0),id_temp,ii)
          if (acclzs(id_temp)+offset_temp < acclzs(id)+offset) then
            id = id_temp
            offset = offset_temp
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine applyuniformfields(np,ez,bz)
      use InGen
      integer(ISZ):: np
      real(kind=8):: ez(np),bz(np)

c Apply uniform fields.

      ez = ez + ez0
      bz = bz + bz0

      return
      end
c============================================================================ 
      subroutine getbend(np,npz,zp,uzp,gaminv,bendres,bendradi,dtl,dtr,lslice)
      use InGen
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: dtl,dtr
      real(kind=8):: zp(npz),uzp(npz),gaminv(npz)
      logical(ISZ):: lslice
      real(kind=8):: bendres(np),bendradi(np)

c  Gets residence factor and radius for bends.
c  For periodic runs, assumes mesh period length = lattice period length

      integer(ISZ):: ip,iz
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac
      real(kind=8):: vz(npz),dzi(npz)

      if (.not. bends) then
        bendres = 0.
        return
      endif

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

c  Extract the local bend properties for each particle off the general lattice

      if (lslice) then
c       --- find z-cell in which particle lies
        iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c       --- Precalculate these quantities. zl is the min of the two and
c       --- zr is the max. This generalizes the routine, allowing left
c       --- moving particles, vz < 0.
        z1 = zp(1) + vz(1)*dtl
        z2 = zp(1) + vz(1)*dtr
c       --- "left" end of velocity advance step
        zl = min(z1,z2)
        fl = 0.
        if (zl >= cbendzs(iz) .and. zl < cbendze(iz)) fl = 1.
c       --- "right" end of velocity advance step
        zr = max(z1,z2)
        fr = 0.
        if (zr >= cbendzs(iz) .and. zr < cbendze(iz)) fr = 1.
c       --- residence fraction
        frac = fl
        if (fl > fr) frac = (cbendze(iz)-zl)*dzi(1)
        if (fr > fl) frac = (zr-cbendzs(iz))*dzi(1)
        bendres(1) = frac
        bendradi(1) = cbendrc(iz)
        return
      endif

      do ip = 1, np

c       --- find z-cell in which particle lies
        iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c       --- Precalculate these quantities. zl is the min of the two and
c       --- zr is the max. This generalizes the routine, allowing left
c       --- moving particles, vz < 0.
        z1 = zp(ip) + vz(ip)*dtl
        z2 = zp(ip) + vz(ip)*dtr
c       --- "left" end of velocity advance step
        zl = min(z1,z2)
        fl = 0.
        if (zl >= cbendzs(iz) .and. zl < cbendze(iz)) fl = 1.
c       --- "right" end of velocity advance step
        zr = max(z1,z2)
        fr = 0.
        if (zr >= cbendzs(iz) .and. zr < cbendze(iz)) fr = 1.
c       --- residence fraction
        frac = fl
        if (fl > fr) frac = (cbendze(iz)-zl)*dzi(ip)
        if (fr > fl) frac = (zr-cbendzs(iz))*dzi(ip)

        bendres(ip) = frac
        bendradi(ip) = cbendrc(iz)
      enddo

      return
      end
c=============================================================================
      subroutine applybend(np,xp,uzp,npz,bendres,bendradi,m,q,lslice,by)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),uzp(np),bendres(npz),bendradi(npz)
      real(kind=8):: m,q
      logical(ISZ):: lslice
      real(kind=8):: by(np)

c Apply transformation for bending field.

      integer(ISZ):: ip
      real(kind=8):: qi,bres,brad

      if (.not. bends) return

      if (lslice) then
        bres = bendres(1)
        brad = bendradi(1)
        if (bres == 0.) return
      endif
      qi = 1./q
      do ip = 1,np
        if (.not. lslice) then
          bres = bendres(ip)
          brad = bendradi(ip)
        endif
        if (bres > 0.) then
          by(ip) = by(ip) - bres*(m*qi)*uzp(ip)/(brad + xp(ip))
        endif
      enddo

      return
      end
c=============================================================================
      subroutine applyquad(np,xp,yp,npz,zp,uzp,gaminv,dtl,dtr,dt,lslice,
     &                     ex,ey,bx,by)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz),uzp(npz),gaminv(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),bx(np),by(np)

c Apply the transverse focusing element quad.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey   quadrupole electric field
c   bx, by   quadrupole magnetic field

      integer(ISZ):: io,ip,iz
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac,xpmqoff,ypmqoff
      real(kind=8):: dodec
      real(kind=8):: vz(npz),dzi(npz)

      if (.not. quads) return

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

      do io=1,nquadol

        if (lslice) then
c         --- find z-cell in which particle lies
          iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c         --- Precalculate these quantities. zl is the min of the two and
c         --- zr is the max. This generalizes the routine, allowing left
c         --- moving particles, vz < 0.
          z1 = zp(1) + vz(1)*dtl
          z2 = zp(1) + vz(1)*dtr
c         --- "left" end of velocity advance step
          zl = min(z1,z2)
          fl = 0.
          if (zl >= cquadzs(iz,io) .and. zl < cquadze(iz,io)) fl = 1.
c         --- "right" end of velocity advance step
          zr = max(z1,z2)
          fr = 0.
          if (zr >= cquadzs(iz,io) .and. zr < cquadze(iz,io)) fr = 1.
c         --- residence fraction
          frac = fl
          if (fl > fr) frac = (cquadze(iz,io)-zl)*dzi(1)
          if (fr > fl) frac = (zr-cquadzs(iz,io))*dzi(1)
c         --- if frac is zero, no field will be applied
          if (frac == 0.) cycle
        endif
        do ip = 1, np

          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c           --- Precalculate these quantities. zl is the min of the two and
c           --- zr is the max. This generalizes the routine, allowing left
c           --- moving particles, vz < 0.
            z1 = zp(ip) + vz(ip)*dtl
            z2 = zp(ip) + vz(ip)*dtr
c           --- "left" end of velocity advance step
            zl = min(z1,z2)
            fl = 0.
            if (zl >= cquadzs(iz,io) .and. zl < cquadze(iz,io)) fl = 1.
c           --- "right" end of velocity advance step
            zr = max(z1,z2)
            fr = 0.
            if (zr >= cquadzs(iz,io) .and. zr < cquadze(iz,io)) fr = 1.
c           --- residence fraction
            frac = fl
            if (fl > fr) frac = (cquadze(iz,io)-zl)*dzi(ip)
            if (fr > fl) frac = (zr-cquadzs(iz,io))*dzi(ip)
          endif
      
c         --- set the field
          if (frac > 0.) then
            xpmqoff = xp(ip) - cqoffx(iz,io)
            ypmqoff = yp(ip) - cqoffy(iz,io)
            bx(ip) = bx(ip) + cquaddb(iz,io)*frac*ypmqoff
            by(ip) = by(ip) + cquaddb(iz,io)*frac*xpmqoff
            ex(ip) = ex(ip) + cquadde(iz,io)*frac*xpmqoff
            ey(ip) = ey(ip) - cquadde(iz,io)*frac*ypmqoff
            dodec = quaddo(cquadid(iz,io))
            if (dodec .ne. 0.) then
              ex(ip) = ex(ip) + dodec*cquadde(iz,io)*frac*
     &                 (6.*xpmqoff**5 - 60.*xpmqoff**3*ypmqoff**2 +
     &                 30.*xpmqoff*ypmqoff**4)
              ey(ip) = ey(ip) - dodec*cquadde(iz,io)*frac*
     &                 (6.*ypmqoff**5 - 60.*ypmqoff**3*xpmqoff**2 +
     &                 30.*ypmqoff*xpmqoff**4)
            endif
          endif
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine applydipo(np,npz,zp,uzp,gaminv,dtl,dtr,dt,lslice,
     &                     ex,ey,bx,by)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: zp(npz),uzp(npz),gaminv(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),bx(np),by(np)

c Apply the transverse bending element dipo.
c Input:
c   np       number of particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey   dipole electric field
c   bx, by   dipole magnetic field

      integer(ISZ):: io,ip,iz
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac
      real(kind=8):: vz(npz),dzi(npz)

      if (.not. dipos) return

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

      do io=1,ndipool

        if (lslice) then
c         --- find z-cell in which particle lies
          iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c         --- Precalculate these quantities. zl is the min of the two and
c         --- zr is the max. This generalizes the routine, allowing left
c         --- moving particles, vz < 0.
          z1 = zp(1) + vz(1)*dtl
          z2 = zp(1) + vz(1)*dtr
c         --- "left" end of velocity advance step
          zl = min(z1,z2)
          fl = 0.
          if (zl >= cdipozs(iz,io) .and. zl < cdipoze(iz,io)) fl = 1.
c         --- "right" end of velocity advance step
          zr = max(z1,z2)
          fr = 0.
          if (zr >= cdipozs(iz,io) .and. zr < cdipoze(iz,io)) fr = 1.
c         --- residence fraction
          frac = fl
          if (fl > fr) frac = (cdipoze(iz,io)-zl)*dzi(1)
          if (fr > fl) frac = (zr-cdipozs(iz,io))*dzi(1)
c         --- if frac is zero, no field will be applied
          if (frac == 0.) cycle
        endif
        do ip = 1, np

          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c           --- Precalculate these quantities. zl is the min of the two and
c           --- zr is the max. This generalizes the routine, allowing left
c           --- moving particles, vz < 0.
            z1 = zp(ip) + vz(ip)*dtl
            z2 = zp(ip) + vz(ip)*dtr
c           --- "left" end of velocity advance step
            zl = min(z1,z2)
            fl = 0.
            if (zl >= cdipozs(iz,io) .and. zl < cdipoze(iz,io)) fl = 1.
c           --- "right" end of velocity advance step
            zr = max(z1,z2)
            fr = 0.
            if (zr >= cdipozs(iz,io) .and. zr < cdipoze(iz,io)) fr = 1.
c           --- residence fraction
            frac = fl
            if (fl > fr) frac = (cdipoze(iz,io)-zl)*dzi(ip)
            if (fr > fl) frac = (zr-cdipozs(iz,io))*dzi(ip)
          endif

c         --- set the field
          if (frac > 0.) then
            by(ip) = by(ip) + cdipoby(iz,io)*frac
            bx(ip) = bx(ip) + cdipobx(iz,io)*frac
            ex(ip) = ex(ip) + cdipoex(iz,io)*frac
            ey(ip) = ey(ip) + cdipoey(iz,io)*frac
          endif
        enddo
      enddo

      return
      end
c=========================================================================== 
      subroutine sledgcor(np,xp,yp,npz,zp,zpo,uxp,uyp,uzp,gaminv,m,q,lslice)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np), yp(np), zp(npz), zpo(npz) 
      real(kind=8):: uxp(np), uyp(np), uzp(npz), gaminv(npz) 
      real(kind=8):: q, m
      logical(ISZ):: lslice

c Applies position and velocity jump corrections on entry/exit of dipoles 
c with slanted faces.  (Must be called after all z advancement is done this 
c timestep, so that entry/exit of the bend is detected exactly once.)    
c Accounts for the fact that the jumped velocities act on the particles for 
c less than the entire andvance step on entry/exit of the dipole.  

      integer(ISZ):: io,ip,iz
      real(kind=8):: vzi(np),qoverm

      if (.not.dipos) return

      qoverm = q/m
      if (npz == 1) then
        vzi = 1./dvnz(uzp(1)*gaminv(1))
      else
        vzi = 1./dvnz(uzp*gaminv)
      endif

      do io=1,ndipool

        if (lslice) then
          iz = 0
          if (.not.((zpo(1) < cdipozs(iz,io) .and. zp(1) >= cdipozs(iz,io)) .or.
     &              (zpo(1) < cdipoze(iz,io) .and. zp(1) >= cdipoze(iz,io))))
     &      cycle
        endif
        do ip = 1,np
          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = min(nzl,int(max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)))
          endif
          if ((zpo(ip) < cdipozs(iz,io)) .and. (zp(ip) >= cdipozs(iz,io))) then 
c           --- entered dipole
            uxp(ip) = uxp(ip) + xp(ip)*cdipoby(iz,io)*cdipota(iz,io)*qoverm
            uyp(ip) = uyp(ip) - yp(ip)*cdipoby(iz,io)*cdipota(iz,io)*qoverm
            xp(ip)  = xp(ip) + (zp(ip)-cdipozs(iz,io))*vzi(ip)*
     &                         xp(ip)*cdipoby(iz,io)*cdipota(iz,io)*qoverm
            yp(ip)  = yp(ip) - (zp(ip)-cdipozs(iz,io))*vzi(ip)*
     &                         yp(ip)*cdipoby(iz,io)*cdipota(iz,io)*qoverm
          elseif ((zpo(ip) <  cdipoze(iz,io)) .and.
     &             (zp(ip) >= cdipoze(iz,io))) then  
c           --- exited dipole
            uxp(ip) = uxp(ip) + xp(ip)*cdipoby(iz,io)*cdipotb(iz,io)*qoverm
            uyp(ip) = uyp(ip) - yp(ip)*cdipoby(iz,io)*cdipotb(iz,io)*qoverm
            xp(ip)  = xp(ip) + (zp(ip)-cdipoze(iz,io))*vzi(ip)*
     &                         xp(ip)*cdipoby(iz,io)*cdipotb(iz,io)*qoverm
            yp(ip)  = yp(ip) - (zp(ip)-cdipoze(iz,io))*vzi(ip)*
     &                         yp(ip)*cdipoby(iz,io)*cdipotb(iz,io)*qoverm
          endif 
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine applysext(np,xp,yp,npz,zp,uzp,gaminv,dtl,dtr,dt,lslice,
     &                     ex,ey,bx,by)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz),uzp(npz),gaminv(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),bx(np),by(np)

c Apply the transverse element sext.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey   sextupole electric field
c   bx, by   sextupole magnetic field

      integer(ISZ):: io,ip,iz
      real(kind=8):: xpmqoff,ypmqoff
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac
      real(kind=8):: vz(npz),dzi(npz)

      if (.not. sexts) return

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

      do io=1,nsextol

        if (lslice) then
c         --- find z-cell in which particle lies
          iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c         --- Precalculate these quantities. zl is the min of the two and
c         --- zr is the max. This generalizes the routine, allowing left
c         --- moving particles, vz < 0.
          z1 = zp(1) + vz(1)*dtl
          z2 = zp(1) + vz(1)*dtr
c         --- "left" end of velocity advance step
          zl = min(z1,z2)
          fl = 0.
          if (zl >= csextzs(iz,io) .and. zl < csextze(iz,io)) fl = 1.
c         --- "right" end of velocity advance step
          zr = max(z1,z2)
          fr = 0.
          if (zr >= csextzs(iz,io) .and. zr < csextze(iz,io)) fr = 1.
c         --- residence fraction
          frac = fl
          if (fl > fr) frac = (csextze(iz,io)-zl)*dzi(1)
          if (fr > fl) frac = (zr-csextzs(iz,io))*dzi(1)
c         --- if frac is zero, no field will be applied
          if (frac == 0.) cycle
        endif
        do ip = 1, np

          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c           --- Precalculate these quantities. zl is the min of the two and
c           --- zr is the max. This generalizes the routine, allowing left
c           --- moving particles, vz < 0.
            z1 = zp(ip) + vz(ip)*dtl
            z2 = zp(ip) + vz(ip)*dtr
c           --- "left" end of velocity advance step
            zl = min(z1,z2)
            fl = 0.
            if (zl >= csextzs(iz,io) .and. zl < csextze(iz,io)) fl = 1.
c           --- "right" end of velocity advance step
            zr = max(z1,z2)
            fr = 0.
            if (zr >= csextzs(iz,io) .and. zr < csextze(iz,io)) fr = 1.
c           --- residence fraction
            frac = fl
            if (fl > fr) frac = (csextze(iz,io)-zl)*dzi(ip)
            if (fr > fl) frac = (zr-csextzs(iz,io))*dzi(ip)
          endif

c         --- set the field
          if (frac > 0.) then
            xpmqoff = xp(ip)
            ypmqoff = yp(ip)
            bx(ip) = bx(ip) + csextdb(iz,io)*frac*3.*(xpmqoff**2 - ypmqoff**2)
            by(ip) = by(ip) + csextdb(iz,io)*frac*(-6.)*xpmqoff*ypmqoff
            ex(ip) = ex(ip) + csextde(iz,io)*frac*3.*(xpmqoff**2 - ypmqoff**2)
            ey(ip) = ey(ip) + csextde(iz,io)*frac*(-6.)*xpmqoff*ypmqoff
          endif
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine applyhele(np,xp,yp,npz,zp,uzp,gaminv,dtl,dtr,dt,lslice,
     &                     ex,ey,ez,bx,by,bz)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz),uzp(npz),gaminv(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),ez(np),bx(np),by(np),bz(np)

c Apply the transverse element hele.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   uzp      massless momentum of the particles
c   gaminv   one over gamma of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey, ez   electric field
c   bx, by, bz   magnetic field

      integer(ISZ):: io,ip,iz,iele,ii
      real(kind=8):: z1,z2,zl,fl,zr,fr,frac
      real(kind=8):: vz(npz),dzi(npz)
      real(kind=8):: n,v,f,rpow,cosnt,sinnt
      real(kind=8):: xpmqoff,ypmqoff
      integer(ISZ):: tiele(np)
      real(kind=8):: tfrac(np),txpmqoff(np),typmqoff(np)
      real(kind=8):: trpmqoff(np),ttpmqoff(np)

      if (.not. heles) return

      vz  = uzp*gaminv
c     --- Note that the absolute value is taken since dzi is used to scale
c     --- the fraction of the step inside the element and the sign is
c     --- not needed. This only matters if vz < 0.
      dzi = abs(1./dvnz(vz*(dtr-dtl)))

      do io=1,nheleol

c       --- calculate coordinates, element indices, and residence
c       --- fractions over particle block and store in temp arrays
c       --- These factors will be used repeatedly when the fields
c       --- are accumulated at each particle resulting from
        if (lslice) then
c         --- find z-cell in which particle lies
          iz = max(0., (zp(1) - zlmin - zlframe)*dzli + 0.5)
c         --- Precalculate these quantities. zl is the min of the two and
c         --- zr is the max. This generalizes the routine, allowing left
c         --- moving particles, vz < 0.
          z1 = zp(1) + vz(1)*dtl
          z2 = zp(1) + vz(1)*dtr
c         --- "left" end of velocity advance step
          zl = min(z1,z2)
          fl = 0.
          if (zl >= chelezs(iz,io) .and. zl < cheleze(iz,io)) fl = 1.
c         --- "right" end of velocity advance step
          zr = max(z1,z2)
          fr = 0.
          if (zr >= chelezs(iz,io) .and. zr < cheleze(iz,io)) fr = 1.
c         --- residence fraction
          frac = fl
          if (fl > fr) frac = (cheleze(iz,io)-zl)*dzi(1)
          if (fr > fl) frac = (zr-chelezs(iz,io))*dzi(1)
c         --- if frac is zero, no field will be applied
          if (frac == 0.) cycle
        endif
        do ip = 1, np

          if (.not. lslice) then
c           --- find z-cell in which particle lies
            iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
c           --- Precalculate these quantities. zl is the min of the two and
c           --- zr is the max. This generalizes the routine, allowing left
c           --- moving particles, vz < 0.
            z1 = zp(ip) + vz(ip)*dtl
            z2 = zp(ip) + vz(ip)*dtr
c           --- "left" end of velocity advance step
            zl = min(z1,z2)
            fl = 0.
            if (zl >= chelezs(iz,io) .and. zl < cheleze(iz,io)) fl = 1.
c           --- "right" end of velocity advance step
            zr = max(z1,z2)
            fr = 0.
            if (zr >= chelezs(iz,io) .and. zr < cheleze(iz,io)) fr = 1.
c           --- residence fraction
            frac = fl
            if (fl > fr) frac = (cheleze(iz,io)-zl)*dzi(ip)
            if (fr > fl) frac = (zr-chelezs(iz,io))*dzi(ip)
          endif

c         --- find element index
          iele = cheleid(iz,io)
          tiele(ip) = iele
          tfrac(ip) = frac
c         --- set the field
          if (frac > 0.) then
c           --- x,y coordinates about multipole center
            xpmqoff = xp(ip) - heleox(iele)
            ypmqoff = yp(ip) - heleoy(iele)
            txpmqoff(ip) = xpmqoff
            typmqoff(ip) = ypmqoff
c           --- r,theta coordinates about multipole center
            trpmqoff(ip) = sqrt(xpmqoff**2 + ypmqoff**2)
            ttpmqoff(ip) = atan2(ypmqoff,dvnz(xpmqoff))
          endif
        enddo
c       --- accumulate field contributions at each particle in
c       --- block looping (outermost) over each multipole component
c       --- for vectorization
        do ii = 1, nhmlt
c         --- loop over particle block using stored coordinates, element
c         --- indices, and residence fractions
          do ip = 1,np
            if (tfrac(ip) == 0.) cycle
            iele = tiele(ip)
            n = hele_n(ii,iele)
            v = hele_v(ii,iele)
c           --- electric multipoles
            if (ii <= helene(iele)) then
              if (n == 0) then
c               --- accelerating field
c               --- Note that this really should not be used since the
c               --- accl element does a better job.  The transverse
c               --- fields should only be kicks at the entrance and exit.
                f = 1./(2.*v+2.)
                rpow = dvnz(trpmqoff(ip))**(2*v)
                ex(ip) = ex(ip)-tfrac(ip)*heleep(ii,iele)*rpow*f*txpmqoff(ip)
                ey(ip) = ey(ip)-tfrac(ip)*heleep(ii,iele)*rpow*f*typmqoff(ip)
                ez(ip) = ez(ip)+tfrac(ip)*heleae(ii,iele)*rpow
              else
c               --- rest of components
                f = 1. + 2.*v/n
                cosnt = cos(n*(ttpmqoff(ip)-helepe(ii,iele)))
                sinnt = sin(n*(ttpmqoff(ip)-helepe(ii,iele)))
                rpow = dvnz(trpmqoff(ip))**(n-2+2*v)
                ex(ip) = ex(ip) + tfrac(ip)*heleae(ii,iele)*rpow*
     &                         (f*txpmqoff(ip)*cosnt + typmqoff(ip)*sinnt)
                ey(ip) = ey(ip) + tfrac(ip)*heleae(ii,iele)*rpow*
     &                         (f*typmqoff(ip)*cosnt - txpmqoff(ip)*sinnt)
              endif
            endif
c           --- magnetic multipoles
            if (ii <= helenm(iele)) then
              if (n == 0) then
c               --- Solenoid field
c               --- The transverse fields should only be kicks at the entrance
c               --- and exit.
                f = 1./(2.*v+2.)
                rpow = dvnz(trpmqoff(ip))**(2*v)
                bx(ip) = bx(ip)-tfrac(ip)*helemp(ii,iele)*rpow*f*txpmqoff(ip)
                by(ip) = by(ip)-tfrac(ip)*helemp(ii,iele)*rpow*f*typmqoff(ip)
                bz(ip) = bz(ip)+tfrac(ip)*heleam(ii,iele)*rpow
              else
c               --- rest of components
                f = 1. + 2.*v/n
                cosnt = cos(n*(ttpmqoff(ip)-helepm(ii,iele)))
                sinnt = sin(n*(ttpmqoff(ip)-helepm(ii,iele)))
                rpow = dvnz(trpmqoff(ip))**(n-2+2*v)
                bx(ip) = bx(ip) + tfrac(ip)*heleam(ii,iele)*rpow*
     &                         (f*txpmqoff(ip)*sinnt - typmqoff(ip)*cosnt)
                by(ip) = by(ip) + tfrac(ip)*heleam(ii,iele)*rpow*
     &                         (f*typmqoff(ip)*sinnt + txpmqoff(ip)*cosnt)
              endif
            endif
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine applyemlt(np,xp,yp,npz,zp,dtl,dtr,dt,lslice,ex,ey,ez)
      use Lattice
      use LatticeInternal
      use Mult_data
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),ez(np)

c Apply the element emlt.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   ex, ey, ez   electric field

c  See HIF-note 96-10, equation D17, with E replacing B and
c  -n*psi replacing psi.
c  Field is of the following form, where n and v are integers, n>0, v>=0
c    Er = - sum_v  [ E'0v/(2*v+2)*r**(2*v+1) ] +
c         sum_nv [ Env*(1+2*v/n)*r**(n-1+2*v)*cos(n*(theta-psi)) ]
c    Et = - sum_nv [ Env*r**(n-1+2*v)*sin(n*(theta-psi)) ]
c    Ez = sum_v  [ E0v*r**(2*v) ] +
c         sum_nv [ E'nv/n*r**(2*v+n)*cos(n*(theta-psi)) ]
c
c    Ex = Er*cos(theta) - Et*sin(theta)
c    Ey = Er*sin(theta) + Et*cos(theta)
c
c  Note that the total phase angle has two parts, one associated with the
c  multipole data, and one associated with the lattice element.  This allows
c  an element with multipoles with different phases to be physically rotated
c  by changing only one variable, the phase associated with the element.
c  This allows the same set of data to be used by lattice elements with
c  different errors in the angle.
c  Additional scale factors are included in the field, allowing different
c  elements to use the same data set scaled to different values.  The scale
c  factors are added together (in mltlocat) and multiply the applied field.
c  The addition is done in mltlocat so that it is only done once per particle
c  and not once per particle per multipole component.

      integer(ISZ):: io,in,ip
      real(kind=8):: n,v,f,fz,coeffp,coeff,rpow,alpha,alphap,cosnt,sinnt
      integer(ISZ):: izm(np),imm(np)
      real(kind=8):: xxx(np),yyy(np),zzz(np),ttt(np)
      real(kind=8):: wzm(np),rr(np),tt(np)
      real(kind=8):: mltph(np),mltsf(np)
      real(kind=8):: tex(np),tey(np),tez(np)

      if (.not. emlts .or. .not. linemlt(0)) return

      do io=1,nemltol
        if (.not. linemlt(io)) cycle

c       --- Get location of particle relative to moment data.
        call mltlocat(np,xp,yp,npz,zp,nemltsets,dzemlt,nzemlt,nzl,nzlmax,
     &                nemltol,io,cemltzs,cemltze,cemltph,cemltsf,cemltsc,
     &                cemltim,cemltox,cemltoy,imm,izm,wzm,xxx,yyy,zzz,ttt,
     &                rr,tt,mltph,mltsf,zlmin,zlframe,dzli,
     &                bends,linbend,cbendzs,cbendze,cbendrc,lslice)

c       --- zero temporary field arrays
        tex = 0.
        tey = 0.
        tez = 0.

c       --- accumulate the E field of the electrostatic multipoles
        do in=1,nesmult
          n = emlt_n(in)
          v = emlt_v(in)
          if (n == 0) then
c           --- Apply accelerating field and it's pseudomultipoles.
            f = 1./(2.*v+2.)
            do ip=1,np
              if (imm(ip) > 0) then
                coeffp = mltsf(ip)*
     &                     (esemltp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      esemltp(izm(ip)+1,in,imm(ip))*wzm(ip))
                coeff = mltsf(ip)*
     &                     (esemlt(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemlt(izm(ip)+1,in,imm(ip))*wzm(ip))
                rpow = dvnz(rr(ip))**(2*v)
                tex(ip) = tex(ip) - coeffp*rpow*f*xxx(ip)
                tey(ip) = tey(ip) - coeffp*rpow*f*yyy(ip)
                tez(ip) = tez(ip) + coeff*rpow
              endif
            enddo
          else
c           --- Apply the rest of the multipoles.
            f = 1. + 2.*v/n
            fz = 1./n
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = tt(ip) - mltph(ip) -
     &                      esemltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      esemltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (esemltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (esemlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      esemlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (esemltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(n*alpha)
                sinnt = sin(n*alpha)
                rpow = dvnz(rr(ip))**(n-2+2*v)
                tex(ip)=tex(ip) + coeff*rpow*(f*xxx(ip)*cosnt + yyy(ip)*sinnt)
                tey(ip)=tey(ip) + coeff*rpow*(f*yyy(ip)*cosnt - xxx(ip)*sinnt)
                tez(ip)=tez(ip) + fz*dvnz(rr(ip))**(n+2*v)*
     &                            (coeffp*cosnt - alphap*coeff*sinnt)
              endif
            enddo
          endif
        enddo

c       --- do coordinate transform on fields back to warped coordinates
        if (bends .and. linbend) then
          do ip=1,np
            ex(ip) = ex(ip) + tex(ip)*cos(ttt(ip)) - tez(ip)*sin(ttt(ip))
            ey(ip) = ey(ip) + tey(ip)
            ez(ip) = ez(ip) + tex(ip)*sin(ttt(ip)) + tez(ip)*cos(ttt(ip))
          enddo
        else
          do ip=1,np
            ex(ip) = ex(ip) + tex(ip)
            ey(ip) = ey(ip) + tey(ip)
            ez(ip) = ez(ip) + tez(ip)
          enddo
        endif

      enddo

      return
      end
c=============================================================================
      subroutine applymmlt(np,xp,yp,npz,zp,dtl,dtr,dt,lslice,bx,by,bz)
      use Lattice
      use LatticeInternal
      use Mult_data
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      real(kind=8):: dtl,dtr,dt
      logical(ISZ):: lslice
      real(kind=8):: bx(np),by(np),bz(np)

c Apply the element mmlt.
c Input:
c   np       number of particles
c   xp       x position of the particles
c   yp       y position of the particles
c   npz      number of z position data values (must be either 1 or == np)
c   zp       z position of the particles
c   dtl      size of left half of time step
c   dtr      size of right half of time step
c   dt       full step size
c Output:
c   bx, by, bz   magnetic field

c  See HIF-note 96-10, equation D17, with -n*psi-pi/2 replacing psi.
c  Field is of the following form, where n and v are integers, n>0, v>=0
c    Br = - sum_v  [ B'0v/(2*v+2)*r**(2*v+1) ] +
c         sum_nv [ Bnv*(1+2*v/n)*r**(n-1+2*v)*sin(n*(theta-psi)) ]
c    Bt = sum_nv [ Bnv*r**(n-1+2*v)*cos(n*(theta-psi)) ]
c    Bz = sum_v  [ B0v*r**(2*v) ] +
c         sum_nv [ B'nv/n*r**(2*v)*sin(n*(theta-psi)) ]
c
c    Bx = Br*cos(theta) - Bt*sin(theta)
c    By = Br*sin(theta) + Bt*cos(theta)
c
c  Note that the total phase angle has two parts, one associated with the
c  multipole data, and one associated with the lattice element.  This allows
c  an element with multipoles with different phases to be physically rotated
c  by changing only one variable, the phase associated with the element.
c  This allows the same set of data to be used by lattice elements with
c  different errors in the angle.
c  Additional scale factors are included in the field, allowing different
c  elements to use the same data set scaled to different values.  The scale
c  factors are added together (in mltlocat) and multiply the applied field.
c  The addition is done in mltlocat so that it is only done once per particle
c  and not once per particle per multipole component.

      integer(ISZ):: io,in,ip
      real(kind=8):: n,v,f,fz,coeffp,coeff,rpow,alpha,alphap,cosnt,sinnt
      integer(ISZ):: izm(np),imm(np)
      real(kind=8):: xxx(np),yyy(np),zzz(np),ttt(np)
      real(kind=8):: wzm(np),rr(np),tt(np)
      real(kind=8):: mltph(np),mltsf(np)
      real(kind=8):: tbx(np),tby(np),tbz(np)

      if (.not. mmlts .or. .not. linmmlt(0)) return

      do io=1,nmmltol
        if (.not. linmmlt(io)) cycle

c       --- Get location of particle relative to moment data.
        call mltlocat(np,xp,yp,npz,zp,nmmltsets,dzmmlt,nzmmlt,nzl,nzlmax,
     &                nmmltol,io,cmmltzs,cmmltze,cmmltph,cmmltsf,cmmltsc,
     &                cmmltim,cmmltox,cmmltoy,imm,izm,wzm,xxx,yyy,zzz,ttt,
     &                rr,tt,mltph,mltsf,zlmin,zlframe,dzli,
     &                bends,linbend,cbendzs,cbendze,cbendrc,lslice)

c       --- zero the temporary field arrays
        tbx = 0.
        tby = 0.
        tbz = 0.

c       --- accumulate the B field of the magnetostatic multipoles
        do in=1,nmsmult
          n = mmlt_n(in)
          v = mmlt_v(in)
          if (n == 0) then
c           --- Apply solenoidal field component and it's pseudomultipoles.
            f = 1./(2.*v+2.)
            do ip=1,np
              if (imm(ip) > 0) then
                coeffp = mltsf(ip)*
     &                    (msmmltp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                     msmmltp(izm(ip)+1,in,imm(ip))*wzm(ip))
                coeff = mltsf(ip)*
     &                    (msmmlt(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                     msmmlt(izm(ip)+1,in,imm(ip))*wzm(ip))
                rpow = dvnz(rr(ip))**(2*v)
                tbx(ip) = tbx(ip) - coeffp*rpow*f*xxx(ip)
                tby(ip) = tby(ip) - coeffp*rpow*f*yyy(ip)
                tbz(ip) = tbz(ip) + coeff*rpow
              endif
            enddo
          else
c           --- Apply the rest of the multipoles.
            f = 1. + 2.*v/n
            fz = 1./n
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = tt(ip) - mltph(ip) -
     &                      msmmltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      msmmltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (msmmltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (msmmlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      msmmlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (msmmltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(n*alpha)
                sinnt = sin(n*alpha)
                rpow = dvnz(rr(ip))**(n-2+2*v)
                tbx(ip)=tbx(ip) + coeff*rpow*(f*xxx(ip)*sinnt - yyy(ip)*cosnt)
                tby(ip)=tby(ip) + coeff*rpow*(f*yyy(ip)*sinnt + xxx(ip)*cosnt)
                tbz(ip)=tbz(ip) + fz*dvnz(rr(ip))**(n+2*v)*
     &                            (coeffp*sinnt + alphap*coeff*cosnt)
              endif
            enddo
          endif
        enddo

c       --- do coordinate transform on fields back to warped coordinates
        if (bends .and. linbend) then
          do ip=1,np
            bx(ip) = bx(ip) + tbx(ip)*cos(ttt(ip)) - tbz(ip)*sin(ttt(ip))
            by(ip) = by(ip) + tby(ip)
            bz(ip) = bz(ip) + tbx(ip)*sin(ttt(ip)) + tbz(ip)*cos(ttt(ip))
          enddo
        else
          do ip=1,np
            bx(ip) = bx(ip) + tbx(ip)
            by(ip) = by(ip) + tby(ip)
            bz(ip) = bz(ip) + tbz(ip)
          enddo
        endif

      enddo

      return
      end
c=============================================================================
      subroutine mltlocat(np,xp,yp,npz,zp,nmltsets,dzmlt,nzmlt,nzl,nzlmax,
     &                    nol,io,cmltzs,cmltze,cmltph,cmltsf,cmltsc,
     &                    cmltim,cmltox,cmltoy,
     &                    imm,izm,wzm,xxx,yyy,zzz,ttt,rr,tt,
     &                    mltph,mltsf,zlmin,zlframe,dzli,bends,linbend,
     &                    cbendzs,cbendze,cbendrc,lslice)
      integer(ISZ):: np,npz,nmltsets,nzl,nzlmax,nol,io
      real(kind=8):: xp(np),yp(np),zp(npz)
      real(kind=8):: dzmlt(nmltsets)
      integer(ISZ):: nzmlt(nmltsets)
      real(kind=8):: cmltzs(0:nzlmax,nol),cmltze(0:nzlmax,nol)
      real(kind=8):: cmltph(0:nzlmax,nol)
      real(kind=8):: cmltsf(0:nzlmax,nol),cmltsc(0:nzlmax,nol)
      integer(ISZ):: cmltim(0:nzlmax,nol)
      real(kind=8):: cmltox(0:nzlmax,nol),cmltoy(0:nzlmax,nol)
      integer(ISZ):: imm(np),izm(np)
      real(kind=8):: wzm(np),xxx(np),yyy(np),zzz(np)
      real(kind=8):: ttt(np),rr(np),tt(np),mltph(np),mltsf(np)
      logical(ISZ):: bends,linbend
      real(kind=8):: cbendzs(0:nzlmax),cbendze(0:nzlmax),cbendrc(0:nzlmax)
      real(kind=8):: zlmin,zlframe,dzli
      logical(ISZ):: lslice

c     --- Calculate grid location and polar coordinates of particles
c     --- for mlt elements.  Includes change of coordinates out of
c     --- warped coordinates (when a bend is between the
c     --- particle and the nearest mult).
c  Assumes that at most one bend will be between particle and the mult center.

      integer(ISZ):: iz,ip
      real(kind=8):: zz,cent,dzi,rrr,bzs,bze,bendrci

      if (lslice) then
        zz = zp(1)
        iz = 0
      endif
      if (.not. bends .or. .not. linbend) then
        do ip=1,np
          if (.not. lslice) then
            zz = zp(ip)
c           --- find z-cell in which particle lies
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
          endif
          imm(ip) = cmltim(iz,io)
          dzi = 1./dzmlt(imm(ip))
          mltph(ip) = cmltph(iz,io)
          mltsf(ip) = cmltsc(iz,io) + cmltsf(iz,io)
c         --- find z location relative to multipole data
          izm(ip) = int((zz-cmltzs(iz,io))*dzi)
          wzm(ip) =     (zz-cmltzs(iz,io))*dzi - izm(ip)
          if (zz < cmltzs(iz,io)) then
            izm(ip) = 0
            wzm(ip) = 0
            imm(ip) = 0
          elseif (zz > cmltze(iz,io)) then
            izm(ip) = nzmlt(imm(ip))
            wzm(ip) = 0
            imm(ip) = 0
          endif
c         --- Set temporary coordinates and add offset if particle within
c         --- an element.
          if (imm(ip) > 0) then
            xxx(ip) = xp(ip) - cmltox(iz,io)
            yyy(ip) = yp(ip) - cmltoy(iz,io)
            zzz(ip) = zz
            rr(ip) = sqrt(xxx(ip)**2 + yyy(ip)**2)
            tt(ip) = atan2(yyy(ip),dvnz(xxx(ip)))
          endif
        enddo
      else
        do ip=1,np
          if (.not. lslice) then
            zz = zp(ip)
c           --- find z-cell in which particle lies
            iz = max(0., (zz - zlmin - zlframe)*dzli)
          endif
c         --- set temporaries
          imm(ip) = cmltim(iz,io)
          dzi = 1./dzmlt(imm(ip))
          cent = 0.5*(cmltzs(iz,io) + cmltze(iz,io))
          xxx(ip) = xp(ip)
          zzz(ip) = zz
          if (cmltzs(iz,io) <= zz .and. zz <= cmltze(iz,io)) then
            mltph(ip) = cmltph(iz,io)
            mltsf(ip) = cmltsc(iz,io) + cmltsf(iz,io)
c           --- apply coordinate change left of mult, using index of (iz)
            if (zz <= cent) then
              bendrci = 1./cbendrc(iz)
c             --- find end of bend: if bend extends beyond center of mult, use
c             --- that as the end of the bend since moments are in coordinate
c             --- system at the center of the mult
              if (cbendze(iz) < cent) then
                bze = cbendze(iz)
              else
                bze = cent
              endif
c             --- if particle is not in bend and full bend is between particle
c             --- and mult
              rrr = xp(ip) + cbendrc(iz)
              if (zz < cbendzs(iz) .and. cbendzs(iz) < cent) then
                ttt(ip) = (bze - cbendzs(iz))*bendrci
                zzz(ip) = bze - rrr*sin(ttt(ip)) -
     &                    (cbendzs(iz) - zz)*cos(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(iz) -
     &                    (cbendzs(iz) - zz)*sin(ttt(ip))
c             --- if particle is in bend
              elseif (zz < bze .and. cbendzs(iz) < cent) then
                ttt(ip) = (bze - zz)*bendrci
                zzz(ip) = bze - rrr*sin(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(iz)
              endif
c           --- apply coordinate change right of mult, using index of (iz+1)
            elseif (zz >= cent) then
              bendrci = 1./cbendrc(iz+1)
c             --- find start of bend: if bend extends beyond center of mult, use
c             --- that as the start of the bend since moments are in coordinate
c             --- system at the center of the mult
              if (cbendzs(iz+1) < cent) then
                bzs = cent
              else
                bzs = cbendzs(iz+1)
              endif
c             --- if particle is not in bend and full bend is between particle
c             --- and mult
              rrr = xp(ip) + cbendrc(iz+1)
              if (zz > cbendze(iz+1) .and. cbendze(iz+1) > cent) then
                ttt(ip) =  - (cbendze(iz+1) - bzs)*bendrci
                zzz(ip) = bzs - rrr*sin(ttt(ip)) +
     &                    (zz - cbendze(iz+1))*cos(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(iz+1) +
     &                    (zz - cbendze(iz+1))*sin(ttt(ip))
c             --- if particle is in bend
              elseif (zz > bzs .and. cbendze(iz+1) > cent) then
                ttt(ip) = - (zz - bzs)*bendrci
                zzz(ip) = bzs - rrr*sin(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(iz+1)
              endif
            endif
          endif
c         --- end of coordinate transformation

          izm(ip) = int((zzz(ip)-cmltzs(iz,io))*dzi)
          wzm(ip) =     (zzz(ip)-cmltzs(iz,io))*dzi - izm(ip)
          if (zzz(ip) < cmltzs(iz,io)) then
            izm(ip) = 0
            wzm(ip) = 0
            imm(ip) = 0
          elseif (zzz(ip) > cmltze(iz,io)) then
            izm(ip) = nzmlt(imm(ip))
            wzm(ip) = 0
            imm(ip) = 0
          endif
c         --- Add offset
          if (imm(ip) > 0) then
            xxx(ip) = xxx(ip) - cmltox(iz,io)
            yyy(ip) = yp(ip) -  cmltoy(iz,io)
            rr(ip) = sqrt(xxx(ip)**2 + yyy(ip)**2)
            tt(ip) = atan2(yyy(ip),dvnz(xxx(ip)))
          endif
        enddo
      endif

      return
      end
c=========================================================================== 
      subroutine applybgrd(np,xp,yp,npz,zp,lslice,bx,by,bz)
      use Lattice
      use LatticeInternal
      use BGRDdata
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      logical(ISZ):: lslice
      real(kind=8):: bx(np),by(np),bz(np)

c  Sets magnetic field for particles from data sets
c  containing Bx, By, and Bz on 3-D grids.
c  Calling arguments:
c     np          number of particles
c     npz         number of z position data values (must be either 1 or == np)
c     xp,yp,zp    coordinates of particles
c  Output:
c     bx,by,bz       magnetic field

c  The field is:
c     Bx = u0 * v0 * w0 * bgrdbx(i  ,j  ,k  ,bgrdid)
c        + u1 * v0 * w0 * bgrdbx(i+1,j  ,k  ,bgrdid)
c        + u0 * v1 * w0 * bgrdbx(i  ,j+1,k  ,bgrdid)
c        + ...

      integer(ISZ):: io,ip,i,j,k,iz,ii,ib
      real(kind=8):: zz,u0,u1,v0,v1,w0,w1,bfac,xsign,ysign,zsign,txp,typ,temp
      real(kind=8):: tbx,tby

c  Evaluation of B, vectorized over particles

      if (.not. bgrds .or. .not. linbgrd(0)) return

      do io=1,nbgrdol
        if (.not. linbgrd(io)) cycle

        if (lslice) then
c         --- All particles are in the same z-cell
          zz = zp(1)
          iz = 0
          ii = cbgrdid(iz,io)
          ib = bgrdid(ii)
          k = (zz - cbgrdzs(iz,io))*bgrddzi(ib)
c         --- skip field accumulation if particles are outside of axial grid
          if (zz < cbgrdzs(iz,io) .or. cbgrdze(iz,io) < zz) cycle
        endif

        do ip = 1, np

          if (.not. lslice) then
            zz = zp(ip)
c           --- find the location of the particle in the internal lattice arrays
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
            ii = cbgrdid(iz,io)
            ib = bgrdid(ii)
            if (zz < cbgrdzs(iz,io) .or. cbgrdze(iz,io) < zz) cycle
          endif

c         --- find transverse particle coordinate in frame of gridded field 
c         --- transverse offsets
          txp = xp(ip) - bgrdox(ii)
          typ = yp(ip) - bgrdoy(ii)
c         --- transverse rotation to take into account an active rotation 
c         --- of the field element.  Particles are rotated in that oposite 
c         --- sense of the element. Later the field components accumulated 
c         --- must be rotated back. 
          if ( bgrdph(ii) .ne. 0. ) then
            temp = txp  
            txp =  temp*bgrdcp(ii) + typ*bgrdsp(ii) 
            typ = -temp*bgrdsp(ii) + typ*bgrdcp(ii) 
          endif 

c         --- Shift coordinates to measure from the edge of the field grid
          txp = txp - bgrdxs(ii)
          typ = typ - bgrdys(ii)

c         --- Set default sign of B field
          xsign = 1.
          ysign = 1.
          zsign = 1.

c         --- If B is quadrupolar symmetric, make transformations.
c         --- When the particle is in one of the even quadrants (either
c         --- x<0 or y<0 but not both), the transformation is done by
c         --- swapping x and y, and by swapping Bx and By (done at the end
c         --- of the loop). The sign of Bz is also changed. In the third
c         --- quadrant, the signs of both Bx and By are changed.
c                 Quadrupole symmetries on field grid
c           
c             Quadrant      B_x            B_y             B_z 
c             --------------------------------------------------------
c                I          B_x( x, y,z)   B_y( x, y,z)    B_z( x, y,z) 
c                II         B_y( y,-x,z)  -B_x( y,-x,z)   -B_z( y,-x,z)  
c                III       -B_x(-x,-y,z)  -B_y(-x,-y,z)    B_z(-x,-y,z) 
c                IV        -B_y(-y, x,z)   B_x(-y, x,z)   -B_z(-y, x,z) 

          if (bgrdsy(ii) == 2) then
c           --- Get quadrant that the particle is in.
            if (txp < 0.) then
              xsign = -1.
              txp = -txp
            endif
            if (typ < 0.) then
              ysign = -1.
              typ = -typ
            endif
c           --- If in even quadrant...
            if (xsign*ysign < 0.) then
c             --- Switch sign of Bz.
              zsign = -1.
c             --- Swap x and y
              temp = txp
              txp = typ
              typ = temp
            endif
          endif

c         --- find location of particle in B field grid
          i =  txp*bgrddxi(ib)
          j =  typ*bgrddyi(ib)
          k = (zz - cbgrdzs(iz,io))*bgrddzi(ib)

c         --- Calculate linear weights
          u1 = txp*bgrddxi(ib) - i
          v1 = typ*bgrddyi(ib) - j
          w1 = (zz - cbgrdzs(iz,io))*bgrddzi(ib) - k
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Only calculate for particles inside the B field grid
          if (0 <= i .and. i < bgrdnx .and.
     &        0 <= j .and. j < bgrdny .and.
     &        0 <= k .and. k < bgrdnz) then

            bfac = bgrdsc(ii) + bgrdsf(ii)

            tbx = xsign*bfac*(u0*v0*w0*bgrdbx(i  ,j  ,k  ,ib) +
     &                        u1*v0*w0*bgrdbx(i+1,j  ,k  ,ib) +
     &                        u0*v1*w0*bgrdbx(i  ,j+1,k  ,ib) +
     &                        u1*v1*w0*bgrdbx(i+1,j+1,k  ,ib) +
     &                        u0*v0*w1*bgrdbx(i  ,j  ,k+1,ib) +
     &                        u1*v0*w1*bgrdbx(i+1,j  ,k+1,ib) +
     &                        u0*v1*w1*bgrdbx(i  ,j+1,k+1,ib) +
     &                        u1*v1*w1*bgrdbx(i+1,j+1,k+1,ib))

            tby = ysign*bfac*(u0*v0*w0*bgrdby(i  ,j  ,k  ,ib) +
     &                        u1*v0*w0*bgrdby(i+1,j  ,k  ,ib) +
     &                        u0*v1*w0*bgrdby(i  ,j+1,k  ,ib) +
     &                        u1*v1*w0*bgrdby(i+1,j+1,k  ,ib) +
     &                        u0*v0*w1*bgrdby(i  ,j  ,k+1,ib) +
     &                        u1*v0*w1*bgrdby(i+1,j  ,k+1,ib) +
     &                        u0*v1*w1*bgrdby(i  ,j+1,k+1,ib) +
     &                        u1*v1*w1*bgrdby(i+1,j+1,k+1,ib))

            if (bgrdsy(ii) == 2 .and. xsign*ysign < 0.) then
              temp = tbx 
              tbx = tby
              tby = temp
            endif 

c           --- rotate transverse field componets for correct lab frame 
c               orientation if the element is rotated.  
            if ( bgrdph(ii) .ne. 0. ) then
              temp = tbx 
              tbx = temp*bgrdcp(ii) - tby*bgrdsp(ii) 
              tby = temp*bgrdsp(ii) + tby*bgrdcp(ii)         
            endif 

c           --- tranverse field 
            bx(ip) = bx(ip) + tbx
            by(ip) = by(ip) + tby

            bz(ip) = bz(ip) + zsign*bfac*(u0*v0*w0*bgrdbz(i  ,j  ,k  ,ib) +
     &                                    u1*v0*w0*bgrdbz(i+1,j  ,k  ,ib) +
     &                                    u0*v1*w0*bgrdbz(i  ,j+1,k  ,ib) +
     &                                    u1*v1*w0*bgrdbz(i+1,j+1,k  ,ib) +
     &                                    u0*v0*w1*bgrdbz(i  ,j  ,k+1,ib) +
     &                                    u1*v0*w1*bgrdbz(i+1,j  ,k+1,ib) +
     &                                    u0*v1*w1*bgrdbz(i  ,j+1,k+1,ib) +
     &                                    u1*v1*w1*bgrdbz(i+1,j+1,k+1,ib))

          endif

        enddo
      enddo

      return
      end
c=========================================================================== 
      subroutine applybgrd2(np,xp,yp,npz,zp,lslice,b)
      use Lattice
      use LatticeInternal
      use BGRDdata
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      logical(ISZ):: lslice
      real(kind=8):: b(bgrdnc,np)

c  Sets magnetic field for particles from data sets
c  containing Bx, By, and Bz on 3-D grids. Grid can also contain
c  additional vector information. It is assumed that the additional
c  information has the same symmetry as the magnetic field.
c  Calling arguments:
c     np          number of particles
c     npz         number of z position data values (must be either 1 or == np)
c     xp,yp,zp    coordinates of particles
c  Output:
c     b           magnetic field and additional information

c  The field is:
c     b(:,ip) = u0 * v0 * w0 * bgrd(:,i  ,j  ,k  ,bgrdid)
c             + u1 * v0 * w0 * bgrd(:,i+1,j  ,k  ,bgrdid)
c             + u0 * v1 * w0 * bgrd(:,i  ,j+1,k  ,bgrdid)
c             + ...

      integer(ISZ):: io,ip,i,j,k,iz,ii,ib,ic
      real(kind=8):: zz,u0,u1,v0,v1,w0,w1,bfac,xsign,ysign,zsign,txp,typ,temp
      real(kind=8):: tb(bgrdnc)

c  Evaluation of B, vectorized over particles

      if (.not. bgrds .or. .not. linbgrd(0)) return

      do io=1,nbgrdol
        if (.not. linbgrd(io)) cycle

        if (lslice) then
c         --- All particles are in the same z-cell
          zz = zp(1)
          iz = 0
          ii = cbgrdid(iz,io)
          ib = bgrdid(ii)
          k = (zz - cbgrdzs(iz,io))*bgrddzi(ib)
c         --- skip field accumulation if particles are outside of axial grid
          if (zz < cbgrdzs(iz,io) .or. cbgrdze(iz,io) < zz) cycle
        endif

        do ip = 1, np

          if (.not. lslice) then
            zz = zp(ip)
c           --- find the location of the particle in the internal lattice arrays
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
            ii = cbgrdid(iz,io)
            ib = bgrdid(ii)
            if (zz < cbgrdzs(iz,io) .or. cbgrdze(iz,io) < zz) cycle
          endif

c         --- find transverse particle coordinate in frame of gridded field 
c         --- transverse offsets
          txp = xp(ip) - bgrdox(ii)
          typ = yp(ip) - bgrdoy(ii)
c         --- transverse rotation to take into account an active rotation 
c         --- of the field element.  Particles are rotated in that oposite 
c         --- sense of the element. Later the field components accumulated 
c         --- must be rotated back. 
          if ( bgrdph(ii) .ne. 0. ) then
            temp = txp  
            txp =  temp*bgrdcp(ii) + typ*bgrdsp(ii) 
            typ = -temp*bgrdsp(ii) + typ*bgrdcp(ii) 
          endif 

c         --- Shift coordinates to measure from the edge of the field grid
          txp = txp - bgrdxs(ii)
          typ = typ - bgrdys(ii)

c         --- Set default sign of B field
          xsign = 1.
          ysign = 1.
          zsign = 1.

c         --- If B is quadrupolar symmetric, make transformations.
c         --- When the particle is in one of the even quadrants (either
c         --- x<0 or y<0 but not both), the transformation is done by
c         --- swapping x and y, and by swapping Bx and By (done at the end
c         --- of the loop). The sign of Bz is also changed. In the third
c         --- quadrant, the signs of both Bx and By are changed.
c                 Quadrupole symmetries on field grid
c           
c             Quadrant      B_x            B_y             B_z 
c             --------------------------------------------------------
c                I          B_x( x, y,z)   B_y( x, y,z)    B_z( x, y,z) 
c                II         B_y( y,-x,z)  -B_x( y,-x,z)   -B_z( y,-x,z)  
c                III       -B_x(-x,-y,z)  -B_y(-x,-y,z)    B_z(-x,-y,z) 
c                IV        -B_y(-y, x,z)   B_x(-y, x,z)   -B_z(-y, x,z) 

          if (bgrdsy(ii) == 2) then
c           --- Get quadrant that the particle is in.
            if (txp < 0.) then
              xsign = -1.
              txp = -txp
            endif
            if (typ < 0.) then
              ysign = -1.
              typ = -typ
            endif
c           --- If in even quadrant...
            if (xsign*ysign < 0.) then
c             --- Switch sign of Bz.
              zsign = -1.
c             --- Swap x and y
              temp = txp
              txp = typ
              typ = temp
            endif
          endif

c         --- find location of particle in B field grid
          i =  txp*bgrddxi(ib)
          j =  typ*bgrddyi(ib)
          k = (zz - cbgrdzs(iz,io))*bgrddzi(ib)

c         --- Calculate linear weights
          u1 = txp*bgrddxi(ib) - i
          v1 = typ*bgrddyi(ib) - j
          w1 = (zz - cbgrdzs(iz,io))*bgrddzi(ib) - k
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Only calculate for particles inside the B field grid
          if (0 <= i .and. i < bgrdnx .and.
     &        0 <= j .and. j < bgrdny .and.
     &        0 <= k .and. k < bgrdnz) then

            bfac = bgrdsc(ii) + bgrdsf(ii)

            tb = bfac*(u0*v0*w0*bgrd(:,i  ,j  ,k  ,ib) +
     &                 u1*v0*w0*bgrd(:,i+1,j  ,k  ,ib) +
     &                 u0*v1*w0*bgrd(:,i  ,j+1,k  ,ib) +
     &                 u1*v1*w0*bgrd(:,i+1,j+1,k  ,ib) +
     &                 u0*v0*w1*bgrd(:,i  ,j  ,k+1,ib) +
     &                 u1*v0*w1*bgrd(:,i+1,j  ,k+1,ib) +
     &                 u0*v1*w1*bgrd(:,i  ,j+1,k+1,ib) +
     &                 u1*v1*w1*bgrd(:,i+1,j+1,k+1,ib))


            if (bgrdsy(ii) == 2 .and. xsign*ysign < 0.) then
              do ic=0,bgrdnc/3-1
                tb(1+ic*3) = xsign*tb(1+ic*3)
                tb(2+ic*3) = ysign*tb(2+ic*3)
                tb(3+ic*3) = zsign*tb(3+ic*3)
              enddo
              if (xsign*ysign < 0.) then
                do ic=0,bgrdnc/3-1
                  temp = tb(1+ic*3)
                  tb(1+ic*3) = tb(2+ic*3)
                  tb(2+ic*3) = temp
                enddo
              endif 
            endif 

c           --- rotate transverse field componets for correct lab frame 
c               orientation if the element is rotated.  
            if ( bgrdph(ii) .ne. 0. ) then
              do ic=0,bgrdnc/3-1
                temp = tb(1+ic*3)
                tb(1+ic*3) = temp*bgrdcp(ii) - tb(2+ic*3)*bgrdsp(ii) 
                tb(2+ic*3) = temp*bgrdsp(ii) + tb(2+ic*3)*bgrdcp(ii)         
              enddo
            endif 

c           --- Accumulate the interpolated field.
            b(:,ip) = b(:,ip) + tb

          endif

        enddo
      enddo

      return
      end
c=========================================================================== 
      subroutine applypgrd(np,xp,yp,npz,zp,lslice,ex,ey,ez)
      use Lattice
      use LatticeInternal
      use PGRDdata
      integer(ISZ):: np,npz
      real(kind=8):: xp(np),yp(np),zp(npz)
      logical(ISZ):: lslice
      real(kind=8):: ex(np),ey(np),ez(np)

c  Sets electric field for particles from data
c  sets containing the potential on a 3-D grid.
c  Calling arguments:
c     np          number of particles
c     npz         number of z position data values (must be either 1 or == np)
c     xp,yp,zp    coordinates of particles
c  Output:
c     ex,ey,ez    electric field

c  The field is:
c     Ex = u0 * v0 * w0 * ex(i  ,j  ,k  ,pgrdid)
c        + u1 * v0 * w0 * ex(i+1,j  ,k  ,pgrdid)
c        + u0 * v1 * w0 * ex(i  ,j+1,k  ,pgrdid)
c        + ...
c
c  Note that this routine is very similar to sete3d routine.  It might be
c  a good idea to combine them to reduce code complexity, but that would make
c  sete3d more complex and possibly slower.  Since the pgrd routine will not
c  be used very often and sete3d is used always, the author felt it was better
c  to have a seperate routine for pgrd.

      integer(ISZ):: io,ip,i,j,k,iz,ii,ib,im1,jm1
      real(kind=8):: zz,u0,u1,v0,v1,w0,w1,xfac,yfac,zfac,txp,typ,temp

c  Evaluation of E, vectorized over particles

      if (.not. pgrds .or. .not. linpgrd(0)) return

      do io=1,npgrdol
        if (.not. linpgrd(io)) cycle

        if (lslice) then
c         --- find the location of the particle in the internal lattice arrays
          zz = zp(1)
          iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
          ii = cpgrdid(iz,io)
          ib = pgrdid(ii)
          k = (zz - cpgrdzs(iz,io))*pgrddzi(ib)
          if (zz < cpgrdzs(iz,io) .or. cpgrdze(iz,io) < zz) cycle
        endif

        do ip = 1, np

          if (.not. lslice) then
c           --- find the location of the particle in the internal lattice arrays
            zz = zp(ip)
            iz = max(0., (zz - zlmin - zlframe)*dzli + 0.5)
            ii = cpgrdid(iz,io)
            ib = pgrdid(ii)
c           --- cycle if lost particle
            if (zz < cpgrdzs(iz,io) .or. cpgrdze(iz,io) < zz) cycle
          endif

c         --- find transverse particle coordinate in frame of gridded field 
c         --- transverse offsets
          txp = xp(ip) - pgrdox(ii)
          typ = yp(ip) - pgrdoy(ii)

c         --- transverse rotation to take into account an active rotation 
c         --- of the field element.  Particles are rotated in that oposite 
c         --- sense of the element. Later the field components accumulated 
c         --- must be rotated back. 
          if ( pgrdph(ii) .ne. 0. ) then
            temp = txp  
            txp =  temp*pgrdcp(ii) + typ*pgrdsp(ii) 
            typ = -temp*pgrdsp(ii) + typ*pgrdcp(ii) 
          endif 

c         --- find location of particle in potential grid
          i =  abs(txp - pgrdxs(ii))*pgrddxi(ib)
          j =  abs(typ - pgrdys(ii))*pgrddyi(ib)
          k = (zz - cpgrdzs(iz,io))*pgrddzi(ib)

c         --- Calculate linear weights
          u1 =  abs(txp - pgrdxs(ii))*pgrddxi(ib) - i
          v1 =  abs(typ - pgrdys(ii))*pgrddyi(ib) - j
          w1 = (zz - cpgrdzs(iz,io))*pgrddzi(ib) - k
          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Only calculate for particles inside the B field grid
          if (0 <= i .and. i < pgrdnx .and.
     &        0 <= j .and. j < pgrdny .and.
     &        0 <= k .and. k < pgrdnz) then

c           --- make sure all indices refer to first quadrant
            im1 = i - 1
            jm1 = j - 1
            if (i == 0) im1 = 1
            if (j == 0) jm1 = 1

c           --- adjust sign of E field for approiate quadrant
            xfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddxi(ib)*0.5
            yfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddyi(ib)*0.5
            zfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddzi(ib)*0.5
            if (xp(ip) < pgrdxs(ii)) xfac = -xfac
            if (yp(ip) < pgrdys(ii)) yfac = -yfac

            ex(ip) = ex(ip) + xfac*
     &            (u0*v0*w0*(pgrd(im1,j  ,k  ,ib) - pgrd(i+1,j  ,k  ,ib))
     &           + u1*v0*w0*(pgrd(i  ,j  ,k  ,ib) - pgrd(i+2,j  ,k  ,ib))
     &           + u0*v1*w0*(pgrd(im1,j+1,k  ,ib) - pgrd(i+1,j+1,k  ,ib))
     &           + u1*v1*w0*(pgrd(i  ,j+1,k  ,ib) - pgrd(i+2,j+1,k  ,ib))
     &           + u0*v0*w1*(pgrd(im1,j  ,k+1,ib) - pgrd(i+1,j  ,k+1,ib))
     &           + u1*v0*w1*(pgrd(i  ,j  ,k+1,ib) - pgrd(i+2,j  ,k+1,ib))
     &           + u0*v1*w1*(pgrd(im1,j+1,k+1,ib) - pgrd(i+1,j+1,k+1,ib))
     &           + u1*v1*w1*(pgrd(i  ,j+1,k+1,ib) - pgrd(i+2,j+1,k+1,ib)))

            ey(ip) = ey(ip) + yfac*
     &            (u0*v0*w0*(pgrd(i  ,jm1,k  ,ib) - pgrd(i  ,j+1,k  ,ib))
     &           + u1*v0*w0*(pgrd(i+1,jm1,k  ,ib) - pgrd(i+1,j+1,k  ,ib))
     &           + u0*v1*w0*(pgrd(i  ,j  ,k  ,ib) - pgrd(i  ,j+2,k  ,ib))
     &           + u1*v1*w0*(pgrd(i+1,j  ,k  ,ib) - pgrd(i+1,j+2,k  ,ib))
     &           + u0*v0*w1*(pgrd(i  ,jm1,k+1,ib) - pgrd(i  ,j+1,k+1,ib))
     &           + u1*v0*w1*(pgrd(i+1,jm1,k+1,ib) - pgrd(i+1,j+1,k+1,ib))
     &           + u0*v1*w1*(pgrd(i  ,j  ,k+1,ib) - pgrd(i  ,j+2,k+1,ib))
     &           + u1*v1*w1*(pgrd(i+1,j  ,k+1,ib) - pgrd(i+1,j+2,k+1,ib)))

            ez(ip) = ez(ip) + zfac*
     &            (u0*v0*w0*(pgrd(i  ,j  ,k-1,ib) - pgrd(i  ,j  ,k+1,ib))
     &           + u1*v0*w0*(pgrd(i+1,j  ,k-1,ib) - pgrd(i+1,j  ,k+1,ib))
     &           + u0*v1*w0*(pgrd(i  ,j+1,k-1,ib) - pgrd(i  ,j+1,k+1,ib))
     &           + u1*v1*w0*(pgrd(i+1,j+1,k-1,ib) - pgrd(i+1,j+1,k+1,ib))
     &           + u0*v0*w1*(pgrd(i  ,j  ,k  ,ib) - pgrd(i  ,j  ,k+2,ib))
     &           + u1*v0*w1*(pgrd(i+1,j  ,k  ,ib) - pgrd(i+1,j  ,k+2,ib))
     &           + u0*v1*w1*(pgrd(i  ,j+1,k  ,ib) - pgrd(i  ,j+1,k+2,ib))
     &           + u1*v1*w1*(pgrd(i+1,j+1,k  ,ib) - pgrd(i+1,j+1,k+2,ib)))

          endif

        enddo
      enddo

      return
      end
c=============================================================================
      subroutine applyaccl(np,xp,zp,uzp,gaminv,dtl,dtr,dt,qoverm,lslice,ez)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np
      real(kind=8):: xp(np),zp(np),uzp(np),gaminv(np)
      real(kind=8):: dtl,dtr,dt,qoverm
      logical(ISZ):: lslice
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

      real(kind=8):: dti,vz,gapez,frac,z1,z2,vn,zl,zr,cacclz,caccll
      real(kind=8):: frac1,frac2
      integer(ISZ):: io,ip,iz

      if (.not. accls) return

      dti = 1./(dtr-dtl)

c     --- Finite length gaps
      if (.not. lacclzl) then

        do io=1,nacclol

c         --- For velocity correction, first calculate the velocity at
c         --- time level n.
          do ip=1,np
            vz = uzp(ip)*gaminv(ip)
            if (vz == 0.) vz = LARGEPOS
c           --- find z-cell in which particle lies
            iz = int(max(0., (zp(ip) - zlmin - zlframe) * dzli + 0.5))
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
            z1 = zp(ip) + dtl*vz - 0.5*gapez*qoverm*dtl**2
            z2 = zp(ip) + dtl*vz
            if (zp(ip) <= cacclzs(iz,io)) then
              vn = vz
            elseif (z1 < cacclzs(iz,io)) then
              vn = sqrt(vz**2 + 2.*gapez*qoverm*(zp(ip) - cacclzs(iz,io)))
            elseif (zp(ip) <= cacclze(iz,io)) then
              vn = vz - gapez*qoverm*dtl
            elseif (z2 < cacclze(iz,io)) then
              vn = 0.5*(vz + 0.5*gapez*qoverm*dt + sqrt((vz +
     &             0.5*gapez*qoverm*dt)**2 -
     &             4.*gapez*qoverm*(zp(ip) - cacclze(iz,io))))
            else
              vn = vz
            endif
c           --- Calculate the fraction of time in the gap.  Cases inside
c           --- and outside are included implicitly in the max and min calls.
c           --- Note that the max's inside the sqrt are for idiot proofing.
            if (zp(ip) <= cacclzs(iz,io)) then
              frac = max((dtr + (zp(ip) - cacclzs(iz,io))/vn)*dti , 0.)
            elseif (zp(ip) <= cacclze(iz,io)) then
              frac1 = min((dtr + 2.*(zp(ip) - cacclzs(iz,io))/
     &                    (vz + vn))*dti, 1.)
              frac2 = min((-dtl + 2.*(cacclze(iz,io) - zp(ip))/
     &                (sqrt(max(0.,vn**2 + 2.*gapez*qoverm*
     &                (cacclze(iz,io) - zp(ip)))) + vn))*dti, 1.)
              frac = min(frac1,frac2)
            else
              frac = max((-dtl - (zp(ip) - cacclze(iz,io))/vn)*dti , 0.)
            endif
c           --- add acceleration field to Ez field
            ez(ip) = ez(ip) + gapez*frac
          enddo
        enddo

      else

        do io=1,nacclol

c         --- Zero length gaps
          do ip = 1, np
            vz = uzp(ip)*gaminv(ip)
            iz = max(0., (zp(ip) - zlmin - zlframe) * dzli + 0.5)
c           --- "left" end of velocity advance step
            zl = zp(ip) + vz*dtl
c           --- "right" end of velocity advance step
            zr = zp(ip) + vz*dtr
c           --- Gap center and length
            cacclz = 0.5*(cacclzs(iz,io) + cacclze(iz,io))
            caccll = cacclze(iz,io) - cacclzs(iz,io)
c           --- Add acceleration to velocity.
c           --- Calculate the change in velocity and convert that into
c           --- and Ez - the expression reduces to q*V/m for small V.
            if (zl <= cacclz .and. cacclz < zr) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              ez(ip) = ez(ip) + (sqrt(vz**2 + 2.*qoverm*gapez*caccll)-vz)/
     &                          qoverm*dti
            endif
          enddo
        enddo

      endif

      return
      end
c=========================================================================== 
      subroutine zgapcorr(np,zp,xp,uzp,gaminv,dtl,dtr,dt,m,q,time)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np
      real(kind=8):: zp(np),xp(np),uzp(np),gaminv(np)
      real(kind=8):: dtl,dtr,dt,m,q,time

      integer(ISZ):: io,ip,iz
      real(kind=8):: gapez,zp2,delta,qoverm,cacclz,caccll,deltav,zr
      real(kind=8):: vz(np),vzi(np)

      qoverm = q/m

c     Add residence correction to particle position when acceleration is done

      if (accls) then

        vz = uzp*gaminv
        vzi = 1./dvnz(vz)

c       --- Finite length gap
        if (.not. lacclzl) then

          do io=1,nacclol

            do ip=1,np
c             --- find index of nearest accl element
              iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
c             --- Approximate position of particle on next step
c             --- with simple residence correction for gaps.
c             --- Calculate time of entrance or exit, delta.
c             --- Space charge forces are ignored.
              zp2 = zp(ip) + vz(ip)*dt
              delta = 0
              if (zp(ip) <= cacclzs(iz,io) .and. cacclzs(iz,io) <= zp2) then
                delta = (cacclzs(iz,io) - zp(ip))*vzi(ip)
              elseif (cacclzs(iz,io) < zp(ip) .and. zp(ip) <= cacclze(iz,io)
     &            .and. cacclze(iz,io) <= zp2 + gapez*qoverm*dt*(dtr-dtl)) then
                delta = 2.*(cacclze(iz,io) - zp(ip))/
     &                  (sqrt((vz(ip) - gapez*qoverm*dtl)**2 +
     &                  2.*(cacclze(iz,io) - zp(ip))*gapez*qoverm) +
     &                  vz(ip) - gapez*qoverm*dtl)
                zp2 = zp2 + gapez*qoverm*dt*(dtr-dtl)
              endif
c             --- calculate correction on particle position
              if (delta > 0.5*dt) delta = dt - delta
              if (zp(ip) <= cacclzs(iz,io) .and. cacclzs(iz,io) <= zp2) then
                zp(ip) = zp(ip) + 0.5*gapez*qoverm*delta**2
              elseif (zp(ip) <= cacclze(iz,io) .and. cacclze(iz,io) <= zp2) then
                zp(ip) = zp(ip) - 0.5*gapez*qoverm*delta**2
              endif
            enddo
          enddo

        else
c         --- Zero length gap

          do io=1,nacclol

            do ip = 1, np
              iz = max(0., (zp(ip) - zlmin - zlframe)*dzli + 0.5)
              cacclz = 0.5*(cacclzs(iz,io) + cacclze(iz,io))
              caccll = cacclze(iz,io) - cacclzs(iz,io)
c             --- "right" end of velocity advance step
              zr = zp(ip) + vz(ip)*dtr
c             --- calculate correction on position
              if (zp(ip) > cacclz .and. cacclz > zr) then
                gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
                deltav = sqrt(vz(ip)**2 + 2.*qoverm*gapez*caccll) - vz(ip)
                zp(ip) = zp(ip) - (cacclz - zp(ip))*vzi(ip)*deltav
              endif
              if (zr > cacclz .and. cacclz > zp(ip) + vz(ip)*dt) then
                gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
                deltav = sqrt(vz(ip)**2 + 2.*qoverm*gapez*caccll) - vz(ip)
                zp(ip) = zp(ip) + (dt - (cacclz - zp(ip))*vzi(ip))*deltav
              endif
            enddo
          enddo

        endif
      endif

      return
      end
c=============================================================================
      subroutine applyacclxy(np,xp,zp,uzp,gaminv,dtp,dtl,dtr,m,q,dt,lslice,ez)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np
      real(kind=8):: dtl,dtr,m,q,dt
      real(kind=8):: xp(np),zp(np)
      real(kind=8):: uzp(np),gaminv(np),dtp(np)
      logical(ISZ):: lslice
      real(kind=8):: ez(np)

c Apply acceleration for finite length gaps.
c This is done seperately from the rest of the elements (which are done
c in extebxy) since the residence corrections depends on dtp which
c changes during the iteration dealing with changing Vz.  It would be much
c less efficient to put the entire extebxy within that iteration loop so
c the accl elements was seperated out.

      real(kind=8):: dti,oneodt,qoverm,moverq,gapez,z1,z2,vn,frac
      real(kind=8):: cacclz,caccll,zl,zr
      real(kind=8):: dtl_a,dtr_a,dti_a,vz_a
      integer(ISZ):: io,iz,ip

c  --- Note that the accl elements in the slice code are slightly different
c  --- than other hard edged elements since the residence correction is
c  --- based on fraction of time spent inside element rather than fraction of
c  --- distance.  So, the particles individual velocity and dt are used.

      if (.not. accls) return

      dti = 1./(dtr-dtl)
      oneodt = 1./dt
      qoverm = q/m

c     --- All particles lie at the same z-cell
      iz = 0

      if (.not. lacclzl) then

        do io=1,nacclol

c         --- For velocity correction, first calculate the velocity at
c         --- time level n.
          do ip=1,np
            dtl_a = dtl*oneodt*dtp(ip)
            dtr_a = dtr*oneodt*dtp(ip)
            dti_a = dti*dt/dtp(ip)
            vz_a = uzp(ip)*gaminv(ip)
            if (vz_a == 0.) vz_a = LARGEPOS
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
            z1 = zlframe + dtl_a*vz_a - 0.5*gapez*qoverm*dtl_a**2
            z2 = zlframe + dtl_a*vz_a
            if (zlframe <= cacclzs(iz,io)) then
              vn = vz_a
            elseif (z1 < cacclzs(iz,io)) then
              vn = sqrt(vz_a**2 + 2.*gapez*qoverm*(zlframe - cacclzs(iz,io)))
            elseif (zlframe <= cacclze(iz,io)) then
              vn = vz_a - gapez*qoverm*dtl_a
            elseif (z2 < cacclze(iz,io)) then
              vn = 0.5*(vz_a + 0.5*gapez*qoverm*dtp(ip) + sqrt((vz_a +
     &             0.5*gapez*qoverm*dtp(ip))**2 -
     &             4.*gapez*qoverm*(zlframe - cacclze(iz,io))))
            else
              vn = vz_a
            endif
c           --- Calculate the fraction of time in the gap.  Cases inside
c           --- and outside are included implicitly in the max and min calls.
c           --- Note that the max's inside the sqrt are for idiot proofing.
            if (zlframe <= cacclzs(iz,io)) then
              frac = max((dtr_a + (zlframe - cacclzs(iz,io))/vn)*dti_a , 0.)
            elseif (zlframe <= (cacclzs(iz,io) + cacclze(iz,io))*0.5) then
              frac = min((dtr_a + 2.*(zlframe - cacclzs(iz,io))/
     &               (sqrt(max(0.,vn**2 - 2.*gapez*qoverm*
     &               (zlframe - cacclzs(iz,io)))) + vn))*dti_a, 1.)
            elseif (zlframe <= cacclze(iz,io)) then
              frac = min((-dtl_a + 2.*(cacclze(iz,io) - zlframe)/
     &               (sqrt(max(0.,vn**2 + 2.*gapez*qoverm*
     &               (cacclze(iz,io) - zlframe))) + vn))*dti_a, 1.)
            else
              frac = max((-dtl_a - (zlframe - cacclze(iz,io))/vn)*dti_a , 0.)
            endif
c           --- add acceleration field to Ez field
            ez(ip) = ez(ip) + gapez*frac
          enddo
        enddo

      else

c       --- Zero-length gaps.
        do io=1,nacclol

          moverq = m/q
          cacclz = 0.5*(cacclzs(iz,io) + cacclze(iz,io))
          caccll = cacclze(iz,io) - cacclzs(iz,io)
c         --- Add acceleration to velocity.
c         --- Calculate the change in velocity and convert that into
c         --- and Ez - the expression reduces to q*V/m for small V.
          do ip = 1, np
            dtl_a = dtl*oneodt*dtp(ip)
            dtr_a = dtr*oneodt*dtp(ip)
            dti_a = dti*dt/dtp(ip)
            vz_a = uzp(ip)*gaminv(ip)
            zl = zlframe + dtl_a*vz_a
            zr = zlframe + dtr_a*vz_a
            if (zl <= cacclz .and. cacclz < zr) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              ez(ip) = ez(ip) + (sqrt(vz_a**2 + 2.*qoverm*gapez*caccll) - vz_a)
     &                 *moverq*dti_a
            endif
          enddo
        enddo

      endif

c  ------------------------------------------------------------
c  End of acceleration
c  ------------------------------------------------------------

      return
      end
c=========================================================================== 
      subroutine zgapcorrxy(np,zp,xp,uzp,gaminv,dtp,dtl,dtr,dt,m,q,time)
      use Lattice
      use LatticeInternal
      integer(ISZ):: np
      real(kind=8):: zp(np), xp(np), uzp(np), gaminv(np), dtp(np)
      real(kind=8):: dtl,dtr,dt,m,q,time

c Add residence correction to particle position when acceleration is done

      integer(ISZ):: io,ip,iz
      real(kind=8):: gapez,zp2,delta,qoverm,vz(np),vzi(np)
      real(kind=8):: dtl_a,dtr_a
      real(kind=8):: cacclz,caccll,zr,deltav

      if (.not. accls) return

      qoverm = q/m
      vz = uzp*gaminv
      vzi = 1./dvnz(vz)

c     --- All particles lie at the same z-cell
      iz = 0

c     --- Finite length gap
      if (.not. lacclzl) then

        do io=1,nacclol

          do ip=1,np
c           --- Scale time step size by particles time step.
            dtl_a = dtl/dt*dtp(ip)
            dtr_a = dtr/dt*dtp(ip)
c           --- Fetch accelerating voltage
            gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
c           --- Approximate position of particle on next step
c           --- with simple residence correction for gaps.
c           --- Calculate time of entrance or exit, delta.
c           --- Space charge forces are ignored.
            zp2 = zlframe + vz(ip)*dtp(ip)
            delta = 0
            if (zlframe <= cacclzs(iz,io) .and. cacclzs(iz,io) <= zp2) then
              delta = (cacclzs(iz,io) - zlframe)*vzi(ip)
            elseif (cacclzs(iz,io) < zlframe .and. zlframe <= cacclze(iz,io)
     &   .and. cacclze(iz,io) <= zp2 + gapez*qoverm*dtp(ip)*(dtr_a-dtl_a)) then
              delta = 2.*(cacclze(iz,io) - zlframe)/
     &                (sqrt((vz(ip) - gapez*qoverm*dtl_a)**2 +
     &                2.*(cacclze(iz,io) - zlframe)*gapez*qoverm) +
     &                vz(ip) - gapez*qoverm*dtl_a)
              zp2 = zp2 + gapez*qoverm*dtp(ip)*(dtr_a-dtl_a)
            endif
c           --- calculate correction on particle position
            if (delta > 0.5*dtp(ip)) delta = dtp(ip) - delta
            if (zlframe <= cacclzs(iz,io) .and. cacclzs(iz,io) <= zp2) then
              zp(ip) = zp(ip) + 0.5*gapez*qoverm*delta**2
            elseif (zlframe <= cacclze(iz,io) .and. cacclze(iz,io) <= zp2) then
              zp(ip) = zp(ip) - 0.5*gapez*qoverm*delta**2
            endif
          enddo
        enddo

      else

c       --- Zero length gap
        do io=1,nacclol

          cacclz = 0.5*(cacclzs(iz,io) + cacclze(iz,io))
          caccll = cacclze(iz,io) - cacclzs(iz,io)
          do ip = 1, np
c           --- "right" end of velocity advance step
            zr = zlframe + vz(ip)*dtr
c           --- calculate correction on position
            if (zlframe < cacclz .and. cacclz < zr) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              deltav = sqrt(vz(ip)**2 + 2.*qoverm*gapez*caccll) - vz(ip)
              zp(ip) = zp(ip) - (cacclz - zlframe)*vzi(ip)*deltav
            endif
            if (zr < cacclz .and. cacclz < zlframe + vz(ip)*dtp(ip)) then
              gapez = cacclez(iz,io) + cacclxw(iz,io)*xp(ip)
              deltav = sqrt(vz(ip)**2 + 2.*qoverm*gapez*caccll) - vz(ip)
              zp(ip) = zp(ip) + (dtp(ip) - (cacclz - zlframe)*vzi(ip))*deltav
            endif
          enddo
        enddo

      endif

      return
      end
c=============================================================================
      subroutine acclbfrm(zcorrection)
      use Beam_acc
      use InGen
      use InPart
      use Lattice
      use LatticeInternal
      use Picglb
      use Particles
      real(kind=8):: zcorrection
c Accelerate the beam frame using the accl lattice elements.
c Beam frame is treated as a particle of species 1 located at the center of
c the beam frame.
c
c Output argument, zcorrection, is the correction to the beam frame's
c longitudinal postion. Note that zbeam can not be directly changed here
c since that may force some particles off of the mesh. Instead, the
c zcorrection is returned so that it can be dealt with in a reasonable way
c by the calling routine.

      integer(ISZ):: io,astat,iz
      real(kind=8):: zs,ze,zc
      real(kind=8):: gapez,z1,z2,vn,frac,qoverm
      real(kind=8):: zz,zl,zr,vbeamfrmo,deltav
      real(kind=8):: zznew

      qoverm = sq(1)/sm(1)
      zcorrection = 0.

      if (accls .and. vbeamfrm .ne. 0.) then

        if (acclbeamframe .ne. 0.) then
          zz = acclbeamframe + zlframe
          iz = nint((acclbeamframe-zlmin)/dzl)
          iz = max(0,min(nzl,iz))
        else
          iz = nzl/2
          zz = zlmin + iz*dzl + zlframe
        endif

        do io=1,nacclol

c         --- Check if this gap should accelerate the beam frame.
          if (cacclsw(iz,io) .ne. 0) cycle

c         --- Set some temporaries
          gapez = cacclez(iz,io)
          zs = cacclzs(iz,io)
          ze = cacclze(iz,io)

c         --- If not using zero length gaps then do residence corrections
c         --- across the gap.
          if (.not. lacclzl) then

c           --- Calculate grid frame velocity at time level n
            z1 = zz - 0.5*dt*vbeamfrm - 0.5*gapez*qoverm*(dt*0.5)**2
            z2 = zz - 0.5*dt*vbeamfrm
            if (zz <= zs) then
              vn = vbeamfrm
            elseif (z1 < zs) then
              vn = sqrt(vbeamfrm**2 + 2.*gapez*qoverm*(zz - zs))
            elseif (zz <= ze) then
              vn = vbeamfrm + 0.5*gapez*qoverm*dt
            elseif (z2 < ze) then
              vn = 0.5*(vbeamfrm + 0.5*gapez*qoverm*dt + sqrt((vbeamfrm +
     &             0.5*gapez*qoverm*dt)**2 -
     &             4.*gapez*qoverm*(zz - ze)))
            else
              vn = vbeamfrm
            endif
c           --- Calculate the fraction of time in the gap.  Cases inside
c           --- and outside are included implicitly in the max and min calls.
c           --- Note that the max's inside the sqrt are for idiot proofing.
            if (zz <= zs) then
              frac = max(0.5 + (zz - zs)/(vn*dt) , 0.)
            elseif (zz <= (zs + ze)*0.5) then
              frac = min(0.5 + 2.*(zz - zs)/
     &               (dt*(sqrt(max(0.,vn**2 - 2.*gapez*qoverm*
     &               (zz - zs))) + vn)), 1.)
            elseif (zz <= ze) then
              frac = min(0.5 + 2.*(ze - zz)/
     &               (dt*(sqrt(max(0.,vn**2 + 2.*gapez*qoverm*
     &               (ze - zz))) + vn)), 1.)
            else
              frac = max(0.5 - (zz - ze)/(vn*dt) , 0.)
            endif
c           --- do correction on position
            zznew = zz
            call zgapcorr(1,zznew,0.,vbeamfrm,1.,-dt*0.5,dt*0.5,dt,
     &                    sm(1),sq(1),time)
            zcorrection = zcorrection + zznew - zz
c           --- add acceleration to beam frame velocity
            vbeamfrm = vbeamfrm + gapez*qoverm*dt*frac
c           --- also change beam frame velocity
            vbeam = vbeamfrm

          else
c           --- Zero-length gap

c           --- location of gap center
            zc = 0.5*(zs + ze)
c           --- "left" end of velocity advance step
            zl = zz - 0.5*vbeamfrm*dt
c           --- "right" end of velocity advance step
            zr = zz + 0.5*vbeamfrm*dt
c           --- Save current beam velocity
            vbeamfrmo = vbeamfrm
c           --- add acceleration to beam velocity
            if (zl .le. zc .and. zc .lt. zr) then
              vbeamfrm = sqrt(vbeamfrm**2 + 2.*qoverm*gapez*(ze - zs))
              vbeam = vbeamfrm
            endif
c           --- calculate correction on position of beam frame
            frac = 0.
            deltav = 0.
            if (zz < zc .and. zc < zr) then
              frac = - (zc - zz)/(vbeamfrmo*dt)
              deltav = vbeamfrm - vbeamfrmo
            endif
            if (zr < zc .and. zc < zz+vbeamfrmo*dt) then
              frac = 1. - (zc - zz)/(vbeamfrmo*dt)
              deltav = sqrt(vbeamfrm**2 + 2.*qoverm*gapez*(ze - zs))
     &                 - vbeamfrm
            endif
c           --- Set the correction to the beam frame position
            zcorrection = zcorrection + frac*deltav*dt

          endif
        enddo
      endif

      return
      end
c=============================================================================
      subroutine getzmmnt(np,xp,yp,zp,uxp,uyp,uzp,gaminv,q,m,w,dt,itask,nplive,
     &                    uxpo,uypo,uzpo,is,ns,maxp,minp,zmmnts0,zmmnts)
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
      integer(ISZ):: np,itask,nplive,is,ns
      real(kind=8):: q,m,w,dt
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np)
      real(kind=8):: maxp(6),minp(6)
      real(kind=8):: zmmnts0(NUMZMMNT)
      real(kind=8):: zmmnts(0:nzmmnt,NUMZMMNT)

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


      real(kind=8):: pp(0:17,np)
      integer(ISZ):: iz(np)
      integer(ISZ):: im,izm,iw,ip,icell,i
      real(kind=8):: dti,dtip1,wz1,wz0,zwin,tg,th
      real(kind=8):: oneondt,clighti,vzi,pnumi
      real(kind=8):: delxsq,delxxp,delxpsq,delysq,delyyp,delypsq
      real(kind=8):: delxy,delxyp,delyxp,delxpyp
      real(kind=8):: delzsq,delvzsq,delzvz
      real(kind=8):: delvxsq,delvysq
      real(kind=8):: gamma
      real(kind=8):: timetemp,wtime
#ifdef J90
      real(kind=8):: pmmnts(NUMZMMNT,np)
#endif

      timetemp = wtime()

      if (ifzmmnt == 0) return

      oneondt = 1./dvnz(dt)

      if (itask == 1 .and. .not. laccumulate_zmoments) then
         maxp = -LARGEPOS
         minp = +LARGEPOS
         zmmnts0 = 0.
         zmmnts = 0.
         if (nepwin > 0) nep = 0
         if(wpid>0) nplive = 0
      endif

      if (itask == 2) then
c       --- Set maximum and minimum particle coordinate variables
c       --- set multiplier so lost particles don't contribute
        do ip = 1, np
          if (uzp(ip) /= 0.) then
            maxp(1) = max(maxp(1), xp(ip))
            minp(1) = min(minp(1), xp(ip))
            maxp(2) = max(maxp(2), yp(ip))
            minp(2) = min(minp(2), yp(ip))
            maxp(3) = max(maxp(3), zp(ip))
            minp(3) = min(minp(3), zp(ip))
            maxp(4) = max(maxp(4), gaminv(ip)*uxp(ip))
            minp(4) = min(minp(4), gaminv(ip)*uxp(ip))
            maxp(5) = max(maxp(5), gaminv(ip)*uyp(ip))
            minp(5) = min(minp(5), gaminv(ip)*uyp(ip))
            maxp(6) = max(maxp(6), gaminv(ip)*uzp(ip))
            minp(6) = min(minp(6), gaminv(ip)*uzp(ip))

            pp(0,ip) = 1.
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

c-------------------------------------------------------------------------
#ifdef J90
c       --- This method is faster for vectorized machines, primarily the J90.
c       --- put particle moments into pmmnts array for total moments
        do ip=1,np
          vzi = 1./(uzp(ip)+SMALLPOS)
          if (pp(0,ip) > 0.) then
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
          if (pp(0,ip) > 0.) then
c           do im=1,NUMZMMNT
c             zmmnts0(im) = zmmnts0(im) + pmmnts(ip,im)
c           enddo
            zmmnts0( 1) = zmmnts0( 1) + pmmnts( 1,ip)
            zmmnts0( 2) = zmmnts0( 2) + pmmnts( 2,ip)
            zmmnts0( 3) = zmmnts0( 3) + pmmnts( 3,ip)
            zmmnts0( 4) = zmmnts0( 4) + pmmnts( 4,ip)
            zmmnts0( 5) = zmmnts0( 5) + pmmnts( 5,ip)
            zmmnts0( 6) = zmmnts0( 6) + pmmnts( 6,ip)
            zmmnts0( 7) = zmmnts0( 7) + pmmnts( 7,ip)
            zmmnts0( 8) = zmmnts0( 8) + pmmnts( 8,ip)
            zmmnts0( 9) = zmmnts0( 9) + pmmnts( 9,ip)
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
          endif
        enddo

#else

c-------------------------------------------------------------------------
c       --- This method is faster for RISC and pentium machines.
c       --- sum particle moments into total moments
        do ip=1,np
          vzi = 1./(uzp(ip)+SMALLPOS)
          if (pp(0,ip) > 0.) then
            zmmnts0( 1) = zmmnts0( 1) + 1.
            zmmnts0( 2) = zmmnts0( 2) + xp(ip)
            zmmnts0( 3) = zmmnts0( 3) + yp(ip)
            zmmnts0( 4) = zmmnts0( 4) + zp(ip)
            zmmnts0( 7) = zmmnts0( 7) + uxp(ip)*gaminv(ip)
            zmmnts0( 8) = zmmnts0( 8) + uyp(ip)*gaminv(ip)
            zmmnts0( 9) = zmmnts0( 9) + uzp(ip)*gaminv(ip)
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
            zmmnts0( 5) = zmmnts0( 5) + uxp(ip)*vzi
            zmmnts0( 6) = zmmnts0( 6) + uyp(ip)*vzi
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
            pp(16,ip) = dti
            pp(17,ip) = dtip1
          enddo

c         --- Add particles to ExtPart arrays (skips if nepwin = 0).
c         --- This allows access to the particle data when it is extrapolated
c         --- to the moments grid cell centers.
          do iw=1,nepwin
            if (izepwin(iw) < 0 .or. izepwin(iw) > nzmmnt) cycle
            do ip=1,np
              if (iz(ip) == izepwin(iw) .and. nep(iw,is) < nepmax .and.
     &            pp(0,ip) > 0.) then
!$OMP ATOMIC
                nep(iw,is) = nep(iw,is) + 1
                tep(nep(iw,is),iw,is) = pp(16,ip) + time
                xep(nep(iw,is),iw,is) = pp(1,ip)
                yep(nep(iw,is),iw,is) = pp(5,ip)
                uxep(nep(iw,is),iw,is) = pp(3,ip)
                uyep(nep(iw,is),iw,is) = pp(7,ip)
                uzep(nep(iw,is),iw,is) = pp(11,ip)
              else if (iz(ip)+1 == izepwin(iw) .and.
     &                 nep(iw,is) < nepmax .and.
     &                 pp(0,ip) > 0.) then
!$OMP ATOMIC
                nep(iw,is) = nep(iw,is) + 1
                tep(nep(iw,is),iw,is) = pp(17,ip) + time
                xep(nep(iw,is),iw,is) = pp(2,ip)
                yep(nep(iw,is),iw,is) = pp(6,ip)
                uxep(nep(iw,is),iw,is) = pp(4,ip)
                uyep(nep(iw,is),iw,is) = pp(8,ip)
                uzep(nep(iw,is),iw,is) = pp(12,ip)
              endif
            enddo
          enddo

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

      if (itask == 3 .and. .not. laccumulate_zmoments) then
        clighti = 1./clight
        gamma = 1./sqrt(1. - (vbeam*clighti)**2)

        xmaxp  = maxp(1)
        xminp  = minp(1)
        ymaxp  = maxp(2)
        yminp  = minp(2)
        zmaxp  = maxp(3)
        zminp  = minp(3)
        vxmaxp = maxp(4)
        vxminp = minp(4)
        vymaxp = maxp(5)
        vyminp = minp(5)
        vzmaxp = maxp(6)
        vzminp = minp(6)

c       --- For slave, call routine which sums moments over processors.
#ifdef MPIPARALLEL
        call sum_mmnts(zmmnts0,zmmnts)
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
        rrms(0) = sqrt(xrms(0)**2 + yrms(0)**2)

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
        if(wpid==0) nplive = pnum(0) + 0.5

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
        if (ifzmmnt == 2) then
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
            rrmsz(izm) = sqrt(xrmsz(izm)**2 + yrmsz(izm)**2)

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
            if (zmmntmin <= zwin .and. zwin < zmmntmax) then
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
              rrms(iw)    = rrmsz(icell)*wz1    + rrmsz(icell+1)*wz0
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

!$OMP ATOMIC
      momentstime = momentstime + (wtime() - timetemp)

      return
      end
c=============================================================================
      subroutine getzmmnt_weights(np,xp,yp,zp,uxp,uyp,uzp,gaminv,wp,q,m,w,dt,itask,nplive,
     &                    uxpo,uypo,uzpo,is,ns,maxp,minp,zmmnts0,zmmnts)
      use Constant
      use Beam_acc
      use InDiag
      use Z_Moments
      use Win_Moments
      use Moments
      use Picglb
      use ExtPart
      use Timers
      integer(ISZ):: np,itask,nplive,is,ns
      real(kind=8):: q,m,w,dt
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np), wp(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np)
      real(kind=8):: maxp(6),minp(6)
      real(kind=8):: zmmnts0(NUMZMMNT)
      real(kind=8):: zmmnts(0:nzmmnt,NUMZMMNT)

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


      real(kind=8):: pp(1:17,np)
      integer(ISZ):: iz(np)
      integer(ISZ):: im,izm,iw,ip,icell,i
      real(kind=8):: dti,dtip1,wz1,wz0,zwin,tg,th
      real(kind=8):: oneondt,clighti,vzi,pnumi
      real(kind=8):: delxsq,delxxp,delxpsq,delysq,delyyp,delypsq
      real(kind=8):: delxy,delxyp,delyxp,delxpyp
      real(kind=8):: delzsq,delvzsq,delzvz
      real(kind=8):: delvxsq,delvysq
      real(kind=8):: gamma
      real(kind=8):: timetemp,wtime
#ifdef J90
      real(kind=8):: pmmnts(NUMZMMNT,np)
#endif

      timetemp = wtime()

      if (ifzmmnt == 0) return

      oneondt = 1./dvnz(dt)

      if (itask == 1 .and. .not. laccumulate_zmoments) then
        write(0,*) 'Error in getzmmnt_weights: call getzmmnt when itask = 1'
        stop
      endif

      if (itask == 2) then
c       --- Set maximum and minimum particle coordinate variables
        do ip = 1, np
          if (uzp(ip) /= 0.) then
            maxp(1) = max(maxp(1), xp(ip))
            minp(1) = min(minp(1), xp(ip))
            maxp(2) = max(maxp(2), yp(ip))
            minp(2) = min(minp(2), yp(ip))
            maxp(3) = max(maxp(3), zp(ip))
            minp(3) = min(minp(3), zp(ip))
            maxp(4) = max(maxp(4), gaminv(ip)*uxp(ip))
            minp(4) = min(minp(4), gaminv(ip)*uxp(ip))
            maxp(5) = max(maxp(5), gaminv(ip)*uyp(ip))
            minp(5) = min(minp(5), gaminv(ip)*uyp(ip))
            maxp(6) = max(maxp(6), gaminv(ip)*uzp(ip))
            minp(6) = min(minp(6), gaminv(ip)*uzp(ip))
            nplive = nplive + 1
          endif
        enddo

c       --- Accumate global moments first.

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
c             zmmnts0(im) = zmmnts0(im) + pmmnts(ip,im)
c           enddo
            zmmnts0( 1) = zmmnts0( 1) + pmmnts( 1,ip)
            zmmnts0( 2) = zmmnts0( 2) + pmmnts( 2,ip)
            zmmnts0( 3) = zmmnts0( 3) + pmmnts( 3,ip)
            zmmnts0( 4) = zmmnts0( 4) + pmmnts( 4,ip)
            zmmnts0( 5) = zmmnts0( 5) + pmmnts( 5,ip)
            zmmnts0( 6) = zmmnts0( 6) + pmmnts( 6,ip)
            zmmnts0( 7) = zmmnts0( 7) + pmmnts( 7,ip)
            zmmnts0( 8) = zmmnts0( 8) + pmmnts( 8,ip)
            zmmnts0( 9) = zmmnts0( 9) + pmmnts( 9,ip)
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
          endif
        enddo

#else

c-------------------------------------------------------------------------
c       --- This method is faster for RISC and pentium machines.
c       --- sum particle moments into total moments
        do ip=1,np
          vzi = 1./(uzp(ip)+SMALLPOS)
          if (wp(ip) > 0.) then
            zmmnts0( 1) = zmmnts0( 1) + wp(ip)*wp(ip)
            zmmnts0( 2) = zmmnts0( 2) + xp(ip)*wp(ip)
            zmmnts0( 3) = zmmnts0( 3) + yp(ip)*wp(ip)
            zmmnts0( 4) = zmmnts0( 4) + zp(ip)*wp(ip)
            zmmnts0( 7) = zmmnts0( 7) + uxp(ip)*gaminv(ip)*wp(ip)
            zmmnts0( 8) = zmmnts0( 8) + uyp(ip)*gaminv(ip)*wp(ip)
            zmmnts0( 9) = zmmnts0( 9) + uzp(ip)*gaminv(ip)*wp(ip)
            zmmnts0(10) = zmmnts0(10) + xp(ip)*yp(ip)*wp(ip)
            zmmnts0(12) = zmmnts0(12) + xp(ip)**2*wp(ip)
            zmmnts0(13) = zmmnts0(13) + yp(ip)**2*wp(ip)
            zmmnts0(14) = zmmnts0(14) + zp(ip)**2*wp(ip)
            zmmnts0(17) = zmmnts0(17) + (gaminv(ip)*uxp(ip))**2*wp(ip)
            zmmnts0(18) = zmmnts0(18) + (gaminv(ip)*uyp(ip))**2*wp(ip)
            zmmnts0(19) = zmmnts0(19) + (gaminv(ip)*uzp(ip))**2*wp(ip)
            zmmnts0(22) = zmmnts0(22) + zp(ip)*gaminv(ip)*uzp(ip)*wp(ip)
            zmmnts0(25) = zmmnts0(25) + xp(ip)*uzp(ip)*gaminv(ip)*wp(ip)
            zmmnts0(26) = zmmnts0(26) + yp(ip)*uzp(ip)*gaminv(ip)*wp(ip)
            zmmnts0(27) = zmmnts0(27) + uxp(ip)*uzp(ip)*gaminv(ip)**2*wp(ip)
            zmmnts0(28) = zmmnts0(28) + uyp(ip)*uzp(ip)*gaminv(ip)**2*wp(ip)
            zmmnts0( 5) = zmmnts0( 5) + uxp(ip)*vzi*wp(ip)
            zmmnts0( 6) = zmmnts0( 6) + uyp(ip)*vzi*wp(ip)
            zmmnts0(11) = zmmnts0(11) + uxp(ip)*uyp(ip)*vzi**2*wp(ip)
            zmmnts0(15) = zmmnts0(15) + (uxp(ip)*vzi)**2*wp(ip)
            zmmnts0(16) = zmmnts0(16) + (uyp(ip)*vzi)**2*wp(ip)
            zmmnts0(20) = zmmnts0(20) + xp(ip)*uxp(ip)*vzi*wp(ip)
            zmmnts0(21) = zmmnts0(21) + yp(ip)*uyp(ip)*vzi*wp(ip)
            zmmnts0(23) = zmmnts0(23) + xp(ip)*uyp(ip)*vzi*wp(ip)
            zmmnts0(24) = zmmnts0(24) + yp(ip)*uxp(ip)*vzi*wp(ip)
          endif
        enddo

c-------------------------------------------------------------------------
#endif


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
            pp(14,ip) = (1. - (zp(ip)-zbeam-zmmntmin)*dzmi + iz(ip))*wp(ip)
            pp(15,ip) = ((zp(ip) - zbeam - zmmntmin)*dzmi - iz(ip))*wp(ip)
            pp(16,ip) = dti
            pp(17,ip) = dtip1
          enddo

c         --- Add particles to ExtPart arrays (skips if nepwin = 0).
c         --- This allows access to the particle data when it is extrapolated
c         --- to the moments grid cell centers.
          do iw=1,nepwin
            if (izepwin(iw) < 0 .or. izepwin(iw) > nzmmnt) cycle
            do ip=1,np
              if (iz(ip) == izepwin(iw) .and. nep(iw,is) < nepmax .and.
     &            wp(ip) > 0.) then
!$OMP ATOMIC
                nep(iw,is) = nep(iw,is) + 1
                tep(nep(iw,is),iw,is) = pp(16,ip) + time
                xep(nep(iw,is),iw,is) = pp(1,ip)
                yep(nep(iw,is),iw,is) = pp(5,ip)
                uxep(nep(iw,is),iw,is) = pp(3,ip)
                uyep(nep(iw,is),iw,is) = pp(7,ip)
                uzep(nep(iw,is),iw,is) = pp(11,ip)
              else if (iz(ip)+1 == izepwin(iw) .and.
     &                 nep(iw,is) < nepmax .and.
     &                 wp(ip) > 0.) then
!$OMP ATOMIC
                nep(iw,is) = nep(iw,is) + 1
                tep(nep(iw,is),iw,is) = pp(17,ip) + time
                xep(nep(iw,is),iw,is) = pp(2,ip)
                yep(nep(iw,is),iw,is) = pp(6,ip)
                uxep(nep(iw,is),iw,is) = pp(4,ip)
                uyep(nep(iw,is),iw,is) = pp(8,ip)
                uzep(nep(iw,is),iw,is) = pp(12,ip)
              endif
            enddo
          enddo

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

      if (itask == 3 .and. .not. laccumulate_zmoments) then
        write(0,*) 'Error in getzmmnt_weights: call getzmmnt when itask = 3'
        stop
      endif

!$OMP ATOMIC
      momentstime = momentstime + (wtime() - timetemp)

      return
      end
c=============================================================================
      subroutine getextrapolatedparticles(np,xp,yp,zp,uxp,uyp,uzp,gaminv,
     &                                    dt,uxpo,uypo,uzpo,is,time,zbeam)
      use Z_Moments
      use ExtPart
      use Timers
      integer(ISZ):: np,is
      real(kind=8):: dt,time,zbeam
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np)

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

      real(kind=8):: oneondt,dtstar,delz,vx,vy,vz,ax,ay,az
      integer(ISZ):: iw,ip,nn

      real(kind=8):: timetemp,wtime

      oneondt = 1./dt

c     --- Loop over windows
      do iw=1,nepwin

c       --- Skip windows with iz defined since those are already dealt with
c       --- in the moments calculation
        if (izepwin(iw) >= 0) cycle

c       --- Skip windows that are outside of the grid frame.
        if (zzepwin(iw)+wzepwin(iw) <= zmmntmin+zbeam .or.
     &      zzepwin(iw)-wzepwin(iw) >= zmmntmax+zbeam) cycle

c       --- Make sure that there is enough room - assuming that all particles
c       --- could be saved.
        if (nep(iw,is)+np > nepmax) then
          nepmax = nep(iw,is)+np
          call gchange("ExtPart",0)
        endif

c       --- extrapolate to window if near it
        do ip=1,np
          if (zp(ip) < zzepwin(iw)-wzepwin(iw) .or.
     &        zp(ip) > zzepwin(iw)+wzepwin(iw)) cycle
c         if (nep(iw,is) == nepmax) cycle

          vx = uxp(ip)*gaminv(ip)
          vy = uyp(ip)*gaminv(ip)
          vz = uzp(ip)*gaminv(ip)
          ax = (uxp(ip) - uxpo(ip))*oneondt
          ay = (uyp(ip) - uypo(ip))*oneondt
          az = (uzp(ip) - uzpo(ip))*oneondt

          delz = (zzepwin(iw) - zp(ip))
          dtstar = 2.*delz/(sqrt(vz**2 - 2.*az*delz) + vz)

          nep(iw,is) = nep(iw,is) + 1
          nn = nep(iw,is)

          tep(nn,iw,is) = time + dtstar
          xep(nn,iw,is) = xp(ip) + vx*dtstar + 0.5*ax*dtstar**2
          yep(nn,iw,is) = yp(ip) + vy*dtstar + 0.5*ay*dtstar**2
          uxep(nn,iw,is) = vx + ax*dtstar
          uyep(nn,iw,is) = vy + ay*dtstar
          uzep(nn,iw,is) = vz + az*dtstar

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
         if (zlw(iwin) /= 0.) nlabwn = iwin
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

      real(kind=8):: w0,w1
      integer(ISZ):: i,ilw,iz

      if (nlabwn == 0) return
      if (mod(it,itlabwn) > 0) return
 
      do i = 1, nlabwn
        if (zlw(i) > (zbeam+zmmntmin) .and.
     &      zlw(i) < (zbeam+zmmntmax)) then
          if (ilabwn(i) == ntlabwn) then
            ntlabwn = ntlabwn + max(10,int(ntlabwn/10))
            call gchange("Lab_Moments",0)
          endif
          ilabwn(i) = ilabwn(i) + 1
          ilw = ilabwn(i)
          iz = (zlw(i) - zbeam - zmmntmin)*dzmi
          w0 = (zlw(i) - zbeam - zmmntmin)*dzmi - iz
          w1 = 1. - w0
          timelw(ilw,i) = time
          pnumlw(ilw,i)   = w1*pnumz(iz)   + w0*pnumz(iz+1)
          xbarlw(ilw,i)   = w1*xbarz(iz)   + w0*xbarz(iz+1)
          ybarlw(ilw,i)   = w1*ybarz(iz)   + w0*ybarz(iz+1)
          xrmslw(ilw,i)   = w1*xrmsz(iz)   + w0*xrmsz(iz+1)
          yrmslw(ilw,i)   = w1*yrmsz(iz)   + w0*yrmsz(iz+1)
          rrmslw(ilw,i)   = w1*rrmsz(iz)   + w0*rrmsz(iz+1)
          xxpbarlw(ilw,i) = w1*xxpbarz(iz) + w0*xxpbarz(iz+1)
          yypbarlw(ilw,i) = w1*yypbarz(iz) + w0*yypbarz(iz+1)
          vzbarlw(ilw,i)  = w1*vzbarz(iz)  + w0*vzbarz(iz+1)
          epsxlw(ilw,i)   = w1*epsxz(iz)   + w0*epsxz(iz+1)
          epsylw(ilw,i)   = w1*epsyz(iz)   + w0*epsyz(iz+1)
          epszlw(ilw,i)   = w1*epszz(iz)   + w0*epszz(iz+1)
          vxrmslw(ilw,i)  = w1*vxrmsz(iz)  + w0*vxrmsz(iz+1)
          vyrmslw(ilw,i)  = w1*vyrmsz(iz)  + w0*vyrmsz(iz+1)
          vzrmslw(ilw,i)  = w1*vzrmsz(iz)  + w0*vzrmsz(iz+1)
          iz = (zlw(i) - zbeam - zzmin)*dzzi
          w0 = (zlw(i) - zbeam - zzmin)*dzzi - iz
          w1 = 1. - w0
          currlw(ilw,i)  = w1*curr(iz)   + w0*curr(iz+1)
          linechglw(ilw,i)  = w1*linechg(iz)   + w0*linechg(iz+1)
          lostparslw(ilw,i)  = w1*lostpars(iz)   + w0*lostpars(iz+1)
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
      subroutine setcurr (curr,np,zp,uzp,gaminv,q,wght,zbeam,dzz,zzmin)
      integer(ISZ):: np
      real(kind=8):: q,wght,zbeam,dzz,zzmin
      real(kind=8):: curr(0:*), zp(np), uzp(np), gaminv(np)

c  Sets 1d beam current as a function of z, directly from particle data

      real(kind=8):: s0(np), s1(np)
      integer(ISZ):: kk(np)
      real(kind=8):: g,w1,w0,dzzi
      integer(ISZ):: ip

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
!$OMP CRITICAL (CRITICAL_SETCURR1)
      do ip = 1,np
         if (uzp(ip) /= 0) then
            curr(kk(ip)) = curr(kk(ip)) + s0(ip)
            curr(kk(ip)+1) = curr(kk(ip)+1) + s1(ip)
         endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETCURR1)

      return
      end
c=============================================================================
      subroutine setcurr_weights (curr,np,zp,uzp,w,gaminv,q,wght,zbeam,dzz,zzmin)
      integer(ISZ):: np
      real(kind=8):: q,wght,zbeam,dzz,zzmin
      real(kind=8):: curr(0:*), zp(np), uzp(np), gaminv(np), w(np)

c  Sets 1d beam current as a function of z, directly from particle data

      real(kind=8):: s0(np), s1(np)
      integer(ISZ):: kk(np)
      real(kind=8):: g,w1,w0,dzzi
      integer(ISZ):: ip

      g = wght * q / dzz
      dzzi = 1. / dzz
c     --- vector loop
      do ip = 1,np
         kk(ip) = (zp(ip) - zbeam - zzmin) * dzzi
         w1 =  (zp(ip) - zbeam - zzmin) * dzzi - kk(ip)
         w0 = 1. - w1
         s0(ip) = w0 * g * gaminv(ip) * uzp(ip) * w(ip)
         s1(ip) = w1 * g * gaminv(ip) * uzp(ip) * w(ip)
      enddo
c     --- scalar loop
!$OMP CRITICAL (CRITICAL_SETCURR_WEIGHTS1)
      do ip = 1,np
         if (uzp(ip) /= 0) then
            curr(kk(ip)) = curr(kk(ip)) + s0(ip)
            curr(kk(ip)+1) = curr(kk(ip)+1) + s1(ip)
         endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETCURR_WEIGHTS1)

      return
      end
c=============================================================================
      subroutine percurr(curr,nzzarr,bound0)
      use GlobalVars
      integer(ISZ):: nzzarr,bound0
      real(kind=8):: curr(0:nzzarr)

c Make current periodic.

#ifdef MPIPARALLEL
c Parallel version sums curr over all processors.
      call parallelsumrealarray(curr,nzzarr+1)
#endif

      if (bound0==periodic) then
        curr(0) = curr(0) + curr(nzzarr)
        curr(nzzarr) = curr(0)
      endif

      return
      end
c=============================================================================
      subroutine getvzofz
      use Z_arrays

c  Computes the mean beam z velocity from the current and line charge density

      integer(ISZ):: i

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

      integer(ISZ):: i
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
      logical(ISZ):: lrelativ
c  Converts v to u, sets inverse gamma factor for all particles


      real(kind=8):: clghtisq,gamma
      integer(ISZ):: is,ip

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
      use GlobalVars
      use InGen
      use InPart
      use Particles
      real(kind=8):: zmmax,zmmin,dz,zgrid

c  Impose boundary conditions on zp

      real(kind=8):: syslen
      integer(ISZ):: is,ip

c Parallel version makes call to special routine zpartbnd_slave
#ifdef MPIPARALLEL
      call zpartbnd_slave(zmmax,zmmin,dz,zgrid)
#else

      if (pbound0==periodic) then
c     --- periodic boundary condition
        syslen = zmmax - zmmin
        do is=1,ns
          do ip=ins(is),ins(is)+nps(is)-1
            if (zp(ip) < zgrid + zmmin)
     &        zp(ip) = zp(ip) + syslen*int((zgrid + zmmax - zp(ip))/syslen)
            if (zp(ip) > zgrid + zmmax)
     &        zp(ip) = zp(ip) - syslen*int((zp(ip) - zgrid - zmmin)/syslen)
          enddo
        enddo
      else
c     --- other boundary condition (particles absorbed or reflected)
        do is=1,ns
          do ip=ins(is),ins(is)+nps(is)-1
            if ((zp(ip)-zgrid) >= zmmax) then
              if(pboundnz==absorb) then
                gaminv(ip) = 0.
              elseif(pboundnz==reflect) then
                zp(ip) = 2.*zgrid + 2.*zmmax - zp(ip)
                uzp(ip) = -uzp(ip)
              end if
            elseif ((zp(ip)-zgrid) <= zmmin) then
              if(pbound0==absorb) then
                gaminv(ip) = 0.
              elseif(pbound0==reflect) then
                zp(ip) = 2.*zgrid + 2.*zmmin - zp(ip)
                uzp(ip) = -uzp(ip)
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
      subroutine reorgparticles()
#ifdef MPIPARALLEL
      call reorgparticles_parallel()
#endif
      return
      end
c=============================================================================
      subroutine semitransparent_disc(dz)
      use Particles
      use InGen
      use Moments
      use InPart
      use SemiTransparentDisc
c Random absorption of particles passing through a semitranparent disc.
c This assumes that no particle travels more than dz in one time step.
c dz is used only for efficiency purpose: only particles within dz 
c of semi-transparent disc location are checked for disc crossing.

      real(kind=8), intent(in) :: dz 
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
      do is=1,ns
c       --- loop on particles
        do ip=ins(is),ins(is)+nps(is)-1
c         --- loop on number of discs
          do idisc = 1, n_STdiscs
c           --- do only the work for particles close to disc location
            if(zp(ip)>=zmind(idisc) .and. zp(ip)<=zmaxd(idisc)) then
              zpold = zp(ip)-dt*uzp(ip)*gaminv(ip)
c             --- test if particle crossed disc location
              if(zpold<z_STdiscs(idisc) .and. zp(ip)>=z_STdiscs(idisc)) then
                dtcross = dt*(zp(ip)-z_STdiscs(idisc))/(zp(ip)-zpold)
                xpatdisc = xp(ip)-dtcross*uxp(ip)*gaminv(ip)
                ypatdisc = yp(ip)-dtcross*uyp(ip)*gaminv(ip)
                r = sqrt(xpatdisc**2+ypatdisc**2)
c               --- test if particle passed through disc
                if(r<=r_STdiscs(idisc)) then
c                  --- randomly absorb particles according to transparency coefficient
                   if(wranf()>t_STdiscs(idisc)) then
                     gaminv(ip) = 0.
                   end if
                end if
              end if
            end if
          end do
        end do
      end do
      
      end
c=============================================================================
c A general particle scraping routine, which allows particles to be scraped
c efficiently using complex geometries.
c=============================================================================
      subroutine griddedparticlescraper(is,distance,nx,ny,nz,dx,dy,dz,
     &                                  xmin,ymin,zmin,zbeam,l2symtry,l4symtry)
      use Particles
      use Scraped_Particles
      use Z_arrays
      use InGen
      integer(ISZ):: is,nx,ny,nz
      real(kind=8):: dx,dy,dz,xmin,ymin,zmin,zbeam
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: distance(0:nx,0:ny,0:nz)

c This provides a method of scraping particles in a complex geometry.
c The 3-d array distances holds the distance of each node to the nearest
c scraping boundary, a negative value being inside.
c This distance is linearly interpolated to the particle
c position. If the resulting value is negative, then the particle is
c considered inside the boundary and is removed.

c Input:
c   is: species number of particles to process
c   distance(0:nx,0:ny,0:nz): distance to boundary surface
c   nx,ny,nz: grid dimensions
c   dx,dy,dz: grid cell sizes
c   xmin,ymin,zmin: coordinates of lower left hand corner of grid
c   zbeam: location of the beam frame (only needed to set lostpars)
c   l2symtry,l4symtry: level of symmetry to use

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz
      real(kind=8):: dd,xpold,ypold,zpold,ddratio
 
c     --- set grid cell inverse sizes
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c     --- initializes number of scraped particles
      if(is==1) scr_np = 0
      scr_ins(is) = scr_np + 1
      scr_nps(is) = scr_np

c     --- loop over particles
      do ip = ins(is), ins(is) + nps(is)

c       --- find location on grid
        if (l4symtry) then
          ix = abs((xp(ip) - xmin)*dxi)
          iy = abs((yp(ip) - ymin)*dyi)
          iz =     (zp(ip) - zmin)*dzi
          wx = abs((xp(ip) - xmin)*dxi) - ix
          wy = abs((yp(ip) - ymin)*dyi) - iy
          wz =     (zp(ip) - zmin)*dzi - iz
        else if (l2symtry) then
          ix =     (xp(ip) - xmin)*dxi
          iy = abs((yp(ip) - ymin)*dyi)
          iz =     (zp(ip) - zmin)*dzi
          wx =     (xp(ip) - xmin)*dxi - ix
          wy = abs((yp(ip) - ymin)*dyi) - iy
          wz =     (zp(ip) - zmin)*dzi - iz
        else
          ix =     (xp(ip) - xmin)*dxi
          iy =     (yp(ip) - ymin)*dyi
          iz =     (zp(ip) - zmin)*dzi
          wx =     (xp(ip) - xmin)*dxi - ix
          wy =     (yp(ip) - ymin)*dyi - iy
          wz =     (zp(ip) - zmin)*dzi - iz
        endif
 
c       --- if within grid, examine the particle, otherwise assume its OK.
        if (ix >= 0 .and. ix < nx .and.
     &      iy >= 0 .and. iy < ny .and.
     &      iz >= 0 .and. iz < nz) then

          dd = distance(ix  ,iy  ,iz  )*(1. - wx)*(1. - wy)*(1. - wz) +
     &         distance(ix+1,iy  ,iz  )*      wx *(1. - wy)*(1. - wz) +
     &         distance(ix  ,iy+1,iz  )*(1. - wx)*      wy *(1. - wz) +
     &         distance(ix+1,iy+1,iz  )*      wx *      wy *(1. - wz) +
     &         distance(ix  ,iy  ,iz+1)*(1. - wx)*(1. - wy)*      wz  +
     &         distance(ix+1,iy  ,iz+1)*      wx *(1. - wy)*      wz  +
     &         distance(ix  ,iy+1,iz+1)*(1. - wx)*      wy *      wz  +
     &         distance(ix+1,iy+1,iz+1)*      wx *      wy *      wz

          if (dd < 0.) then
c           --- Scrape the particle
            if (uzp(ip) /= 0.) lostpars(iz) = lostpars(iz) + 1
            scr_np = scr_np + 1
            scr_nps(is) = scr_nps(is) + 1
            if(scr_npmax<scr_np) then
              scr_npmax = scr_npmax + scr_npbunch
              call gchange("Scraped_Particles",0)
            end if
            xpold = xp(ip) - uxp(ip) * gaminv(ip) * dt
            ypold = yp(ip) - uyp(ip) * gaminv(ip) * dt
            zpold = zp(ip) - uzp(ip) * gaminv(ip) * dt
            ddratio = min(1.,1.00001*abs(dd)/sqrt((xp(ip)-xpold)**2
     &                                           +(yp(ip)-ypold)**2
     &                                           +(zp(ip)-zpold)**2))
            scr_xp(scr_np) = xp(ip) - ddratio * uxp(ip) * gaminv(ip)*dt
            scr_yp(scr_np) = yp(ip) - ddratio * uyp(ip) * gaminv(ip)*dt
            scr_zp(scr_np) = zp(ip) - ddratio * uzp(ip) * gaminv(ip)*dt
            scr_uxp(scr_np) = uxp(ip)
            scr_uyp(scr_np) = uyp(ip)
            scr_uzp(scr_np) = uzp(ip)
            iz = max(0,min(nzzarr,int((zp(ip) - zzmin - zbeam)*dzzi)))
            if (uzp(ip) /= 0.) then
              iz = max(0,min(nzzarr,int((zp(ip) - zzmin - zbeam)*dzzi)))
              lostpars(iz) = lostpars(iz) + 1
            endif
            gaminv(ip) = 0.
          endif

        endif
 
      enddo

      return
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
      subroutine zbeamcom(zbeam)
      use InPart
      use Particles
      real(kind=8):: zbeam

c Calculate the center of mass of the beam along the z axis from the particles


      integer(ISZ):: is,ip
      real(kind=8):: totz,totmass

      totz = 0.
      totmass = 0.
      do is=1,ns
        do ip=ins(is),ins(is)+nps(is)-1
          if (uzp(ip) /= 0.) then
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

      integer(ISZ):: is

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
        if (np_s(is) == 0) np_s(is) = npmax*sp_fract(is)
      enddo
#ifndef MPIPARALLEL
      npmax = 0
      do is=1,ns
        npmax = npmax + np_s(is)
      enddo
      np = npmax
#endif
      npmaxb = npmax
      if (npid > 0) then
        npmaxi = npmax
      else
        npmaxi = 1
      endif
      npidmax = max(1,npid)

c     --- Allocate space. A gchange is done primarily so that any data
c     --- in npmax_s is preserved.
      call gchange("Particles",0)

c     --- setup npmax_s array
      npmax_s(0) = 0
      do is=1,ns
        if (npmax_s(is) == 0) npmax_s(is) = npmax_s(is-1) + np_s(is)
        ins(is) = npmax_s(is-1) + 1
        nps(is) = 0
      enddo

      return
      end
c=============================================================================
      subroutine chckpart(is,nlower,nhigher,lfullshft)
      use InPart
      use Particles
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

      integer(ISZ):: nadd,i,ishft

c     --- If there is already enough space, then return.
      if (ins(is)-npmax_s(is-1)-1 >= nlower .and.
     &    npmax_s(is)-ins(is)-nps(is)+1 >= nhigher) return

c     --- If there is not enough space, then allocate new space.
      if (ins(is)-npmax_s(is-1)-1 + npmax_s(is)-ins(is)-nps(is)+1 <
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
        if (npid > 0) then
          npmaxi = npmax
        else
          npmaxi = 1
        endif
        npidmax = max(1,npid)
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
      integer(ISZ):: is,ishft
      call shftpartwork(is,ishft,npmax,npmaxi,ns,npmax_s,ins,nps,npid,
     &                  xp,yp,zp,uxp,uyp,uzp,gaminv,pid)
      return
      end
c=============================================================================
      subroutine shftpartwork(is,ishft,npmax,npmaxi,ns,npmax_s,ins,nps,npid,
     &                        xp,yp,zp,uxp,uyp,uzp,gaminv,pid)
      integer(ISZ):: is,ishft
      integer(ISZ):: npmax,npmaxi,ns,npmax_s(0:ns),ins(ns),nps(ns),npid
      real(kind=8):: xp(npmax),yp(npmax),zp(npmax)
      real(kind=8):: uxp(npmax),uyp(npmax),uzp(npmax)
      real(kind=8):: gaminv(npmax),pid(npmaxi,npid)

c Shift particles by 'ishft'.  All of the data from the low end of the
c species block to the top of the live particles are moved.  The area
c below the live particles may include newly injected particles which
c have not been included into the simulation yet.
c Assumes that the is enought space in the particle arrays to make the
c shift without clobbering other particles or going past the ends of
c the arrays.

      integer(ISZ):: ip

c     --- return if shift is zero
      if (ishft == 0) return

c     --- Positive and negative shifts are treated differently since
c     --- a negative shift may result in lost data.
      if (ishft > 0) then

c       --- Loop from top of live particles to bottom of species block.
        do ip=ins(is)+nps(is)-1, npmax_s(is-1)+1, -1
          xp(ip+ishft) = xp(ip)
          yp(ip+ishft) = yp(ip)
          zp(ip+ishft) = zp(ip)
          uxp(ip+ishft) = uxp(ip)
          uyp(ip+ishft) = uyp(ip)
          uzp(ip+ishft) = uzp(ip)
          gaminv(ip+ishft) = gaminv(ip)
          if (npid > 0) then
            pid(ip+ishft,:) = pid(ip,:)
            pid(ip,:) = 0
          endif
          uzp(ip) = 0.
        enddo

c       --- change lower bound of live particles
        ins(is) = ins(is) + ishft

      else
c       --- Loop from -ishft above bottom of species block to top of live
c       --- particles.  Note that particles will be clobbered if
c       --- (ins(is) - npmax_s(is-1)) < -ishft
        if (ins(is) - npmax_s(is-1) < -ishft) then
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
          if (npid > 0) then
            pid(ip+ishft,:) = pid(ip,:)
            pid(ip,:) = 0
          endif
          uzp(ip) = 0.
        enddo

c       --- change number of particles if data is lost
        if (ins(is) - npmax_s(is-1) < -ishft) then
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
      integer(ISZ):: it,nn,istart
      integer(ISZ):: ii(0:nn-1)

c Copies particle data from locations given by ii to the locations starting
c at 'it', or if istart > 0, copy sequential particles from locations
c starting at 'istart' to locations starting at 'it'.

      integer(ISZ):: i,in,ij
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
          if (npid > 0) then
            pid(in,:) = pid(ij,:)
          endif
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
          if (npid > 0) then
            pid(in,:) = pid(ij,:)
          endif
        enddo
      endif
      return
      end
c=============================================================================
      subroutine addpart(nn,nnpid,x,y,z,vx,vy,vz,gi,ppid,is,lallindomain,zmmin,zmmax,
     &                   lmomentum)
      use Constant
      use Particles
      use Beam_acc
      integer(ISZ):: is,nn,nnpid
      real(kind=8):: x(nn),y(nn),z(nn),vx(nn),vy(nn),vz(nn),gi(nn),ppid(nn,nnpid)
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

      integer(ISZ):: ip,n1,i1,i2

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
      call chckpart(is,0,n1,.false.)

c     --- Now, copy data into arrays
      i1 = ins(is) + nps(is)
      i2 = i1 + n1 - 1
      xp(i1:i2) = x(1:n1)
      yp(i1:i2) = y(1:n1)
      zp(i1:i2) = z(1:n1)
      uxp(i1:i2) = vx(1:n1)
      uyp(i1:i2) = vy(1:n1)
      uzp(i1:i2) = vz(1:n1)
      gaminv(i1:i2) = gi(1:n1)
      if (nnpid > 0) pid(i1:i2,1:nnpid) = ppid(1:n1,:)
      nps(is) = nps(is) + n1

      return
      end
c=============================================================================
      subroutine clearpart(js,fillmethod)
      use InPart
      use Particles
      integer(ISZ):: js,fillmethod

c Clears out all of the lost particles. If js < 0, the clears out all
c species, otherwise only species requested.
c The empty spaces can be fill in multiple ways, based on value of fillmethod
c   0 : don't do anything
c   1 : fill spaces with particles from the end of the list
c   2 : shift particles above down (they keep the same order)

      integer(ISZ):: js1,js2,j,ip
      integer(ISZ):: ips,ds,ipslimit

      if (js < 0) then
        js1 = 1
        js2 = ns
      else
        js1 = js
        js2 = js
      endif

      if (fillmethod == 1) then

c       --- Live particles from the end of the species block are moved
c       --- to fill in the empty spaces.
        do j=js1,js2
          ips = ins(j) + nps(j) - 1
          ip = ins(j)
          nps(j) = 0
          do while (ip <= ips)
            if (uzp(ip) == 0.) then
              do while (uzp(ips) == 0. .and. ips > ip)
                ips = ips - 1
              enddo
              if (ips > ip) then
                xp(ip) = xp(ips)
                yp(ip) = yp(ips)
                zp(ip) = zp(ips)
                uxp(ip) = uxp(ips)
                uyp(ip) = uyp(ips)
                uzp(ip) = uzp(ips)
                gaminv(ip) = gaminv(ips)
                if (npid > 0) then
                  pid(ip,:) = pid(ips,:)
                  pid(ips,:) = 0
                endif
                uzp(ips) = 0.
                gaminv(ips) = 0.
                if(wpid>0) pid(ips,wpid) = 0.
                ips = ips - 1
                nps(j) = nps(j) + 1
              endif
            else
              nps(j) = nps(j) + 1
            endif
            ip = ip + 1
          enddo
        enddo

      else if (fillmethod == 2) then

c       --- Live particles are shifted downward to fill in the empty spaces
        do j=js1,js2
          ipslimit = ins(j)+nps(j)-1
          nps(j) = 0
          ip = ins(j)
          do ips=ip,ipslimit
            if (uzp(ips) /= 0.) then
              if (ips > ip) then
                xp(ip) = xp(ips)
                yp(ip) = yp(ips)
                zp(ip) = zp(ips)
                uxp(ip) = uxp(ips)
                uyp(ip) = uyp(ips)
                uzp(ip) = uzp(ips)
                gaminv(ip) = gaminv(ips)
                if (npid > 0) then
                  pid(ip,:) = pid(ips,:)
                  pid(ips,:) = 0
                endif
                uzp(ips) = 0.
                gaminv(ips) = 0.
                if(wpid>0) pid(ips,wpid) = 0.
              endif
              ip = ip + 1
              nps(j) = nps(j) + 1
            endif
          enddo
        enddo

      endif

      return
      end
c=============================================================================
      subroutine processlostpart(is,clearlostpart,time,zbeam)
      use Constant
      use Beam_acc
      use Particles
      use LostParticles
      use Z_arrays
      integer(ISZ):: is,clearlostpart
      real(kind=8):: time,zbeam

c Processes lost particles, which are flagged by having gaminv set to zero.
c The particle data is optionally saved. The velocities are set to zero and
c the routine which clears out the lost particles is called.

      integer(ISZ):: nlost,i1,i2
      integer(ISZ):: ip,iz
      real(kind=8):: usq

      if (lsavelostpart) then

c       --- Make sure that if npidmax has changed, then npidlostmax is also
c       --- changed.
        if (npidmax > npidlostmax) then
          npidlostmax = npidmax
          call gchange("LostParticles",0)
        endif

c       --- Get the number of lost particles and make sure there is
c       --- enough space to save them in the LostParticles arrays.
        i1 = ins(is)
        i2 = ins(is)+nps(is)-1
        nlost = count((gaminv(i1:i2)==0 .and. uzp(i1:i2)/=0.))
        call chcklostpart(is,0,nlost)

      endif

c     --- Loop over the particles arrays and process the ones with gaminv
c     --- set to zero, the flag for scraped particles.
      do ip=ins(is),ins(is)+nps(is)-1
        if (gaminv(ip) == 0.) then

c         --- Accumulate the array holding the location where particles
c         --- are lost.
          if (uzp(ip) /= 0.) then
            iz = max(0,min(nzzarr,int((zp(ip) - zzmin - zbeam)*dzzi)))
            lostpars(iz) = lostpars(iz) + 1
          endif

c         --- If requested to save lost particles, do it. Particle data is
c         --- only save if uzp /= 0, since if uzp is 0, the particle was
c         --- lost on a previous step but wasn't cleared out of the array.
          if (lsavelostpart .and. uzp(ip) /= 0.) then
            i1 = inslost(is) + npslost(is)
            xplost(i1) = xp(ip)
            yplost(i1) = yp(ip)
            zplost(i1) = zp(ip)
            uxplost(i1) = uxp(ip)
            uyplost(i1) = uyp(ip)
            uzplost(i1) = uzp(ip)
            if (lrelativ) then
              usq = uxp(ip)**2 + uyp(ip)**2 + uzp(ip)**2
              gaminvlost(i1) = 1./sqrt(1. + usq/clight**2)
            else
              gaminvlost(i1) = 1.
            endif
            tplost(i1) = time
            if (npid > 0) pidlost(i1,1:npidmax) = pid(ip,:)
            npslost(is) = npslost(is) + 1
          endif

c         --- Zero out the velocity, the flag for a lost particle.
          uxp(ip) = 0.0
          uyp(ip) = 0.0
          uzp(ip) = 0.0
          if(wpid>0) pid(ip,wpid) = 0.

        endif
      enddo

      call clearpart(is,clearlostpart)

      return
      end
c=============================================================================
      subroutine alotlostpart
      use InPart
      use Particles
      use LostParticles

c Allocate space for lost particles and set up associated arrays.

      integer(ISZ):: is

      npidlostmax = max(npidmax,npidlostmax)
      call gchange("LostParticles",0)

      npmaxlost_s(0) = 0
      do is=1,ns
        npmaxlost_s(is) = npmaxlost_s(is-1) + npmaxlost*sp_fract(is)
        inslost(is) = npmaxlost_s(is-1) + 1
        npslost(is) = 0
      enddo

      return
      end
c=============================================================================
      subroutine chcklostpart(is,nlower,nhigher)
      use InPart
      use LostParticles
      integer(ISZ):: is,nlower,nhigher

c Make sure that there is enough space in the lost particle arrays for nlower
c new particles below and nhigher above the lost particles.  Returns if
c there is already enough space above and below.  If there is enough total
c space but not enough room above or below, the lost particles are shifted
c appropriately. If there is not enough space, add more to the arrays.
c Particle data is shifted appropriately.

      integer(ISZ):: nadd,i,ishft

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
      ishft = max(0,nlower - (inslost(is)-npmaxlost_s(is-1)-1)) +
     &        min(0,(npmaxlost_s(is)-inslost(is)-npslost(is)+1) - nhigher)
      call shftlostpart(is,ishft)

      return
      end
c=============================================================================
      subroutine shftlostpart(is,ishft)
      use InPart
      use LostParticles
      use Particles
      integer(ISZ):: is,ishft
      call shftpartwork(is,ishft,npmaxlost,npmaxlost,ns,npmaxlost_s,
     &                  inslost,npslost,npid,
     &                  xplost,yplost,zplost,uxplost,uyplost,uzplost,
     &                  gaminvlost,pidlost)
      return
      end
c=============================================================================
c  VARIOUS UTILITIES
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

      real(kind=8):: dyi,ni,nip1,ninti,sumn,sumnx,wy
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
      subroutine gettemperature(np,xp,yp,zp,uxp,uyp,uzp,gaminv,q,m,w,dt,
     &                          itask,uxpo,uypo,uzpo,is,wp,lw,laxisymmetric)
c     computes temperature in slices.
      use Constant
      use Beam_acc
      use InDiag
      use Picglb
      use ExtPart
      use Temperatures
      use Timers
      integer(ISZ):: np,itask,is
      real(kind=8):: q,m,w,dt
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np), wp(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np)
      logical(ISZ):: lw,laxisymmetric

c     Three parts
c     When itask=1  zeros out all moments
c          itask=2  sums moments from particles
c          itask=3  divides by number of particles, calculates emittances and rms

      integer(ISZ):: ip,its,izl,ixt,iyt,izt
      real(kind=8):: dti,wt,ddx,ddy,oddx,oddy,xt,yt,tfact,wpp
      real(kind=8):: oneondt,clighti,vzi,pnumi,zc,tottmp
      real(kind=8):: xpt,ypt,zpt,vxpt,vypt,vzpt,cost,sint,rpt
      real(kind=8):: delxsq,delvxsq,delxvx
      real(kind=8):: delysq,delvysq,delyvy
      real(kind=8):: delzsq,delvzsq,delzvz
      real(kind=8):: timetemp,wtime

      timetemp = wtime()

      if (.not.l_temp .or. np==0) return

      wpp = 1.

      oneondt = 1./dvnz(dt)

      if (itask == 1 .or. itask == -1 .and. .not. l_accumulate_temperatures) then
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
      endif

      if (itask == 2 .or. itask == -1) then
c       --- loop over particles
        do ip=1,np
          izl  = 1+int((zp(ip) - zbeam - tloc_zmin)*tloc_dzi)

c         --- cycle if particle out of zone of calculation
          if(izl<1 .or. izl>nztlocator) cycle

c         --- loop over temperature slices
          do its = 1, ntl(izl)
            izt = tslice_locator(izl,its)

c           --- cycle if particle not in slice
            if(xp(ip)<=tslicexmin(izt) .or. xp(ip)>tslicexmax(izt) .or.
     &         yp(ip)<=tsliceymin(izt) .or. yp(ip)>tsliceymax(izt) .or.
     &         zp(ip)<=tslicezmin(izt) .or. zp(ip)>tslicezmax(izt)) cycle

            if(l_temp_collapseinz) then
c             --- collapse slice in z
              zc = 0.5*(tslicezmin(izt)+tslicezmax(izt))
              vzi = 1./(uzp(ip)*gaminv(ip)+SMALLPOS)
              dti  = (zbeam+zc-zp(ip))*vzi
              xpt  = xp(ip) + uxp(ip)*dti*gaminv(ip)
              ypt  = yp(ip) + uyp(ip)*dti*gaminv(ip)
              zpt  = zp(ip) + uzp(ip)*dti*gaminv(ip)
              vxpt = (uxp(ip)*(1. + dti*oneondt) - uxpo(ip)*dti*oneondt) * gaminv(ip)
              vypt = (uyp(ip)*(1. + dti*oneondt) - uypo(ip)*dti*oneondt) * gaminv(ip)
              vzpt = (uzp(ip)*(1. + dti*oneondt) - uzpo(ip)*dti*oneondt) * gaminv(ip)
            else
              xpt  = xp(ip)
              ypt  = yp(ip)
              zpt  = zp(ip)
              vxpt = uxp(ip) * gaminv(ip)
              vypt = uyp(ip) * gaminv(ip)
              vzpt = uzp(ip) * gaminv(ip)
            end if
            if(laxisymmetric) then
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
            end if

c           --- compute coefficients for transverse linear deposition
            xt = (xp(ip)-tslicexmin(izt))*dxti(izt)
            yt = (yp(ip)-tsliceymin(izt))*dyti(izt)
            ixt = max(0,min(nxtslices-1,int(xt)))
            iyt = max(0,min(nytslices-1,int(yt)))
            ddx = xt-ixt
            ddy = yt-iyt
            oddx = 1.-ddx
            oddy = 1.-ddy

            if(lw) wpp = wp(ip)

c           --- perform transverse linear deposition
c           --- node x, y
            wt = oddx*oddy*wpp
            pnumt  (ixt,iyt,izt) = pnumt  (ixt,iyt,izt) + oddx*oddy
            pnumtw (ixt,iyt,izt) = pnumtw (ixt,iyt,izt) + wt

            vxbart  (ixt,iyt,izt) = vxbart  (ixt,iyt,izt) + wt * vxpt
            vybart  (ixt,iyt,izt) = vybart  (ixt,iyt,izt) + wt * vypt
            vzbart  (ixt,iyt,izt) = vzbart  (ixt,iyt,izt) + wt * vzpt

            vxsqbart(ixt,iyt,izt) = vxsqbart(ixt,iyt,izt) + wt * vxpt**2
            vysqbart(ixt,iyt,izt) = vysqbart(ixt,iyt,izt) + wt * vypt**2
            vzsqbart(ixt,iyt,izt) = vzsqbart(ixt,iyt,izt) + wt * vzpt**2

c           --- node x+1, y
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

c           --- node x, y+1
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

c           --- node x+1, y+1
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
c             --- node x, y
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

c             --- node x+1, y
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

c             --- node x, y+1
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

c             --- node x+1, y+1
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
      endif

      if (itask == 3 .or. itask == -1 .and. .not. l_accumulate_temperatures) then

c       --- For slave, call routine which sums moments over processors.
#ifdef MPIPARALLEL
c        call sum_mmnts
       write(0,*) 'Error in gettemperature. This is not yet implemented in parallel.'
       stop
#endif

c       --- Complete the calculation of moments: divide by particle number

c       --- Now scale the full Z moments by particle accumulated weights
        tfact = m/echarge
        do izt = 1, nztslices
          do iyt = 0, nytslices
            do ixt = 0, nxtslices

              pnumi = 1./(pnumtw(ixt,iyt,izt)+SMALLPOS)

c             --- Compute averages
              vxbart  (ixt,iyt,izt) = vxbart  (ixt,iyt,izt) * pnumi
              vybart  (ixt,iyt,izt) = vybart  (ixt,iyt,izt) * pnumi
              vzbart  (ixt,iyt,izt) = vzbart  (ixt,iyt,izt) * pnumi

              vxsqbart(ixt,iyt,izt) = vxsqbart(ixt,iyt,izt) * pnumi
              vysqbart(ixt,iyt,izt) = vysqbart(ixt,iyt,izt) * pnumi
              vzsqbart(ixt,iyt,izt) = vzsqbart(ixt,iyt,izt) * pnumi

c             --- Compute second order moments with averages subtracted
              delvxsq = vxsqbart(ixt,iyt,izt) - vxbart(ixt,iyt,izt)**2
              delvysq = vysqbart(ixt,iyt,izt) - vybart(ixt,iyt,izt)**2
              delvzsq = vzsqbart(ixt,iyt,izt) - vzbart(ixt,iyt,izt)**2

c             --- Compute linear correlations
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

c             --- Compute temperatures
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
      endif

      temperaturestime = temperaturestime + (wtime() - timetemp)
      return
      end
c=============================================================================
      subroutine gett(is,itask,laxisymmetric)
        use InGen
        use Particles
        integer(ISZ):: is, ipmin, itask
        logical(ISZ):: laxisymmetric

        ipmin = ins(is)
        if(wpid>0) then
          call gettemperature(nps(is),xp(ipmin),yp(ipmin),zp(ipmin),
     &                           uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                           gaminv(ipmin),sq(is),sm(is),sw(is),dt,itask,
     &                           uxp(ipmin),uyp(ipmin),uzp(ipmin),is,
     7                           pid(ipmin,wpid),.true.,laxisymmetric)
        else
          call gettemperature(nps(is),xp(ipmin),yp(ipmin),zp(ipmin),
     &                           uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                           gaminv(ipmin),sq(is),sm(is),sw(is),dt,itask,
     &                           uxp(ipmin),uyp(ipmin),uzp(ipmin),is,
     7                           pid(ipmin,1),.false.,laxisymmetric)
        end if
        return
      end
c=============================================================================
      subroutine setregulartgrid(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,dz,nzloc,lcollapse,lcorrel)
        use Temperatures
        use InPart

        integer(ISZ) :: nx, ny, nz, nzloc
        logical(ISZ) :: lcollapse, lcorrel
        real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax, dz

        integer(ISZ) :: i

        nstmp = ns

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
        do i = 1, nztslices
          tslicezmin(i) = (i-1)*(zmax-zmin)/(nztslices-1)-0.5*dz
          tslicezmax(i) = tslicezmin(i) + dz
        end do
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

        allocate(ntsloc(nztlocator,nztslices))

        write(0,*) nztlocator,nztslices

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
        write(0,*) 'deallocate'!,allocated(ntsloc)
        deallocate(ntsloc)
        write(0,*) 'return'

        return
      end
c=============================================================================
      subroutine impact_ion(is1,is2,nbp,wp,shiftx,shifty,shiftz,deltax,deltay,deltaz,condid)
      USE Constant
      use Particles
      USE InPart
      use LostParticles
      INTEGER, INTENT(IN) :: is1, is2, nbp, condid
      REAL(8), INTENT(IN) :: wp,shiftx,shifty,shiftz,deltax,deltay,deltaz

c Create particles of specie is2 created by impact of particles of species is1
c - nbp: number of new particles are generated for each lost partcle
c - wp: energy in volts of new particles
c - shiftz,y,z: the location where the new particles are created relative
c               to where the original particles are lost
c - deltax,y,z: the length over which the new particles are distributed
c               relative to the location where the original particles are lost.
c               The distribution is uniform between -delta/2 and +delta/2.
c - condid: id of conductor - only particles that are lost on the specified
c           conductor generate new particles.

      INTEGER :: ip,i1,i2,ir
      REAL(8) :: gamma, uz_tmp
      REAL(8), EXTERNAL :: wranf

        IF(nbp==0) return

        gamma = 1.+wp*echarge/(sm(is2)*clight**2)
        IF(gamma-1.<1.e-6) then
          uz_tmp = gamma*SQRT(1.-1./(gamma**2))*clight
        else
          uz_tmp = SQRT(2.*wp*echarge/sm(is2))
        END if
        call chckpart(is2,0,nbp*npslost(is1),.false.)
        do ip = inslost(is1), inslost(is1)+npslost(is1)-1
          if (condid > 0 .and. pidlost(ip,npidlostmax) .ne. condid) cycle
          i1 = ins(is2) + nps(is2)
          i2 = i1 + nbp - 1
          do ir = 1, nbp
            xp(i1+ir-1) = xplost(ip)+shiftx+(wranf()-0.5)*deltax
            yp(i1+ir-1) = yplost(ip)+shifty+(wranf()-0.5)*deltay
            zp(i1+ir-1) = zplost(ip)+shiftz+(wranf()-0.5)*deltaz
          end do
          uxp(i1:i2) = 0.
          uyp(i1:i2) = 0.
          uzp(i1:i2) = -uz_tmp
          if(wpid>0) pid(i1:i2,wpid) = pidlost(ip,wpid)
          nps(is2) = nps(is2) + nbp
        end do

      return
      end subroutine impact_ion
c=============================================================================


