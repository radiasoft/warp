#include "top.h"
c=======================================================================
c@(#) File FRZ.M, version $Revision: 3.13 $, $Date: 2002/04/23 14:15:47 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  Written by Debbie Callahan
c  It is based on the package SRC.F3D by Alex Friedman.
c  It is the source file for the package FRZ of the PIC code WARP6,
c  but it may be useful by itself.  It contains:
c
c     1) Dale Nielsen's vectorized version of Langdon's cpft, VCPFT
c     2) Vectorized two at a time real(kind=8):: transforms: VRPFT2, VRPFTI2
c     4) A vectorized periodic (z) and tridiagonal (r) Poisson solver: VPOISRZ
c     5) A BASIS test driver (invoke via PACKAGE FRZ;GENERATE)
c        The arrays and scalars used are local to package FRZ; some have the
c        same names as those in the associated particle code package W6,
c        but they are <<not>> the same variables.  After calling this driver
c        and examining the results, the user should type FINISH to deallocate.
c
c  Call chain for the vector rz fieldsolve is (roughly):
c  VPOISRZ - VCPFT,VRPFT2,TRIDI,VRPFTI2
c
c=======================================================================
      subroutine frzinit
      use FRZversion

c  Called at first reference to package (not nec. a "run" etc.).


      call frzvers (STDOUT)

      return
      end
c=============================================================================
      subroutine frzvers (iout)
      use FRZversion
      integer(ISZ):: iout
c  Echoes code version,etc. to output files when they're created
      call printpkgversion(iout,"Fieldsolver FRZ",versfrz)
      return
      end
c=============================================================================
      subroutine frzgen()
      use Constant
      use FRZvars
      real(kind=8):: time1,time2
      real(kind=8):: lr4

c  Invoked by the GENERATE command, it tests the fieldsolver.
c  (The only BASIS/MPPL dependent coding is use of STDOUT, return (value),
c  and the next few lines).

c  Bring in the common block; allocate dynamic arrays ("group allot")

      integer(ISZ):: ir,iz
      real(kind=8):: dr,dz,sz,cz,fact1,fact2,rjph,rjmh,dr2j,dr20
      real(kind=8):: errmax,error,dz2,delsq,dt,a0,temp
      real(kind=8):: wtimeoff,wranf

      call gallot("FRZvars",0)

c  Check for improper input; we're using power-of-2 fft's

      if  ( nz .ne. 2**int( .5 + log(real(nz)) / log(2.) ) ) then
         write (STDOUT,17) nz
   17    format ( " ---> Bad nz ... must be power of 2",/,
     &            " nz =",i4)
         return
      endif

c  Echo system size to the standard output (generally the user's terminal)

      write (STDOUT,19) nr, nz
   19 format( /," ***  testing rz solver  ***",/,
     &          "      nr =",i4,"  nz =",i4)

c  Set constants: pi, grid spacing

      dr = lr / nr
      dz = lz / nz

c **** Test on a known function ****
c      phi(r,z) = (lr**4 - r**4)*(sin(6*pi*iz/nz)+cos(4.*pi*iz/nz))

      lr4 = lr**4
      do 99 ir = 0,nr
         r(ir) = dr*float(ir)
 99   continue

      sz = 6.*pi/nz
      cz = 4.*pi/nz
      fact1 = (2./dz*sin(sz/2.))**2
      fact2 = (2./dz*sin(cz/2.))**2
      do iz=0,nz
        do ir=0,nr
          bsav(ir,iz) = (lr4 - r(ir)**4)
        enddo
      enddo
      do ir = 1, nr-1
         rjph = 0.5*(r(ir+1) + r(ir))
         rjmh = 0.5*(r(ir) + r(ir-1))
         dr2j = rjph**2 - rjmh**2
         do iz = 0, nz
            b(ir,iz) =(-2.*rjph/(dr2j*dr)*(bsav(ir+1,iz) - bsav(ir,iz))
     &                +2.*rjmh/(dr2j*dr)*(bsav(ir,iz) - bsav(ir-1,iz)))
     &                *(sin(sz*iz)+cos(cz*iz))
     &                +bsav(ir,iz)*(fact1*sin(sz*iz)+fact2*cos(cz*iz))
            b(ir,iz) = b(ir,iz)*eps0
         enddo
      enddo
      do 101 iz = 0,nz
         dr20 = (0.5*r(1))**2
         b(0,iz) = (-r(1)/(dr20*dr)*(bsav(1,iz) - bsav(0,iz)))
     &             *(sin(sz*iz)+cos(cz*iz))
     &             +bsav(0,iz)*(fact1*sin(sz*iz)+fact2*cos(cz*iz))
         b(0,iz) = eps0*b(0,iz)
         b(nr,iz) = 0.0
  101 continue
c     --- No resistivity
      eta = 0.
      taurc = 0.

c     --- Call the solver
      call wtimeon
      call vpoisrz (0,b,b,kzsq,schrg,eta,phikold,taurc,attz,filt,
     &              0.,0.,lr,lz,nr,nz,rfsmat,scrtch,scrtch2,0)
      time1 = wtimeoff()
c     --- Subtract the analytic solution to get the error
      errmax = 1.0e-30
      do ir = 0, nr
        do iz = 0, nz
          error = b(ir,iz) - bsav(ir,iz)*(sin(sz*iz)+cos(cz*iz))
          errmax = max(errmax,abs(error))
        enddo
      enddo

      write(STDOUT,26)
26    format("Testing analytic solution -- ")
      write(STDOUT,27)errmax
27    format("Maximum Error = ",1pe12.5)
      write(STDOUT,29)
      write(6,'(" time =",f9.3," milliseconds")')time1
29    format("Calling VPOISRZ with iwhich = 0")
      write(STDOUT,*)" "

c **** Test on random data ****

      do ir = 0,nr-1
        do iz = 0,nz-1
          bsav(ir,iz) = wranf()
          b(ir,iz) = bsav(ir,iz)*eps0
        enddo
      enddo
c     --- Apply the conducting wall boundary condition at r = lr
      do iz = 0,nz-1
          b(nr,iz) = 0.0
          bsav(nr,iz) = b(nr,iz)
      enddo
c     --- Fix the periodic boundary conditions in z
      do ir = 0,nr
          b(ir,nz) = b(ir,0)
          bsav(ir,nz) = b(ir,nz)
      enddo
c     --- Call the solver in two parts--first with iwhich = 1 to
c         set up the kzsq and attz, then with iwhich = -1 to do the
c         solve.
      call wtimeon
      call vpoisrz (1,b,b,kzsq,schrg,eta,phikold,taurc,attz,
     &              filt,0.,0.,
     &              lr,lz,nr,nz,rfsmat,scrtch,scrtch2,0)
      time1 = wtimeoff()
      call vpoisrz (-1,b,b,kzsq,schrg,eta,phikold,taurc,attz,
     &              filt,0.,0.,
     &              lr,lz,nr,nz,rfsmat,scrtch,scrtch2,0)
      time2 = wtimeoff()

c     --- check using 5 point finite difference form of Poisson's equation
      errmax = 1.0e-30
      dz2 = dz*dz
      do ir = 1, nr-1
         rjph = 0.5*(r(ir+1) + r(ir))
         rjmh = 0.5*(r(ir) + r(ir-1))
         dr2j = rjph**2 - rjmh**2
        do iz = 1, nz-1
           delsq =(-2.*rjph/(dr2j*dr)*(b(ir+1,iz) - b(ir,iz))
     &               +2.*rjmh/(dr2j*dr)*(b(ir,iz) - b(ir-1,iz)))
     &             - 1./dz2*(b(ir,iz+1)-2.*b(ir,iz)+b(ir,iz-1))
           error = delsq - bsav(ir,iz)
           errmax = max(errmax,abs(error))
        enddo
      enddo
      do 401 iz = 1,nz-1
         dr20 = (0.5*r(1))**2
         delsq = (-r(1)/(dr20*dr)*(b(1,iz) - b(0,iz)))
     &           - 1./dz2*(b(0,iz+1)-2.*b(0,iz)+b(0,iz-1))
         error = delsq - bsav(0,iz)
         errmax = max(errmax,abs(error))
  401 continue

      write(STDOUT,399)
 399  format("Test on random data -- ")
      write(STDOUT,411)errmax
  411 format("Error checked by 5 point Poisson solve ",1pe12.5)
      write(STDOUT,410)time1
  410 format("Set up time (iwhich = 1) = ",1pe12.5," milliseconds")
      write(STDOUT,412)time2
 412  format("Solve time (iwhich = -1) = ",1pe12.5," milliseconds")

c **** Test with resistivity ****

c   From Ed Lee's talk, resistivity is ~100 ohms/meter
      eta = 100. 
      dt = 1.0e-09
c   Calculate the analytic solution for rho = 2. for r < a = r_wall/2

      a0 = r(nr/2)
      do iz = 0,nz-1
      do ir = 0,nr/2
          bsav(ir,iz) = 2./(4.*eps0)*(a0**2 - r(ir)**2)
     &                  +2.*a0**2/(2.*eps0)*log((lr)/a0)
          b(ir,iz) = bsav(ir,iz)
      enddo
      do ir = nr/2+1,nr
          bsav(ir,iz) = 2.*a0**2/(eps0*2.)*log((lr)/r(ir))
          b(ir,iz) = bsav(ir,iz)
      enddo
      enddo

c  Set the periodic boundary conditions

      do ir = 0,nr
         bsav(ir,nz) = bsav(ir,0)
         b(ir,nz) = bsav(ir,nz)
      enddo
c  Plug the potential into the difference formula to get the 
c  rho.

      dz2 = dz*dz
      do ir = 1, nr-1
         rjph = 0.5*(r(ir+1) + r(ir))
         rjmh = 0.5*(r(ir) + r(ir-1))
         dr2j = rjph**2 - rjmh**2
      do iz = 1, nz-1
         bsav(ir,iz) =(-2.*rjph/(dr2j*dr)*( b(ir+1,iz) - b(ir,iz))
     &             +2.*rjmh/(dr2j*dr)*(b(ir,iz) - b(ir-1,iz)))
     &           - 1./dz2*(b(ir,iz+1)-2.*b(ir,iz)+b(ir,iz-1))
      enddo
      enddo
      do iz = 1,nz-1
         dr20 = (0.5*r(1))**2
         bsav(0,iz) = (-r(1)/(dr20*dr)*(b(1,iz) - b(0,iz)))
     &           - 1./dz2*(b(0,iz+1)-2.*b(0,iz)+b(0,iz-1))
      enddo
      rjph = 0.5*(2.*r(nr) + dr)
      rjmh = 0.5*(2.*r(nr) - dr)
      dr2j = 2.*r(nr)*dr
      do iz = 1, nz-1
         bsav(nr,iz) =(-2.*rjph/(dr2j*dr)*(- b(nr,iz))
     &             +2.*rjmh/(dr2j*dr)*(b(nr,iz) - b(nr-1,iz)))
     &           - 1./dz2*(b(nr,iz+1)-2.*b(nr,iz)+b(nr,iz-1))
     &           - dt/(pi*dr2j*dz2*eta*eps0)*(b(nr,iz+1) - 
     &             2.*b(nr,iz) + b(nr,iz-1))
      enddo

      do ir = 0,nr
         bsav(ir,0) = bsav(ir,1)
         bsav(ir,nz) = bsav(ir,0)
      enddo

      do ir = 0,nr
      do iz = 0,nz
         temp = bsav(ir,iz)
         bsav(ir,iz) = b(ir,iz)
         b(ir,iz) = temp * eps0
      enddo
      enddo
c   Set the surface charge
      do iz = 0,nz
         schrg(iz) = b(ir,iz) * 2.*pi*lr*dr*dz
         b(nr,iz) = 0.
      enddo

c     --- Call the solver in two parts--first with iwhich = 1 to
c         set up the kzsq and attz, then with iwhich = -1 to do the
c         solve.
      call wtimeon
      call vpoisrz (1,b,b,kzsq,schrg,eta,phikold,taurc,attz,
     &              filt,dt,0.,
     &              lr,lz,nr,nz,rfsmat,scrtch,scrtch2,0)
      time1 = wtimeoff()
      call vpoisrz (-1,b,b,kzsq,schrg,eta,phikold,taurc,attz,
     &              filt,dt,0.,
     &              lr,lz,nr,nz,rfsmat,scrtch,scrtch2,0)
      time2 = wtimeoff()

c  Adjust the potential at the wall
      do ir = 0,nr
      do iz = 0,nz
         bsav(ir,iz) = bsav(ir,iz) + b(nr,iz)
      enddo
      enddo

      errmax = 1.0e-30

      do iz = 0,nz
      do ir = 0,nr
         error = abs((bsav(ir,iz) - b(ir,iz))/b(ir,iz))
         errmax = max(errmax,abs(error))
      enddo
      enddo

      write(STDOUT,*)" "
      write(STDOUT,699)
 699  format("Test on analytic solution with resistive wall-- ")
      write(STDOUT,611)errmax
  611 format("Error checked by 5 point Poisson solve ",1pe12.5)
      write(STDOUT,610)time1
  610 format("Set up time (iwhich = 1) = ",1pe12.5," milliseconds")
      write(STDOUT,612)time2
 612  format("Solve time (iwhich = -1) = ",1pe12.5," milliseconds")

c  Test resistive wall with random data ...
      
      do ir = 0,nr-1
      do iz = 0,nz-1
          bsav(ir,iz) = wranf()
          b(ir,iz) = bsav(ir,iz)
      enddo
      enddo

      do iz = 0,nz-1
          bsav(nr,iz) = 0.
          b(nr,iz) = 0.
      enddo

      do iz = 0,nz-1
         schrg(iz) = 0.
         do ir = 0,nr
            schrg(iz) = schrg(iz) + b(ir,iz)*r(ir)
         enddo
         schrg(iz) = -schrg(iz)*dr*2.*pi*dz
      enddo

      do iz = 0,nz-1
         bsav(nr,iz) = bsav(nr,iz) + 0.5*schrg(iz)/(pi*dr*lr*dz)
      enddo

      do ir = 0,nr
         bsav(ir,nz) = bsav(ir,0)
         b(ir,nz) = b(ir,0)
      enddo
      schrg(nz) = schrg(0)

      call wtimeon
      call vpoisrz (1,b,b,kzsq,schrg,eta,phikold,taurc,attz,
     &              filt,dt,0.,
     &              lr,lz,nr,nz,rfsmat,scrtch,scrtch2,0)
      time1 = wtimeoff()

      call vpoisrz (-1,b,b,kzsq,schrg,eta,phikold,taurc,attz,
     &              filt,dt,0.,
     &              lr,lz,nr,nz,rfsmat,scrtch,scrtch2,0)
      time2 = wtimeoff()

c  Check the solution using a 5 point Poisson formula
      dz2 = dz*dz
      errmax = 1.0e-30
      do ir = 1, nr-1
         rjph = 0.5*(r(ir+1) + r(ir))
         rjmh = 0.5*(r(ir) + r(ir-1))
         dr2j = rjph**2 - rjmh**2
      do iz = 1, nz-1
         delsq =(-2.*rjph/(dr2j*dr)*( b(ir+1,iz) - b(ir,iz))
     &             +2.*rjmh/(dr2j*dr)*(b(ir,iz) - b(ir-1,iz)))
     &           - 1./dz2*(b(ir,iz+1)-2.*b(ir,iz)+b(ir,iz-1))
         delsq = delsq*eps0
         err1(ir,iz) = delsq
         error = delsq - bsav(ir,iz)
         errmax = max(errmax,abs(error))
      enddo
      enddo
      do iz = 1,nz-1
         dr20 = (0.5*r(1))**2
         delsq= (-r(1)/(dr20*dr)*(b(1,iz) - b(0,iz)))
     &           - 1./dz2*(b(0,iz+1)-2.*b(0,iz)+b(0,iz-1))
         delsq = delsq*eps0
         err1(0,iz) = delsq
         error = delsq - bsav(0,iz)
         errmax = max(errmax,abs(error))
      enddo
      rjph = 0.5*(2.*r(nr) + dr)
      rjmh = 0.5*(2.*r(nr) - dr)
      dr2j = 2.*r(nr)*dr
      do iz = 1, nz-1
c        delsq =(-2.*rjph/(dr2j*dr)*(- b(nr,iz))
c        delsq =(
c    &             +2.*rjmh/(dr2j*dr)*(b(nr,iz) - b(nr-1,iz)))
c    &           - 1./dz2*(b(nr,iz+1)-2.*b(nr,iz)+b(nr,iz-1))
c    &           - dt/(pi*dr2j*dz2*eta*eps0)*(b(nr,iz+1) - 
c    &             2.*b(nr,iz) + b(nr,iz-1))
         delsq = 1./(dr*dr)*(1. - 1./(2.*nr))*(b(nr,iz) - b(nr-1,iz))
     &          -1./(dz*dz)*(1.+0.5*dt/(eps0*eta*pi*lr*dr))*
     &           (b(nr,iz+1) - 2.*b(nr,iz) + b(nr,iz-1))
         delsq = delsq*eps0
         err1(nr,iz) = delsq
         error = delsq - bsav(nr,iz)
         errmax = max(errmax,abs(error))
      enddo
          

      write(STDOUT,*)" "
      write(STDOUT,799)
 799  format("Test on random data with resistive wall-- ")
      write(STDOUT,711)errmax
  711 format("Error checked by 5 point Poisson solve ",1pe12.5)
      write(STDOUT,710)time1
  710 format("Set up time (iwhich = 1) = ",1pe12.5," milliseconds")
      write(STDOUT,712)time2
 712  format("Solve time (iwhich = -1) = ",1pe12.5," milliseconds")

      return
      end
c=======================================================================
      subroutine frzfin()
      use FRZvars

c  Deallocates dynamic storage used by test driver


      call gfree ("FRZvars")

      return
      end
c=======================================================================
      subroutine vprzx (iwhich,dt)
      use FRZvars
      use Constant
      integer(ISZ):: iwhich
      real(kind=8):: dt

c  BASIS interface to VPOISRZ, using variables from the database for
c  package FRZ.  The user's application package should contain a
c  similar subroutine for convenience, using its own database variables.
c  To avoid confusion, it should have a different name, say "vprz"


      call vpoisrz (iwhich,b,b,kzsq,schrg,eta,phikold,taurc,attz,filt,
     &              dt,vbeam,lr,lz,nr,nz,rfsmat,scrtch,scrtch2,ibc)

      return
      end
c=======================================================================
      subroutine vpoisrz (iwhich, a, ak, kzsq, schrg, eta, phikold,
     & taurc, attz, filt, dt, vbeam, lr, lz, nr, nz, rfsmat, 
     & scrtch, scrtch2, ibc)
      use Constant
      integer(ISZ):: iwhich,nr,nz,ibc
      real(kind=8):: eta,taurc,dt,vbeam
      real(kind=8)::  a(0:nr, 0:nz)
      real(kind=8):: ak(0:nr,  -nz/2:nz/2)
      real(kind=8):: schrg(0:nz),phikold(0:nz)
      real(kind=8):: kzsq(0:nz), lr, lz
      real(kind=8):: attz(0:nz/2)
      real(kind=8):: filt(5)
      real(kind=8):: rfsmat(0:nr,3,0:nz),scrtch(0:nz),scrtch2(0:nz)

c-----------------------------------------------------------------------------
c               Vectorized RZ Poisson solver:
c               Periodic FFT in z, tridiagonal solve in r
c               Debbie Callahan, LLNL, October 1989.
c
c          The periodic FFT in z for this routine is based on the
c          routines in VPOIS3D by Alex Friedman.  The difference
c          scheme in the radial direction is based on the scheme
c          presented in Birdsall and Langdon, "Plasma Physics via
c          Computer Simulation," pp.333-334.
c
c          Returned value: 0 if ok, 1 if error
c
c          IWHICH = ... -1: full fieldsolve; assumes kzsq and attz
c                           have been set already.  This is equivalent to a
c                           sequence of calls with iwhich = 2,3,4,5.
c                        0: full fieldsolve; kzsq will be set herein.
c                           This is equivalent to a sequence of calls with
c                           iwhich = 1,2,3,4,5. (THIS IS THE SIMPLEST USAGE.)
c                        1: set kzsq etc. and attz etc., return.
c                        2: forward transform in z, return.
c                        3: apply k-space filter using attz, return.
c                        4: solve tridiagonal system in r, return.
c                        5: inverse transform, return.
c                        8: apply reciprocal of filter in k-space, return.
c                        9: apply k^2 factor (turns phi into rho), return.
c              A, AK ... Array to transform (pass it to this routine twice)
c                        See array declarations below these comments.
c               KZSQ ... Discrete analog to dr*dr*kz^2
c                        See array declarations below these comments.
c              SCHRG ... Array containing the surface charge 
c                ETA ... Surface resistivity  (in Ohms/m)
c            PHIKOLD ... Transform of Phi at r_wall from previous step
c                        Used when capacitance is added
c              TAURC ... RC time
c               ATTZ ... Attenuation factor for each mode number in z.
c                        See array declarations below these comments.
c            FILT(i) ... Spatial filtering coefficients; array should be
c                        dimensioned (5) in calling routine.
c                        These refer to filtering in the z direction
c                        and there are four coefficients:
c                        filt(1): coefficient a1, of sin^2 term (boosting)
c                        filt(2): coefficient a2, of tan^4 term (smoothing)
c                        filt(3): cutoff k dx, k dy, or k dz, e.g. pi/2
c                        filt(4): exponent N (sharpness of cutoff) e.g. 8.
c                        filt(5): W factor in binomial filter, e.g. 0.5
c                        (A value of .5 gives the familiar 1/4:1/2:1/4 filter)
c                        See Birdsall and Langdon, "Plasma Physics via
c                        Computer Simulation," Appendices B and C.
c               EPS0 ... Epsilon-0; use 1/4Pi for CGS, 1 for rationalized.
c                 DT ... Time step -- only used if eta .ne. 0
c              VBEAM ... Beam velocity -- used to advect the surface charge
c                        if eta .ne. 0
c            LR,  LZ ... System lengths
c            NR,  NZ ... Number of mesh points (mesh runs from 0 to NR, etc.)
c                        At present these must be powers of two.
c             RFSMAT ... Array of dimension (0:nr,3,0:nz) -- workspace for
c                        storing the radial fieldsolve tridiagonal system
c             SCRTCH ... Workspace used for resistive wall calculations--
c                        not used if eta .eq. 0
c                IBC ... Boundary condition switch, char*8 (for future use)
c
c  For the (periodic) z direction, k dz = 2 pi j/nz, j = 0, ..., nz/2.
c        The longest wavelength (excepting the k=0 mode) is nz dz = lz.
c
c  The periodic FFT in z for this routine is based on the routines in
c  VPOIS3D by Alex Friedman.  The difference scheme in the radial direction
c  is based on the scheme presented in Birdsall and Langdon,
c  "Plasma Physics via Computer Simulation," pp.333-334.
c
c  This routine was inspired by Langdon's scalar sine-periodic solver "sppois"
c-----------------------------------------------------------------------------


c  Use of both a and ak simplifies subscripting of the Fourier coefficients:

      integer(ISZ):: ikz,ikr,iz,ir,nrlast,nz2,nr2
      real(kind=8):: klast,kdzb2,kdz
      real(kind=8):: dr,dz,norm,dr2
c  Error checking; so far, only possible error is non-power-of-2 nz

      if  ( nz .ne. 2**int( .5 + log(real(nz)) / log(2.) ) ) then
         return
      endif

c  Set useful constants

      nz2 = nz / 2
      nr2 = nr / 2
      dr = lr / nr
      dr2 = dr*dr
      dz = lz / nz

c ----------------------------------------------------------------------------
c  For poisson equation, kzsq(ikz) is a discrete
c  analog to dr*dr*kz^2, etc.  If the user has requested it,
c  we set up ksq arrays for a seven point scheme now; the coding for kzsq
c  is arranged so that it will vectorize on "most" compilers.
c  Also, compute attenuation factors as functions of mode numbers in z.
c ----------------------------------------------------------------------------

      if (iwhich .eq. 0 .or. iwhich .eq. 1) then

c        --- compute z direction coefficients---
         do ikz = 0, nz
            kzsq(ikz) = ikz
         enddo
         do ikz = 0, nz2
            kzsq(ikz) = dr2*(2./dz*sin(pi*(nz2-kzsq(ikz))/nz))**2
         enddo
         do ikz = (nz2+1),nz
            kzsq(ikz) = dr2*(2./dz*sin(pi*(kzsq(ikz)-nz2)/nz))**2
         enddo


c        --- spatial filtering in z; binomial filter, or unity
         do ikz = 0, nz2
            kdz = 2. * pi * ikz / nz
            attz(ikz) = (1.+2.*filt(5)*cos(kdz)) / (1.+2.*filt(5))
         enddo
         if (filt(1) .ne. 0. .or. filt(2) .ne. 0.) then
c           --- compute first form of spatial filtering
            do ikz = 0, nz2
               kdzb2 = pi * ikz / nz
               attz(ikz) = attz(ikz) * (exp( filt(1)*(sin(kdzb2))**2
     &          - filt(2)*(tan(kdzb2))**4 ))**2
         enddo
         endif
         if (filt(3) .ne. 0.) then
c           --- compute second form of spatial filtering
            klast = filt(3) * nz / (2. * pi)
            do ikz = 0, nz2
               attz(ikz) = attz(ikz) * exp(-(ikz/klast)**filt(4))
           enddo
         endif

      endif

c ----------------------------------------------------------------------------
c  Do the forward transform
c ----------------------------------------------------------------------------

      if (iwhich.le.0 .or. iwhich.eq.2) then

      norm = .5 * dz
c        --- Normalize the array
      do iz = 0, nz-1
        do ir = 0, nr
          a(ir,iz) = a(ir,iz) * norm
        enddo
      enddo
      if (eta .ne. 0.) then
c        --- Normalize the surface charge too.
         do iz = 0,nz-1
            schrg(iz) = schrg(iz) * norm
            scrtch(iz) = 0.
         enddo
      endif

c  Here we perform the periodic transforms in z, two r's at a time.

c        --- shift k-space origin to middle.
         do iz = 1, nz-1, 2
           do ir = 0, nr
             a(ir,iz) = -a(ir,iz)
           enddo
         enddo
         if (eta .ne. 0.) then
            do iz = 1, nz-1, 2
               schrg(iz) = - schrg(iz)
            enddo
         endif
c        --- do the transforms
         call vcpft (a(0,0), a(1,0), nz, (nr+1),-1,nr2,2)
         call vrpft2(a(0,0), a(1,0), nz, (nr+1),   nr2,2)

         if (eta .ne. 0.) then

c        --- First, tranform the surface charge
            do iz = 0,nz
               scrtch(iz) = schrg(iz)
            enddo
            
c        --- do the transform of the surface charge 
            call vcpft (schrg(0),scrtch(0),nz,1,-1,1,1)
            call vrpft2(schrg(0),scrtch(0),nz,1,   1,1)

c        --- then transform the nr'th row of the a matrix
            do iz= 0,nz
               scrtch(iz) = a(nr,iz)
               scrtch2(iz) = a(nr,iz)
            enddo

c        --- do the transform of the nr'th row of a
            call vcpft (scrtch(0),scrtch2(0),nz,1,-1,1,1)
            call vrpft2(scrtch(0),scrtch2(0),nz,1,   1,1)
            do iz = 0,nz
                a(nr,iz) = scrtch(iz)
            enddo

         endif


      endif

c ----------------------------------------------------------------------------
c  Apply k-space filter; attenuate mode-by-mode
c ----------------------------------------------------------------------------

      if (iwhich.le.0 .or. iwhich.eq.3) then

         do ikz = 1, nz2-1
           do ikr = 1, nr
             ak(ikr,-ikz) = ak(ikr,-ikz)*attz(ikz)
             ak(ikr, ikz) = ak(ikr, ikz)*attz(ikz)
           enddo
         enddo

      endif

c ----------------------------------------------------------------------------
c  Convert rhok to phik. This is done by solving a tri-diagonal system
c  for each value of kzsq
c ----------------------------------------------------------------------------

      if (iwhich.le.0 .or. iwhich.eq.4) then

c  Fill the tridiagonal matrix, stored in band storage mode.  For each
c  inversion, we need to subtract kzsq from the diagonal (rfsmat(ikr,2))
c  term.  The first and last rows of the matrix are special cases.
c  The solve is done so that it is vectorized in z.

c  Fill the matrix...
c   If we have a conducting wall, we have eqns 0 - (nr-1)  (phi(nr) = 0)
c   but if we have resistive walls then we have 0 - (nr) eqns
c   (phi(nr+1)=0), so set nrlast.
      nrlast = nr-1
      if (eta .ne. 0.) nrlast = nr

      do ikz = 0,nz
         rfsmat(0,1,ikz) =  0.
         rfsmat(0,2,ikz) =  4. + kzsq(ikz)
         rfsmat(0,3,ikz) = -4.
      enddo
      do ikr = 1, nrlast-1
        do ikz = 0,nz
            rfsmat(ikr,1,ikz) = -1. + 1./(2.*ikr)
            rfsmat(ikr,2,ikz) =  2. + kzsq(ikz)
            rfsmat(ikr,3,ikz) = -1. - 1./(2.*ikr)
        enddo
      enddo
      if (eta .eq. 0.) then
         do ikz = 0,nz
            rfsmat(nrlast,1,ikz) = -1. + 1./(2.*nrlast)
            rfsmat(nrlast,2,ikz) =  2. + kzsq(ikz)
            rfsmat(nrlast,3,ikz) =  0.
         enddo
      else
         do ikz = 0,nz
c           rfsmat(nrlast,1,ikz) = -1. + 1./(2.*nrlast)
c           rfsmat(nrlast,2,ikz) =  2. + kzsq(ikz)*
c    &                           (1. + 0.5*dt/(eps0*eta*pi*lr*dr)*
c    &                            (1. + taurc/dt))
            rfsmat(nrlast,1,ikz) = -1. + 1./(2.*nrlast)
            rfsmat(nrlast,2,ikz) =  1. - 1./(2.*nrlast) +
     &                       kzsq(ikz)*(1.+0.5*dt/(eps0*eta*pi*lr*dr)
     &                           * (1. + taurc/dt))
            rfsmat(nrlast,3,ikz) = 0.
         enddo
      endif

c   Add the surface charge to the right hand side of the 
c   equation
      
      if (eta .ne. 0.) then
         do ikz = 0,nz
            a(nr,ikz) = a(nr,ikz) + 0.5*schrg(ikz)/(pi*lr*dr*dz)
     &                  + kzsq(ikz)/(dr*dr) * taurc * phikold(ikz)
     &                    / (eta*2.*pi*lr*dr)
         enddo
      endif



c   Multiply rho by dr*dr and divide by epsilon_0...

      do ikr=0,nrlast
        do ikz = 0,nz
           a(ikr,ikz) = a(ikr,ikz)*dr2/eps0
        enddo
      enddo

c   Eliminate the lower diagonal from the matrix... 
      do ikr=1,nrlast
        do ikz = 0,nz
          rfsmat(ikr,2,ikz) = rfsmat(ikr,2,ikz)
     &                        - rfsmat(ikr-1,3,ikz)*rfsmat(ikr,1,ikz
     &                        )/rfsmat(ikr-1,2,ikz)
          a(ikr,ikz) = a(ikr,ikz)
     &               - a(ikr-1,ikz)*rfsmat(ikr,1,ikz)/rfsmat(ikr-1,2,ikz)
          rfsmat(ikr,1,ikz) = 0.
        enddo
      enddo


c   Divide each row by the diagonal...
c   If we're using resistivity, take care of the last radial cell
c   separately
      if (eta .eq. 0.) then
      do ikr=0,nrlast
        do ikz = 0,nz
          rfsmat(ikr,3,ikz) = rfsmat(ikr,3,ikz)/rfsmat(ikr,2,ikz)
          a(ikr,ikz) = a(ikr,ikz)/rfsmat(ikr,2,ikz)
          rfsmat(ikr,2,ikz) = 1.
        enddo
      enddo
      else
      do ikr=0,nrlast-1
        do ikz = 0,nz
          rfsmat(ikr,3,ikz) = rfsmat(ikr,3,ikz)/rfsmat(ikr,2,ikz)
          a(ikr,ikz) = a(ikr,ikz)/rfsmat(ikr,2,ikz)
          rfsmat(ikr,2,ikz) = 1.
        enddo
      enddo
      do ikz = 0,nz/2 -1
         rfsmat(nrlast,3,ikz) = rfsmat(nrlast,3,ikz)/rfsmat(nrlast,2,ikz)
         a(nrlast,ikz) = a(nrlast,ikz)/rfsmat(nrlast,2,ikz)
         rfsmat(nrlast,2,ikz) = 1.
      enddo
      do ikz = nz/2 +1,nz
         rfsmat(nrlast,3,ikz) = rfsmat(nrlast,3,ikz)/rfsmat(nrlast,2,ikz)
         a(nrlast,ikz) = a(nrlast,ikz)/rfsmat(nrlast,2,ikz)
         rfsmat(nrlast,2,ikz) = 1.
      enddo

c   If we're using resistivity, get rid of the large constant part of the
c   potential at the wall
      
      if (eta .ne. 0.) a(nrlast,nz/2) = 0.

      endif

c   Back substitute to get the solution
      do ikr=(nrlast-1),0,-1
        do ikz = 0,nz
           a(ikr,ikz) = a(ikr,ikz) - a(ikr+1,ikz)*rfsmat(ikr,3,ikz)
           rfsmat(ikr,3,ikz) = 0.
        enddo
      enddo

c   Set boundary condition -- phi = 0 at r = r_wall if eta = 0.
c   Otherwise we've used d phi/dr (r_wall+dr) = 0.

      if (eta .eq. 0.) then
        do ikz = 0,nz
          a(nr,ikz) = 0.
        enddo
      endif


c   Advance the surface charge to it's value at the new time
      if (eta .ne. 0.) then
      call advsc(schrg,eta,a,kzsq,dt,dr,dz,nr,nz,phikold,taurc)
      endif

c *** This is now done in advect ***
c   Save the transform of phi at the wall for the next time step
c     if (taurc .ne. 0.) then
c     do iz = 0,nz
c         phikold(iz) = a(nr,iz)
c     enddo
c     endif

      endif

c ----------------------------------------------------------------------------
c  Do the inverse transform
c ----------------------------------------------------------------------------

      if (iwhich.le.0 .or. iwhich.eq.5) then

c  Inverse transform and shift in z.

c        --- do the transforms
         call vrpfti2(a(0,0), a(1,0), nz, (nr+1),   nr2,2)
         call vcpft  (a(0,0), a(1,0), nz, (nr+1),+1,nr2,2)

         if (eta .ne. 0.) then
            do ikz = 0,nz
               scrtch(ikz) = schrg(ikz)
            enddo
c        --- do the transform of the surface charge
            call vrpfti2(schrg(0),scrtch(0),nz,1,   1,1)
            call vcpft  (schrg(0),scrtch(0),nz,1,+1,1,1)

            do ikz = 0,nz
               scrtch(ikz) = a(nr,ikz)
               scrtch2(ikz) = a(nr,ikz)
            enddo

c        --- do the transform of the nr'th row of a
            call vrpfti2(scrtch(0),scrtch2(0),nz,1,   1,1)
            call vcpft  (scrtch(0),scrtch2(0),nz,1,+1,1,1)
            do iz = 0,nz
               a(nr,iz) = scrtch(iz)
            enddo

         endif

c        --- shift k-space origin back.
         do iz = 1, nz-1, 2
           do ir = 0, nr
            a(ir,iz) = -a(ir,iz)
           enddo
         enddo
         if (eta .ne. 0.) then
            do iz = 1, nz-1, 2
               schrg(iz) = -schrg(iz)
            enddo
         endif

c  Re-normalize the array

      norm = 1. / lz
         do ir = 0, nr
           do iz = 0, nz-1
               a(ir,iz) = a(ir,iz) * norm
           enddo
         enddo
         if (eta .ne. 0.) then
            do iz = 0, nz-1
               schrg(iz) = schrg(iz) * norm
            enddo
         endif

c  Enforce periodic bc's in z

      do ir = 0, nr
         a(ir,nz) = a(ir,0)
      enddo

      schrg(nz) = schrg(0)
c  Advect the surface charge with the moving window.
      if (eta .ne. 0.) then
         norm = .5 * dz
         call advect(schrg,a,phikold,vbeam,taurc,dt,nr,nz,dz,scrtch)

c        --- Transform the old potential at the wall here since it's
c            nice to have this available in the interpretor

         do iz = 0,nz-1
            phikold(iz) = phikold(iz) * norm
            scrtch(iz) = 0.
         enddo
         do iz = 1, nz-1, 2
            phikold(iz) = - phikold(iz)
         enddo
         call vcpft (phikold(0),scrtch(0),nz,1,-1,1,1)
         call vrpft2(phikold(0),scrtch(0),nz,1,   1,1)
      endif

      endif

c ----------------------------------------------------------------------------
c  The usual fieldsolve has been completed.  Utility operations follow.
c ----------------------------------------------------------------------------

c ----------------------------------------------------------------------------
c  Un-filter in k space
c ----------------------------------------------------------------------------

c     if (iwhich.eq.8) then

c     do 1740 iky = 1, ny-1
c        do 1720 ikx = 1, nx-1
c           ak(ikx,iky,0   ) = ak(ikx,iky,0   )
c    &       / (attx(ikx) * atty(iky) * attz(0))
c           ak(ikx,iky,-nz2) = ak(ikx,iky,-nz2)
c    &       / (attx(ikx) * atty(iky) * attz(nz2))
c1720    continue
c        do 1740 ikz = 1, nz2-1
c        do 1740 ikx = 1, nx-1
c           ak(ikx,iky,-ikz) = ak(ikx,iky,-ikz)
c    &       / (attx(ikx) * atty(iky) * attz(ikz))
c           ak(ikx,iky, ikz) = ak(ikx,iky, ikz)
c    &       / (attx(ikx) * atty(iky) * attz(ikz))
c1740 continue

c     endif

c ----------------------------------------------------------------------------
c  Multiply by k^2 factor, in contrast to the usual division.
c  This is useful if one wants to recreate rho from phi.
c ----------------------------------------------------------------------------

c     if (iwhich.eq.9) then

c     do 1850 iky = 1, ny-1
c        do 1800 ikx = 1, nx-1
c           ak(ikx,iky,0   ) = ak(ikx,iky,0   )
c    &                            * ( kxsq(ikx) + kysq(iky) )
c           ak(ikx,iky,-nz2) = ak(ikx,iky,-nz2)
c    &                            * ( kxsq(ikx) + kysq(iky) + kzsq(nz2) )
c1800    continue
c        do 1850 ikz = 1, nz2-1
c        do 1850 ikx = 1, nx-1
c           ak(ikx,iky,-ikz) = ak(ikx,iky,-ikz)
c    &                         * ( kxsq(ikx) + kysq(iky) + kzsq(ikz) )
c           ak(ikx,iky, ikz) = ak(ikx,iky, ikz)
c    &                         * ( kxsq(ikx) + kysq(iky) + kzsq(ikz) )
c1850 continue

c     endif

c ----------------------------------------------------------------------------
c  End of VPOIS3D
c ----------------------------------------------------------------------------

      return
      end
c======================================================================
      subroutine tridi(amat,rhs,n,ibn,dr2,nr,nz,ierr)
      use Constant
      integer(ISZ):: n,ibn,nz,nr,ierr
      real(kind=8):: dr2
      real(kind=8):: amat(0:nr,3),rhs(0:nr,0:nz)
      integer(ISZ):: i

      do i=0,n
         rhs(i,ibn) = rhs(i,ibn)*dr2/eps0
      enddo

      do i=1,n
         amat(i,2) = amat(i,2) - amat(i-1,3)*amat(i,1)/amat(i-1,2)
         rhs(i,ibn) = rhs(i,ibn) - rhs(i-1,ibn)*amat(i,1)/amat(i-1,2)
         amat(i,1) = 0.
      enddo

      do i=0,n
         amat(i,3) = amat(i,3)/amat(i,2)
         rhs(i,ibn) = rhs(i,ibn)/amat(i,2)
         amat(i,2) = 1.
      enddo

      do i=n-1,0,-1
         rhs(i,ibn) = rhs(i,ibn) - rhs(i+1,ibn)*amat(i,3)
         amat(i,3) = 0.
      enddo

      return
      end

c======================================================================
      subroutine advsc(schrg,eta,a,kzsq,dt,dr,dz,nr,nz,phikold,taurc)
      integer(ISZ):: nr,nz
      real(kind=8):: eta,dt,dr,dz,taurc
      real(kind=8):: schrg(0:nz), kzsq(0:nz), a(0:nr,0:nz)
      real(kind=8):: phikold(0:nz)
c  Advance the surface charge using the Fourier transform of 
c  d Q/dt = (dz/eta) * (d^2 Phi/dz^2 - eta*C d/dt( d^2 Phi/dz^2))
      integer(ISZ):: iz

      do iz = 0,nz
         schrg(iz) = schrg(iz) - (dt*dz)/eta * kzsq(iz)/(dr*dr) * 
     &               ((1.+taurc/dt)*a(nr,iz)-taurc/dt*phikold(iz))
      enddo

      return
      end

c======================================================================
      subroutine advect(schrg,a,phikold,vbeam,taurc,dt,nr,nz,dz,scrtch)
      integer(ISZ):: nr,nz
      real(kind=8):: vbeam,taurc,dt,dz
      real(kind=8):: schrg(0:nz),scrtch(0:nz),a(0:nr,0:nz),phikold(0:nz)
c  Maps the surface charge along as the mesh moves
      integer(ISZ):: nl,iz
      
      nl = int(vbeam*dt/dz)

      do iz = 0, (nz-nl)
         scrtch(iz) = schrg(iz + nl)
      enddo

      do iz = (nz-nl+1),nz
         scrtch(iz) = schrg(iz-nz+nl-1)
      enddo

      do iz = 0,nz
         schrg(iz) = scrtch(iz)
      enddo
c  Also map phi at the wall for the next time step if using capacitance
      do iz = 0, (nz-nl)
         scrtch(iz) = a(nr,iz + nl)
      enddo

      do iz = (nz-nl+1),nz
         scrtch(iz) = a(nr,iz-nz+nl-1)
      enddo

      do iz = 0,nz
         phikold(iz) = scrtch(iz)
      enddo

      return
      end

