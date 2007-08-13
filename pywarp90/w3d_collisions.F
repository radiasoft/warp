#include "top.h"
c=============================================================================
c@(#) File w3d_collision.F, version $Revision: 1.1 $, $Date: 2007/08/13 22:05:08 $
c# Copyright (c) 2007-2007, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This has routines handling particle collisions
c  Alex Friedman, LLNL, (925)422-0827
c  David P. Grote, LLNL, (925)423-7194
c  Bruce Cohen, LLNL, (925)422-9823
c=============================================================================
      subroutine langevincollisions2d(lself,np,ux1,uy1,uz1,
     &                                vxmi,vymi,vzmi,density2,vthsqi2,
     &                                q1,q2,m1,m2,vthsqinit1,vthsqinit2,
     &                                dt,loglambda,epvth)
      use Constant
      logical(ISZ):: lself
      integer(ISZ):: np
      real(kind=8):: ux1(np),uy1(np),uz1(np)
      real(kind=8):: vxmi(np),vymi(np),vzmi(np),density2(np),vthsqi2(np)
      real(kind=8):: q1,q2,m1,m2,vthsqinit1,vthsqinit2
      real(kind=8):: dt,loglambda,epvth

cc CIC collision operator, 2 space, 3 velocity dimensions
cc Implements Manheimer, Lampe, Joyce, JCP 138, 563 (1997)
cc test particle species i1 scatters of field particle species i2

cc nudt input to this routine corresponds to:
cc     nudt = dt * (2^1.5 /3pi^.5) nu_NRL (i1,i1) for species i1
cc where  nu_NRL = (4pi Z_1^4 e^4 n_1 ln Lambda)/(m_1^2 v_th_1^3)

cc For collisions of i1 on i2, this routine internally converts
cc nudt to nudt = dt * (2^1.5 /3pi^.5) nu_NRL (i1,i2) for species i1 on i2
cc where  nu_NRL(1,2) = (4pi Z_1^2 Z_2^2 e^4 n_2 ln Lambda)/(m_1^2 v_th_2^3)

cc This routine is to be called every ncint timesteps.
c ncint = (1/nudt)*(1/ncoll) where nudt=nu_i *dt and ncoll=no. of
c of times we wish to apply collision operator over the period 1/nudt
c nudt and ncoll are input by user.

c Note:
c v_thermal near edges of plasma is prone to big errors that
c can cause a runaway of particle energy.--DON'T NEED TO WORRY THIS.

      integer(ISZ):: i
      real(kind=8):: sqrt2,twopi
      real(kind=8):: nudtcoeff,nudt,eta,vthsq2
      real(kind=8):: etotf1,etotf2,escale
      real(kind=8):: uxf,uyf,uzf,uxf0
      real(kind=8):: vv,vv1,vv2,fdv,d11,d33
      real(kind=8):: xrand1
      real(kind=8):: wranf

c     --- First, calculate the part of the collision frequency that
c     --- depends on constants or global quantities.
      nudtcoeff = dt*2**1.5/(3.*sqrt(pi))*
     &            4.*pi*q1**2*q2**2*loglambda/m1**2


c     --- Weight factor to the thermal velocity incorporated to deal with
c     --- a cooling instability.
c     --- See Cohen et al PoP 13, 22705 (2006) after eq.(4).
      if(epvth == 0.) epvth=0.95

c     --- Useful constants
      sqrt2 = sqrt(2.)
      twopi = 2.*pi

c     --- Mass scaling factor. Note that for self collisions, eta=1
      eta = 0.5*(1.+m1/m2)

c     --- Initialize the total energy. Note that these are only used
c     --- to scale the final velocity for energy conservation for
c     --- self collisions.
      etotf1 = 0.
      etotf2 = 0.

c     --- Now loop over the test particles.
      do i=1,np

        vthsq2 = epvth*vthsqi2(i) + (1. - epvth)*vthsqinit2
        nudt = nudtcoeff*density2(i)/vthsq2**1.5
        nudt = abs(nudt)

        if(nudt >= 0.99) then
c         --- This should print some kind of warning perhaps?
          write(*,*) "nudt=",nudt
        endif

c Compute particle velocities in local drift frame

        uxf = ux1(i) - vxmi(i)
        uyf = uy1(i) - vymi(i)
        uzf = uz1(i) - vzmi(i)

c compute the temperature before the collision

        etotf1 = (uxf**2 + uyf**2 + uzf**2)/3. + etotf1

c Compute 3D speed and some other speed normalizations for relative
c angles in coordinate frame oriented parallel to local velocity

        vv  =  sqrt(uxf**2+uyf**2+uzf**2) + 1.e-06*sqrt(vthsqinit1)
        vv1 =  sqrt(uxf**2+uyf**2)+ 1.e-06*sqrt(vthsqinit1)
        vv2 =  1./(vv*vv1 +1.e-06*vthsqinit1)

c calculate diffusion and drag coefficients from Pade approximations
c for isotropic F(v) and resulting kicks

        fdv = abs( eta*nudt/(1.+1.08*(vv/(sqrt2*sqrt(vthsqi2(i))))**2.76) )

        d33 = sqrt(abs(2.*vthsqi2(i)*(fdv/eta)*abs(-log(1.-wranf()))))
     &      *cos(twopi*wranf())

        d11 = sqrt(nudt*2.*vthsqi2(i)/(1.+0.1856*
     &      (vv/(sqrt2*sqrt(vthsqi2(i))))**1.848)*abs(-log(1.-wranf())))

c need to multiply d11 by cos and sin of ranf*twopi for kicks

        fdv = -vv*fdv*(1.+0.5*fdv)      !! Includes energy cons. correction ~dt**2

c Do the Monte Carlo collisional kicks in the particle velocity local frame
c and transform back to x-y-z velocity coordinates in drift frame

        xrand1 = wranf()

c Later on we need a duplicate copy of uxf (after uxf is
c overwritten with a new value after the collision) for computing a
c cosine = uxf/vv1.

        uxf0 = uxf

c reduce number of divides by replacing vv with 1./vv

        vv = 1./vv

        uxf = uxf + (fdv+d33)*uxf*vv + d11*cos(twopi*xrand1)*uyf/vv1 
     &      + d11*sin(twopi*xrand1)*uzf*uxf*vv2

        uyf = uyf + (fdv+d33)*uyf*vv + d11*cos(twopi*xrand1)*(-uxf0)/vv1 
     &      + d11*sin(twopi*xrand1)*uzf*uyf*vv2

        uzf = uzf + (fdv+d33)*uzf*vv + d11*sin(twopi*xrand1)*(-vv1**2)*vv2

c Galilean tranformation to lab velocity coordinates from drift frame

        ux1(i) = uxf
        uy1(i) = uyf
        uz1(i) = uzf

c compute temperature after the collision
        etotf2 = (uxf**2+uyf**2+uzf**2)/3. + etotf2

      enddo   !!loop over particle blocks

      escale=sqrt(etotf1/etotf2)

      do i = 1,np

        uxf = ux1(i)
        uyf = uy1(i)
        uzf = uz1(i)

c Scale the particle velocities to conserve energy after the
c collisions.  Manheimer et al. suggest doing this.
c This collision scheme does NOT conserve energy like Jones
c algorithm does, except statistically for an ideal Maxwellian
c distribution of test particles, and erf function drag-diff coeffs,
c none of which we satisfy.

        if(lself) then
          uxf = uxf*escale
          uyf = uyf*escale
          uzf = uzf*escale
        endif

c Galilean tranformation to lab velocity coordinates from drift frame

        ux1(i) = uxf + vxmi(i)
        uy1(i) = uyf + vymi(i)
        uz1(i) = uzf + vzmi(i)

      enddo

      return
      end  !!langevincollisions2d
c=============================================================================
