#include "top.h"
c=============================================================================
c@(#) File w3d_load.F, version $Revision: 1.8 $, $Date: 2004/12/10 23:18:37 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is a source code file of package W3D of the WARP 
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c=============================================================================
c=============================================================================
      subroutine perp_cfe_den() 

      use Constant
      use Beam_acc
      use InPart3d
      use W3Dload

c Calculate the density profile of a continuous focused f(H) equilibrium 
c distribution that will be used to load an rms equivalent beam.  
c
c The routine works with the following distributions:
c   distrbtn = WB  or Waterbag
c            = PA  or Parabolic              
c            = TE  or ThermalEquilibrium 
c
c The axisymmetric equivalent beam density profile is returned via a grid 
c defined by:
c
c   hl_nrdist = number of radial data points of CF Hamiltonian eq 
c                 ranging from [0:hl_nrdist] -> [0,hl_rmax].   
c                 This variable is user set.  
c   hl_rmax   = Code set: max radius that the density profile is 
c                 defined by (must contain all non-negligible density) 
c   hl_rdist(0:hl_nrdist) = Code set: particle density in r 
c
c Auxillary variables used in loading the density profile are also set:
c
c   hl_nrmrdist(0:hl_nrdist) = Code set: normalized hl_rdist with 
c                                sum(hl_nrmrdist)=1. 
c   hl_intrdist(0:hl_nrdist) = Code set: integral of r*rdist from 
c                                r=0 to r where r is a scaled radius.

c     --- beam rms descriptions 
      real(kind=8):: r_x, r_y, emit_x, emit_y
      integer(ISZ):: ir
c     --- waterbag equilibrium 
      real(kind=8):: sc_param, k1re1, k1re2, rho
      real(kind=8):: zbrent, bessi0, bessi
      external wbagconstr
      common /wbagconstrcom/ sc_param
c     --- distribution loading 
      real(kind=8):: hl_sumrdist

   
c     Define local rms equivalent beam parameters to use in the subroutine

      r_x = a0 
      r_y = b0 
      emit_x = emitx
      emit_y = emity 
      q_perv = zion*echarge*ibeam/
     &         dvnz(2.*Pi*eps0*aion*amu*(vbeam*gammabar)**3)

c     Continuous focusing, matched rms equivlaent beam parameters 

      r_b = sqrt( r_x*r_y )
      emit_b = sqrt( emit_x*emit_y ) 
      k_beta0 = sqrt( q_perv/r_b**2 + emit_b**2/r_b**4 )  

c     Allocate dynamic variables used to define the gridded density 
c     and check that the density grid has been defined.   
      call gchange("W3Dload",0)
      if (hl_nrdist .lt. 2) then
        call remark("routine perp_cf_eq: hl_nrdist >= 2 must be set") 
        call kaboom(0)  
      endif 

c     Calculate gridded density for classes of equilbria 
c     --- waterbag equilibrium 
      if (distrbtn == "WB" .or. distrbtn == "Waterbag") then 

c       --- calculate k1re = k1*r_e using bracketed brent numerical root 
c           find method in zbrent  
c             k1re1 = left  bracket established from small argument exapansion 
c             k1re2 = right bracket set to a point where the bessel function 
c                     part of waterbagconstr() changes sign for sc_param = 0.

        sc_param = (r_b*k_beta0)**2/q_perv
        k1re1 = (64./sc_param)**0.25
        k1re2 = 2. + 4./(sc_param-1.)  

        k1re = zbrent(wbagconstr,k1re1,k1re2,1.e-6,1000)

c       --- calculate properties of the distribution 
c             r_e = edge radius where the density falls to zero 
c             f_2 = distribution normalization coefficient 

        r_e = r_b*( bessi0(k1re)/bessi(2,k1re) - 
     &              4./(k1re**2)*( 2. + k1re*bessi(3,k1re)/bessi(2,k1re) ) 
     &            )**(-0.5) 

        f_2 = aion*amu*eps0*vbeam**2*gammabar**3/
     &        (2.*Pi*(zion*echarge)**2)*(k1re/r_e)**2 

c       ---- calculate radial density profile
        hl_rmax = r_e 
        do ir = 0,hl_nrdist
c         --- radial coordinate in distribution 
          rho = float(ir)/float(hl_nrdist) 
c         --- density 
          hl_rdist(ir) = 4.*Pi*f_2*(k_beta0*r_e/k1re)**2 * 
     &                   ( 1. - bessi0(k1re*rho)/bessi0(k1re) ) 
        enddo 

c     --- parabolic equilibrium 
      elseif(distrbtn == "PA"  .or. distrbtn == "Parabolic") then 

        call remark("routine perp_cf_eq: Parabolic Eq option not implemented") 
        call kaboom(0) 

c     --- thermal equilibrium 
      elseif(distrbtn == "TE"  .or. distrbtn == "ThermalEquilibrium") then 

        call remark("routine perp_cf_eq: Thermal Equilibrium not implemented") 
        call kaboom(0) 

c       --- undefined option trap 
      else 

        call remark("routine perp_cf_eq: distrbtn option not implemented") 
        call kaboom(0) 

      endif 

c     Prepare for loading a specified particle density distribution in r,
c     based on hl_rdist by calculating needed auxillary variables that will 
c     employed in the loading process.    

c     --- Multiply hl_rdist by r.
      do ir=0,hl_nrdist
        hl_nrmrdist(ir) = hl_rdist(ir)*float(ir)/float(hl_nrdist)
      enddo
c     --- Integrate r*hl_rdist so it can be normalized.  
c     --- Note that hl_nrmrdist(0) is ignored since it is always zero.
      hl_sumrdist = 0.
      do ir=1,hl_nrdist-1
        hl_sumrdist = hl_sumrdist + hl_nrmrdist(ir)
      enddo
      hl_sumrdist = hl_sumrdist + 0.5*hl_nrmrdist(hl_nrdist)
c     --- Calculate integral of hl_rdist from 0 to ir
c     --- Integral = 0.5*hl_rdist(0) + sum(hl_rdist(1:ir-1)) + 0.5*hl_rdist(ir)
c     --- Also, normalize rdist.
      hl_intrdist(0) = 0.
      do ir=1,hl_nrdist
        hl_nrmrdist(ir) = hl_nrmrdist(ir)/hl_sumrdist
        hl_intrdist(ir) = hl_intrdist(ir-1) + 
     &                    0.5*( hl_nrmrdist(ir-1) + hl_nrmrdist(ir) )
      enddo

      return 
      end 
c=============================================================================
      real(kind=8) function wbagconstr(k1re)

      real(kind=8):: k1re

c watbagconstr returns a scaled constraint equation assoicated with 
c a matched waterbag equilbrium distribution function in a continuous 
c focusing channel.  
c
c    k1re = k_1 * r_e       k_1 = wavenumber parameter 
c                           r_e = beam edge radius
c

      real(kind=8):: bessi0,bessi
      real(kind=8):: sc_param
      common /wbagconstrcom/ sc_param

      wbagconstr = sc_param - bessi0(k1re)**2/bessi(2,k1re)**2 + 
     &  4./(k1re**2)*( 2.*bessi(2,k1re)/bessi0(k1re) +
     &                 k1re*bessi0(k1re)*bessi(3,k1re)/(bessi(2,k1re))**2
     &               )

      return 
      end 
c=============================================================================
      subroutine loadperpdist(kioff,np,x,y,r,t,xp,yp,rx,ry,rxp,ryp,
     &                        epsx,epsy)
      use Constant
      use InPart3d
      use W3Dload
      integer(ISZ):: kioff, np
      real(kind=8):: x(np), y(np), xp(np), yp(np), r(np), t(np) 
      real(kind=8):: rx(np), ry(np), rxp(np), ryp(np)
      real(kind=8):: epsx(np), epsy(np) 
c
c Load a transverse distribution of particles based on transformation of 
c an "equivalent" Hamiltonian defined beam equilibrium in a continuous 
c focusing channel.  Each particle can have 
c separate envelope coordinates and emittances so the script can be used 
c in 2D or 3D loads.  The continuous focused equilibrium is tranformed 
c consistent with KV like elliptical beam cross-sections and flutter to 
c better match a focusing channel with axially varying focusing forces.  
c These loads are not strict equilibria of such focusing channels but in 
c many situaitons can provide reasonable approximate initial beam 
c distributions.  
c   
c Particle coordinates loaded:
c   np = number of particles to load 
c   x,y    = x, y  coordinates of particle 
c   r,t    = radius, theta polar coordinates of particle 
c   xp,yp  = x',y' coordinates of particle
c
c Each particle can be in different parts of the beam envelope and the 
c beam envelope in the axial slice of the each particle (note the 
c z-coordinates of the particles are unknown within this routine) are 
c defined by the arrays:
c
c   rx,   ry    = x- and y-envelope coordinates of beam: r_x,   r_y    
c   rxp,  ryp   = x- and y-envelope angles      of beam: r_x',  r_y' 
c   emitx,emity = x- and y-emittances           of beam: eps_x, eps_y 
c
c The method can be applied to an infinity of distributions.  Forms 
c implemented include:
c   distrbtn = "WB" or "Waterbag":  waterbag form (theta function) 
c            = "PA" or "Parabolic": parabolic form 
c            = "GA" or "Gaussian":  Gaussian form 
c            = "TE" or "ThermalEquilibrium": Thermal Equilibrium form 
c
      integer(ISZ):: ip, irr
      real(kind=8):: wrandom, wrandomgauss
      real(kind=8):: u_max, ang_mag, ang_phase  
      real(kind=8):: bessi0
 
c     Load axisymmetric particle coordinates x,y based on the equivalent 
c     continuous focused equilibrium profile calculated in 
c     perp_cfe_den() and returned via the hl_rdist(0:hl_nrdist) and 
c     related variables.  This is the same (simplified) procedure as setup in 
c     stptcl3d() for use in loading "arbitrary" radial density profiles 
c     under the semiGaussian load option.  

      do ip = 1,np  

c       --- radius = r in [0,1) and theta = t in [-Pi,Pi) variables 
c       --- used to construct needed probability transforms
        r(ip) = wrandom(xrandom,ip+kioff-1,dig1,1,1)
        t(ip) = Pi*( 2.*wrandom(xrandom,ip+kioff-1,dig2,1,1) - 1. )

c       --- Transform the radial positions to match the specified
c       --- distribution.

c       --- Find what grid cell the particle is in.
        irr = 0
        do while (r(ip) > hl_intrdist(irr+1))
          irr = irr + 1
        enddo
c       --- Do the transformation.  The transformation assumes a
c       --- linearly varying distribution between grid locations.
c       --- The transformation is derived by equating an integral
c       --- over the uniform distribution to an integral over the
c       --- specified distribution.  The particle remains within
c       --- the range [0,1].
        if (hl_nrmrdist(irr+1) == hl_nrmrdist(irr)) then
          r(ip) = ( r(ip) - hl_intrdist(irr) )/hl_nrmrdist(irr)
        else
          r(ip) = ( ( sqrt( hl_nrmrdist(irr)**2 -
     &                      2.*( hl_nrmrdist(irr+1) - hl_nrmrdist(irr) )*
     &                         ( hl_intrdist(irr) - r(ip) )
     &                    ) - hl_nrmrdist(irr) 
     &               )/( hl_nrmrdist(irr+1) - hl_nrmrdist(irr) ) + 
     &               irr 
     &            )/hl_nrdist
        endif

c       --- rescale radius r to physical units and 
c       --- convert radius r and theta t to x and y
        r(ip) = hl_rmax*r(ip)
        x(ip) = r(ip)*cos(t(ip))
        y(ip) = r(ip)*sin(t(ip))

      enddo

c     Load axisymmetric particle coordinates x',y' based on the class of 
c     distribution loaded and the particle x,y coordinates loaded

c     --- waterbag equilibrium 
      if (distrbtn == "WB" .or. distrbtn == "Waterbag") then 

        do ip = 1,np

c         --- calculate maximum value of u = (1/2)*( x'^2 + y'^2 ) 
c         --- acceptable at the particle radius r = r(ip) 
          u_max = 2.*(k_beta0*r_e/k1re)**2*
     &               ( 1. - bessi0(k1re*r(ip)/r_e)/bessi0(k1re) )  

c         --- set particle x' and y' to be uniformly distributed within 
c         --- a unit disk of radius u_max 
          ang_mag   = sqrt(abs(2.*u_max))*
     &                  wrandom(xrandom,ip+kioff-1,dig3,1,1) 
          ang_phase = Pi*( 2.*wrandom(xrandom,ip+kioff-1,dig4,1,1) - 1. ) 
          xp(ip) = ang_mag*cos(ang_phase)  
          yp(ip) = ang_mag*sin(ang_phase) 
          
        enddo  

c     --- parabolic equilibrium 
      elseif(distrbtn == "PA"  .or. distrbtn == "Parabolic") then 

        call remark("routine loadperpdist: distrbtn = PA not implemented") 
        call kaboom(0) 

c     --- thermal equilibrium 
      elseif(distrbtn == "TE"  .or. distrbtn == "ThermalEquilibrium") then 

        call remark("routine loadperpdist: distrbtn = TE not implemented") 
        call kaboom(0)

c       --- undefined option trap 
      else 

        call remark("routine loadperpdist: distrbtn option not implemented") 
        call kaboom(0) 
 
      endif 

c     Transform the axisymmetric particle coordinates for elliptical rms 
c     equivalency.   

      do ip = 1,np 
 
c       --- tranform particle coordinates x,y 
        x(ip) = rx(ip)/r_b * x(ip) 
        y(ip) = ry(ip)/r_b * y(ip) 
 
c       --- tranform particle angles x',y'
        xp(ip) = (r_b/rx(ip))*(epsx(ip)/emit_b)*xp(ip) 
     &            + (rxp(ip)/rx(ip))*x(ip)  
        yp(ip) = (r_b/ry(ip))*(epsy(ip)/emit_b)*yp(ip)
     &            + (ryp(ip)/ry(ip))*y(ip)  

      enddo 

      return 
      end 
c=============================================================================
      subroutine loadperpdist0(kioff,np,x,y,xp,yp,rx,ry,rxp,ryp,
     &                         epsx,epsy)
      use Constant
      use InPart3d
      integer(ISZ):: kioff, np
      real(kind=8):: x(np), y(np), xp(np), yp(np)
      real(kind=8):: rx(np), ry(np), rxp(np), ryp(np)
      real(kind=8):: epsx(np), epsy(np)
c
c Load a transverse distribution of particles based Courant-Snyder invariants 
c in applied fields.  Each particle can have 
c separate envelope coordinates and emittances so the script can be used 
c in 2D or 3D loads.  The method employed here is based in part on 
c T. Kikuchi's implementation of a procedure first developed by:
c   Y.R. Batygin, "Particle Distribution generator in 4D Phase-Space"
c   AIP Conf. Proc. 297 (1994) pp. 419-426  
c The method uses local Courant-Snyder invariant type forms to define 
c distributions that should be direct analogs to the equilibrium forms in 
c continuous focusing channels. These distributions are equilibria in 
c periodic focusing lattices for beams of zero current.  Because the
c equilibria are defined in terms of functions of Courant Snyder invariants 
c of the linear forcusing lattice with axially varying applied fields, the 
c forms automatically hit the target envelope coordinates r_x, r_y, 
c r_x', r_y' and emittances emitx, emity.   These forms are consitent only 
c for linear total forces within the ellptic symmetry beam core.  Consequently 
c the forms only reflect rms equivalent beam compensations and can be far 
c from equilibrium form for strong space-charge depressions where the internal 
c beam space-charge structure is Debye screened.  The case 
c KV0 is correct for finite beam current because the internal forces are linear
c in this case.  This loading method thereby provides an alternative 
c loading method for a KV beam of arbitrary current/space-charge strength to 
c the option "K-V" or "KV" which is implemented differently.  
c   
c Particle coordinates loaded:
c   np = number of particles to load 
c   x,y    = x, y  coordinates of particle 
c   xp,yp  = x',y' coordinates of particle
c
c Each particle can be in different parts of the beam envelope and the 
c beam envelope in the axial slice of the each particle (note the 
c z-coordinates of the particles are unknown within this routine) are 
c defined by the arrays:
c
c   rx,  ry   =   x- and y-envelope coordinates of beam: r_x,   r_y   
c   rxp, ryp  =   x- and y-envelope angles      of beam: r_x',  r_y' 
c   epsx,epsy =   x- and y-emittances           of beam: eps_x, eps_y
c
c The method can be applied to an infinity of distributions.  Forms 
c presently implemented include:
c   distrbtn = "KV0":  KV or Kapchinskij-Vladimirskij (delta function) 
c            = "WB0" or "Waterbag0":  waterbag form (theta function) 
c            = "PA0" or "Parabolic0": parabolic form 
c            = "GA0" or "Gaussian0":  Gaussian form 
c
      integer(ISZ):: ip
      real(kind=8):: alphax, alphay, alphaxp, alphayp
      real(kind=8):: F0, F, G
      real(kind=8):: Ax, Ay, betax, betay 
      real(kind=8):: wrandom, wrandomgauss
      real(kind=8):: alpha3,r
      character(72):: errline  

      
      if (distrbtn == "KV0" .or. 
     &    distrbtn == "WB0" .or. distrbtn == "Waterbag0" .or. 
     &    distrbtn == "PA0" .or. distrbtn == "Parabolic0") then 

c       --- Loop over particles to load 
        do ip = 1,np

c         --- phase ellipse variables (see Batygin papers and Fortran 
c         --- subroutine PDS1 for initializing particles) used in defining 
c         --- local KV like linear force psudo-invariants
          alphax = rx(ip)/sqrt(epsx(ip))
          alphay = ry(ip)/sqrt(epsy(ip))
          alphaxp = rxp(ip)/sqrt(epsx(ip))
          alphayp = ryp(ip)/sqrt(epsy(ip))

c         --- calculate value of F0, F 
          if (distrbtn == "KV0") then
            F0 = epsx(ip)
            F  = F0 

          elseif (distrbtn == "WB0" .or. distrbtn == "Waterbag0") then
            F0 = 1.5*epsx(ip)
            F  = F0*sqrt( wrandom(xrandom,ip+kioff-1,dig1,1,1) ) 

          elseif (distrbtn == "PA0" .or. distrbtn == "Parabolic0") then 

            F0 = 2.0*epsx(ip)
          
            G = wrandom(xrandom,ip+kioff-1,dig1,1,1)

            alpha3 = ( acos( 1.-2.*G ) )/3.
            F = (F0/2.)*( 1.-2.*cos( alpha3 - 2.*Pi/3. ) )
            if (G .eq. 0.) F = 0. 
            if (G .eq. 1.) F = F0
          endif 

c         --- Set particle coordinates 
          Ax = sqrt( F*wrandom(xrandom,ip+kioff-1,dig2,1,1) )
          Ay = sqrt( (F - Ax*Ax)*epsy(ip)/epsx(ip) )

          betax = 2.0*Pi*wrandom(xrandom,ip+kioff-1,dig3,1,1)
          betay = 2.0*Pi*wrandom(xrandom,ip+kioff-1,dig4,1,1)

          x(ip) = Ax * alphax * cos( betax )
          y(ip) = Ay * alphay * cos( betay )

          xp(ip) = Ax * ( alphaxp * cos( betax ) - sin( betax ) / alphax )
          yp(ip) = Ay * ( alphayp * cos( betay ) - sin( betay ) / alphay )

        enddo 
      
c     --- Implement Gaussian as a special case using a Gaussian distributed 
c     --- random number generator with unit sigma.  This alternative 
c     --- implementation is equivalent to loading the Courant-Snyder type 
c     --- psudoinvariants but saves numerical work in avoiding the need to 
c     --- solve an axillary equation numerically to solve for "F".  
      elseif (distrbtn == "GA0" .or. distrbtn == "Gaussian0") then

c       --- Loop over particles to load 
c       --- Because of the way wrandomgauss is written, the code runs somwhat
c       --- more efficiently when the loops are broken up. Note, very
c       --- important, if the loops are not broken up, then the last argument
c       --- to wrandomgauss MUST be false in order to the code to work
c       --- properly.
        do ip = 1,np
          r = wrandomgauss(xrandom,ip+kioff-1,dig1,dig2,1,1,.true.)
          x(ip) = 0.5*rx(ip)*r
        enddo 
        do ip = 1,np
          r = wrandomgauss(xrandom,ip+kioff-1,dig3,dig4,1,1,.true.)
          y(ip) = 0.5*ry(ip)*r
        enddo 

        do ip = 1,np
          r = wrandomgauss(xrandom,ip+kioff-1,dig5,dig6,1,1,.true.)
          xp(ip) = 0.5*(epsx(ip)/rx(ip))*r + (rxp(ip)/rx(ip))*x(ip) 
        enddo 
        do ip = 1,np
c         --- Note that dig7 and dig8 are already used for the longitudinal
c         --- velocity spread.
          r = wrandomgauss(xrandom,ip+kioff-1,dig9,dig10,1,1,.true.)
          yp(ip) = 0.5*(epsy(ip)/ry(ip))*r + (ryp(ip)/ry(ip))*y(ip)
        enddo 

      else 
        write (errline,'("ERROR: distrbtn has an improper value = ",a8)')
     &        distrbtn
        call remark (errline)
        call kaboom (0)
      endif

      return 
      end
c=============================================================================
      real(kind=8) function bessi0(x)
      real(kind=8):: x

c bessi0(x) calculates the Bessel function I_0(x) for real x using 
c approximate polynomial coefficients.  Because I_0 is smooth, this 
c provides a good, efficient approximation using only a few terms and
c together with recusion relations can be used to calculate higher order 
c I_n(x) with integer n > 0.    
c
c The routine is slightly modified from the one presented in:
c    Numerical Recipes The Art of Scientific Computing, 1st Edition 
c    Press, Flannery, Teukolsky, and Vetterling (1986)
 
      real(kind=8):: P1 =  1.0e0
      real(kind=8):: P2 =  3.5156229e0
      real(kind=8):: P3 =  3.0899424e0
      real(kind=8):: P4 =  1.2067492e0
      real(kind=8):: P5 =  0.2659732e0
      real(kind=8):: P6 =  0.360768e-1
      real(kind=8):: P7 =  0.45813e-2
      real(kind=8):: Q1 =  0.39894228e0
      real(kind=8):: Q2 =  0.1328592e-1
      real(kind=8):: Q3 =  0.225319e-2
      real(kind=8):: Q4 = -0.157565e-2 
      real(kind=8):: Q5 =  0.916281e-2
      real(kind=8):: Q6 = -0.2057706e-1
      real(kind=8):: Q7 =  0.2635537e-1
      real(kind=8):: Q8 = -0.1647633e-1
      real(kind=8):: Q9 =  0.392377e-2
      
      real(kind=8):: Y, AX

      if ( abs(x) .lt. 3.75) then
        Y = (x/3.75)**2
        bessi0 = P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      else
        AX=abs(x)
        Y=3.75/AX
        bessi0 = (exp(AX)/sqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     &      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      endif
 
      return
      end
c=============================================================================
      real(kind=8) function bessi1(x)
      real(kind=8):: x

c bessi1(x) calculates the Bessel function I_1(x) for real x using 
c approximate polynomial coefficients.  Because I_1 is smooth, this 
c provides a good, efficient approximation using only a few terms and
c together with recusion relations can be used to calculate higher order 
c I_n(x) with integer n > 0.    
c
c The routine is slightly modified from the one presented in:
c    Numerical Recipes The Art of Scientific Computing, 1st Edition 
c    Press, Flannery, Teukolsky, and Vetterling (1986)

      real(kind=8):: P1 =  0.5e0
      real(kind=8):: P2 =  0.87890594e0
      real(kind=8):: P3 =  0.51498869e0
      real(kind=8):: P4 =  0.15084934e0
      real(kind=8):: P5 =  0.2658733e-1
      real(kind=8):: P6 =  0.301532e-2
      real(kind=8):: P7 =  0.32411e-3
      real(kind=8):: Q1 =  0.39894228e0 
      real(kind=8):: Q2 = -0.3988024e-1 
      real(kind=8):: Q3 = -0.362018e-2
      real(kind=8):: Q4 =  0.163801e-2 
      real(kind=8):: Q5 = -0.1031555e-1
      real(kind=8):: Q6 =  0.2282967e-1
      real(kind=8):: Q7 = -0.2895312e-1 
      real(kind=8):: Q8 =  0.1787654e-1
      real(kind=8):: Q9 = -0.420059e-2
      
      real(kind=8):: Y, AX

      if (abs(x) .lt. 3.75) then
        Y=(x/3.75)**2
        bessi1 = x*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      else
        AX=abs(x)
        Y=3.75/AX
        bessi1 = (exp(AX)/sqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     &      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
        if (x .lt. 0.) bessi1 = -bessi1
      endif

      return
      end
c=============================================================================
      real(kind=8) function bessi(n,x)
      integer(ISZ):: n 
      real(kind=8):: x

c bessi(n,x) calculates the Bessel function I_n(x) for real x 
c and integer n >= 2 using downward recursion operators seeded using 
c calls to bessi0(x) for I_0(x).  
c
c The routine is slightly modified from the one presented in:
c    Numerical Recipes The Art of Scientific Computing, 1st Edition 
c    Press, Flannery, Teukolsky, and Vetterling (1986)

      integer(ISZ):: J,M
      real(kind=8):: BI,BIM,BIP,TOX,bessi0

      integer(ISZ):: IACC=40
      real(kind=8):: BIGNO = 1.0e10
      real(kind=8):: BIGNI = 1.0e-10

      if ((n .lt. 2) .or. mod(n,1) .ne. 0.) then 
        call remark('bad argument n in bessi(n,x)')
        call kaboom(1) 
      endif 
      
      if (x .eq. 0.) then 
        bessi = 0. 
      else
        TOX=2.0/x
        BIP=0.0
        BI=1.0
        bessi = 0.
        M=2*((N+INT(sqrt(float(IACC*N)))))
        do J=M,1,-1
          BIM=BIP+float(J)*TOX*BI
          BIP=BI
          BI=BIM
          if (abs(BI) .gt. BIGNO) then
            bessi=bessi*BIGNI
            BI=BI*BIGNI
            BIP=BIP*BIGNI
          endif 
          if (J .eq. N) bessi = BIP
        enddo
        bessi = bessi*bessi0(x)/BI
        if (x .lt. 0. .and. mod(n,2) .eq. 1) bessi = -bessi 
      endif 

      return
      end
c=============================================================================
      real(kind=8) function zbrent(func,x1,x2,tol,itmax)

      integer(ISZ):: itmax 
      real(kind=8):: func,x1,x2,tol

c zbrent caluclates the root of a func(x) that is bracketed between 
c x1 and x2.  The root is found to accuracy tol using up to itmax 
c iterations.  
c
c The routine is slightly modified from the one presented in:
c    Numerical Recipes The Art of Scientific Computing, 1st Edition 
c    Press, Flannery, Teukolsky, and Vetterling (1986)

  
      integer(ISZ):: ITER
      real(kind=8):: A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM

      real(kind=8):: EPS = 3.e-8

      A=x1
      B=x2
      FA=func(A)
      FB=func(B)
      if(FB*FA .gt. 0.) then
        call remark('Root must be bracketed for zbrent.') 
        call kaboom(1)  
      endif 
      FC=FB
      do ITER=1,itmax
        if(FB*FC .gt. 0.) then
          C=A
          FC=FA
          D=B-A
          E=D
        endif
        if(abs(FC) .lt. abs(FB)) then
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        endif 
        TOL1=2.*EPS*abs(B)+0.5*tol
        XM=.5*(C-B)
        if(abs(XM) .le. TOL1 .or. FB .eq. 0.) then
          zbrent = B
          return
        endif
        if(abs(E) .ge. TOL1 .and. abs(FA) .gt. abs(FB)) then
          S=FB/FA
          if(A .eq. C) then
            P=2.*XM*S
            Q=1.-S
          else
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          endif
          if(P .gt. 0.) Q=-Q
          P=ABS(P)
          if(2.*P .lt. min(3.*XM*Q-abs(TOL1*Q),abs(E*Q))) then
            E=D
            D=P/Q
          else
            D=XM
            E=D
          endif
        else
          D=XM
          E=D
        endif
        A=B
        FA=FB
        if(abs(D) .gt. TOL1) then
          B=B+D
        else
          B=B+sign(TOL1,XM)
        endif
        FB=func(B)
      enddo
      zbrent = B
      call remark('zbrent exceeding maximum iterations.')
      call kaboom(1)

      return
      end
c=============================================================================
      subroutine rk4(y,dydx,n,x,h,yout,derivs)

      integer(ISZ):: n 
      real(kind=8):: x,h
      real(kind=8):: y(n),dydx(n),yout(n) 
      external derivs 

c 4th order Runge-Kutta advance 
c   x          initial coordinate 
c   y(1:n)     initial y values at x 
c   dydx(1:n)  initial dy/dx derivaties at x 
c   h          step 
c   yout(1:n)  final y values at x + h
c                y and yout need not be distinct  
c
c   external user supplied subroutine derivs(x,y,dydx) 
c   returns the derivatives dydx evaluated at x, y(x) 
c   where y and dydx have dimension nvar.     
c
c The routine is slightly modified from the one presented in:
c    Numerical Recipes The Art of Scientific Computing, 1st Edition 
c    Press, Flannery, Teukolsky, and Vetterling (1986)

      integer(ISZ),parameter:: NMAX=10
      real(kind=8):: YT(NMAX),DYT(NMAX),DYM(NMAX)
      integer(ISZ):: I
      real(kind=8):: HH,H6,XH

      HH=H*0.5
      H6=H/6.
      XH=X+HH

c     First step
      do I=1,n
        YT(I)=y(I)+HH*dydx(I)
      enddo
c     Second step
      call derivs(XH,n,YT,DYT)
      do I=1,N
        YT(I)=Y(I)+HH*DYT(I)
      enddo 
c     Third step 
      call derivs(XH,n,YT,DYM)
      do I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
      enddo
c     Fourth step 
      call derivs(X+H,n,YT,DYT)
      do I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
      enddo

      return
      end
c=============================================================================
      subroutine rkdumb(y,x,nvar,nstep,x1,x2,y_work,dy_work,derivs)

      integer(ISZ):: nvar,nstep 
      real(kind=8):: x1,x2
      real(kind=8):: y(nvar,nstep),x(nstep),y_work(nvar),dy_work(nvar)
      external derivs 

c 4th order Runge-Kutta advance with uniform step size  
c
c   external user supplied subroutine derivs(x,y,dydx) 
c   returns the derivatives dydx evaluated at x, y(x) 
c   where y and dydx have dimension nvar.   
c
c The routine is slightly modified from the one presented in:
c    Numerical Recipes The Art of Scientific Computing, 1st Edition 
c    Press, Flannery, Teukolsky, and Vetterling (1986)

      integer(ISZ):: I,K
      real(kind=8):: H,XX 

c     load starting values 
      do I=1,nvar
        y_work(I)=y(I,1)
      enddo 
      X(1)=X1
      XX=X1
      H=(X2-X1)/float(nstep)
c     take nstep 4th order Runge-Kutta steps using rk4
      do K=1,NSTEP
        call derivs(XX,nvar,y_work,dy_work)
        call rk4(y_work,dy_work,nvar,XX,H,y_work,DERIVS)
        if(XX+H .eq. XX) then
          call remark('Stepsize not significant in rkdumb') 
          call kaboom(1) 
        endif 
c       store intermediate steps 
        XX=XX+H
        x(K+1)=XX
        do I=1,NVAR
          y(I,K+1)=y_work(I)
        enddo
      enddo 

      return
      end
c=============================================================================
      subroutine derivs_cf(x,nvar,y,dydx)
      
      integer(ISZ):: nvar
      real(kind=8):: x 
      real(kind=8):: y(nvar),dydx(nvar) 
      
c calculates derivatives dydx at x,y
c    nvar       = number equations 
c    x          = independent coordinate 
c    y(nvar)    = dependant variable at x 
c    dydx(nvar) = derivatives at x,y 
c

      dydx(1) =  y(2) 
      dydx(2) = -y(1) 
      
      return
      end
c=============================================================================
      subroutine integrate_test(y,x,nvar,nstep,x1,x2,y_work,dy_work)

      integer(ISZ):: nvar,nstep 
      real(kind=8):: x1,x2
      real(kind=8):: y(nvar,nstep),x(nstep),y_work(nvar),dy_work(nvar)
      external derivs_cf 

      call rkdumb(y,x,nvar,nstep,x1,x2,y_work,dy_work,derivs_cf)
      
      return
      end
c=============================================================================
