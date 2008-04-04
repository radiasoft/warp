#include "top.h"
c=============================================================================
c@(#) File w3d_collision.F, version $Revision: 1.5 $, $Date: 2008/04/04 21:36:07 $
c# Copyright (c) 2007-2007, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This has routines handling particle collisions
c  Alex Friedman, LLNL, (925)422-0827
c  David P. Grote, LLNL, (925)423-7194
c  Bruce Cohen, LLNL, (925)422-9823
c=============================================================================
      subroutine langevincollisions3d(lself,np,ux1,uy1,uz1,
     &                                vxmi,vymi,vzmi,density2,vthsqi2,
     &                                q1,q2,m1,m2,vthsqinit1,vthsqinit2,
     &                                dt,loglambda,epvth)
      use Constant
      logical(ISZ):: lself
      integer(ISZ):: np
      real(kind=8):: ux1(np),uy1(np),uz1(np)
      real(kind=8):: vxmi(np),vymi(np),vzmi(np),density2(np),vthsqi2(np)
      real(kind=8):: q1,q2,m1,m2,vthsqinit1,vthsqinit2
      real(kind=8):: dt,loglambda(np),epvth

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
     &            q1**2*q2**2/m1**2/(4.*pi*eps0**2)

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
        nudt = nudtcoeff*density2(i)*loglambda(i)/vthsq2**1.5
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

c       fdv = abs( eta*nudt/(1.+1.08*(vv/(sqrt2*sqrt(vthsqi2(i))))**2.76) )
        fdv = abs( eta*nudt/(1.+1.06483*(vv/(sqrt2*sqrt(vthsqi2(i))))**2.488) )

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
      end  !!langevincollisions3d
c=============================================================================

c=============================================================================
c=============================================================================
c=============================================================================

c=============================================================================

      subroutine ijcoll2d(i1,i2,np1,np2,nxp,nyp,ux1,ux2,uy1,uy2,
     &   uz1,uz2,epvth,x1,y1,x2,y2,lx,ly,aion1,zion1,aion2,zion2,nion1,nion2,
     &   vthmisq1,vthmisq2,nudt,ncint,ixbc,izbc)
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

cc ixbc=0,1,2  izbc=0,1,2  bounded,bounded,periodic

cc Add whatever "Use" statements as needed here

      Implicit None

      integer i,ip,ix,iy,ix1,iy1,ncint,nx1ny1,ixbc,izbc
      integer i1,i2,np1,np2,nxp,nyp

      integer :: istat1, istat2, istat3, istat4, istat
      integer :: istat5, istat6, istat7, istat8

      real escale,eta,vthmisqa,vthmisqb,aiona,ziona,aionb,zionb
      real vthmisqhyd,vtha,vthb,na,nb,vth1,lx,ly
      real ke1tot,ke2tot,delke12,twopi,wxa,wxb,wya,wyb,s0i

      real ux1(np1),ux2(np2),uy1(np1),uy2(np2),uz1(np1),uz2(np2)
      real x1(np1),y1(np1),x2(np2),y2(np2),xf,yf

      real nuadt,s00
      real ddxp,ddyp,epvth	
      real etotf1,etotf2,vv,vv1,vv2,fdv,d33,d11,xrand1,uxf,uxf0,uyf,uzf
      real sqrt2
      real aion1,zion1,aion2,zion2,nion1,nion2,vthmisq1,vthmisq2
      real nudt,vthsqi
      real sddx2d,wranf



cc Declare Allocatable new variables

      real, allocatable, dimension(:) :: pdmy2          ! dummy array
      real, allocatable, dimension(:) :: pdmy1          ! dummy array
      real, allocatable, dimension(:) :: vxmi           ! local mean x velocity of background
      real, allocatable, dimension(:) :: vymi           ! local mean y velocity of background
      real, allocatable, dimension(:) :: vzmi           ! local mean z velocity of background


      real, allocatable, dimension(:,:) :: s0             ! density array
      real, allocatable, dimension(:,:) :: svx             ! x-momentum array
      real, allocatable, dimension(:,:) :: svy             ! y-momentum array
      real, allocatable, dimension(:,:) :: svz             ! z-momentum array
      real, allocatable, dimension(:,:) :: se             ! energy array



      allocate(pdmy1(np1), stat=istat1)
      allocate(pdmy2(np2), stat=istat2)
      allocate(vxmi(np1), stat=istat2)
      allocate(vymi(np1), stat=istat2)
      allocate(vzmi(np1), stat=istat2)

      allocate(s0(0:nxp,0:nyp), stat=istat3)
      allocate(se(0:nxp,0:nyp), stat=istat4)
      allocate(svx(0:nxp,0:nyp), stat=istat5)
      allocate(svy(0:nxp,0:nyp), stat=istat6)
      allocate(svz(0:nxp,0:nyp), stat=istat7)




cc The strategy to use here is to pass the velocity and position
cc arrays for both species needed in the collisions

cc The test particle is species 1 and the field particle is species 2



      ddxp=sddx2d(lx,nxp)
      ddyp=sddx2d(ly,nyp)
      sqrt2 = sqrt(2.)

      nx1ny1 = (nxp+1)*(nyp+1)

cncoll ; ncint; nudt ; step ; i1; i2
c call this routine every ncint timesteps

c ncint = (1/nudt)*(1/ncoll) where nudt=nu_i *dt and ncoll=no. of
c of times we wish to apply collision operator over the period 1/nudt
c nudt and ncoll are input by user. Then ncint is calculated in ib_deck
c before calling this routine.

      twopi=8.*atan(1.)

      ke1tot=sum(ux1**2+uy1**2+uz1**2)/(3.*np1*vthmisq1)  ! ke1tot
      ke2tot=sum(ux2**2+uy2**2+uz2**2)/(3.*np2*vthmisq2); ! ke2tot


c call pdiag to compute moments of 2D particle array on grid

      pdmy1=0.
      pdmy2=0.



      call pdiag2d(x2,y2,ux2,uy2,uz2,pdmy2,np2,lx,ly,nxp,nyp,s0,
     &        svx,svy,svz,se,ixbc,izbc)

cc      write(*,*) "sum(s0)/((nxp+1)*(nyp+1))=",sum(s0)/((nxp+1)*(nyp+1))

cc      write(*,*) "x2(1:np2)=",x2(1:np2)
cc      write(*,*) "y2(1:np2)=",y2(1:np2)
cc      write(*,*) "ux2(1:np2)=",ux2(1:np2)
cc      write(*,*) "uy2(1:np2)=",uy2(1:np2)
cc      write(*,*) "uz2(1:np2)=",x2(1:np2)
cc      write(*,*) "s0(0:nxp,0:nyp)=",s0(0:nxp,0:nyp)
cc      write(*,*) "se(0:nxp,0:nyp)=",se(0:nxp,0:nyp)


      vth1 = sqrt(vthmisq1)

cc Need to define charges and masses of species 1 and 2, depends
cc on the particular code, nion1,2 & vthmisq1,2 are species 1,2 (absolute index)
cc while i1 and i2 are relative indices for collisions

      vtha = sqrt(vthmisq1)
      na=nion1
      ziona=zion1
      aiona=aion1

      vthb = sqrt(vthmisq2)
      nb=nion2
      zionb=zion2
      aionb=aion2


      if(i1.eq.i2) then
         eta = 1.
      else
         eta = 0.5*(1.+aion1/aion2)
       endif

c pdiag returns <(v_x-<v_x>)^2+"y & z">/3.=v_th^2


c v_thermal near edges of plasma is prone to big errors that
c can cause a runaway of particle energy.--DON'T NEED TO WORRY THIS.

      vthmisqb=vthmisq2
      vthmisqa=vthmisq1

cCollision rate nudt is input = nu_eff (1st species, whatever it is)
c and frquency ncint is computed elsewhere


      etotf1 = 0.
      etotf2 = 0.

cFor i1 on i2, with nudt defined relative to first species:

      do i=1,np1

         nuadt = nudt*(nb/nion1)*float(ncint)*(ziona/zion1)**2*
     &         (zionb/zion1)**2*(vth1/vthb)**3*(aion1/aiona)**2


         xf = x1(i)
         yf = y1(i)

         If(i.eq.1) write(*,*) "nudt,nuadt=",nudt,nuadt
 
c Compute grid indices near particle from (x,y), indices span 0 to nx, 0 to ny

         ix = mod(int(xf*ddxp),nxp+1)
         ix1 = mod(int(xf*ddxp)+1,nxp+1)

         iy = mod(int(yf*ddyp),nyp+1)
         iy1 = mod(int(yf*ddyp)+1,nyp+1)

c Compute linear interpolation factors

         wxb = xf*ddxp-aint(xf*ddxp)
         wxa = 1.-wxb
         wyb = yf*ddyp-aint(yf*ddyp)
         wya = 1.-wyb



c Compute local average velocities and vx**2+vy**2 using linear interpolation

         vxmi(i) = wxa*wya*svx(ix,iy)+wxb*wya*svx(ix1,iy)
     &          + wxa*wyb*svx(ix,iy1) + wxb*wyb*svx(ix1,iy1)

         vymi(i) = wxa*wya*svy(ix,iy)+wxb*wya*svy(ix1,iy)
     &          + wxa*wyb*svy(ix,iy1) + wxb*wyb*svy(ix1,iy1)
         vzmi(i) = wxa*wya*svz(ix,iy)+wxb*wya*svz(ix1,iy)
     &          + wxa*wyb*svz(ix,iy1) + wxb*wyb*svz(ix1,iy1)

         vthsqi = wxa*wya*se(ix,iy)+wxb*wya*se(ix1,iy)+1.e-6*vthmisqb
     &          + wxa*wyb*se(ix,iy1) + wxb*wyb*se(ix1,iy1)
         vthsqi = abs(vthsqi)

         s0i = wxa*wya*s0(ix,iy)+wxb*wya*s0(ix1,iy)
     &          + wxa*wyb*s0(ix,iy1)+ wxb*wyb*s0(ix1,iy1)

 
         s00 = sum(s0)/nx1ny1 

cc         if(i.eq.10) write(*,*) "s0i,s0(ix,iy),s00,nx1ny1=",
cc     &          s0i,s0(ix,iy),s00,nx1ny1
cc         if(i.eq.10) write(*,*) "s0(1:nxp,1:nyp)=",s0(1:nxp,1:nyp)
cc         if(i.eq.10) write(*,*) "se(1:nxp,1:nyp)=",se(1:nxp,1:nyp)
   
cassumes that the average density is const and moats are negligible


         if(epvth.eq.0.) epvth=0.95		!!epvth=0.95 works fine
c Apply collisions with local correction to the temperature and density
c The nu_0 NRL = 2**1.5 * nu_alpha Jones et al.

         nuadt = nuadt*(s0i/s00)*(vthmisqb/(epvth*vthsqi
     &        +(1.-epvth)*vthmisqb))**1.5
         nuadt= abs(nuadt)

         if(i.eq.1) write(*,*) "s0i,s00,vthsqi,vthmisqb,nuadt=",
     &              s0i,s00,vthsqi,vthmisqb,nuadt

cc         if(i.eq.10) write(*,*) "vthmisqb,vthsqi,nuadt=",vthmisqb,vthsqi,nuadt
 
         if(nuadt.ge.0.99) then
            write(*,*) "nuadt=",nuadt
         endif

c Compute particle velocities in local drift frame

         uxf = ux1(i)-vxmi(i)
         uyf = uy1(i)-vymi(i)
         uzf = uz1(i)-vzmi(i)

c compute the temperature before the collision

         etotf1=(uxf**2+uyf**2+uzf**2)/(3.*np1*vthmisqa)+etotf1


c Compute 3D speed and some other speed normalizations for relative
c angles in coordinate frame oriented parallel to local velocity

         vv  =  sqrt(uxf**2+uyf**2+uzf**2) + 1.e-06*sqrt(vthmisqa)
         vv1 =  sqrt(uxf**2+uyf**2)+ 1.e-06*sqrt(vthmisqa)
         vv2 =  1./(vv*vv1 +1.e-06*vthmisqa)

c calculate diffusion and drag coefficients from Pade approximations
c for isotropic F(v) and resulting kicks

c        fdv = abs( eta*nuadt/(1.+1.08*(vv/(sqrt2*sqrt(vthsqi)))**2.76) )
         fdv = abs( eta*nuadt/(1.+1.06483*(vv/(sqrt2*sqrt(vthsqi)))**2.488) )

         d33 = sqrt(abs(2.*vthsqi*(fdv/eta)*abs(-log(1.-wranf()))))
     &       *cos(twopi*wranf())

         d11 = sqrt(nuadt*2.*vthsqi/(1.+0.1856*
     &       (vv/(sqrt2*sqrt(vthsqi)))**1.848)*abs(-log(1.-wranf())))
c need to multiply d11 by cos and sin of ranf*twopi for kicks

         fdv = -vv*fdv*(1.+0.5*fdv)      !! Includes energy cons. correction ~dt**2

c Do the Monte Carlo collisional kicks in the particle velocity local frame
c and transform back to x-y-z velocity coordinates in drift frame

         xrand1=wranf()

c Later on we need a duplicate copy of uxf (after uxf is
c overwritten with a new value after the collision) for computing a
c cosine = uxf/vv1, so stash a copy of uxf in s0i



         uxf0=uxf

c reduce number of divides by replacing vv with 1./vv

         vv=1./vv

         uxf = uxf + (fdv+d33)*uxf*vv + d11*cos(twopi*xrand1)*uyf/vv1 
     &       + d11*sin(twopi*xrand1)*uzf*uxf*vv2

         uyf = uyf + (fdv+d33)*uyf*vv + d11*cos(twopi*xrand1)*(-uxf0)/vv1 
     &       + d11*sin(twopi*xrand1)*uzf*uyf*vv2

         uzf = uzf + (fdv+d33)*uzf*vv + d11*sin(twopi*xrand1)*(-vv1**2)*vv2


c Galilean tranformation to lab velocity coordinates from drift frame

         ux1(i) = uxf
         uy1(i) = uyf
         uz1(i) = uzf

c compute temperature after the collision
         etotf2=(uxf**2+uyf**2+uzf**2)/(3.*np1*vthmisqa)+etotf2

      enddo   !!loop over particle blocks

      escale=sqrt(etotf1/etotf2)

      write(*,*) "escale=",escale

      do i = 1,np1

         uxf = ux1(i)
         uyf = uy1(i)
         uzf = uz1(i)



cScale the particle velocities to conserve energy after the
ccollisions.  Manheimer et al. suggest doing this.
cThis collision scheme does NOT conserve energy like Jones
calgorithm does, except statistically for an ideal Maxwellian
cdistribution of test particles, and erf function drag-diff coeffs,
cnone of which we satisfy.


         if(i1.eq.i2) then
            uxf=uxf*escale
            uyf=uyf*escale
            uzf=uzf*escale
         endif

c Galilean tranformation to lab velocity coordinates from drift frame

         ux1(i) = uxf + vxmi(i)
         uy1(i) = uyf + vymi(i)
         uz1(i) = uzf + vzmi(i)

      enddo



!! optional writes of key collision quantities
      write(*,*)  "nuadt=",nuadt
      write(*,*)   "etotf1=",etotf1
      write(*,*)   "etotf2=",etotf2
      write(*,*)   "escale=",escale


      ke1tot=sum(ux1**2+uy1**2+uz1**2)/(3.*np1*vthmisq1)
      write(*,*) "ke1tot=",ke1tot

      ke2tot=sum(ux2**2+uy2**2+uz2**2)/(3.*np2*vthmisq2); 
      write(*,*) "ke2tot=",ke2tot
!!      delke12 = ke1tot-(mion2/mion1)*ke2tot*(vthmisq2/vthmisq1) !! delke12


ccpaws

cc Deallocate arrays

      deallocate(pdmy1, stat=istat1)
      deallocate(pdmy2, stat=istat2)
      deallocate(vxmi, stat=istat1)
      deallocate(vymi, stat=istat1)
      deallocate(vzmi, stat=istat1)

      deallocate(s0, stat=istat3)
      deallocate(se, stat=istat4)
      deallocate(svx, stat=istat5)
      deallocate(svy, stat=istat6)
      deallocate(svz, stat=istat7)



      end  !!ijcoll2d


c pdiag.m

c   nxp,nzp here can differ from field arrays nx,nz
c    need to be able to do only a part of the system?

      subroutine pdiag2d(x,z,ux,uy,uz,gammai,np,lxp,lzp,nxp,nzp,
     & s0,sx,sy,sz,se,ixbc,izbc)
      
      Implicit None

      integer np, nxp, nzp
      real x(np), z(np), ux(np), uy(np), uz(np), gammai(np)
c  only need for nx,nz

      real s0(0:nxp,0:nzp), sx(0:nxp,0:nzp), sy(0:nxp,0:nzp), sz(0:nxp,0:nzp), se(0:nxp,0:nzp)
      real lxp, lzp

      integer i, ix, iz,ixbc,izbc
      real ddxp, ddzp, rx, rz, wxa, wxb, wza, wzb
      real wxaza, wxazb, wxbza, wxbzb
      real qvc, tx, ty, tz, tem
      real sddx2d

      ddxp = sddx2d(lxp,nxp)
      ddzp = sddx2d(lzp,nzp)

cc      write(*,*) "ddxp=","ddzp=",ddxp,ddzp
cc      write(*,*) "lxp,lzp,nxp,nzp=",lxp,lzp,nxp,nzp
      qvc = ddxp*ddzp

      do ix = 0,nxp
      do iz = 0,nzp
       s0(ix,iz) = 0.
       sx(ix,iz) = 0.
       sy(ix,iz) = 0.
       sz(ix,iz) = 0.
       se(ix,iz) = 0.
      enddo
      enddo

      do i = 1,np
       rx = x(i)*ddxp
       rz = z(i)*ddzp
       ix = rx
       iz = rz
c      Bilinear (area) weighting.
       wxb = rx - ix
       wzb = rz - iz
       wxa = 1. - wxb
       wza = 1. - wzb
       wxaza = wxa*wza
       wxazb = wxa*wzb
       wxbza = wxb*wza
       wxbzb = wxb*wzb
       s0(ix  ,iz  ) = s0(ix  ,iz  ) + qvc*wxaza
       s0(ix  ,iz+1) = s0(ix  ,iz+1) + qvc*wxazb
       s0(ix+1,iz  ) = s0(ix+1,iz  ) + qvc*wxbza
       s0(ix+1,iz+1) = s0(ix+1,iz+1) + qvc*wxbzb
       tx = qvc*ux(i)
       sx(ix  ,iz  ) = sx(ix  ,iz  ) + tx*wxaza
       sx(ix  ,iz+1) = sx(ix  ,iz+1) + tx*wxazb
       sx(ix+1,iz  ) = sx(ix+1,iz  ) + tx*wxbza
       sx(ix+1,iz+1) = sx(ix+1,iz+1) + tx*wxbzb
       ty = qvc*uy(i)
       sy(ix  ,iz  ) = sy(ix  ,iz  ) + ty*wxaza
       sy(ix  ,iz+1) = sy(ix  ,iz+1) + ty*wxazb
       sy(ix+1,iz  ) = sy(ix+1,iz  ) + ty*wxbza
       sy(ix+1,iz+1) = sy(ix+1,iz+1) + ty*wxbzb
       tz = qvc*uz(i)
       sz(ix  ,iz  ) = sz(ix  ,iz  ) + tz*wxaza
       sz(ix  ,iz+1) = sz(ix  ,iz+1) + tz*wxazb
       sz(ix+1,iz  ) = sz(ix+1,iz  ) + tz*wxbza
       sz(ix+1,iz+1) = sz(ix+1,iz+1) + tz*wxbzb
      enddo

c  Establish periodicity in z.-- assume periodic in z
cc      if(nyp = ny) then

      if(izbc.eq.2) then

         do ix = 0,nxp
            s0(ix, 0) = s0(ix,0) + s0(ix,nzp)
            s0(ix,nzp) = s0(ix,0)
            sx(ix, 0) = sx(ix,0) + sx(ix,nzp)
            sx(ix,nzp) = sx(ix,0)
            sy(ix, 0) = sy(ix,0) + sy(ix,nzp)
            sy(ix,nzp) = sy(ix,0)
            sz(ix, 0) = sz(ix,0) + sz(ix,nzp)
            sz(ix,nzp) = sz(ix,0)
         enddo
      endif

c  Establish periodicity in x.
c***not if open sided bc
       if(ixbc.eq.2) then
          do iz = 0,nzp
             s0(0, iz) = s0(0,iz) + s0(nxp,iz)
             s0(nxp,iz) = s0(0,iz)
             sx(0, iz) = sx(0,iz) + sx(nxp,iz)
             sx(nxp,iz) = sx(0,iz)
             sy(0, iz) = sy(0,iz) + sy(nxp,iz)
             sy(nxp,iz) = sy(0,iz)
             sz(0, iz) = sz(0,iz) + sz(nxp,iz)
             sz(nxp,iz) = sz(0,iz)
          enddo
       endif

c  For now, just divide by s0 to get mean velocities.
      do ix = 0,nxp
      do iz = 0,nzp
       if(s0(ix,iz).ne. 0.) then
        sx(ix,iz) = sx(ix,iz)/s0(ix,iz)
        sy(ix,iz) = sy(ix,iz)/s0(ix,iz)
        sz(ix,iz) = sz(ix,iz)/s0(ix,iz)
       endif
      enddo
      enddo

      do i = 1,np
       rx = x(i)*ddxp
       rz = z(i)*ddzp
       ix = rx
       iz = rz
c      Bilinear (area) weighting.
       wxb = rx - ix
       wzb = rz - iz
       wxa = 1. - wxb
       wza = 1. - wzb
       wxaza = wxa*wza
       wxazb = wxa*wzb
       wxbza = wxb*wza
       wxbzb = wxb*wzb
       tx = wxaza*sx(ix,iz  ) + wxbza*sx(ix+1,iz) +
     &            wxazb*sx(ix,iz+1) + wxbzb*sx(ix+1,iz+1)
       ty = wxaza*sy(ix,iz  ) + wxbza*sy(ix+1,iz) +
     &            wxazb*sy(ix,iz+1) + wxbzb*sy(ix+1,iz+1)
       tz = wxaza*sz(ix,iz  ) + wxbza*sz(ix+1,iz) +
     &            wxazb*sz(ix,iz+1) + wxbzb*sz(ix+1,iz+1)
       tem = ((ux(i) - tx)**2 + (uy(i) - ty)**2 + (uz(i) - tz)**2)/3.
       se(ix  ,iz  ) = se(ix  ,iz  ) + tem*wxaza*qvc
       se(ix  ,iz+1) = se(ix  ,iz+1) + tem*wxazb*qvc
       se(ix+1,iz  ) = se(ix+1,iz  ) + tem*wxbza*qvc
       se(ix+1,iz+1) = se(ix+1,iz+1) + tem*wxbzb*qvc
cc       write(*,*) "i,x(i),z(i)=",i,x(i),z(i)
cc       write(*,*) "ix,iz=",ix,iz
cc       write(*,*) "i,rx,rz,wxb,wzb,se(ix,iz),se(ix,iz+1),se(ix+1,iz),
cc     &     se(ix+1,iz+1)="
cc       write(*,*) i,rx,rz,wxb,wzb,se(ix,iz),se(ix,iz+1),se(ix+1,iz),
cc     &     se(ix+1,iz+1)
  
 
      enddo

c  Establish periodicity in z.  assume periodic in z
cc      if(nyp = ny) then
      if(izbc.eq.2) then
         do ix = 0,nxp
            se(ix, 0) = se(ix,0) + se(ix,nzp)
            se(ix,nzp) = se(ix,0)
         enddo
      endif

c  Establish periodicity in x.
c***not if open sided bc
      if(ixbc.eq.2) then
         do iz = 0,nzp
            se(0, iz) = se(0,iz) + se(nxp,iz)
            se(nxp,iz) = se(0,iz)
         enddo
      endif

c  Now have thermal energy, except for missing mass factor

c  For now, just divide by s0 to get mean particle thermal energy.
      do ix = 0,nxp
      do iz = 0,nzp
       if(s0(ix,iz).ne. 0.) se(ix,iz) = se(ix,iz)/s0(ix,iz)
      enddo
      enddo

cc      do ix = 0,nxp
cc         do iz = 0,nzp
cc            write(*,*) "ix,iz,s0,se,svx,svy,svz=",
cc     &              ix,iz,s0(ix,iz),se(ix,iz),sx(ix,iz)
cc     &              ,sy(ix,iz),sz(ix,iz)
cc            write(*,*)

cc         enddo
cc      enddo

      return
      end

      real function sddx2d(lx,nx)

      Implicit None
      real lx
      integer nx
      real ddx, eps, loddx, hiddx
c  sddx = nx/lx, or is reduced to guarantee lx*ddx <= nx.

      sddx2d = nx/lx
cc      write(*,*) "sddx2d=",sddx2d
      eps = lx*sddx2d - nx
cc      write(*,*) "eps=",eps
      if(eps .lt. 0.) return

c  Sigh. Bisection method to find a suitable ddx.
      ddx = sddx2d
      loddx = ddx*(1-2.e-7)
      hiddx = ddx
      do
        ddx = 0.5*(loddx + hiddx)
        eps = lx*ddx - nx
        if(ddx .eq. loddx)goto 666
        if(ddx .eq. hiddx) then
          ddx = loddx
          eps = lx*ddx - nx
          goto 666
        endif
        if(eps < 0) then
          loddx = ddx
        else
          hiddx = ddx
          ddx = loddx
        endif
      enddo
 666  continue
      sddx2d = ddx
cc      write(*,*) "sddx2d=",sddx2d

      if(lx*ddx .ge. nx) then
        write(*,*) 'sddx2d error'
      endif
      return
      end




