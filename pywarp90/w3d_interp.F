#include "top.h"
c============================================================================ 
c     ROUTINES TO INTERPOLATE BETWEEN DRIFT KINETICS AND FULL ION ORBIT
c     These are:
c     setptrs, sets temporary E and B arrays in database equal to the E's and
c        B's being passed as arguments elsewhere in w3d.
c        Thus bxd is the database equivalent of bx, etc.
c        Presently not using this; had trouble making it work.
c     mugrdbpush, does the mu grad B parallel acceleration and corresponding
c       change to vperp
c     xpush3dintrp, does the interpolated x push
c     getvperpparsq, finds v_perp^2, v_parallel^2, vparallel/B,v^2
c         of particles
c     getveff, gets components of interpolated velocity used in x push
c       Also calls getvperpparsq, and calls setfields on corrector step
c     getvdrift, calculates vdrift from ExB and gradB
c     setfields, sets E,B at particle arrays
c     getgradbsq, calculates or fetches grad B^2 components, and 
c       interpolation parameter alpha and its complement alphabar
c     NOTE CURRENTLY DOES NOT UPDATE GAMINV, i.e. nonrelativistic
c     
c============================================================================ 
c      subroutine setptrs(bx,by,bz,ex,ey,ez)
c      use DKInterptmp
c      real(kind=8), target:: bx(:),by(:),bz(:),ex(:),ey(:),ez(:)
c      bxd => bx
c      byd => by
c      bzd => bz
c      exd => ex
c      eyd => ey
c      ezd => ez
c      return
c      end
c============================================================================
      subroutine mugrdbpush(npd,is,ipmin,dtb,needcalcgradb)
c     Determines effective B to use in an extra bpush to implement
c     the mu dB/ds correction for the parallel velocity and the
c     corresponding change in vperp.  This change is of the form
c     of a rotation of v about an axis in the direction of
c     v x B, thorugh an angle alphabar mu dB/ds delta t/m vperp.
c     The required Beff is thus
c      Beff = alphabar*((m dB/ds)/2 q B^2) v x B
c           = 0.25 alphabar (m/qB^4)(B dot grad B^2)(v x B)
c     In all of above alphabar=1-alpha, complement of interp param.
c      q is electron charge.
      use Particles
      use DKInterp
      use DKInterptmp
      use InGen
      integer(ISZ):: npd,is,ipmin,needcalcgradb
      real(kind=8):: dtb
c
      integer(ISZ):: ip,ipmin1,i
      real(kind=8):: bdbsqds,coef,bxeff(npd),byeff(npd),bzeff(npd)
      real(kind=8):: vxbx,vxby,vxbz

c     Note: fields are set in padvnc3d before calling mugrdbpush.
c     So we don't need to calculate them here, but we do need
c     to calculate grad B^2 and the interpolation parameter
c     if this is the grad B push before the B push.
c     Note grdbsq(i,) is dbsqdz, dbsqdx, dbsqdy for i=1,2,3.

      if (needcalcgradb .ne. 0) call getgradbsq(npd,is,ipmin,xp(ipmin),
     & yp(ipmin),zp(ipmin))
c     print*,"start mugrdb, uxp,dtb = ", uxp(ipmin),needcalcgradb,dtb

      ipmin1=ipmin-1
      do i = 1,npd
        ip = i+ipmin1
        vxbx=uyp(ip)*bz(i)-uzp(ip)*by(i)
        vxby=uzp(ip)*bx(i)-uxp(ip)*bz(i)
        vxbz=uxp(ip)*by(i)-uyp(ip)*bx(i)
        bdbsqds = bx(i)*grdbsq(2,i)+by(i)*grdbsq(3,i)+bz(i)*grdbsq(1,i)
        coef=0.25*alphabar(i)*m_over_q(is)*bsqi(i)**2*bdbsqds
        bxeff(i)=coef*vxbx
        byeff(i)=coef*vxby
        bzeff(i)=coef*vxbz
      enddo
      call bpush3d (npd,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &           bxeff,byeff,bzeff,sq(is),sm(is),dtb,ibpush)
c     print*,"end mugrdb, uxp = ", uxp(ipmin)
      return
      end
c=========================================================================
      subroutine xpush3dintrp(npd,is,ipmin)
c     3D xpush with interpolation between drift kinetic and full PIC
      use Particles
      use DKInterp
      use DKInterptmp
      use InGen
      use InGen3d
      use Picglb
      use InMesh3d
      use Picglb3d
      integer(ISZ):: npd,is,ipmin
      integer(ISZ):: ip,i,ipmin1,ipc
      real(kind=8):: xpr(npd),ypr(npd),zpr(npd),xbar(npd)
      real(kind=8):: ybar(npd),zbar(npd)
      real(kind=8):: acntrbar
      real(kind=8):: uxscratch(npd),uyscratch(npd),uzscratch(npd)
      real(kind=8):: gaminvscratch(npd)
c     get effective velocities for old positions
      call getveff(npd,is,ipmin,xp(ipmin),yp(ipmin),zp(ipmin),"predictor")
c     push to get x': (predictor)
c     If acntr is 0, skip to push to get final x
      if (acntr(is) > 1.e-20) then

c      BEGIN predictor-corrector loop for v_eff
       do ipc = 1,npredcor
        ipmin1=ipmin-1
        do i=1,npd      
          ip=i+ipmin1
          xpr(i)=xp(ip)
          ypr(i)=yp(ip)
          zpr(i)=zp(ip)
        enddo
c       print*,"xpush1",ipmin,dt,alpha(1)
c       print*,npd,xpr(1),ypr(1),zpr(1),uxeff(1),uyeff(1),uzeff(1),gaminv(1)
        call xpush3d (npd,xpr,ypr,zpr,uxeff,uyeff,uzeff,gaminv(ipmin),dt)
c       Now construct averaged positions; self.actr is centering param
        acntrbar=1.-acntr(is)
        do i=1,npd      
          ip=i+ipmin1
          xbar(i)=acntr(is)*xpr(i)+acntrbar*xp(ip)
          ybar(i)=acntr(is)*ypr(i)+acntrbar*yp(ip)
          zbar(i)=acntr(is)*zpr(i)+acntrbar*zp(ip)
        enddo
c       Note need to make sure the trial position doesn't move a 
c       particle out of bounds of the field array.  IF it does,
c       apply particle b.c.s.  Don't change the v's so use
c       scratch arrays.  Also this won't do the right
c       thing for absorbing b.c.'s.
        call stckxy3d(npd,xbar,xmmax,xmmin,dx,ybar,ymmax,
     &                   ymmin,dy,zbar,zmmin,dz,uxscratch,
     &                   uyscratch,uzscratch,gaminvscratch,
     &                   zgrid,zbeam,l2symtry,l4symtry,pboundxy,.true.)
        call zpartbndwithdata(npd,zbar,uzscratch,gaminvscratch,
     &                        zmmax,zmmin,dz,zgrid)

c       Now get average velocity using these positions
        call getveff(npd,is,ipmin,xbar,ybar,zbar,"corrector")
       enddo
c      END predictor-corrector loop for v_eff

      endif

c     Now do the final particle push (corrector) using these veffs.
c     This time call xpush3d with xp, etc, so that the
c     real particle arrays  will be updated to the new values
c     note if interpolation is such that veff is the pure v,
c     (that is if alpha = 1), then
c     the result of the final push will be identical to
c     a full conventional push.
        
c     print*,"xpush2",dt,alpha(1)
c     print*,npd,xp(1),yp(1),zp(1),uxeff(1),uyeff(1),uzeff(1),gaminv(1)
      call xpush3d (npd,xp(ipmin),yp(ipmin),zp(ipmin),uxeff,
     &               uyeff,uzeff,gaminv(ipmin),dt)

      return
      end

c=========================================================================
      subroutine getvperpparsq(npd,ipmin)
c     find vperp^2 and vpar^2
c     calculate vpar dot B
      use Particles
      use DKInterptmp
      integer(ISZ):: npd,ipmin

      real(kind=8):: vparB
      integer(ISZ):: i,ip,ipmin1

c     print*,"*******SETTING VPERPPAR"
c     print*,uxp(ipmin),uzp(ipmin),bx(1),bz(1)
      ipmin1=ipmin-1
      do i=1,npd
         ip = i+ipmin1
         vparB = uxp(ip)*bx(i) + uyp(ip)*by(i) + uzp(ip)*bz(i)
         vparoverB(i) = vparB*bsqi(i)
         vparsq(i) = vparB**2*bsqi(i)
c        total velocity squared:
         vsq(i) = uxp(ip)**2 + uyp(ip)**2 + uzp(ip)**2
c         vperpsq(i)=vsq(i)-vparsq(i)
       enddo
       return
       end

c==========================================================================
      subroutine getveff(npd,is,ipmin,x,y,z,predcor)
c     find the effective velocity to use in the xpush.
      use DKInterptmp
      use Particles
      integer(ISZ):: npd,is,ipmin
      character(*):: predcor
      real(kind=8):: x(npd),y(npd),z(npd)

      real(kind=8):: vdkx,vdky,vdkz
      integer(ISZ):: i,ip,ipmin1

c     On the predictor step, fields have been pre-computed.
c     On the corrector step, need to recompute, as the
c      positions have changed.  
      if (predcor .eq. "corrector") call setfields(npd,is,ipmin,x,y,z)
c     in either case need proper vperp, vpar:
      call getvperpparsq(npd,ipmin)

c     First, get the drift velocity
      call getvdrift(npd,is,ipmin,x,y,z)

      ipmin1=ipmin-1
      do i=1,npd
         ip = i+ipmin1
c     add in vparallel to get the drift-kinetic velocity
         vdkx = vdx(i) + vparoverB(i)*bx(i)
         vdky = vdy(i) + vparoverB(i)*by(i)
         vdkz = vdz(i) + vparoverB(i)*bz(i)

c         if (i==1 .and. ipmin==1) print*, "in getveff",
c     1      vdx(i),vparoverB(i),bx(i),vdkx,alpha(1),uxp(1)
c     calculate interpolated velocity
         uxeff(i)=alpha(i)*uxp(ip) + alphabar(i)*vdkx
         uyeff(i)=alpha(i)*uyp(ip) + alphabar(i)*vdky
         uzeff(i)=alpha(i)*uzp(ip) + alphabar(i)*vdkz
       enddo
       return
       end

c==========================================================================
      subroutine getvdrift(np,is,ipmin,x,y,z)
c     Assumes E, B, grad B^2 are pre-computed
      use DKInterp
      use DKInterptmp
      integer(ISZ):: np,is,ipmin
      real(kind=8):: x(np),y(np),z(np)

      integer(ISZ):: i
      real(kind=8):: bgradbbx,bgradbby,bgradbbz,coeff

c     first, vgradb:
c     calculate sum of curvature and grad-B drifts, assuming vacuum B field
c     formula: vd = [m/(2 q B^4)](vpar^2 + vperp^2/2)(Bvec x grad B^2)
c     first get gradb; this also sets up B_j and |B| in the database.

c    call getvperpparsq() -- no longer needed
c      call getgradbsq(np,is,ipmin,x,y,z) -- only do it when fields
c      are recomputed, so moved into setfields

c     comps of B x grad B**2:
      do i=1,np
        bgradbbx = by(i)*grdbsq(1,i) - bz(i)*grdbsq(3,i)
        bgradbby = bz(i)*grdbsq(2,i) - bx(i)*grdbsq(1,i)
        bgradbbz = bx(i)*grdbsq(3,i) - by(i)*grdbsq(2,i)

c      common coeff:
        coeff = 0.25*m_over_q(is)*bsqi(i)**2*(vsq(i) + vparsq(i))

c      and now the components:
        vbx(i) = coeff*bgradbbx
        vby(i) = coeff*bgradbby
        vbz(i) = coeff*bgradbbz

c      now get v_ExB:
        vex(i) = (ey(i)*bz(i)-ez(i)*by(i))*bsqi(i)
        vey(i) = (ez(i)*bx(i)-ex(i)*bz(i))*bsqi(i)
        vez(i) = (ex(i)*by(i)-ey(i)*bx(i))*bsqi(i)

c      If needed, add computation of vpolarization here.

c      Construct vdrift as sums of individual drifts.
c        vdx(i) = vex(i) + vbx(i) + vpolx(i)
c        vdy(i) = vey(i) + vby(i) + vpoly(i)
c        vdz(i) = vez(i) + vbz(i) + vpolz(i)
        vdx(i) = vex(i) + vbx(i)
        vdy(i) = vey(i) + vby(i)
        vdz(i) = vez(i) + vbz(i)
      enddo
c      if (ipmin == 1) print*, "in getvdrift",
c    1      vex(1),ey(1),bz(1),ez(1),by(1),bsqi(1),vdx(1)
c     print*,"drifts, x,ex, bx,bsqi ", x(1),ex(1),bx(1),bsqi(1)
c     print*,"vex,vbx,vdx",vex(1),vbx(1),vdx(1)
      return
      end

c=========================================================================
c      subroutine getvpol(is,np,dt)
c
c      use DKInterptmep
c      integer(ISZ):: is,np,dt
c      real(kind=8):: dvexdt,dveydt,dvezdg,qbsqi
c      integer(ISZ):: i
c      dti=1./dt
c      do i=1,np
c        dvexdt=(vex(i)-vexold(i))*dti
c        dveydt=(vey(i)-veyold(i))*dti
c        dvezdt=(vez(i)-vezold(i))*dti
c        qbsqi=m_over_q(is)*bsqi(i)
c        vpolx(i) = (by(i)*dvezdt-bz(i)*dveydt)*qbsqi
c        vpoly(i) = (bz(i)*dvexdt-bx(i)*dvezdt)*qbsqi
c        vpolz(i) = (bx(i)*dveydt-by(i)*dvexdt)*qbsqi
c      enddo
c      return
c      end  

c=========================================================================
      subroutine setfields(npd,is,ipmin,x,y,z)
c     sets ej and bj.
      use DKInterp
      use DKInterptmp
      use InGen
      use InPart
      use Beam_acc
      use Z_arrays
      use Picglb
      use Particles
      use GlobalVars
      integer(ISZ):: npd,is,ipmin
      real(kind=8):: x(npd),y(npd),z(npd)
      integer(ISZ):: i
      real(kind=8):: dtr
      real(kind=8):: bendres(nparpgrp), bendradi(nparpgrp)
c     print*,"*****SETTING FIELDS"
      call fetche3d1(npd,x,y,z)
c      call fetche3d(ipmin,npd,is)
c     add in ears and uniform focusing E field pieces
      call othere3d(npd,x,y,z,zbeam,zimax,zimin,straight,ifeears,eears,
     &                     eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,
     &                     ex,ey,ez,bx,by,bz)
c      call geteb(npd,is,ipmin,x,y,z)
      dtr=.5*dt
      bx(:)=0.;by(:)=0.;bz(:)=0.
      call exteb3d(npd,x,y,z,uzp(ipmin),gaminv(ipmin),-dtr,dtr,
     &         bx,by,bz,ex,ey,ez,sm(is),
     &         sq(is),bendres,bendradi,gammabar,dt)
c     Now get the gradients of B**2
      call getgradbsq(npd,is,ipmin,x,y,z)
      return
      end

c=========================================================================
      subroutine getgradbsq(np,is,ipmin,x,y,z)
c      calculate grad B**2;  also calculates B^2, 1/B^2, and 
c      interpolation parameter alpha and its complement
c      Note shape of gradbsq is (3,np), and the components
c       are grdbsq(1,:) = dbsq/dz; grdbsq(2,:) = dbsq/dx;
c       grdbsq(3,:)=dbsq/dy
c      How grad B^2 is gotten depends on igradb:
c       igradb=1, interpolate from precalculated lookup table
c       igradb=2, calculate assuming pure quadrupole field,
c         ignoring dB^2/dz
c       igradb=3, calculate transverse components from pure
c          quadrupole; dB^2/dz from interpolation of lookup table.
c       To use igradb=3 efficiently, should invoke setup of top.bsqgrad
c         with "zonly" set to true, so that only the dB^2/dz data
c         is stored and fetched.  In this case also w3d.ngrdb should
c         be set to 1 and BSQGRADdata gchanged.
      use DKInterp
      use DKInterptmp
      use InGen
      use BSQGRADdata
      integer(ISZ):: np,is,ipmin
      real(kind=8):: x(np),y(np),z(np)
      integer(ISZ):: i,ip
      real(kind=8):: twobbrsqi,bsq_safe
      real(kind=8):: omegacesq,omegadtsq
c
c     Assumes  E and B are already set, either via setfields
c      or via calls in padvnc3d.
c     Calculate B**2 and interpolation parameter
      if (ipalpha == 1) then
       do i=1,np
          bsq(i)=bx(i)**2 + by(i)**2 + bz(i)**2
          bsq_safe=bsq(i)+dksmall
          bsqi(i) = 1./bsq_safe
c         calculate interpolation parameter, 0 if magnetized, 1 if not
          omegacesq = qovermsq(is)*bsq(i)
          omegadtsq = omegacesq*dt**2
          alpha(i) = usealphacalc(is)/sqrt(1.+alphcoef*omegadtsq)  
     &       + notusealphcalc(is)
c         alpha(i) = 1./sqrt(1.+omegadtsq)
c         alpha = 1./sqrt(1.+4.*omegadtsq)
c         alpha=exp(-omegadtsq)
          alphabar(i)=1.-alpha(i)
       enddo
      elseif (ipalpha == 2) then
       do i=1,np
          bsq(i)=bx(i)**2 + by(i)**2 + bz(i)**2
          bsq_safe=bsq(i)+dksmall
          bsqi(i) = 1./bsq_safe
c         calculate interpolation parameter, 0 if magnetized, 1 if not
          omegacesq = qovermsq(is)*bsq(i)
          omegadtsq = omegacesq*dt**2
          alpha(i) = usealphacalc(is)/(1.+alphcoef*omegadtsq)  
     &       + notusealphcalc(is)
          alphabar(i)=1.-alpha(i)
       enddo
      else
       do i=1,np
          bsq(i)=bx(i)**2 + by(i)**2 + bz(i)**2
          bsq_safe=bsq(i)+dksmall
          bsqi(i) = 1./bsq_safe
c         calculate interpolation parameter, 0 if magnetized, 1 if not
          omegacesq = qovermsq(is)*bsq(i)
          omegadtsq = omegacesq*dt**2
          alpha(i) = usealphacalc(is)/(1.+alphcoef*omegadtsq)**palpha  
     &       + notusealphcalc(is)
          alphabar(i)=1.-alpha(i)
       enddo
      endif

      if (igradb == 1) then
c       igradb == 1, get data from lookup table.
c       Must zero the array first, as applybsqgrad accumulates
        grdbsq(:,1:np)=0.
        call applybsqgrad(np,x,y,np,z,.false.,grdbsq)
      elseif (igradb == 2) then
        do i=1,np
          twobbrsqi = 2.*bsq(i)/(x(i)**2+y(i)**2+dksmall)
c          dbsqdx(i) = twobbrsqi*x(i)
c          dbsqdy(i) = twobbrsqi*y(i)
c          dbsqdz(i) = 0.
          grdbsq(1,i)=0.
          grdbsq(2,i)=twobbrsqi*x(i)
          grdbsq(3,i)=twobbrsqi*y(i)
        enddo
      elseif (igradb == 3) then
c       fetch dB^2/dz.  We need only this from the bsqgrad array.
c       Verified, this works OK whether bsqgradnc is 1 or 3; no
c       harm in dimensioning grdbsq(3,..) even if bsqgradnc is 1.
c       Warning if bsqgradnc isn't 1 or 3, dbsq/dz undefined.
        call applybsqgrad(np,x,y,np,z,.false.,grdbsq)
        do i=1,np
c         Note grdbsq(1,:) directly set by applybsqgrad
          twobbrsqi = 2.*bsq(i)/(x(i)**2+y(i)**2+dksmall)
          grdbsq(2,i)=twobbrsqi*x(i)
          grdbsq(3,i)=twobbrsqi*y(i)
        enddo
      endif
      return
      end

c=========================================================================
      subroutine geteb(npd,is,ipmin,x,y,z)
c      clean fetch of electric and magnetic fields.
c      bx,by,bz written in bx(i), by(i), bz(i).  Returns ex, ey, ez
      use DKInterptmp
      use DKInterp
      use Particles
      use InGen
      use Beam_acc
      use GlobalVars
      integer(ISZ):: npd,ipmin,is
      real(kind=8):: x(npd),y(npd),z(npd)
      integer(ISZ):: i,ip
      real(kind=8):: dtr,bsq_safe
      real(kind=8):: bendres(nparpgrp), bendradi(nparpgrp)
      dtr=.5*dt
      bx(:)=0.;by(:)=0.;bz(:)=0.
      call exteb3d(npd,x,y,z,uzp(ipmin),gaminv(ipmin),-dtr,dtr,
     &         bx,by,bz,ex,ey,ez,sm(is),
     &         sq(is),bendres,bendradi,gammabar,dt)
      do i=1,npd
        bsq(i)=bx(i)**2 + by(i)**2 + bz(i)**2
        bsq_safe=bsq(i)+dksmall
        bsqi(i) = 1./bsq_safe
      enddo
      return
      end
    

c============================================================================ 
c END INTERPOLATION ROUTINES
c============================================================================
c    Temporarily sets top.xp=x, etc, calls fetch3d, and resets top.xp
      subroutine fetche3d1(npd,x,y,z)
      use Particles
      integer(ISZ):: npd,i
      real(kind=8):: x(npd),y(npd),z(npd),xpd(npd),ypd(npd),zpd(npd)
c     save values of xp, yp, zp and set xp,yp,zp to the x,y,z arrays
      do i=1,npd
         xpd(i)=xp(i)
         ypd(i)=yp(i)
         zpd(i)=zp(i)
         xp(i)=x(i)
         yp(i)=y(i)
         zp(i)=z(i)
      enddo
      call fetche3d(1,npd,1)
c     reset xp,yp,zp
      do i=1,npd
         xp(i)=xpd(i)
         yp(i)=ypd(i)
         zp(i)=zpd(i)
      enddo
      return
      end
