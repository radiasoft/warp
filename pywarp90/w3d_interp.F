#include "top.h"
c============================================================================ 
c     ROUTINES TO INTERPOLATE BETWEEN DRIFT KINETICS AND FULL ION ORBIT
c     These are:
c     setptrs, sets temporary E and B arrays in database equal to the E's and
c        B's being passed as arguments elsewhere in w3d.
c        Thus bxd is the database equivalent of bx, etc.
c        Presently not using this; had trouble making it work.
c     mugrdbpush, does the mu grad B parallel acceleration and corresponding
c       change to uperp
c     xpush3dintrp, does the interpolated x push
c     getvperpparsq, finds u_perp^2, u_parallel^2, uparallel/B,v^2
c         of particles
c     getveff, gets components of interpolated velocity used in x push
c       Also calls getvperpparsq, and calls setfields on corrector step
c     getvdrift, calculates vdrift from ExB and gradB
c     setfields, sets E,B at particle arrays
c     getgradbsq, calculates or fetches grad B^2 components, and 
c       interpolation parameter alpha and its complement alphabar
c     This version has relativistic corrections.  Basic approach:
c       gamma should be updated from full Boris velocity push,
c       as it is only the Epush that affects gamma.
c     
c============================================================================
      subroutine mugrdbpush(pgroup,npd,is,ipmin,dtb,needcalcgradb)
c     Determines effective B to use in an extra bpush to implement
c     the mu dB/ds correction for the parallel velocity and the
c     corresponding change in uperp.  This change is of the form
c     of a rotation of v about an axis in the direction of
c     v x B, thorugh an angle alphabar mu dB/ds delta t/m uperp.
c     The required Beff is thus
c      Beff = alphabar*((gamma*m dB/ds)/2 q B^2) v x B
c           = 0.25 alphabar (m/qB^4)(B dot grad B^2)(v x B)
c     In all of above alphabar=1-alpha, complement of interp param.
c      q is electron charge.
      use ParticleGroupmodule
      use DKInterp
      use DKInterptmp
      use InGen, only: ibpush
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,is,isid,ipmin,needcalcgradb
      real(kind=8):: dtb
c
      integer(ISZ):: ip,ipmin1,i
      real(kind=8):: bdbsqds,coef,bxeff(npd),byeff(npd),bzeff(npd)
      real(kind=8):: vxbx,vxby,vxbz

      isid = pgroup%sid(is-1) + 1

c     Note: fields are set in padvnc3d before calling mugrdbpush.
c     So we don't need to calculate them here, but we do need
c     to calculate grad B^2 and the interpolation parameter
c     if this is the grad B push before the B push.
c     Note grdbsq(i,) is dbsqdz, dbsqdx, dbsqdy for i=1,2,3.

      if (needcalcgradb .ne. 0) call getgradbsq(npd,is,pgroup%xp(ipmin),
     & pgroup%yp(ipmin),pgroup%zp(ipmin),
     & pgroup%bx(ipmin),pgroup%by(ipmin),pgroup%bz(ipmin))
      if (needcalcgradb .ne. 0) call getalphas(pgroup,npd,isid,ipmin,
     & 2.*dtb)
c      print*,"start mugrdb, uxp,dtb = ", pgroup%uxp(ipmin),needcalcgradb,dtb

      ipmin1=ipmin-1
      do i = 1,npd
        ip = i+ipmin1
c       The following are gamma*(v x B) so no further gammas needed
        vxbx=pgroup%uyp(ip)*pgroup%bz(ip)-pgroup%uzp(ip)*pgroup%by(ip)
        vxby=pgroup%uzp(ip)*pgroup%bx(ip)-pgroup%uxp(ip)*pgroup%bz(ip)
        vxbz=pgroup%uxp(ip)*pgroup%by(ip)-pgroup%uyp(ip)*pgroup%bx(ip)
        bdbsqds = pgroup%bx(ip)*grdbsq(2,i) +
     &            pgroup%by(ip)*grdbsq(3,i) +
     &            pgroup%bz(ip)*grdbsq(1,i)
        coef=0.25*alphabar(i)*m_over_q(is)*bsqi(i)**2*bdbsqds
        bxeff(i)=coef*vxbx
        byeff(i)=coef*vxby
        bzeff(i)=coef*vxbz
      enddo
      call bpush3d(npd,pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     &             pgroup%gaminv(ipmin),
     &             bxeff,byeff,bzeff,pgroup%sq(is),pgroup%sm(is),dtb,ibpush)
c      print*,"end mugrdb, uxp,uyp,bx = ", pgroup%uxp(ipmin),pgroup%uyp(ipmin),
c     &  pgroup%bx(ipmin)
      return
      end
c=========================================================================
      subroutine xpush3dintrp(pgroup,npd,is,ipmin,dt)
c     3D xpush with interpolation between drift kinetic and full PIC
      use ParticleGroupmodule
      use DKInterp
      use DKInterptmp
      use InGen, only: pboundxy
      use FieldSolveAPI, only: ipminfsapi,npfsapi,jsfsapi
      use InGen3d
      use Picglb
      use InMesh3d
      use Picglb3d
c      use Particles, only: chdtspid
      use Particles,Only: xoldpid,yoldpid,zoldpid,
     &                    bxoldpid,byoldpid,bzoldpid
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,is,isid,ipmin
      integer(ISZ):: ip,i,ipmin1,ipc
      real(kind=8):: xpr(npd),ypr(npd),zpr(npd),xbar(npd),dt
      real(kind=8):: ybar(npd),zbar(npd)
      real(kind=8):: acntrbar,dtuse
      real(kind=8):: uxscratch(npd),uyscratch(npd),uzscratch(npd)
      real(kind=8):: gaminvscratch(npd)
      real(kind=8),pointer:: xold(:),yold(:),zold(:)
      real(kind=8),pointer:: bxold(:),byold(:),bzold(:)
c     declare some pointers
      if (impinterp == 0) then
         xold => pgroup%pid(:,xoldpid)
         yold => pgroup%pid(:,yoldpid)
         zold => pgroup%pid(:,zoldpid)
         bxold => pgroup%pid(:,bxoldpid)
         byold => pgroup%pid(:,byoldpid)
         bzold => pgroup%pid(:,bzoldpid)
      endif

      isid = pgroup%sid(is-1) + 1

c      print*, "ENTERING INTERPOLATION ROUTINES, ey,ez=",pgroup%ey(1),pgroup%ez(1)

c     For implicit interpolated mover, the beginning of the step is
c     the corrector for the previous timestep's predictor; we need to
c     retrieve the E and B at the predicted positions and assign them to
c     the current particle fields.  This is because, for the standard
c     implicit stepping, at this point the particle E's and B's contain
c     the lagged (step n) values, versus the predicted (step n+1) values.
cc      if (impinterp==1 .and. ipredcor > 0) call restoreEBFromPred(pgroup)
c    11/5/07 NO!  E and B should be E_n and B_n as they are.  That gives
c     us vdrift at n, which we want to average with the predicted vdrift
c      which we will calculate from the saved Preds in getvdrift.
c      If this is a predictor step and implicit, we need to gather the current
c      field (saved from end of previous predictor) to the current corrected
c      particle positions.  Half of it is used to advance the particles.
      if (impinterp == 1 .and. ipredcor == 0) 
     &     call setTotalE(pgroup,npd,is,ipmin,dt)
c     For non-implicit interpolated mover with leapfrog predictor,
c      time step for xpush is 2 dt for ipredcor = 0, dt for ipredcor = 1
c      but on first timestep there is no old value, so just iterate
c      from level 0 to level 1 twice
c      print*
c      print*,"start xpushintrp,xold,yold,zold,x,y,z = ",ipredcor,xold(1),yold(1),
c     & zold(1),pgroup%xp(1),pgroup%yp(1),pgroup%zp(1)
c      print*, "and uxp,uyp,uzp =", pgroup%uxp(1),pgroup%uyp(1),pgroup%uzp(1)
      dtuse = dt
c      print*, "ipredcor = ", ipredcor
      if (impinterp == 0 .and. npcmax .ge. 0 .and. ipredcor == 0 .and. 
     &  it > 1) dtuse = 2.*dt
c      print*, "it = ", it
c     get effective velocities for old positions
c     We must calculate grad B and B^2 now for ipredcor = 1; 
c     for ipredcor = 0 it is done by the calls to mugrdbpush in padvnc3d.
      if (ipredcor > 0)
     & call getgradbsq(npd,isid,pgroup%xp(ipmin),
     & pgroup%yp(ipmin),pgroup%zp(ipmin),pgroup%bx(ipmin),pgroup%by(ipmin),
     & pgroup%bz(ipmin))
c     If npcmax >= 0 and we are not using the implicit mover, 
c     we need to calculate the interpolation factors 
c     separately for ipredcor = 0 (timestep different than for mugrdbush).
      if (impinterp == 0 .and. npcmax .ge. 0) 
     &      call getalphas(pgroup,npd,isid,ipmin,dtuse)
      call getveff(pgroup,npd,is,ipmin,pgroup%xp(ipmin),pgroup%yp(ipmin),
     &             pgroup%zp(ipmin),"predictor",dt)
c      print *,"xp,vxeff,vbx,vex,alpha = ",pgroup%xp(ipmin),uxeff(1),vbx(1),vex(1),alpha(1)

c
c     If we are doing correction with average drifts,
c.....save the old B's (prev timestep) for use in constructing 
c......averages for corrector.  This is available during the predictor step.
c     We don't need this for the implicit version.
      if (impinterp == 0 .and. ipredcor == 0) then
         ipmin1=ipmin-1
         do i = 1, npd
            ip=i+ipmin1
            bxold(ip) = pgroup%bx(ip)
            byold(ip) = pgroup%by(ip)
            bzold(ip) = pgroup%bz(ip)
         enddo
      endif
c     push to get x': (predictor)
c     If acntr is 0, skip to push to get final x
      if (acntr(isid) > 1.e-20) then

c      The following is a LOCAL predictor-corrector loop to get an
c      approximation to v_d at half-integer steps.  This will not be used
c      when we do the global predictor-corrector of step3d which entails
c      predicting from n-1 to n+1, field solving, and then correcting
c      n to n+1.   To exercise the latter be sure to set 
c      npredcorv to zero.  It is also not used in base operation of
c      implicit mover.  We use the implicit approx to advanced fields
c      to do the centering.
c
c      BEGIN predictor-corrector loop for v_eff
       do ipc = 1,npredcorv
        ipmin1=ipmin-1
        do i=1,npd      
          ip=i+ipmin1
          xpr(i)=pgroup%xp(ip)
          ypr(i)=pgroup%yp(ip)
          zpr(i)=pgroup%zp(ip)
        enddo
c       print*,"xpush1",ipmin,dt,alpha(1)
c       print*,npd,xpr(1),ypr(1),zpr(1),uxeff(1),uyeff(1),uzeff(1),gaminv(1)
        call xpush3d (npd,xpr,ypr,zpr,uxeff,uyeff,uzeff,pgroup%gaminv(ipmin),
     &                dt)
c       Now construct averaged positions; self.actr is centering param
        acntrbar=1.-acntr(isid)
        do i=1,npd      
          ip=i+ipmin1
          xbar(i)=acntr(isid)*xpr(i)+acntrbar*pgroup%xp(ip)
          ybar(i)=acntr(isid)*ypr(i)+acntrbar*pgroup%yp(ip)
          zbar(i)=acntr(isid)*zpr(i)+acntrbar*pgroup%zp(ip)
        enddo
c       Note need to make sure the trial position doesn't move a 
c       particle out of bounds of the field array.  IF it does,
c       apply particle b.c.s.  Don't change the v's so use
c       scratch arrays.  Also this won't do the right
c       thing for absorbing b.c.'s.
        gaminvscratch = 1.
        call stckxy3d(npd,xbar,xmmax,xmmin,dx,ybar,ymmax,
     &                   ymmin,dy,zbar,zmminlocal,dz,uxscratch,
     &                   uyscratch,uzscratch,gaminvscratch,
     &                   zgrid,zbeam,l2symtry,l4symtry,pboundxy,.true.)
        call zpartbndwithdata(npd,zbar,uzscratch,gaminvscratch,
     &                        zmmaxlocal,zmminlocal,dz,zgrid)

        do i=1,npd      
          if (gaminvscratch(i)==0.) then
            ip=i+ipmin1
            xbar(i) = pgroup%xp(ip)
            ybar(i) = pgroup%yp(ip)
            zbar(i) = pgroup%zp(ip)
          endif
        enddo
c       Now get average velocity using these positions
        call getveff(pgroup,npd,is,ipmin,xbar,ybar,zbar,
     &               "corrector",dt)
c        print *,"xbar,npd,is,vxeff,vex,alpha = ",xbar(ipmin),npd,is,uxeff(1),vex(1),alpha(1)
       enddo
c      END predictor-corrector loop for v_eff

      endif

c     option to add extra drift from a python script called "LeeTangDrift"
c     Note the python script needs to know about ipmin and npd
      ipminfsapi = ipmin
      npfsapi = npd
      jsfsapi = isid
      if (igetextdrift > 0) call execuser("LeeTangDrift")

c     Now do the final particle push (corrector) using these veffs.
c     This time call xpush3d with xp, etc, so that the
c     real particle arrays  will be updated to the new values
c     note if interpolation is such that veff is the pure v,
c     (that is if alpha = 1), then
c     the result of the final push will be identical to
c     a full conventional push.
        
c     If this is a (global) predictor step we need to swap the x arrays
c      so we are pushing from n-1 and also defining a new set of xold's
c      for the next time step.  If it is a corrector step then copy
c      x arrays so that we are pushing from n to n+1, step n in both
c      xp and xold.  For predictor at step 1,just copy xp to xpold.
c      Only do predictor swap if it>1 so there is something to swap.
c      call swapxarrays(isid)
c      Note for implicit version all swapping is handled by the stepper, so
c      don't do it here.
c     temporary option
      if (impinterp == 0 .and. npcmax .ge. 0) call swapxarrays(pgroup,npd,is,ipmin)
      if (impinterp == 0 .and. npcmax < 0) call swapxarr2(pgroup,npd,is,ipmin)
        
c      print*,"xpush2",dt,alpha(1)
c      print*,npd,pgroup%xp(1),pgroup%yp(1),pgroup%zp(1),uxeff(1),uyeff(1),uzeff(1),pgroup%gaminv(1)
      call xpush3d (npd,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &              uxeff,uyeff,uzeff,pgroup%gaminv(ipmin),dtuse)

c      print*,"pushed",ipmin,pgroup%xp(1),pgroup%yp(1),pgroup%zp(1),
c     &    uxeff(1),uyeff(1),uzeff(1),gaminv(1)
c      print*,"pushed, xpushintrp,xold,yold,zold,x,y,z = ",ipredcor,xold(1),yold(1),
c     & zold(1),pgroup%xp(1),pgroup%yp(1),pgroup%zp(1)

      return
      end

c=========================================================================
      subroutine getvperpparsq(pgroup,npd,ipmin)
      use ParticleGroupmodule
      use DKInterptmp
c     find uperp^2 and upar^2  (momentum/mass **2)
c     calculate upar dot B
      use ParticleGroupmodule
      use Particles,Only: uxoldpid,uyoldpid,uzoldpid,bxoldpid,
     &   byoldpid,bzoldpid,uparBopid,
     &  uparoBpredpid,bxpredpid,bypredpid,bzpredpid

      use Picglb
      use DKInterp
      use DKInterptmp
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,ipmin

      real(kind=8):: uparB,uparoldB,uparBnew,Bnewsq
      integer(ISZ):: i,ip,ipmin1,ipmax
      real(kind=8),pointer:: uxold(:),uyold(:),uzold(:),uxp(:),uyp(:),uzp(:)
      real(kind=8),pointer:: bxold(:),byold(:),bzold(:),uparBo(:)
      real(kind=8),pointer:: bx(:),by(:),bz(:)
      real(kind=8),pointer:: bxpred(:),bypred(:),bzpred(:),uparoBpred(:)

      ipmin1=ipmin-1
      ipmax=ipmin1+npd
c     declare some pointers
      uxp => pgroup%uxp(ipmin:ipmax)
      uyp => pgroup%uyp(ipmin:ipmax)
      uzp => pgroup%uzp(ipmin:ipmax)
      if (impinterp == 0) then
         uxold => pgroup%pid(ipmin:ipmax,uxoldpid)
         uyold => pgroup%pid(ipmin:ipmax,uyoldpid)
         uzold => pgroup%pid(ipmin:ipmax,uzoldpid)
         uparBo => pgroup%pid(ipmin:ipmax,uparBopid)
         bxold => pgroup%pid(ipmin:ipmax,bxoldpid)
         byold => pgroup%pid(ipmin:ipmax,byoldpid)
         bzold => pgroup%pid(ipmin:ipmax,bzoldpid)
      else
         uparoBpred => pgroup%pid(ipmin:ipmax,uparoBpredpid)
         bxpred => pgroup%pid(ipmin:ipmax,bxpredpid)
         bypred => pgroup%pid(ipmin:ipmax,bypredpid)
         bzpred => pgroup%pid(ipmin:ipmax,bzpredpid)
      endif
      bx => pgroup%bx(ipmin:ipmax)
      by => pgroup%by(ipmin:ipmax)
      bz => pgroup%bz(ipmin:ipmax)

c     print*,"*******SETTING UPERPPAR"
c     print*,uxp(ipmin),uzp(ipmin),bx(ipmin),bz(ipmin)
      ipmin1=ipmin-1
c     Following coding calculates the v's for the advanced u's.
c     For global corrector this is all we need.
c     average with corresponding quantities from old v's if this is
c     predictor step.  The old v's are stored in uxpbar, ... at this time.
c.... Standard explicit interpolated mover
      if (impinterp == 0) then
        if (ipredcor == 0 .and. npcmax .ge. 0 .and. it > 1) then
          do i=1, npd
c           First generate upar with updated velocities and current-time
c           fields
c           Store it as global array for use on corrector
            uparB = uxp(i)*bx(i) + uyp(i)*by(i) + uzp(i)*bz(i)
            uparBo(i) = uparB*bsqi(i)
c           Calculate upar**2 from new velocities and current-time fields;
            uparsq_new(i) = uparB**2*bsqi(i)
c           now average these with corresponding quantities using old
c           v's for predictor
            uparoldB = uxold(i)*bx(i) + uyold(i)*by(i) + uzold(i)*bz(i)
            uparoverB(i) = 0.5*(uparBo(i)+uparoldB*bsqi(i))
            uparsq(i) = 0.5*(uparsq_new(i)+uparoldB**2*bsqi(i))
c           total momentum per mass squared:
            usq_new(i) = uxp(i)**2+uyp(i)**2+uzp(i)**2
            usq(i) = .5*(usq_new(i) + uxold(i)**2 + 
     &           uyold(i)**2 + uzold(i)**2)
c         if (i == 1) then
c            print*, " second pass:"
c            print*, ip,ipmin1
c            print*, uparB,uparoverB(i)
c         endif
         enddo
        elseif (ipredcor == 0 .and. (npcmax < 0 .or. it == 1)) then
c       Predictor based on simple projection from step 0 to step 1
         do i=1, npd
c           First generate upar with updated velocities and current-time
c           fields
c           Store it as global array for use on corrector
            uparB = uxp(i)*bx(i) + uyp(i)*by(i) + uzp(i)*bz(i)
            uparBo(i) = uparB*bsqi(i)
c           Calculate upar**2 from new velocities and current-time fields;
            uparsq_new(i) = uparB**2*bsqi(i)
c           For npcmax < 0 or it = 1 we are taking a step dt on predictor
c            and using new fields with current unadvanced positions
c            so uparoverB is the same as uparBo, the thing we save
c            for the corrector
            uparoverB(i) = uparBo(i)
            uparsq(i) = uparsq_new(i)
c           total momentum per mass squared:
            usq_new(i) = uxp(i)**2+uyp(i)**2+uzp(i)**2
            usq(i) = usq_new(i)
c         if (i == 1) then
c            print*, " second pass:"
c            print*, ip,ipmin1
c            print*, uparB,uparoverB(i)
c         endif
         enddo
        else
c.....  Corrector
        do i=1,npd
c.......   u dot (B_n + B_predicted)
           uparBnew = uxp(i)*bx(i) + uyp(i)*by(i) + uzp(i)*bz(i)
c           uparB = uparBo(i) + uparBnew
c           bsumsq = (bx(i)+bxold(i))**2 + (by(i)+byold(i))**2
c     &        + (bz(i)+bzold(i))**2
c           uparoverB(i) = 0.5*(uparB*bsqi(i)+uparnoverBo(i))
           uparoverB(i) = uparBnew*bsqi(i)
c          above is uparallel(updated)/B(predicted); will get multiplied
c          by B_j,predicted and combined with corresponding terms using
c          B at current timestep in getveff.
c          This is opposed to the following, which is upar squared
c          using the new B's, for coefficients un vdrift(x_predicted)
           usq(i) = uxp(i)**2+uyp(i)**2+uzp(i)**2
           uparsq(i) = uparBnew**2*bsqi(i)
c           uperpsq(i)=usq(i)-uparsq(i)
c           if (i == 1) then
c              print*, " first pass:"
c              print*, ip,ipmin1
c              print*, uparB,uparoverB(i)
c           endif
         enddo
        endif
c.....IMPLICIT version
      else
        if (ipredcor == 0) then
c       Predictor based on simple projection from step 0 to step 1
         do i=1, npd
c           Generate upar with velocities upated from corrector
c           and current-time fields
c           Store it as global array for use on corrector
            uparB = uxp(i)*bx(i) + uyp(i)*by(i) + uzp(i)*bz(i)
c           upar/B
            uparoverB(i) = uparB*bsqi(i)
c           Calculate upar**2
            uparsq(i) = uparB**2*bsqi(i)
c           total momentum per mass squared, usq, is calculated in the
c           corrector step, which comes first.
          enddo
        else
c.......Corrector
c.......For the implicit scheme this is the first part of the
c.......step, correcting predictor from previous timestep.
c.......We need quantities at both the current and predicted future 
c.......timesteps
          do i=1,npd
            uparB = uxp(i)*bx(i) + uyp(i)*by(i) + uzp(i)*bz(i)
            uparBnew = uxp(i)*bxpred(i)+uyp(i)*bypred(i)+uzp(i)*bzpred(i)
            Bnewsq = bxpred(i)**2+bypred(i)**2+bzpred(i)**2
            uparoverB(i) = uparB*bsqi(i)
            uparoBpred(i) = uparBnew/(Bnewsq+dksmall)
            usq(i) = uxp(i)**2+uyp(i)**2+uzp(i)**2
          enddo
        endif   
      endif
      return
      end

c==========================================================================
      subroutine getveff(pgroup,npd,is,ipmin,x,y,z,predcor,dt)
c     find the effective velocity to use in the xpush.
      use ParticleGroupmodule
      use Particles,Only: uxoldpid,uyoldpid,uzoldpid,
     &  vdxoldpid,vdyoldpid,vdzoldpid,uparBopid,bxoldpid,byoldpid,bzoldpid,
     &  uparoBpredpid,bxpredpid,bypredpid,bzpredpid
      use DKInterp
      use DKInterptmp
      use Picglb
      type(ParticleGroup):: pgroup

      integer(ISZ):: npd,is,isid,ipmin
      character(*):: predcor
      real(kind=8):: dt
      real(kind=8):: x(npd),y(npd),z(npd)

      real(kind=8):: udkx,udky,udkz,gamasafe,acntrbar
      integer(ISZ):: i,ip,ipmin1,ipmax
      real(kind=8),pointer:: uxp(:),uyp(:),uzp(:)
      real(kind=8),pointer:: uxold(:),uyold(:),uzold(:),gaminv(:)
      real(kind=8),pointer:: vdxold(:),vdyold(:),vdzold(:)
      real(kind=8),pointer:: uparBo(:),bxold(:),byold(:),bzold(:)
      real(kind=8),pointer:: ex(:),ey(:),ez(:)
      real(kind=8),pointer:: bx(:),by(:),bz(:)
      real(kind=8),pointer:: bxpred(:),bypred(:),bzpred(:),uparoBpred(:)

      ipmin1=ipmin-1
      ipmax=ipmin1+npd
c     declare some pointers
      uxp => pgroup%uxp(ipmin:ipmax)
      uyp => pgroup%uyp(ipmin:ipmax)
      uzp => pgroup%uzp(ipmin:ipmax)
      gaminv => pgroup%gaminv(ipmin:ipmax)
      vdxold => pgroup%pid(ipmin:ipmax,vdxoldpid)
      vdyold => pgroup%pid(ipmin:ipmax,vdyoldpid)
      vdzold => pgroup%pid(ipmin:ipmax,vdzoldpid)
      if (impinterp == 0) then
         uxold => pgroup%pid(ipmin:ipmax,uxoldpid)
         uyold => pgroup%pid(ipmin:ipmax,uyoldpid)
         uzold => pgroup%pid(ipmin:ipmax,uzoldpid)
         uparBo => pgroup%pid(ipmin:ipmax,uparBopid)
         bxold => pgroup%pid(ipmin:ipmax,bxoldpid)
         byold => pgroup%pid(ipmin:ipmax,byoldpid)
         bzold => pgroup%pid(ipmin:ipmax,bzoldpid)
      else
         uparoBpred => pgroup%pid(ipmin:ipmax,uparoBpredpid)
         bxpred => pgroup%pid(ipmin:ipmax,bxpredpid)
         bypred => pgroup%pid(ipmin:ipmax,bypredpid)
         bzpred => pgroup%pid(ipmin:ipmax,bzpredpid)
      endif
      bx => pgroup%bx(ipmin:ipmax)
      by => pgroup%by(ipmin:ipmax)
      bz => pgroup%bz(ipmin:ipmax)
      ex => pgroup%ex(ipmin:ipmax)
      ey => pgroup%ey(ipmin:ipmax)
      ez => pgroup%ez(ipmin:ipmax)

      isid = pgroup%sid(is-1) + 1

c     On the v-predictor step, fields have been pre-computed.
c     On the v-corrector step, need to recompute, as the
c      positions have changed.  
c     Note at least for now we don't do these correctors for
c      either implicit or the scheme with predictor-corrector in
c      step3d.
      if (predcor .eq. "corrector") then
         call setfields(pgroup,npd,is,ipmin,x,y,z,ex,ey,ez,bx,by,bz,dt)
         call getalphas(pgroup,npd,isid,ipmin,dt)
c         print*,"called setfields, x,ey,ez",x(1),ey(1),ez(1)
      endif

c     in either case need proper uperp, upar:
      call getvperpparsq(pgroup,npd,ipmin)

c     First, get the drift velocity
      call getvdrift(pgroup,npd,is,ipmin,x,y,z)

      acntrbar = 1.-acntr(isid)
c.....Begin coding for explicit interpolated mover
      if (impinterp == 0) then
       if (ipredcor > 0) then
        do i=1,npd
c       add in uparallel to get the drift-kinetic velocity
c       recall that uparoverB is u_||/B momentum/mass/B
          gamasafe = 1./(gaminv(i)+SMALLPOS)
          udkx = (gamasafe*(acntr(isid)*vdx(i)+acntrbar*vdxold(i))
     &        + acntr(isid)*uparoverB(i)*bx(i)+acntrbar*uparBo(i)*bxold(i))
          udky = (gamasafe*(acntr(isid)*vdy(i)+acntrbar*vdyold(i))
     &        + acntr(isid)*uparoverB(i)*by(i)+acntrbar*uparBo(i)*byold(i))
          udkz = (gamasafe*(acntr(isid)*vdz(i)+acntrbar*vdzold(i))
     &        + acntr(isid)*uparoverB(i)*bz(i)+acntrbar*uparBo(i)*bzold(i))

c          if (i==1 .and. ipmin==1) print*, "in getveff",
c     1       vdx(i),uparoverB(i),bx(i),udkx,alpha(1),uxp(1)
c        calculate interpolated velocity
          uxeff(i)=alpha(i)*uxp(i) + alphabar(i)*udkx
          uyeff(i)=alpha(i)*uyp(i) + alphabar(i)*udky
          uzeff(i)=alpha(i)*uzp(i) + alphabar(i)*udkz
        enddo
       endif

       if (ipredcor == 0 .and. npcmax .ge. 0 .and. it > 1) then
        do i=1,npd
c       add in uparallel to get the drift-kinetic velocity
c       recall that uparoverB is u_||/B momentum/mass/B
          gamasafe = 1./(gaminv(i)+SMALLPOS)
          udkx = gamasafe*vdx(i) + uparoverB(i)*bx(i)
          udky = gamasafe*vdy(i) + uparoverB(i)*by(i)
          udkz = gamasafe*vdz(i) + uparoverB(i)*bz(i)

c          if (i==1 .and. ipmin==1) print*, "in getveff",
c     1      vdx(i),uparoverB(i),bx(i),vdkx,alpha(1),uxp(1)
c       calculate interpolated velocity
          uxeff(i)=.5*alpha(i)*(uxp(i)+uxold(i)) + alphabar(i)*udkx
          uyeff(i)=.5*alpha(i)*(uyp(i)+uyold(i)) + alphabar(i)*udky
          uzeff(i)=.5*alpha(i)*(uzp(i)+uzold(i)) + alphabar(i)*udkz
        enddo
       endif
       if (ipredcor == 0 .and. (npcmax < 0 .or. it == 1)) then
c      Predictor based on simple projection from step 0 to step 1
        do i=1,npd
c       add in uparallel to get the drift-kinetic velocity
c       recall that uparoverB is u_||/B momentum/mass/B
          gamasafe = 1./(gaminv(i)+SMALLPOS)
          udkx = gamasafe*vdx(i) + uparoverB(i)*bx(i)
          udky = gamasafe*vdy(i) + uparoverB(i)*by(i)
          udkz = gamasafe*vdz(i) + uparoverB(i)*bz(i)

c          if (i==1 .and. ipmin==1) print*, "in getveff",
c     1      vdx(i),uparoverB(i),bx(i),udkx,alpha(1),uxp(1)
c       calculate interpolated velocity
          uxeff(i)=alpha(i)*uxp(i) + alphabar(i)*udkx
          uyeff(i)=alpha(i)*uyp(i) + alphabar(i)*udky
          uzeff(i)=alpha(i)*uzp(i) + alphabar(i)*udkz
        enddo
       endif
c       print*, "in getveff,uxp, uxeff,vex = ",ipredcor,uxp(1),uxeff(1),vex(1)
c       print*, "ipredcor, uzp,uzold,uzeff",ipredcor,uzp(1),uzold(1),uzeff(1)
      else
c.....IMPLICIT interpolated mover
       if (ipredcor == 0) then
c      Predictor based on simple projection from current to advanced step
c      Recall this comes "after corrector" (from previous step) in
c      WARP's implicit procedure
        do i=1,npd
c       add uparallel to drifts to get the drift-kinetic velocity
c       recall that uparoverB is u_||/B momentum/mass/B
          gamasafe = 1./(gaminv(i)+SMALLPOS)
          udkx = gamasafe*vdx(i) + uparoverB(i)*bx(i)
          udky = gamasafe*vdy(i) + uparoverB(i)*by(i)
          udkz = gamasafe*vdz(i) + uparoverB(i)*bz(i)

c          if (i==1 .and. ipmin==1) print*, "in getveff",
c     1      vdx(i),uparoverB(i),bx(i),udkx,alpha(1),uxp(1)
c     calculate interpolated velocity
          uxeff(i)=alpha(i)*uxp(i) + alphabar(i)*udkx
          uyeff(i)=alpha(i)*uyp(i) + alphabar(i)*udky
          uzeff(i)=alpha(i)*uzp(i) + alphabar(i)*udkz
        enddo
       else
c      CORRECTOR for implicit, corrects previous timestep's predictor
        do i=1,npd
c       add uparallel to drifts to get the drift-kinetic velocity;
c       for the corrector we need to average current and predicted
c       timestep versions.   The drifts for the old step are stored during the 
c       previous timestep's predictor and so are vdxold, etc;
c       The drifts for the current timestep are calculated using the
c       fields calculated at the end of the last timestep (predictor),
c       and saved as Epred, Bpred.
c       The vparallel's are calculated in getvperpvparsq (current timestep)
c       using the fields at the current and advanced timestep's positions.
c       recall that uparoverB is u_||/B momentum/mass/B
          gamasafe = 1./(gaminv(i)+SMALLPOS)
          udkx = (gamasafe*(acntr(isid)*vdx(i)+acntrbar*vdxold(i))
     &        + acntrbar*uparoverB(i)*bx(i)
     &        + acntr(isid)*uparoBpred(i)*bxpred(i))
          udky = (gamasafe*(acntr(isid)*vdy(i)+acntrbar*vdyold(i))
     &        + acntrbar*uparoverB(i)*by(i)
     &        + acntr(isid)*uparoBpred(i)*bypred(i))
          udkz = (gamasafe*(acntr(isid)*vdz(i)+acntrbar*vdzold(i))
     &        + acntrbar*uparoverB(i)*bz(i)
     &        + acntr(isid)*uparoBpred(i)*bzpred(i))

c          if (i==1 .and. ipmin==1) print*, "in getveff",
c     1       vdx(i),uparoverB(i),bx(i),udkx,alpha(1),uxp(1)
c        calculate interpolated velocity
          uxeff(i)=alpha(i)*uxp(i) + alphabar(i)*udkx
          uyeff(i)=alpha(i)*uyp(i) + alphabar(i)*udky
          uzeff(i)=alpha(i)*uzp(i) + alphabar(i)*udkz
        enddo
       endif

      endif
      return
      end

c==========================================================================
      subroutine getvdrift(pgroup,npd,is,ipmin,x,y,z)
c     Assumes E, B, grad B^2 are pre-computed
      use DKInterp
      use DKInterptmp
      use ParticleGroupmodule
      use Particles,Only: vdxoldpid,vdyoldpid,vdzoldpid
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,is,isid,ip,ipmin,ipmin1,ipmax
      real(kind=8):: x(npd),y(npd),z(npd)

      integer(ISZ):: i
      real(kind=8):: bgradbbx,bgradbby,bgradbbz,coeff
      real(kind=8):: newcoef,vbxupdated,vbyupdated,vbzupdated

      real(kind=8),pointer:: vdxold(:),vdyold(:),vdzold(:)
      real(kind=8),pointer:: gaminv(:),bx(:),by(:),bz(:),ex(:),ey(:),ez(:)

      ipmin1=ipmin-1
      ipmax=ipmin1+npd
      vdxold => pgroup%pid(ipmin:ipmax,vdxoldpid)
      vdyold => pgroup%pid(ipmin:ipmax,vdyoldpid)
      vdzold => pgroup%pid(ipmin:ipmax,vdzoldpid)
      gaminv => pgroup%gaminv(ipmin:ipmax)
      bx => pgroup%bx(ipmin:ipmax)
      by => pgroup%by(ipmin:ipmax)
      bz => pgroup%bz(ipmin:ipmax)
      ex => pgroup%ex(ipmin:ipmax)
      ey => pgroup%ey(ipmin:ipmax)
      ez => pgroup%ez(ipmin:ipmax)

      isid = pgroup%sid(is-1) + 1

c     first, vgradb:
c     calculate sum of curvature and grad-B drifts, assuming vacuum B field
c     formula: vd = [m/(2 q B^4)](upar^2 + uperp^2/2)(Bvec x grad B^2)
c     first get gradb; this also sets up B_j and |B| in the database.

c    call getvperpparsq() -- no longer needed
c      call getgradbsq(npd,is,ipmin,x,y,z) -- only do it when fields
c      are recomputed, so moved into setfields

c     comps of B x grad B**2:
c      print*, "IN GETVDRIFT, is,isid,ipmin1", is, isid, ipmin1
c      print*, "ey,ez,by,bz,bsqi",ey(1),ez(1),by(1),bz(1),bsqi(1)
      do i=1,npd
        bgradbbx = by(i)*grdbsq(1,i) - bz(i)*grdbsq(3,i)
        bgradbby = bz(i)*grdbsq(2,i) - bx(i)*grdbsq(1,i)
        bgradbbz = bx(i)*grdbsq(3,i) - by(i)*grdbsq(2,i)

c      common coeff:
c      Note should be gamma*m*v**2 ~ m u**2/gamma where u is momentum/mass
        coeff = 0.25*m_over_q(isid)*bsqi(i)**2*(usq(i) + uparsq(i))*gaminv(i)

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
        if (impinterp == 0 .or. ipredcor == 1) then
           vdx(i) = vex(i) + vbx(i)
           vdy(i) = vey(i) + vby(i)
           vdz(i) = vez(i) + vbz(i)
        else
           vdx(i) = 0.5*vex(i) + vbx(i)
           vdy(i) = 0.5*vey(i) + vby(i)
           vdz(i) = 0.5*vez(i) + vbz(i)
           vdxold(i) = vex(i) + vbx(i)
           vdyold(i) = vey(i) + vby(i)
           vdzold(i) = vez(i) + vbz(i)
        endif

c       Need to store old vdrifts if averaging drifts
        if (ipredcor == 0 .and. impinterp == 0) then
c       components of magnetic drift using timestep-n fields and
c       updated velocities. This will determine vdrift_old for corrector
           newcoef = 0.25*m_over_q(isid)*bsqi(i)**2*(usq_new(i) + uparsq_new(i))
     &          *gaminv(i)
           vbxupdated = newcoef*bgradbbx
           vbyupdated = newcoef*bgradbby
           vbzupdated = newcoef*bgradbbz
           vdxold(i) = vex(i) + vbxupdated
           vdyold(i) = vey(i) + vbyupdated
           vdzold(i) = vez(i) + vbzupdated
         endif
      enddo

c      if (ipmin == 1) print*, "in getvdrift",
c    1      vex(1),ey(1),bz(1),ez(1),by(1),bsqi(1),vdx(1)
c      print*,"drifts, x,y,z,ex, bx,bsqi ", x(1),y(1),z(1),ex(1),bx(1),bsqi(1)
c      print*,"vex,vbx,vdx,vdy,vdz",vex(1),vbx(1),vdx(1),vdy(1),vdz(1)
c      print*, "ey,ez,by,bz,bsqi,vex",ey(5),ez(5),by(5),bz(5),bsqi(5),vex(5)
      return
      end

c=========================================================================
c      subroutine getvpol(is,np,dt)
c
c      use DKInterptmep
c      integer(ISZ):: is,,isid,np,dt
c      real(kind=8):: dvexdt,dveydt,dvezdg,qbsqi
c      integer(ISZ):: i
c      dti=1./dt
c      isid = pgroup%sid(is-1) + 1
c      do i=1,np
c        dvexdt=(vex(i)-vexold(i))*dti
c        dveydt=(vey(i)-veyold(i))*dti
c        dvezdt=(vez(i)-vezold(i))*dti
c        qbsqi=m_over_q(isid)*bsqi(i)
c        vpolx(i) = (by(i)*dvezdt-bz(i)*dveydt)*qbsqi
c        vpoly(i) = (bz(i)*dvexdt-bx(i)*dvezdt)*qbsqi
c        vpolz(i) = (bx(i)*dveydt-by(i)*dvexdt)*qbsqi
c      enddo
c      return
c      end  

c=========================================================================
      subroutine setfields(pgroup,npd,is,ipmin,x,y,z,ex,ey,ez,bx,by,bz,dt)
c     sets ej and bj.
      use ParticleGroupmodule
      use DKInterp
      use DKInterptmp
      use InGen, only: ifeears
      use InPart
      use Beam_acc
      use Z_arrays
      use Picglb
      use GlobalVars
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,is,ipmin
      real(kind=8):: x(npd),y(npd),z(npd),dt
      real(kind=8):: ex(npd),ey(npd),ez(npd)
      real(kind=8):: bx(npd),by(npd),bz(npd)
      integer(ISZ):: i,jsid
      real(kind=8):: dtr
      real(kind=8):: bendres(nparpgrp), bendradi(nparpgrp)
      jsid = pgroup%sid(is-1)
c     print*,"*****SETTING FIELDS"
c      print*,"pre-fetche3d1,x,ey,ez = ",x(1),ey(1),ez(1)
      ex(:) = 0.
      ey(:) = 0.
      ez(:) = 0.
      call fetche3dfrompositions(jsid,pgroup%ndts(is-1),npd,x,y,z,ex,ey,ez)
c      call fetche3d1(pgroup,npd,x,y,z)
c      print*,"after fetche3d,x,ey,ez = ",x(1),ey(1),ez(1)
c      call fetche3d(pgroup,ipmin,npd,is)
c     add in ears and uniform focusing E field pieces
      call othere3d(npd,x,y,z,zbeam,zimax,zimin,straight,ifeears,eears,
     &                     eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,dbxdy,dbydx,
     &                     ex,ey,ez,bx,by,bz)
c      call geteb(pgroup,npd,is,ipmin,x,y,z,ex,ey,ez,bx,by,bz)
      dtr=.5*dt
      bx(:) = 0.
      by(:) = 0.
      bz(:) = 0.
      call exteb3d(npd,x,y,z,pgroup%uzp(ipmin),pgroup%gaminv(ipmin),-dtr,dtr,
     &             bx,by,bz,ex,ey,ez,pgroup%sm(is),
     &             pgroup%sq(is),bendres,bendradi,gammabar,dt)
c     Now get the gradients of B**2
      call getgradbsq(npd,is,x,y,z,bx,by,bz)
      return
      end

c=========================================================================
      subroutine getgradbsq(npd,is,x,y,z,bx,by,bz)
c      calculate grad B**2;  also calculates B^2, 1/B^2, and 
c      interpolation parameter alpha and its complement
c      Note shape of gradbsq is (3,npd), and the components
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
c      use InGen
      use BSQGRADdata
      integer(ISZ):: npd,is
      real(kind=8):: x(npd),y(npd),z(npd)
      real(kind=8):: bx(npd),by(npd),bz(npd)
      integer(ISZ):: i,ip
      real(kind=8):: twobbrsqi,bsq_safe
c      print*,"entering getalpha from getgradbsq"
c
c     Assumes  E and B are already set, either via setfields
c      or via calls in padvnc3d.
c     Calculate B**2 
      do i=1,npd
          bsq(i)=bx(i)**2 + by(i)**2 + bz(i)**2
          bsq_safe=bsq(i)+dksmall
          bsqi(i) = 1./bsq_safe
      enddo

c     Now can calculate the interpolation parameters alpha
c     Do this outside this routine.

      if (igradb == 1) then
c       igradb == 1, get data from lookup table.
c       Must zero the array first, as applybsqgrad accumulates
        grdbsq(:,1:npd)=0.
        call applybsqgrad(npd,x,y,npd,z,.false.,grdbsq)
      elseif (igradb == 2) then
        do i=1,npd
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
        call applybsqgrad(npd,x,y,npd,z,.false.,grdbsq)
        do i=1,npd
c         Note grdbsq(1,:) directly set by applybsqgrad
          twobbrsqi = 2.*bsq(i)/(x(i)**2+y(i)**2+dksmall)
          grdbsq(2,i)=twobbrsqi*x(i)
          grdbsq(3,i)=twobbrsqi*y(i)
        enddo
      endif
      return
      end

c========================================================================
      subroutine getalphas(pgroup,npd,isid,ipmin,dtuse)
      use DKInterp
      use DKInterptmp
      use ParticleGroupmodule
c      use InGen
c      use Picglb
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,is,isid,ipmin,ipmin1,ipmax
      integer(ISZ):: i,ip
      real(kind=8):: omegacesq,omegadtsq,dtuse
      real(kind=8),pointer:: gaminv(:)

      ipmin1=ipmin-1
      ipmax=ipmin1+npd
      gaminv => pgroup%gaminv(ipmin:ipmax)

      if (abs(ipalpha) == 1) then
       do i=1,npd
c         calculate interpolation parameter, 0 if magnetized, 1 if not
          omegacesq = gaminv(i)*gaminv(i)*qovermsq(isid)*bsq(i)
          omegadtsq = omegacesq*dtuse**2
          alpha(i) = usealphacalc(isid)/sqrt(1.+alphcoef*omegadtsq)  
     &       + notusealphcalc(isid)
          alphabar(i)=1.-alpha(i)
       enddo
      elseif (abs(ipalpha) == 2) then
       do i=1,npd
c         calculate interpolation parameter, 0 if magnetized, 1 if not
          omegacesq = gaminv(i)*gaminv(i)*qovermsq(isid)*bsq(i)
          omegadtsq = omegacesq*dtuse**2
          alpha(i) = usealphacalc(isid)/(1.+alphcoef*omegadtsq)  
     &       + notusealphcalc(isid)
          alphabar(i)=1.-alpha(i)
       enddo
      else
       do i=1,npd
c         calculate interpolation parameter, 0 if magnetized, 1 if not
          omegacesq = gaminv(i)*gaminv(i)*qovermsq(isid)*bsq(i)
          omegadtsq = omegacesq*dtuse**2
          alpha(i) = usealphacalc(isid)/(1.+alphcoef*omegadtsq)**palpha  
     &       + notusealphcalc(isid)
          alphabar(i)=1.-alpha(i)
       enddo
      endif
      return
      end

c=========================================================================
      subroutine geteb(pgroup,npd,is,ipmin,x,y,z,dt,ex,ey,ez,bx,by,bz)
c      clean fetch of electric and magnetic fields.
c      bx,by,bz written in bx(i), by(i), bz(i).  Returns ex, ey, ez
      use ParticleGroupmodule
      use DKInterptmp
      use DKInterp
c      use InGen
      use Beam_acc
      use GlobalVars
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,ipmin,is
      real(kind=8):: x(npd),y(npd),z(npd),dt
      real(kind=8):: ex(npd),ey(npd),ez(npd)
      real(kind=8):: bx(npd),by(npd),bz(npd)
      integer(ISZ):: i,ip
      real(kind=8):: dtr,bsq_safe
      real(kind=8):: bendres(nparpgrp), bendradi(nparpgrp)
      dtr=.5*dt
      bx(:)=0.;by(:)=0.;bz(:)=0.
      call exteb3d(npd,x,y,z,pgroup%uzp(ipmin),pgroup%gaminv(ipmin),-dtr,dtr,
     &         bx,by,bz,ex,ey,ez,pgroup%sm(is),
     &         pgroup%sq(is),bendres,bendradi,gammabar,dt)
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
c    Temporarily sets top.pgroup.xp=x, etc, calls fetch3d, and resets
c    top.pgroup.xp
      subroutine fetche3d1(pgroup,npd,x,y,z)
      use ParticleGroupmodule
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,i
      real(kind=8):: x(npd),y(npd),z(npd),xpd(npd),ypd(npd),zpd(npd)
c     save values of xp, yp, zp and set xp,yp,zp to the x,y,z arrays
      do i=1,npd
         xpd(i)=pgroup%xp(i)
         ypd(i)=pgroup%yp(i)
         zpd(i)=pgroup%zp(i)
         pgroup%xp(i)=x(i)
         pgroup%yp(i)=y(i)
         pgroup%zp(i)=z(i)
      enddo
      call fetche3d(pgroup,1,npd,1)
c     reset xp,yp,zp
      do i=1,npd
         pgroup%xp(i)=xpd(i)
         pgroup%yp(i)=ypd(i)
         pgroup%zp(i)=zpd(i)
      enddo
      return
      end

c===========================================================================
c OLDSETUP: sets up storage to store prior-step quantities
c===========================================================================
      subroutine oldsetup()
c     Sets up storage in the pid array and makes pointers to its rows for
c     xold, uxold, etc.
      use Particles,Only: xoldpid,yoldpid,zoldpid,uxoldpid,uyoldpid,uzoldpid,
     &  vdxoldpid,vdyoldpid,vdzoldpid,uparBopid,bxoldpid,byoldpid,bzoldpid,
     &  uparoBpredpid,expredpid,eypredpid,ezpredpid,bxpredpid,bypredpid,
     &  bzpredpid
      use ParticleGroupmodule
      use DKInterp
      integer(ISZ):: nextpid
c      call gchange("Particles",0)
      print*,"executing oldsetup with impinterp = ",impinterp
      if (vdxoldpid == 0) vdxoldpid = nextpid()
      if (vdyoldpid == 0) vdyoldpid = nextpid()
      if (vdzoldpid == 0) vdzoldpid = nextpid()
      if (impinterp == 0) then
         if (xoldpid == 0) xoldpid = nextpid()
         if (yoldpid == 0) yoldpid = nextpid()
         if (zoldpid == 0) zoldpid = nextpid()
         if (uxoldpid == 0) uxoldpid = nextpid()
         if (uyoldpid == 0) uyoldpid = nextpid()
         if (uzoldpid == 0) uzoldpid = nextpid()
         if (bxoldpid == 0) bxoldpid = nextpid()
         if (byoldpid == 0) byoldpid = nextpid()
         if (bzoldpid == 0) bzoldpid = nextpid()
         if (uparBopid == 0) uparBopid = nextpid()
      else
         if (bxpredpid == 0) bxpredpid = nextpid()
         if (bypredpid == 0) bypredpid = nextpid()
         if (bzpredpid == 0) bzpredpid = nextpid()
         if (expredpid == 0) expredpid = nextpid()
         if (eypredpid == 0) eypredpid = nextpid()
         if (ezpredpid == 0) ezpredpid = nextpid()
         if (uparoBpredpid == 0) uparoBpredpid = nextpid()
      endif
      return
      end
c=============================================================================
c ARRAY UTILITIES
c==========================================================================
      subroutine swapxarrays(pgroup,npd,is,ipmin)
c     swapxarrays, swaps old, new particle arrays for predictor-corrector
c     Specifically, on predictor step,
c     this routine interchanges the contents of xp and
c     xold, so that xp can be set to the previous (n-1) timestep x's
c     while xold is set to the current timestep's x's for use in the
c     next timestep.
c     On corrector step it just copies xold into xp
c     At t=0, it just copies xp into xold.
      use DKInterp
      use ParticleGroupmodule
      use Particles,Only: xoldpid,yoldpid,zoldpid
      use Picglb
      type(ParticleGroup):: pgroup
      real(kind=8):: xsave,ysave,zsave
      integer(ISZ):: npd,ipmin,is,i
      if (ipredcor == 0 .and. it > 1) then
       do i = ipmin,ipmin+npd-1
         xsave=pgroup%xp(i)
         ysave=pgroup%yp(i)
         zsave=pgroup%zp(i)
         pgroup%xp(i)=pgroup%pid(i,xoldpid)
         pgroup%yp(i)=pgroup%pid(i,yoldpid)
         pgroup%zp(i)=pgroup%pid(i,zoldpid)
         pgroup%pid(i,xoldpid)=xsave
         pgroup%pid(i,yoldpid)=ysave
         pgroup%pid(i,zoldpid)=zsave
       enddo
       return
      endif
      if (ipredcor > 0) then
         do i = ipmin,ipmin+npd-1
            pgroup%xp(i) =pgroup%pid(i,xoldpid)
            pgroup%yp(i)=pgroup%pid(i,yoldpid)
            pgroup%zp(i)=pgroup%pid(i,zoldpid)
         enddo
         return
       endif
       if (it == 1) then
         do i = ipmin,ipmin+npd-1
            pgroup%pid(i,xoldpid) =pgroup%xp(i)
            pgroup%pid(i,yoldpid)=pgroup%yp(i)
            pgroup%pid(i,zoldpid)=pgroup%zp(i)
         enddo
       endif
      return
      end
c===========================================================================
      subroutine storeoldu(pgroup,is)
c     stores old uxp array in uxold
c     This routine does nothing if implicit solver
      use ParticleGroupmodule
      use Particles,Only: uxoldpid,uyoldpid,uzoldpid
      use DKInterp,Only: impinterp
      type(ParticleGroup):: pgroup
      integer (ISZ):: i,is
      if (impinterp == 0) then
         do i = pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
            pgroup%pid(i,uxoldpid) = pgroup%uxp(i)
            pgroup%pid(i,uyoldpid) = pgroup%uyp(i)
            pgroup%pid(i,uzoldpid) = pgroup%uzp(i)
         enddo
      endif
      return
      end
c==========================================================================
      subroutine swapxarr2(pgroup,npd,is,ipmin)
c     saves current x as xold
      use Particles,Only: xoldpid,yoldpid,zoldpid
      use ParticleGroupmodule
      use DKInterp
      use Picglb
      type(ParticleGroup):: pgroup
      real(kind=8):: xsave,ysave,zsave
      integer(ISZ):: npd,ipmin,is,i
      if (ipredcor == 0) then
        do i = ipmin,ipmin+npd-1
          pgroup%pid(i,xoldpid)=pgroup%xp(i)
          pgroup%pid(i,yoldpid)=pgroup%yp(i)
          pgroup%pid(i,zoldpid)=pgroup%zp(i)
        enddo
      else
        do i = ipmin,ipmin+npd-1
          pgroup%xp(i)=pgroup%pid(i,xoldpid)
          pgroup%yp(i)=pgroup%pid(i,yoldpid)
          pgroup%zp(i)=pgroup%pid(i,zoldpid)
        enddo
      endif
      return
      end
c=========================================================================
      subroutine set_polarization(rho,nx,nz,dx,dz,xmin,zmin)
c     Sets epsilon and epszfac if there is an intpolated species.
c     For now this only works with rz solver
      use InPart
      use InGen
      use DKInterp
      use Constant
      integer(ISZ):: nx,nz,i,ngrd,is,isid
      real(kind=8):: dx,dz,xmin,zmin
      real(kind=8):: xgrid(0:nx,0:nz),zgrid(0:nx,0:nz),ygrid(0:nx,0:nz)
      real(kind=8):: Bx(0:nx,0:nz),By(0:nx,0:nz),Bz(0:nx,0:nz)
      real(kind=8):: Bsq(0:nx,0:nz),alphabar(0:nx,0:nz)
      real(kind=8):: omegapsq(0:nx,0:nz),omegacsq(0:nx,0:nz),rho(0:nx,0:nz)
      logical(ISZ):: check
c     Construct a grid of x and z values, so we can get
c     values of B on it.
c      print*, "IN SET_POL"
 
c     create grid for epsilon and epszfac
      epsilon2d%xmin=0.
      epsilon2d%ymin=0.
      epsilon2d%dx=dx
      epsilon2d%dy=dz
      epsilon2d%nx=nx
      epsilon2d%ny=nz
c     shifted by floors following convention of Dave for density
      call setupgrid2dtype(epsilon2d,check)
      if (.not. check) then
          call kaboom("set_polarization: ERROR: inconsistent data for setting up epsilon")
          return
      endif
      epszfac2d%xmin=0.
      epszfac2d%ymin=0.
      epszfac2d%dx=dx
      epszfac2d%dy=dz
      epszfac2d%nx=nx
      epszfac2d%ny=nz
      call setupgrid2dtype(epszfac2d,check)
      if (.not. check) then
          call kaboom("set_polarization: ERROR: inconsistent data for setting up epszfac")
          return
      endif

      ngrd = (nx+1)*(nz+1)
      do i=0,nz
         zgrid(:,i) = zmin + i*dz/nz
      enddo
      do i = 0,nx
         xgrid(i,:) = xmin + i*dx/nx
      enddo
      ygrid(0:nx,0:nz) = 0.
      Bx = 0.
      By = 0.
      Bz = 0.
c     get B components on this grid
      call applybgrd(ngrd,xgrid,ygrid,ngrd,zgrid,.false.,Bx,By,Bz)
      Bsq=Bx*Bx+By*By+Bz*Bz
c     calculate interpolation parameter, 0 if magnetized, 1 if not
c     For now, nonrelativistic.  Need to think about appropriate
c     generalization for fields
c     For now check to see if ANY species is interpolated; if so
c     calculate epsilon for the last such species.  Need to 
c     generalize.
      do is = 1,ns
       if (interpdk(is) == 1) then
        omegacsq = qovermsq(is)*Bsq
c        print*, "qovermsq,Bsq,dt,omegacsq", qovermsq(is),Bsq(0,0),dt,
c     &       omegacsq(0,0)
        alphabar = 1. - usealphacalc(is)/sqrt(1.+alphcoef*omegacsq*dt**2)  
     &       - notusealphcalc(is)
c       We should introduce partial density fractions.  But for our test
c       just take the total density if a species is interpolated
        if (lnonlinpol) then
          omegapsq = rho/(m_over_q(is)*eps0)
c         temporarily borrow epsilon for polarization piece
c         print*, "omegacsq, omegapsq ",omegacsq(0,0),omegapsq(0,0)
          epsilon2d%grid = omegapsq*alphabar/omegacsq
        else
          epsilon2d%grid = omegpomegcsq(is)*alphabar
        endif
        epszfac2d%grid = 1.+epsilon2d%grid*(1.-Bz*Bz/Bsq)
        epsilon2d%grid = 1. + epsilon2d%grid*(1.-Bx*Bx/Bsq)
        epszfac2d%grid = epszfac2d%grid/(epsilon2d%grid+SMALLPOS)
       endif
c       print*, "epsilon2d = ", epsilon2d(1,1)
      enddo
         
      return
      end
c
c===================================================================
c      subroutine restoreEBFromPred(pgroup)
c      use ParticleGroupmodule
c      use Particles,Only: expredpid,eypredpid,ezpredpid
c      type(ParticleGroup):: pgroup
c      integer(ISZ):: js,i1,i2
c      do js = 1,pgroup%ns
c         i1 = pgroup%ins(js)
c         i2 = pgroup%ins(js)+pgroup%nps(js)-1
c         print*,"js,i1,i2,epredpid's = ",js,i1,i2,expredpid,eypredpid,ezpredpid
c         pgroup%ex(i1:i2) = pgroup%pid(i1:i2,expredpid)
c         pgroup%ey(i1:i2) = pgroup%pid(i1:i2,eypredpid)
c         pgroup%ez(i1:i2) = pgroup%pid(i1:i2,ezpredpid)
c      enddo
c      return
c      end
c===================================================================
      subroutine setTotalE(pgroup,npd,is,ipmin,dtloc)
c     sets total E at particle positions
      use ParticleGroupmodule
      use picglb,Only: zbeam,lresetparticlee
      use Z_arrays
      use Inpart
      use InGen
      use Beam_acc
      type(ParticleGroup):: pgroup
      integer(ISZ):: npd,is,isid,ipmin
      real(kind=8):: dtloc,halfdt
      real(kind=8),allocatable:: bx(:),by(:),bz(:),bendres(:),bendradi(:)
      real(kind=8),pointer:: xp(:),yp(:),zp(:),uzp(:),gaminv(:)
      real(kind=8),pointer:: ex(:),ey(:),ez(:)
      real(kind=8),pointer:: sm(:),sq(:),sw(:)

      xp => pgroup%xp
      yp => pgroup%yp
      zp => pgroup%zp
      uzp => pgroup%uzp
      gaminv => pgroup%gaminv
      ex => pgroup%ex
      ey => pgroup%ey
      ez => pgroup%ez

      sm => pgroup%sm
      sq => pgroup%sq
      sw => pgroup%sw

      halfdt = 0.5*dtloc
      allocate(bendres(npd), bendradi(npd),bx(npd),by(npd),bz(npd))
      bendres = 0.
      bendradi = 0.

      lresetparticlee = .true.
      call fetche3d(pgroup,ipmin,npd,is)
c     add external fields.  Unfortunately this adds B fields, which we
c     don't need to do again.  So give it dummy B fields to set.
      call othere3d(npd,xp(ipmin),yp(ipmin),zp(ipmin),
     &              zbeam,zimax,zimin,straight,ifeears,eears,
     &              eearsofz,dzzi,nzzarr,zzmin,dedr,dexdx,deydy,dbdr,
     &              dbxdy,dbydx,ex(ipmin),ey(ipmin),ez(ipmin),
     &              bx,by,bz)
c    --- Set quad, dipole E and B
      call exteb3d(npd,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &              gaminv(ipmin),-halfdt,halfdt,
     &              bx,by,bz,
     &              ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),
     &              bendres,bendradi,gammabar,dt)
      
      lresetparticlee = .false.
      deallocate(bendres, bendradi,bx,by,bz)
      return
      end


