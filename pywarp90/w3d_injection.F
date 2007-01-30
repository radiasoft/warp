#include "top.h"
c=============================================================================
c@(#) File w3d_injection.F, version $Revision: 1.86 $, $Date: 2007/01/30 01:37:51 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package W3D of code WARP
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
      subroutine inject3d(itask,pgroup)
      use ParticleGroupmodule
      use GlobalVars
      use Subtimers3d
      use Constant
      use Beam_acc
      use InGen
      use InGen3d
      use InPart
      use InPart3d
      use InMesh3d
      use Particles,Only: spid,wpid,tpid,rpid,ssn
      use Picglb
      use InjectVars
      use InjectVars3d
      use Picglb3d
      use Setpwork3d
      use Fields3dParticles
#ifdef MPIPARALLEL
      use Parallel
#endif
      integer(ISZ):: itask
      type(ParticleGroup):: pgroup

c This routine does the injection of particles.  It has three main parts, the
c first cleans up the particle arrays and makes sure that there is enough room
c for the new particles.  The second is the coding for constant current
c injection.  The third part is for space-charge limited injection.

      integer(ISZ):: is,jsid,ij,ii,ip,nn,ix,iy,iz,ipmin,ixx,iyy,i1,i2
      integer(ISZ):: inj_is
      real(kind=8):: rnpinjct,rnn,t
      real(kind=8):: rr,vavez,vnorm,zs,wxx,wyy,aa,az,wz,xm,ym,zm,rsign,dr,dri
      real(kind=8):: qoverm,const,dxi,dyi,dzi,vznorm,rp,ai,bi
      real(kind=8):: clightsqi
      real(kind=8):: zinj,ainj,ainji,ainjmini,binj,binji,binjmini
      real(kind=8):: apinj,bpinj,rinj,rinji,xpinj,ypinj,vzinj,ainjp1,binjp1
      real(kind=8):: sphericalcorrection
      real(kind=8):: vtx,vty,larmorfreq
      real(kind=8):: fulldt_s
      real(kind=8):: zmid,zleni,ztilt,gamma
      real(kind=8):: p1x,p1y,p1z,p2x,p2y,p2z,w1x,w1y,w1z,area,circum
      real(kind=8):: ellipseperimeter
      integer(ISZ):: i1x,i1y,i1z
      real(kind=8):: rnrev,wranf,wrandom,wrandomgauss
      integer(ISZ):: ith,nti
      real(kind=8):: cl_const(0:pgroup%ns-1),jton(0:pgroup%ns-1),vzconst(0:pgroup%ns-1)
      real(kind=8):: zmin_tmp, zmax_tmp, r_tmp, theta_tmp, xt_tmp, yt_tmp
      real(kind=8),pointer:: ppvz(:)
      integer(ISZ):: spreadx, spready
      integer(ISZ):: ixmin,ixmax,iymin,iymax
      logical(ISZ):: linj_cell
      real(kind=8),pointer:: xx(:),yy(:),zz(:),ux(:),uy(:),uz(:),id(:),gg(:)
      real(kind=8),pointer:: ex(:),ey(:),ez(:),bx(:),by(:),bz(:)
      real(kind=8),pointer:: bendres(:),bendradi(:),mask(:)
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (inject <= 0) return

      call setuppgroup(pgroup,pgroup%ns)

      dxi = 1./inj_dx
      dyi = 1./inj_dy
      dzi = 1./inj_dz
      if (lrelativ) clightsqi = 1./clight**2

c Print warning if vbeamfrm is not zero
      if (abs(vbeamfrm) > 0.) then
       call remark("WARNING: vbeamfrm is not zero and injection is being used.")
       call remark("A core dump will result if the emitting surface moves off of the field grid.")
      endif

c Clear out the lost particles.
      if (itask == 1) call clearpart(pgroup,-1,2)

c Set inject z min and max.
      call updateinj_zmminmax()

c     --- Set the potential near the emitting surface.
c     --- Note that these routines are now called twice (once for each itask) which
c     --- may not be necessary in most cases, but it does ensure that the phi arrays
c     --- are setup properly if the fields are change by the user.
      call getinj_phi()
      call gettinj_phi()

c     --- Set factor so correct number of particles are loaded when
c     --- symmetry is used.  With 2-fold, particles need to be loaded
c     --- on both sides of the symmetry axis.  With 4-fold, particles
c     --- need to be loaded in each of the four quadrants.
      if(solvergeom==XYZgeom) then
        spreadx = 1
        spready = 1
      elseif(solvergeom==RZgeom) then
        spreadx = 1
        if(l_inj_rz) then
          spready = 0
        else
          spready = 1
        end if
      elseif(solvergeom==XZgeom) then
        spreadx = 1
        spready = 0
      elseif(solvergeom==Zgeom) then
        spreadx = 0
        spready = 0
      endif

c If itask = 1, zero npinject. This is done here since npinject is used both
c by normal and transverse space-charge limited injection. Both sections
c add onto npinject. Only do this if actually injecting particles.
      if (itask == 1) npinject = 0

c The check of whether or not there is enough space for the injected particles
c is now done in each injection method.

      if (inject == 1) then
c Injects particles with constant (or user varied) current.  It does a split
c leap frog advance over partial time steps, starting with all of the particles
c in the plane at zinject.  The particles are advanced uniformly distribuated
c fractions of a timestep to fill the injection region.  The E and B fields are
c first gathered (assuming that injection is not being done inside a bend) then
c a half advance is done on the velocity over the partial timestep.  The
c positions are then advanced for the full amount of the partial timestep and
c the particles are loaded onto the rho grid.  The routine then exits so the
c fields can be solved for.  The routine is called again to complete the
c velocity advance using the new fields including the injected particles.  The
c velocity ends up at time level it-1/2 to match the time level of the rest of
c the particles.
c
c One caveat is the gathering of the E fields (in the call to sete3d).  If the
c injection plane is at the left hand of the grid, the E gather will use the
c plane iz = -1 which is not necessarily correct.  For example, in injection
c off of a conducting plane, phi at iz=-1 is the same as phi at iz=0, this
c gives a E field that is too small (by about a factor of 2).

        if (itask == 1) then
c Start injected particles advance.  Particles end up with position
c at current time level and velocity at a time level of one half of the time
c each particle was advanced from the plane of injection (i.e. all different).

c         --- Zero out temp array.
          npinjtmp = 0
          npinje_s = 0
          inj_prev = 0.

c         --- loop over injection sources
          do ij=1,ninject
            ainj = ainject(ij)
            binj = binject(ij)
            ainji = 1./ainject(ij)
            binji = 1./binject(ij)
            zinj = zinject(ij)
            apinj = apinject(ij)
            bpinj = bpinject(ij)
            xpinj = xpinject(ij)
            ypinj = ypinject(ij)

c         --- Loop over species.
          do is=1,pgroup%ns
            if (pgroup%sid(is-1) == -1) pgroup%sid(is-1) = is -1 
            jsid = pgroup%sid(is-1)
            inj_is = min(jsid+1,inj_ns)
            vzinj = vzinject(ij,jsid+1)
            if (.not. pgroup%ldts(is-1)) cycle
            fulldt_s    = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)

            if (l_inj_user_particles) then

              call execuser("generateuserparticlesforinjection")
              rnn = npgrp
              nn = npgrp
              if (nn == 0) cycle

c             --- Check if there is enough room for the new injected particles.
              call chckpart(pgroup,is,nn,0,.true.)

            else

              rnn = rnpinje_s(jsid+1)*finject(ij,jsid+1)
c             --- Add a random number to the number of particles injected
c             --- so that the average number of particles injected is
c             --- correct.  For example, if rnn < 1., without the
c             --- addition of the random number, no particles would ever
c             --- be injected.  With the random number, particles will be
c             --- injected and but the average number will be less than 1.
              rnn = rnn + wranf()
              nn = int(rnn)
              if(nn==0) cycle

c             --- Check if there is enough room for the new injected particles.
              call chckpart(pgroup,is,nn,0,.true.)

c             --- Make sure there is enough temp space
              if (nn > npgrp) then
                npgrp = nn
                call gallot("Setpwork3d",0)
              endif

c             --- set random numbers
c             --- rt(ip) used to hold f(r**2)=constant
c             --- tt(ip) used to hold f(theta)=constant if injection
c             --- is over an elliptical surface
              do ip = 1,nn
                rt(ip) = wrandom(xrandom,injctcnt+ip-1,dig1,1,1)
                rt(ip) = rt(ip)*(1.-(ainjmin(ij)*ainji)**2) +
     &                              (ainjmin(ij)*ainji)**2
                tt(ip) = wrandom(xrandom,injctcnt+ip-1,dig2,1,1)
              enddo

              if (l_inj_regular) then
                do ip=1,nn
                  rt(ip) = (ip-1.)/(nn-1.)
                enddo
              endif

c             --- Transform the positions into a Gaussian if distrbtn is GA0
              if (distrbtn == "GA0") then
                do ip = 1,nn
                  rt(ip) = -0.5*log(rt(ip))
                enddo
              endif

c             --- or rt(ip) used to hold f(r)~(h+(1-h)r^2)
              if (hollow == 2) then
                do ip = 1,nn
                  rt(ip) = (1 + hollow_h)*rt(ip)/
     &                 (hollow_h + sqrt(hollow_h**2 + (1. - hollow_h**2)*rt(ip)))
                enddo
              endif

c             --- Get velocity distribution.
              if (distrbtn == "K-V") then
c               --- For the K-V load, need another random number.
                do ip = 1,nn
                  t = wrandom(vtrandom,injctcnt+ip-1,dig3,1,1)
                  rr = sqrt(1. - rt(ip))
                  uxt(ip) = rr*cos(2.*pi*t)
                  uyt(ip) = rr*sin(2.*pi*t)
                enddo
              else
c               --- Default is to use a Semi-Gaussian distribution
c               --- Note that the loops must be broken up in this way since the
c               --- routine wrandomgauss must be called in the correct order.
                do ip = 1,nn
                  uxt(ip) = wrandomgauss(vtrandom,injctcnt+ip-1,dig3,dig4,
     &                                   1,1,.true.)
                enddo
                do ip = 1,nn
                  uyt(ip) = wrandomgauss(vtrandom,injctcnt+ip-1,dig5,dig6,
     &                                   1,1,.true.)
                enddo
              endif

c             --- Axial thermal spread
              do ip = 1,nn
                uzt(ip) = wrandomgauss(vzrandom,injctcnt+ip-1,dig7,dig8,
     &                                   1,1,.true.)
              enddo

c             --- Initialize positions as rectangular or elliptical depending on linj_rectangle
              if(linj_rectangle) then
                do ip=1,nn
                  xt(ip) = ainj*(2.*rt(ip)-1.0)
                  yt(ip) = binj*(2.*tt(ip)-1.0)
                enddo
              else if (l_inj_rz) then
                do ip=1,nn
                  tt(ip) = 2.*pi*tt(ip)
                  xt(ip) = ainj*rt(ip)*cos(tt(ip))
                  yt(ip) = binj*rt(ip)*sin(tt(ip))
                enddo
              else
                do ip=1,nn
c                 --- Convert r**2 into r and tt into a random angle
                  rt(ip) = sqrt(rt(ip))
                  tt(ip) = 2.*pi*tt(ip)
                  xt(ip) = ainj*rt(ip)*cos(tt(ip))
                  yt(ip) = binj*rt(ip)*sin(tt(ip))
                enddo
              endif
            endif

c           --- Precalculate larmor frequency if needed. The code assumes
c           --- that the beam was injected in a field free region.
            if (linj_includebz0) then
              larmorfreq = -0.5*pgroup%sq(is)*bz0/pgroup%sm(is)/gammabar
            endif

            spreadx = 1
            spready = 1
            if (l_inj_rz) then
              spreadx = 1
              spready = 0
            endif

c           --- Increment the particle counter for the digit reversed seed
            injctcnt = injctcnt + nn

            do ip=1,nn
c             --- find injection surface axial location via interpolation
c             --- The absolute values are needed for the cases of symmetry.
              ixx = spreadx*(abs(xt(ip) - inj_xmmin(ij))*dxi)
              iyy = spready*(abs(yt(ip) - inj_ymmin(ij))*dyi)
              wxx = spreadx*(abs(xt(ip) - inj_xmmin(ij))*dxi - ixx)
              wyy = spready*(abs(yt(ip) - inj_ymmin(ij))*dyi - iyy)
              zs = inj_grid(ixx        ,iyy        ,ij)*(1. - wxx)*(1. - wyy) +
     &             inj_grid(ixx+spreadx,iyy        ,ij)*      wxx *(1. - wyy) +
     &             inj_grid(ixx        ,iyy+spready,ij)*(1. - wxx)*      wyy  +
     &             inj_grid(ixx+spreadx,iyy+spready,ij)*      wxx *      wyy
c             --- Only save particle if it is within the injection region.
              if (inj_zmmin <= zs+zinj .and. zs+zinj < inj_zmmax) then
                ii = pgroup%ins(is) - 1 - npinjtmp(is,ij)
c               --- Clear out any old data from pid
                pgroup%pid(ii,:) = 0.
c               --- Save the particle data
                pgroup%xp(ii) = xt(ip)
                pgroup%yp(ii) = yt(ip)
                pgroup%zp(ii) = zs
                vtx = vzinj*0.5*emitx_s(jsid+1)*ainji + vthperp_s(jsid+1)
                vty = vzinj*0.5*emity_s(jsid+1)*binji + vthperp_s(jsid+1)
                pgroup%uxp(ii) = vzinj*(apinj*pgroup%xp(ii)*ainji + xpinj) + vtx*uxt(ip)
                pgroup%uyp(ii) = vzinj*(bpinj*pgroup%yp(ii)*binji + ypinj) + vty*uyt(ip)
                if(l_inj_addtempz_abs) then
                  pgroup%uzp(ii) = vzinj + abs(0.5*vthz_s(jsid+1)*uzt(ip))
                else
                  pgroup%uzp(ii) = vzinj + 0.5*vthz_s(jsid+1)*uzt(ip)
                endif
c               --- Save the particle weight
                if (l_inj_rz) pgroup%pid(ii,wpid) = 2.*rt(ip)
                if (linj_includebz0) then
                  rt(ip) = sqrt((pgroup%xp(ii)+xinject(ij))**2 + (pgroup%yp(ii)+yinject(ij))**2)
                  tt(ip) = atan2(pgroup%yp(ii)+yinject(ij),pgroup%xp(ii)+xinject(ij))
                  pgroup%uxp(ii) = pgroup%uxp(ii) - larmorfreq*rt(ip)*sin(tt(ip))
                  pgroup%uyp(ii) = pgroup%uyp(ii) + larmorfreq*rt(ip)*cos(tt(ip))
                endif
                if (lrelativ) then
                  pgroup%gaminv(ii) = sqrt(1.-
     &     (pgroup%uxp(ii)**2+pgroup%uyp(ii)**2+pgroup%uzp(ii)**2)*clightsqi)
                  gamma = 1./pgroup%gaminv(ii)
                  pgroup%uxp(ii) = pgroup%uxp(ii)*gamma
                  pgroup%uyp(ii) = pgroup%uyp(ii)*gamma
                  pgroup%uzp(ii) = pgroup%uzp(ii)*gamma
                else
                  pgroup%gaminv(ii) = 1.
                endif
c               --- Save the ID of the injection source.
                pgroup%pid(ii,injpid) = ij
c               --- Save SSN 
                if (spid>0) then
                  pgroup%pid(ii,spid) = ssn
                  ssn = ssn+1
                end if
c               --- save time of creation if tpid > 0
                if(tpid>0) pgroup%pid(ii,tpid) = time + fulldt_s
c               --- save radius of creation if rpid > 0
                if(rpid>0) then
                  pgroup%pid(ii,rpid) = sqrt(pgroup%xp(ii)**2+pgroup%yp(ii)**2)
                endif
c               --- increment particle counter
                npinjtmp(jsid+1,ij) = npinjtmp(jsid+1,ij) + 1
                inj_prev(ixx        ,iyy        ,ij,inj_is) =
     &          inj_prev(ixx        ,iyy        ,ij,inj_is) + (1.-wxx)*(1.-wyy)
                inj_prev(ixx+spreadx,iyy        ,ij,inj_is) =
     &          inj_prev(ixx+spreadx,iyy        ,ij,inj_is) + (   wxx)*(1.-wyy)
                inj_prev(ixx        ,iyy+spready,ij,inj_is) =
     &          inj_prev(ixx        ,iyy+spready,ij,inj_is) + (1.-wxx)*(   wyy)
                inj_prev(ixx+spreadx,iyy+spready,ij,inj_is) =
     &          inj_prev(ixx+spreadx,iyy+spready,ij,inj_is) + (   wxx)*(   wyy)
              endif
            enddo

c           --- Get actual number of particles injected and ipmin
            nn = npinjtmp(jsid+1,ij)

c           --- Sum npinje_s to get the total number of particles injected.
            npinje_s(jsid+1) = npinje_s(jsid+1) + nn
            npinject     = npinject     + nn

c           --- make injected particles live particles
            pgroup%ins(is) = pgroup%ins(is) - nn
            pgroup%nps(is) = pgroup%nps(is) + nn

c           --- Transform particles to lab frame
            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is) + nn - 1
            if (nn > 0) then
c             --- The check if nn > 0 is needed since if npmax was zero
c             --- and no particles were injected then npmax stays zero
c             --- and the particle arrays are unallocated.
              call inj_transform(nn,pgroup%xp(i1:i2),pgroup%yp(i1:i2),
     &                           pgroup%zp(i1:i2),1,ij,1,1)
              call inj_transform(nn,pgroup%uxp(i1:i2),pgroup%uyp(i1:i2),
     &                           pgroup%uzp(i1:i2),1,ij,1,0)
            endif

c           --- get E self-field at initial positions
            if (linj_efromgrid) then
              call fetche3dfrompositions(jsid,pgroup%ndts(is-1),nn,
     &                                   pgroup%xp(i1:i2),
     &                                   pgroup%yp(i1:i2),
     &                                   pgroup%zp(i1:i2),at,bt,apt)
              call fetchb3dfrompositions(jsid,nn,pgroup%xp(i1:i2),
     &                                   pgroup%yp(i1:i2),
     &                                   pgroup%zp(i1:i2),xt,yt,zt)
            else
              at = 0.
              bt = 0.
              apt = 0.
              xt = 0.
              yt = 0.
              zt = 0.
            endif
            call inj_sete3d(nn,pgroup%xp(i1:i2),pgroup%yp(i1:i2),
     &                      pgroup%zp(i1:i2),
     &                      pgroup%pid(i1:i2,injpid),at,bt,apt)

c           --- get external fields at initial positions
c           --- (Assumes that injection will never be done inside a bend.)
c           --- The reusing of variable makes for very confusing naming!
            do ip=1,nn
              uxt(ip) = 0.
              uyt(ip) = 1.
            enddo
            do ipmin = pgroup%ins(is),pgroup%ins(is)+nn-1,nparpgrp
              ip = min(nparpgrp, pgroup%ins(is)+nn-ipmin)
              ii = ipmin - pgroup%ins(is) + 1
              call exteb3d(ip,pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                     pgroup%zp(ipmin),pgroup%uzp(ipmin),
     &                     pgroup%gaminv(ipmin),0.,fulldt_s*0.5,
     &                     xt(ii),yt(ii),zt(ii),at(ii),bt(ii),apt(ii),
     &                     pgroup%sm(is),pgroup%sq(is),uxt(ii),uyt(ii),
     &                     1.,fulldt_s)
            enddo

c           --- Set fractional time to advance particles
            rnpinjct = 1./dvnz(real(nn,kind=8))
            do ip=1,nn
              bpt(ip) = (ip - .5)*rnpinjct*fulldt_s
            enddo

c           --- Add in a velocity tilt (this is dependent on the fractional
c           --- time step and so must be done here).
            zmid = 0.5*(zimax_s(jsid+1) - zimin_s(jsid+1)) - vzinj*fulldt_s*(it+1)
            zleni = 1./dvnz(zimax_s(jsid+1) - zimin_s(jsid+1))
            do ip=1,nn
              ii = pgroup%ins(is) - 1 + ip
              ztilt = bpt(ip)*vzinj
              pgroup%uzp(ii) = pgroup%uzp(ii) -
     &                         vzinj*vtilt_s(jsid+1)*(zmid+ztilt)*zleni
            enddo

c           --- do half velocity advance with E fields
            call epusht3d(nn,pgroup%uxp(i1:i2),pgroup%uyp(i1:i2),
     &                    pgroup%uzp(i1:i2),
     &             at(1),bt(1),apt(1),pgroup%sq(is),pgroup%sm(is),bpt(1),0.5)

c           --- Advance relativistic Gamma factor
            call gammaadv(nn,pgroup%gaminv(i1:i2),pgroup%uxp(i1:i2),
     &                    pgroup%uyp(i1:i2),
     &                    pgroup%uzp(i1:i2),
     &                    gamadv,lrelativ)

c           --- do half velocity advance with B fields
            call bpusht3d(nn,pgroup%uxp(i1:i2),pgroup%uyp(i1:i2),
     &                    pgroup%uzp(i1:i2),
     &                    pgroup%gaminv(i1:i2),xt(1),yt(1),zt(1),
     &                    pgroup%sq(is),pgroup%sm(is),
     &                    bpt(1),0.5,ibpush)

c           --- do full position advance
            call xpusht3d(nn,pgroup%xp(i1:i2),pgroup%yp(i1:i2),
     &                    pgroup%zp(i1:i2),
     &                    pgroup%uxp(i1:i2),pgroup%uyp(i1:i2),
     &                    pgroup%uzp(i1:i2),
     &                    pgroup%gaminv(i1:i2),bpt(1))

c           --- Calculate Gamma inverse
c           if (lrelativ) then
c             do ip=pgroup%ins(is),pgroup%ins(is)+nn-1
c               pgroup%gaminv(ip) = 1./sqrt(1. +
c    & (pgroup%uxp(ip)**2+pgroup%uyp(ip)**2+pgroup%uzp(ip)**2)*clightsqi)
c             enddo
c           else
c             do ip=ins(is),ins(is)+nn-1
c               gaminv(ip) = 1.
c             enddo
c           endif

c           --- setrho is now done by the setrho call in padvnc3d.

c         --- end of loop over species
          enddo

c         --- end of loop over injection sources
          enddo

        elseif (itask == 2) then
c Do second part of constant current injection: get new E fields (this is
c done after field solve including injected particles), synchronize velocity
c with position, and then move velocity one half timestep back to match time
c level of rest of particles.

c         loop over species and injection sources
          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
            if (.not. pgroup%ldts(is-1)) cycle
            fulldt_s    = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)

            ipmin = pgroup%ins(is)
            nn = 0
            do ij=1,ninject
              ipmin = ipmin + nn
c             --- Get number of particles injected
              nn = npinjtmp(jsid+1,ij)
              if(nn==0) cycle

c             --- Make sure there is enough temp space
              if (nn > npgrp) then
                npgrp = nn
                call gallot("Setpwork3d",0)
              endif

c             --- calculate new E self-fields
              if (linj_efromgrid) then
                call fetche3dfrompositions(jsid,pgroup%ndts(is-1),nn,
     &                                     pgroup%xp(ipmin),
     &                                     pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                                     at,bt,apt)
                call fetchb3dfrompositions(jsid,nn,pgroup%xp(ipmin),
     &                                     pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                                     xt,yt,zt)
              else
                at = 0.
                bt = 0.
                apt = 0.
                xt = 0.
                yt = 0.
                zt = 0.
              endif
              call inj_sete3d(nn,pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                        pgroup%zp(ipmin),pgroup%pid(ipmin,injpid),
     &                        at,bt,apt)

c             --- get external fields at current positions
              do ip=1,nn
                uxt(ip) = 0.
                uyt(ip) = 1.
              enddo
              do ip = 1,nn,nparpgrp
                ii = min(nparpgrp, nn+1-ip)
                call exteb3d(ii,pgroup%xp(ip+ipmin-1),pgroup%yp(ip+ipmin-1),
     &                       pgroup%zp(ip+ipmin-1),
     &                       pgroup%uzp(ip+ipmin-1),pgroup%gaminv(ip+ipmin-1),
     &                       -fulldt_s*0.5,0.,
     &                       xt(ip),yt(ip),zt(ip),at(ip),bt(ip),apt(ip),
     &                       pgroup%sm(is),pgroup%sq(is),uxt(ip),uyt(ip),
     &                       1.,fulldt_s)
              enddo

c             --- Set fractional time to advance particles
              rnpinjct = 1./dvnz(real(nn,kind=8))
              do ip=1,nn
                bpt(ip) = (ip - .5)*rnpinjct*fulldt_s
              enddo

c             --- complete B advance
              call bpusht3d(nn,pgroup%uxp(ipmin),pgroup%uyp(ipmin),
     &                      pgroup%uzp(ipmin),pgroup%gaminv(ipmin),
     &                      xt(1),yt(1),zt(1),pgroup%sq(is),pgroup%sm(is),
     &                      bpt(1),0.5,ibpush)

c             --- complete the E advance
              call epusht3d(nn,pgroup%uxp(ipmin),pgroup%uyp(ipmin),
     &                      pgroup%uzp(ipmin),
     &                      at(1),bt(1),apt(1),pgroup%sq(is),pgroup%sm(is),
     &                      bpt(1),0.5)

c             --- Advance relativistic Gamma factor
              call gammaadv(nn,pgroup%gaminv(ipmin),pgroup%uxp(ipmin),
     &                      pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     &                      gamadv,lrelativ)

c             --- Now, move velocites back one half a step
c             --- first half of a backward B advance
              call bpush3d(nn,pgroup%uxp(ipmin),pgroup%uyp(ipmin),
     &                     pgroup%uzp(ipmin),
     &                     pgroup%gaminv(ipmin),xt(1),yt(1),zt(1),
     &                     pgroup%sq(is),pgroup%sm(is),
     &                     -0.5*fulldt_s,ibpush)

c             --- then half of a backward E advance
              call epush3d(nn,pgroup%uxp(ipmin),pgroup%uyp(ipmin),
     &                     pgroup%uzp(ipmin),
     &                     at(1),bt(1),apt(1),pgroup%sq(is),pgroup%sm(is),
     &                     -0.5*fulldt_s)

c             --- Advance relativistic Gamma factor
              call gammaadv(nn,pgroup%gaminv(ipmin),pgroup%uxp(ipmin),
     &                      pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     &                      gamadv,lrelativ)

c           --- end of loop over species and injection sources
            enddo
          enddo

c Constant current injection is now complete for this time step.
        endif

c------------------------------------------------------------------------
c Space-Charge Limited injection

c Particles are initially placed on the emitting surface and given a
c time of emission uniformly distributed between 0 and dt. They are given
c an initial normal velocity equal to the Child-Langmuir velocity at dt
c minus the time of emission. This initial velocity gives the particles
c a small kick to get them moving which is needed since the normal field
c is zero at the surface, the particle starting location.  They are then
c advanced a partial isochronous leapfrog timestep, advancing dt minus the
c time of emission.  After the particles are injected, the position and
c velocity are at the same time level so, when this injection is used,
c every step is done with the split-leap frog algorithm.
c
c The voltage drop is calculated between the emitting surface and a
c secondary surface which is parallel too (is concentric too in the
c spherical case) the emitting surface.
c
c For inject == 2
c Space-charge limited injection, version 1. The number of particles injected is
c calculated from the field near the source using the Child-Langmuir relation.
c The number of particles loaded is J*dt*dx*dy/echarge/sw where J is the
c current density as calculated from the Child-Langmuir result.
c This algorithm works well for steady state flow, but does not give a
c good result for time-dependent flow.
c
c For inject == 3
c Space charge limited injection based off of the Gauss's Law.  This is
c the scheme to use with time-dependent injection.

      elseif (inject == 2 .or. inject == 3) then

        if (itask == 1) then
c         --- Zero npinje_s.
          npinje_s = 0
          npinjtmp = 0

c         --- Copy inj_np to inj_prev, and zero inj_np
          inj_prev = inj_np
          inj_np = 0.
          inj_npactual = 0.

c         --- Calculate the charge density on the surface of the emitter.
          if (inject == 3) then
            if (solvergeom==XYZgeom .or. solvergeom==AMRgeom) then
              call inj_setrho3d(pgroup,inj_dz,l2symtry,l4symtry)
            elseif(solvergeom==RZgeom) then
c             --- When using the RZ solver, inj_rho is forced to be
c             --- four-fold symmetric.
              call inj_setrho3d(pgroup,inj_dz,.false.,.true.)
            elseif(solvergeom==Zgeom) then
              call inj_setrho3d_z(pgroup,inj_dz,nz)
            endif
          endif

c         --- loop over injection sources
          do ij=1,ninject

c           --- Set some temporaries.
            ainj = ainject(ij)
            binj = binject(ij)
            ainji = 1./ainject(ij)
            binji = 1./binject(ij)
            ainjmini = 1./dvnz(ainjmin(ij))
            binjmini = 1./dvnz(binjmin(ij))
            zinj = zinject(ij)
            apinj = apinject(ij)
            bpinj = bpinject(ij)
            xpinj = xpinject(ij)
            ypinj = ypinject(ij)
            rinji = 1./rinject(ij)

c           --- The potential drop in front of the surface is now
c           --- calculated by a call to getinj_phi from padvnc3d.

            do is=1,pgroup%ns
              jsid = pgroup%sid(is-1)
              vzinj = vzinject(ij,jsid+1)

c             --- Don't bother doing anything if finject is 0.
              if (.not. pgroup%ldts(is-1) .or. finject(ij,jsid+1) == 0.) cycle
              fulldt_s    = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)

c             --- The factor to convert current density to particle number.
c             --- The absolute value is taken of sq so the factor is always
c             --- positive.
              if(l_inj_rz) then
                jton(jsid) = pi*ainj*inj_dx*fulldt_s/
     &                     abs(pgroup%sq(is))/pgroup%sw(is)
              else
                jton(jsid) = inj_dx*inj_dy*fulldt_s/
     &                     abs(pgroup%sq(is))/pgroup%sw(is)
              end if

c             --- This constant terms includes corrections for injection off
c             --- of a concentric spheres.  See I. Langmuir,
c             --- K. Blodgett, "Currents Limited by Space Charge Between
c             --- Concentric Spheres", PhysRev, 1924.
              if (inject == 2) then
                if (linj_spherical) then
                  zs = (inj_dz*inj_d(ij))*rinji
                  sphericalcorrection = (1. + 1.6*zs + 2.06*zs**2)
                else
                  sphericalcorrection = 1.
                endif
                cl_const(jsid) = 4./9.*eps0*
     &                   (2.*abs(pgroup%sq(is))/pgroup%sm(is))**.5/
     &                   ((inj_dz*inj_d(ij))**2*sphericalcorrection)*jton(jsid)
              endif

c             --- The constant terms in the expression for the normal
c             --- velocity from the Child-Langmuir solution. Note that this
c             --- is only strictly correct for planar emission.
              vzconst(jsid) = 2./9.*sqrt(2.)*
     &                      (abs(pgroup%sq(is))/pgroup%sm(is))**1.5/
     &                      (inj_dz*inj_d(ij))**2

            enddo

c           --- Load particles one grid cell at a time.
            do iy=0,inj_ny
              do ix=0,inj_nx

c               --- Mesh location relative to source center
                xm = inj_xmmin(ij) + ix*inj_dx
                ym = inj_ymmin(ij) + iy*inj_dy

c               --- Check whether the cell is within the injection area.
c               --- This depends on whether the source is elliptical or
c               --- rectangular.
                ai = (ainj+2.*inj_dx)
                bi = (binj+2.*inj_dy)
                if(linj_rectangle) then
                  linj_cell = (abs(xm) <= ai) .and. (abs(ym) <= bi)
                else
                  linj_cell = ((xm*bi)**2 + (ym*ai)**2) <= (ai*bi)**2
                endif

                if (.not. linj_cell) cycle

                zmin_tmp = minval(inj_grid(max(0,ix-1):min(inj_nx,ix+1),
     &                                     max(0,iy-1):min(inj_ny,iy+1),ij)
     &                                     + zinj)
                zmax_tmp = maxval(inj_grid(max(0,ix-1):min(inj_nx,ix+1),
     &                                     max(0,iy-1):min(inj_ny,iy+1),ij)
     &                                     + zinj)
                if(zmin_tmp>=inj_zmmax.or.inj_zmmin>zmax_tmp) cycle
c               --- Loop over the species
                do is=1,pgroup%ns
                  jsid = pgroup%sid(is-1)
                  inj_is = min(jsid+1,inj_ns)

c                 --- Skip this species if the fraction is zero.
c                 if (finject(ij,jsid+1) == 0.) cycle
                  if (.not. pgroup%ldts(is-1) .or. finject(ij,jsid+1) == 0.) cycle

c                 --- number of particles injected in grid cell
                  if (inject == 2) then
c                   --- Child-Langmuir
                    rnn = cl_const(jsid)*abs(inj_phi(ix,iy,ij))**1.5/
     &                    cos(inj_angl(ix,iy,ij))
                    rnn = sign(rnn,pgroup%sq(is)*inj_phi(ix,iy,ij))
                  elseif (inject == 3) then
c                   --- Gauss's law
                    rnn = (eps0*inj_phi(ix,iy,ij)/
     &     (inj_dz*abs(inj_d(ij)))-0.5*inj_rho(ix,iy,ij)*inj_dz)*
     &     inj_dx*inj_dy/pgroup%sq(is)/pgroup%sw(is)/cos(inj_angl(ix,iy,ij))
                  endif

c                 --- Make sure that rnn is >= 0.
c                 --- Don't cycle at this point even if rnn is zero since
c                 --- averaging with the number from the previous step may
c                 --- give a positive value.
                  if (rnn < 0.) rnn = 0.

c                 --- Apply ad-hoc scaling factor
                  rnn = rnn*inj_f(ij)

c                 --- Average the number of particles for this step with
c                 --- that of the previous step. This helps the relaxation
c                 --- toward a steady state, especially in the Egun style
c                 --- iterative mode.
                  rnn = (inj_param*rnn +
     &                  (1. - inj_param)*inj_prev(ix,iy,ij,inj_is))

c                 --- Force injected current to be between the range
c                 --- jmininj to jmaxinj.  The absolute values are used
c                 --- so that this works with either positive ions
c                 --- and electron.
                  rnn = min(rnn,jmaxinj(ij)*jton(jsid))
                  rnn = max(rnn,jmininj(ij)*jton(jsid))

c                 --- Scale the number of particles by the fraction of the
c                 --- total for this species.
                  rnn = rnn*finject(ij,jsid+1)

c                 --- Save the number of particles injected
                  inj_np(ix,iy,ij,inj_is) = inj_np(ix,iy,ij,inj_is) + rnn

c                 --- Add a random number to the number of particles injected
c                 --- so that the average number of particles injected is
c                 --- correct.  For example, if rnn < 1., without the
c                 --- addition of the random number, no particles would ever
c                 --- be injected.  With the random number, particles will be
c                 --- injected and but the average number will be less than 1.
                 if(.not. l_inj_regular) then
                   rnn = rnn + wranf()
c                  --- divide number of particles injected by two on axis
c                  --- when using radial injection
                   if(l_inj_rz .and. ix==0) rnn = 0.5*rnn
                 else
                   rnn = rnn * inj_area(ix,iy,ij)
                 end if

c                 --- Convert to integer.
                  if(.not.l_inj_regular) then
                    nn = int(rnn)
                  else
                    nn = 1
                  end if

c                 --- Only do following work if particles are to be injected
                  if (nn <= 0 .or. rnn<=0.) cycle

c                 --- Make sure there is room for more particles.  Estimate
c                 --- total number of particles to be injected from the
c                 --- source. The second call makes it absolutely sure
c                 --- that there is enough room available in case the
c                 --- first estimate is wrong.
                  if(.not.l_inj_regular) then
                    if(.not.l_inj_rz) then
                      call chckpart(pgroup,is,max(nn+npinje_s(jsid+1),
     &                              int(nn*ainj*binj*pi*dxi*dyi)),0,.true.)
                    else
                      call chckpart(pgroup,is,max(nn+npinje_s(jsid+1),
     &                              int(nn*ainj*dxi)),0,.true.)
                    end if
                  else
                    call chckpart(pgroup,is,max(nn+npinje_s(jsid+1),
     &                            int((inj_nx+1)*(inj_ny+1))),0,.true.)
                  end if

c                 --- loop to load new particles
                  do ip=1,nn
                    ii = pgroup%ins(is) - npinje_s(jsid+1) - 1
c                   --- increment random number counter
                    injctcnt = injctcnt + 1

c                   --- calculate x, y and z of new particle
c                   --- if within injection source load it, else skip it
                    if(.not.l_inj_regular) then
                      pgroup%xp(ii) = wrandom(xrandom,injctcnt,dig1,1,1)
                      pgroup%xp(ii) = xm +
     &                           spreadx*(pgroup%xp(ii) - .5)*inj_dx*inj_xwide
                      if(.not. l_inj_rz) then
                        pgroup%yp(ii) = wrandom(xrandom,injctcnt,dig2,1,1)
                        pgroup%yp(ii) = ym +
     &                           spready*(pgroup%yp(ii) - .5)*inj_dy*inj_ywide
                      else
                        pgroup%yp(ii)=0.
                      endif
                    else
                      pgroup%xp(ii)=xm
                      pgroup%yp(ii)=ym
                    end if
c                   --- Check if particle is within the rectangular or
c                   --- elliptical annulus with outer major and minor
c                   --- size of ainject and binject and inner size of
c                   --- ainjmin and binjmin.
                    if(linj_rectangle) then
                      if (abs(pgroup%xp(ii)) > ainj .or.
     &                    abs(pgroup%yp(ii)) > binj) cycle
                      if (abs(pgroup%xp(ii)) < ainjmin(ij) .or.
     &                    abs(pgroup%yp(ii)) < binjmin(ij)) cycle
                    else
                      if ((pgroup%xp(ii)*ainji)**2 +
     &                    (pgroup%yp(ii)*binji)**2 > 1.) cycle
                      if ((pgroup%xp(ii)*binjmin(ij))**2 +
     &                    (pgroup%yp(ii)*ainjmin(ij))**2 <
     &                    ainjmin(ij)*binjmin(ij)) cycle
                    endif

c                   --- Clear out any old data from pid. This is done here
c                   --- before pid is used anywhere.
                    pgroup%pid(ii,:) = 0.

                    if(l_inj_regular) then
                      if(l_inj_rz) then
                        r_tmp=abs(pgroup%xp(ii))
                        theta_tmp = wrandom(xrandom,injctcnt,dig2,1,1)
                        theta_tmp = 2.*pi*theta_tmp !rnrev(injctcnt,dig2)
                        pgroup%xp(ii) = r_tmp*cos(theta_tmp)
                        pgroup%yp(ii) = r_tmp*sin(theta_tmp)
                        if(ix==0) then
                          pgroup%pid(ii,wpid) = 0.25*inj_dx*ainji
                        else
                          pgroup%pid(ii,wpid) = 2.*r_tmp*ainji
                        end if
                        pgroup%pid(ii,wpid) = pgroup%pid(ii,wpid)*rnn
                      else
                        pgroup%pid(ii,wpid) = rnn
                      endif
                    else
c                     --- if inject only radially, use variable weights and
c                     --- assign pid
                      if(l_inj_rz) then
                        r_tmp=abs(pgroup%xp(ii))
                        theta_tmp = wrandom(xrandom,injctcnt,dig2,1,1)
                        theta_tmp = 2.*pi*theta_tmp !rnrev(injctcnt,dig2)
                        pgroup%xp(ii) = r_tmp*cos(theta_tmp)
                        pgroup%yp(ii) = r_tmp*sin(theta_tmp)
                        pgroup%pid(ii,wpid) = 2.*r_tmp*ainji
                      end if
                    end if

c                   --- Find location of particle relative to source
c                   --- grid. This is needed for zp and az, as well as
c                   --- for vnorm below.
                    ixx = spreadx*(abs(pgroup%xp(ii) - inj_xmmin(ij))*dxi)
                    iyy = spready*(abs(pgroup%yp(ii) - inj_ymmin(ij))*dyi)
                    wxx = spreadx*(abs(pgroup%xp(ii) - inj_xmmin(ij))*dxi - ixx)
                    wyy = spready*(abs(pgroup%yp(ii) - inj_ymmin(ij))*dyi - iyy)

                    if(l_inj_exact .or. l_inj_rz) then
c                     --- calculate injection surface axial location exactly
                      zm = (pgroup%xp(ii)**2+pgroup%yp(ii)**2)/
     &                     (abs(rinject(ij)) +
     &    sqrt(max(0.,rinject(ij)**2 - pgroup%xp(ii)**2 - pgroup%yp(ii)**2)))
                      if (rinject(ij) < 0.) zm = -zm
                      if ((pgroup%xp(ii)**2+pgroup%yp(ii)**2) > rinject(ij)**2)
     &                  zm = rinject(ij)
                      pgroup%zp(ii) = zm + inj_addfdz*inj_dz
                      az = asin(sqrt(pgroup%xp(ii)**2 + pgroup%yp(ii)**2)/rinject(ij))
                    else
c                     --- find injection surface axial location by interpolation
                      pgroup%zp(ii)=inj_grid(ixx        ,iyy        ,ij)*(1.-wxx)*(1.-wyy) +
     &                              inj_grid(ixx+spreadx,iyy        ,ij)*    wxx *(1.-wyy) +
     &                              inj_grid(ixx        ,iyy+spready,ij)*(1.-wxx)*    wyy  +
     &                              inj_grid(ixx+spreadx,iyy+spready,ij)*    wxx *    wyy  +
     &                              inj_zstart(ij)

c                     --- find injection angle due to curvature by interpolation
                      az = inj_angl(ixx        ,iyy        ,ij)*(1. - wxx)*(1. - wyy) +
     &                     inj_angl(ixx+spreadx,iyy        ,ij)*      wxx *(1. - wyy) +
     &                     inj_angl(ixx        ,iyy+spready,ij)*(1. - wxx)*      wyy  +
     &                     inj_angl(ixx+spreadx,iyy+spready,ij)*      wxx *      wyy
                    end if

c                   --- Only inject particles within the injection region.
                    if (pgroup%zp(ii)+zinj < inj_zmmin .or. pgroup%zp(ii)+zinj >= inj_zmmax) cycle

c                   --- Calculate transverse angle
                    aa = atan2(pgroup%yp(ii),pgroup%xp(ii))

c                   --- Calculate and save the fraction of time step this
c                   --- particle is to be advanced. Also equal to dt minus
c                   --- the time of emission, all divided by dt.
c                   --- Note that the number is stored as a fraction between
c                   --- zero and one since pid does double duty, storing
c                   --- also the number of the injection source. Multiplying
c                   --- by dt would lose too many digits when added to an
c                   --- integer.
c                   --- The fractions are chosen using the digit reversed
c                   --- number generator and not particle number since they
c                   --- are correlated with transverse position.
                    if(.not. l_inj_regular) then
                      pgroup%pid(ii,injpid) = rnrev(injctcnt,dig5)
                    else
                      pgroup%pid(ii,injpid) = rnrev(injctcnt,dig1)
                    end if

c                   --- save time of creation if tpid > 0
                    if(tpid>0) pgroup%pid(ii,tpid) = time + fulldt_s

c                   --- reverse sign of pid(ii,tpid) if temperature added
c                   --- after delay.
c                   --- negative pid is then used to track particles on which
c                   --- to add temperature.
                    if(l_inj_delay_temp) then
                      pgroup%pid(ii,tpid) = -pgroup%pid(ii,tpid)
                    endif

c                   --- save radius of creation if rpid > 0
                    if(rpid>0) then
                      pgroup%pid(ii,rpid) =
     &                            sqrt(pgroup%xp(ii)**2+pgroup%yp(ii)**2)
                    endif

c                   --- Get the normal velocity for this particle. This is
c                   --- based on the Child-Langmuir solution, given a
c                   --- voltage drop, vnorm, and a time since emission,
c                   --- pid(ii,injpid)*dt.
                    vnorm = inj_phi(ixx        ,iyy        ,ij)*(1. - wxx)*(1. - wyy) +
     &                      inj_phi(ixx+spreadx,iyy        ,ij)*      wxx *(1. - wyy) +
     &                      inj_phi(ixx        ,iyy+spready,ij)*(1. - wxx)*      wyy  +
     &                      inj_phi(ixx+spreadx,iyy+spready,ij)*      wxx *      wyy
                    vznorm = vzconst(jsid)*abs(vnorm)**1.5*
     &                       (pgroup%pid(ii,injpid)*fulldt_s)**2
                    if (inj_d(ij) < 0.) vznorm = -vznorm

c                   --- Set transverse coordinates. The emittance term
c                   --- is probably not correct since it depends on vbeam.
c                   --- The correct way to add a transverse thermal spread
c                   --- is via vthperp_s. Also add in the transverse
c                   --- component of the normal velocity and the specified
c                   --- vbeam (which is also taken to be normal).
c                   --- Like the emittance term, the term with ap (and bp)
c                   --- are probably not correct.
                    pgroup%uxp(ii) = spreadx*(vzinj*
     &                               (apinj*pgroup%xp(ii)*ainji + xpinj)
     &                      - (vznorm + vzinj)*sin(az)*cos(aa))
                    pgroup%uyp(ii) = spready*(vzinj*
     &                               (bpinj*pgroup%yp(ii)*binji + ypinj)
     &                      - (vznorm + vzinj)*sin(az)*sin(aa))
                    if(.not. l_inj_delay_temp) then
                      vtx = vzinj*0.5*emitx_s(jsid+1)*ainji + vthperp_s(jsid+1)
                      vty = vzinj*0.5*emity_s(jsid+1)*binji + vthperp_s(jsid+1)
                      pgroup%uxp(ii) = pgroup%uxp(ii) + spreadx*vtx*
     &                    wrandomgauss(vtrandom,injctcnt,dig3,dig4,1,1,.false.)
                      pgroup%uyp(ii) = pgroup%uyp(ii) + spready*vty*
     &                    wrandomgauss(vtrandom,injctcnt,dig5,dig6,1,1,.false.)
                    end if
c                   --- Include the user specified vbeam (assumed to be a
c                   --- normal), the thermal velocity, and the term from
c                   --- the Child-Langmuir solution.
c                   --- SMALLPOS is also added so that the field gathering
c                   --- routines don't treat the newly injected particles
c                   --- as lost particles when the other terms are zero.
                    pgroup%uzp(ii) = vzinj*cos(az) + vznorm*cos(az) + SMALLPOS
                    if(.not. l_inj_delay_temp) then
                      if(l_inj_addtempz_abs) then
                        pgroup%uzp(ii) = pgroup%uzp(ii) + abs(vthz_s(jsid+1)*
     &                   wrandomgauss(vzrandom,injctcnt,dig7,dig8,1,1,.false.))
                      else
                        pgroup%uzp(ii) = pgroup%uzp(ii) + vthz_s(jsid+1)*
     &                    wrandomgauss(vzrandom,injctcnt,dig7,dig8,1,1,.false.)
                      end if
                    end if

c                   --- Calculate gamma inverse and reset the particle
c                   --- massless momentums properly.
                    if (lrelativ) then
                      pgroup%gaminv(ii) = sqrt(1.-
     &      (pgroup%uxp(ii)**2+pgroup%uyp(ii)**2+pgroup%uzp(ii)**2)*clightsqi)
                      gamma = 1./pgroup%gaminv(ii)
                      pgroup%uxp(ii) = pgroup%uxp(ii)*gamma
                      pgroup%uyp(ii) = pgroup%uyp(ii)*gamma
                      pgroup%uzp(ii) = pgroup%uzp(ii)*gamma
                    else
                      pgroup%gaminv(ii) = 1.
                    endif

c                   --- Save the ID of the injection source. Note that pid
c                   --- is doing double duty, saving the injection source
c                   --- number and the fraction of the timestep the
c                   --- particle is advanced.
                    pgroup%pid(ii,injpid) = pgroup%pid(ii,injpid) + ij

c                   --- Save SSN 
                    if (spid>0) then
                      pgroup%pid(ii,spid) = ssn
                      ssn = ssn+1
                    end if

c                   --- increment particle counter
                    npinje_s(jsid+1) = npinje_s(jsid+1) + 1
                    npinjtmp(jsid+1,ij) = npinjtmp(jsid+1,ij) + 1
                    inj_npactual(ix,iy,ij,inj_is) = inj_npactual(ix,iy,ij,inj_is) + 1.

                  enddo
c                 --- end of loop over particles
                enddo
c               --- end loop over species
              enddo
c             --- end loop over ix
            enddo
c           --- end loop over iy

c           --- Transform particles to lab frame
            do is=1,pgroup%ns
              jsid = pgroup%sid(is-1)
              ii = pgroup%ins(is) - npinje_s(is)
              nn = npinjtmp(jsid+1,ij)
              if (nn > 0) then
c               --- The check if nn > 0 is needed since if npmax was zero
c               --- and no particles were injected then npmax stays zero
c               --- and the particle arrays are unallocated.
                call inj_transform(nn,pgroup%xp(ii),pgroup%yp(ii),
     &                             pgroup%zp(ii),1,ij,1,1)
                call inj_transform(nn,pgroup%uxp(ii),pgroup%uyp(ii),
     &                             pgroup%uzp(ii),1,ij,1,0)
              endif
            enddo

          enddo
c         --- end of loop over injection sources

c         --- Call transverse particle scraping routine to force
c         --- removal of particles outside the grid.
          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
            if (.not. pgroup%ldts(is-1)) cycle
            ii = pgroup%ins(is) - npinje_s(jsid+1)
            if (npinje_s(jsid+1) > 0) then
              call stckxy3d(npinje_s(jsid+1),pgroup%xp(ii),xmmax,xmmin,dx,
     &                      pgroup%yp(ii),ymmax,
     &                      ymmin,dy,pgroup%zp(ii),zmmin,dz,pgroup%uxp(ii),
     &                      pgroup%uyp(ii),pgroup%uzp(ii),pgroup%gaminv(ii),
     &                      zgrid,zbeam,l2symtry,l4symtry,absorb,.false.)
              call processlostpart(pgroup,is,clearlostpart,time,zbeam)
            endif
          enddo

          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
            if (.not. pgroup%ldts(is-1)) cycle
c           --- Sum npinje_s to get the total number of particles injected.
            npinject = npinject + npinje_s(jsid+1)

c           --- make injected particles live particles
            pgroup%ins(is) = pgroup%ins(is) - npinje_s(jsid+1)
            pgroup%nps(is) = pgroup%nps(is) + npinje_s(jsid+1)
          enddo

c         --- Now, advance the particles off of the emitting surface using
c         --- isochronous leapfrog.  Each particle has it's own time step
c         --- size, uniformly distributed between 0 and dt.
          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
c           --- Skip this if no particles of this species were injected.
            if (.not. pgroup%ldts(is-1)) cycle
            fulldt_s    = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)
            if(npinje_s(jsid+1)==0) cycle

c           --- Make sure that there is enough room in the temporary arrays.
            if (npinje_s(jsid+1) > npgrp) then
              npgrp = npinje_s(jsid+1)
              call gchange("Setpwork3d",0)
            endif

c           --- Make equivalences so names make sense below
            ex => at
            ey => bt
            ez => apt
            bx => bpt
            by => xct
            bz => xpct
            bendres => yct
            bendradi => ypct

c           --- Zero the B field which is accumulated below
            ex = 0.
            ey = 0.
            ez = 0.
            bx = 0.
            by = 0.
            bz = 0.

c           --- Assume that injection will never be done inside a bend.
            bendres = 0.
            bendradi = 1.

c           --- Create equivalences
c           --- In this section of code, this is purely for convenience
c           --- and to make the following code look cleaner.
            nn = npinje_s(jsid+1)
            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is)+npinje_s(jsid+1)-1
            xx => pgroup%xp(i1:i2)
            yy => pgroup%yp(i1:i2)
            zz => pgroup%zp(i1:i2)
            ux => pgroup%uxp(i1:i2)
            uy => pgroup%uyp(i1:i2)
            uz => pgroup%uzp(i1:i2)
            gg => pgroup%gaminv(i1:i2)
            id => pgroup%pid(i1:i2,injpid)

c           --- get E self-field at initial positions
            if (linj_efromgrid) then
              call fetche3dfrompositions(jsid,pgroup%ndts(is-1),nn,xx,yy,zz,
     &                                   ex,ey,ez)
              call fetchb3dfrompositions(jsid,nn,xx,yy,zz,bx,by,bz)
            endif
            call inj_sete3d(nn,xx,yy,zz,id,ex,ey,ez)

c           --- get external fields at initial positions
            call exteb3d(nn,xx,yy,zz,uz,gg,0.,fulldt_s*0.5,bx,by,bz,ex,ey,ez,
     &                   pgroup%sm(is),pgroup%sq(is),bendres,bendradi,1.,fulldt_s)

c           --- Fetch the time step size
            tt(1:nn) = fulldt_s*(id - int(id))

c           --- do half velocity advance with E fields
            call epusht3d(nn,ux,uy,uz,ex,ey,ez,pgroup%sq(is),pgroup%sm(is),tt,0.5)

c           --- Advance relativistic Gamma factor
            call gammaadv(nn,gg,ux,uy,uz,gamadv,lrelativ)

c           --- do half velocity advance with B fields
            call bpusht3d(nn,ux,uy,uz,gg,bx,by,bz,pgroup%sq(is),pgroup%sm(is),
     &                    tt,0.5,ibpush)

c           --- do full position advance
            call xpusht3d(nn,xx,yy,zz,ux,uy,uz,gg,tt)

c           --- setrho is now done by the setrho call in padvnc3d.
          enddo

c         --- print warning if no particles were injected
          if (npinject == 0) then
            call remark("No particles were injected.")
          endif

        elseif (itask == 2) then
c Do second part of constant current injection: get new E fields (this is
c done after field solve including injected particles), synchronize velocity
c with position, and then move velocity one half timestep back to match time
c level of rest of particles.
c         --- loop over species
          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
            if (.not. pgroup%ldts(is-1)) cycle
            fulldt_s = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)

#ifndef MPIPARALLEL
c           --- Skip this if no particles of this species were injected.
            if(npinje_s(jsid+1)==0) cycle
            nn = npinje_s(jsid+1)
#else
c           --- Find particles that were freshly injected. These
c           --- are particles which still have noninteger pid(:,injpid).
            allocate(mask(pgroup%nps(is)),stat=allocerror)
            if (allocerror /= 0) then
              print*,"inject3d: allocation error ",allocerror,
     &               ": could not allocate mask to shape ",pgroup%nps(is)
              stop
            endif

            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is)+pgroup%nps(is)-1
            where (pgroup%pid(i1:i2,injpid) < ninject+1)
              mask = ((pgroup%pid(i1:i2,injpid) -
     &             int(pgroup%pid(i1:i2,injpid))))
            elsewhere
              mask = 0.
            endwhere
            nn = COUNT(mask > 0.)
            if (nn == 0) then
              deallocate(mask)
              cycle
            endif
#endif

c           --- Make sure that there is enough room in the temporary arrays.
            if (nn > npgrp) then
              npgrp = nn
              call gchange("Setpwork3d",0)
            endif

c           --- Make equivalences so names make sense below
            ex => at
            ey => bt
            ez => apt
            bx => bpt
            by => xct
            bz => xpct
            bendres => yct
            bendradi => ypct

c           --- Zero the B field which is accumulated below
            ex = 0.
            ey = 0.
            ez = 0.
            bx = 0.
            by = 0.
            bz = 0.

c           --- Assume that injection will never be done inside a bend.
            bendres = 0.
            bendradi = 1.

#ifndef MPIPARALLEL
c           --- Fetch time step size
            i1 = pgroup%ins(is)+ntinject(jsid+1)
            i2 = pgroup%ins(is)+ntinject(jsid+1)+npinje_s(jsid+1)-1
            tt(1:nn) = fulldt_s*(pgroup%pid(i1:i2,injpid) -
     &                       int(pgroup%pid(i1:i2,injpid)))
c           --- Create equivalences
            xx => pgroup%xp(i1:i2)
            yy => pgroup%yp(i1:i2)
            zz => pgroup%zp(i1:i2)
            ux => pgroup%uxp(i1:i2)
            uy => pgroup%uyp(i1:i2)
            uz => pgroup%uzp(i1:i2)
            gg => pgroup%gaminv(i1:i2)
            id => pgroup%pid(i1:i2,injpid)
#else
c           --- Create equivalences
            xx => xt(1:nn)
            yy => yt(1:nn)
            zz => zt(1:nn)
            ux => uxt(1:nn)
            uy => uyt(1:nn)
            uz => uzt(1:nn)
            gg => rt(1:nn)
            id => perpscal(1:nn)
c           --- Fetch time step size
            tt(1:nn) = fulldt_s*PACK(mask,mask > 0.)
c           --- Copy particle data into temporary arrays, collecting
c           --- only particles that have been injected this time step.
            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is)+pgroup%nps(is)-1
            xx = PACK(pgroup%xp(i1:i2),mask > 0.)
            yy = PACK(pgroup%yp(i1:i2),mask > 0.)
            zz = PACK(pgroup%zp(i1:i2),mask > 0.)
            ux = PACK(pgroup%uxp(i1:i2),mask > 0.)
            uy = PACK(pgroup%uyp(i1:i2),mask > 0.)
            uz = PACK(pgroup%uzp(i1:i2),mask > 0.)
            gg = PACK(pgroup%gaminv(i1:i2),mask > 0.)
            id = PACK(pgroup%pid(i1:i2,injpid),mask > 0.)
#endif

c           --- calculate new E self-fields
            if (linj_efromgrid) then
              call fetche3dfrompositions(jsid,pgroup%ndts(is-1),nn,xx,yy,zz,
     &                                   ex,ey,ez)
              call fetchb3dfrompositions(jsid,nn,xx,yy,zz,bx,by,bz)
            endif
            call inj_sete3d(nn,xx,yy,zz,id,ex,ey,ez)

c           --- Get external fields at current positions.
            call exteb3d(nn,xx,yy,zz,uz,gg,-fulldt_s*0.5,0.,bx,by,bz,ex,ey,ez,
     &                pgroup%sm(is),pgroup%sq(is),bendres,bendradi,1.,fulldt_s)

c           --- complete B advance
            call bpusht3d(nn,ux,uy,uz,gg,bx,by,bz,pgroup%sq(is),pgroup%sm(is),
     &                    tt,0.5,ibpush)

c           --- complete the E advance
            call epusht3d(nn,ux,uy,uz,ex,ey,ez,pgroup%sq(is),pgroup%sm(is),
     &                    tt,0.5)

c           --- Advance relativistic Gamma factor
            call gammaadv(nn,gg,ux,uy,uz,gamadv,lrelativ)

c           --- Now, move velocites back one half a step
c           --- first half of a backward B advance
            call bpush3d(nn,ux,uy,uz,gg,bx,by,bz,pgroup%sq(is),pgroup%sm(is),
     &                    -0.5*fulldt_s,ibpush)

c           --- then half of a backward E advance
            call epush3d(nn,ux,uy,uz,ex,ey,ez,pgroup%sq(is),pgroup%sm(is),
     &                    -0.5*fulldt_s)

c           --- Advance relativistic Gamma factor
            call gammaadv(nn,gg,ux,uy,uz,gamadv,lrelativ)

#ifdef MPIPARALLEL
c           --- Copy data back into the particle arrays
            i1 = pgroup%ins(is)+ntinject(jsid+1)
            i2 = pgroup%ins(is)+ntinject(jsid+1)+pgroup%nps(is)-1
            pgroup%xp(i1:i2) = UNPACK(xx,mask>0.,pgroup%xp(i1:i2))
            pgroup%yp(i1:i2) = UNPACK(yy,mask>0.,pgroup%yp(i1:i2))
            pgroup%zp(i1:i2) = UNPACK(zz,mask>0.,pgroup%zp(i1:i2))
            pgroup%uxp(i1:i2) = UNPACK(ux,mask>0.,pgroup%uxp(i1:i2))
            pgroup%uyp(i1:i2) = UNPACK(uy,mask>0.,pgroup%uyp(i1:i2))
            pgroup%uzp(i1:i2) = UNPACK(uz,mask>0.,pgroup%uzp(i1:i2))
            pgroup%pid(i1:i2,injpid) = int(pgroup%pid(i1:i2,injpid))
            deallocate(mask)
#endif

c        --- end of loop over species
         enddo

        endif

c Axially directed space-charge limited injection is now complete
c for this time step.
      endif

c------------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c Injection off of a tranverse facing surface.
c Uses any of the three methods implemented, based off the value of
c the variable 'inject'.

      if (ntinj > 0 .and.
     &    1 <= inject .and. inject <= 3) then

        if (itask == 1) then

c         --- zero ntinject for inject = 2 or 3
          if (inject == 2 .or. inject == 3) then
            do is=1,pgroup%ns
              jsid = pgroup%sid(is-1)
              if (.not. pgroup%ldts(is-1)) cycle
              ntinject(jsid+1) = 0
            enddo
          endif

c         --- Loop over injection sources
          do ij=1,ntinj

c           --- loop over species
            do is=1,pgroup%ns
              jsid = pgroup%sid(is-1)
              if (.not. pgroup%ldts(is-1)) cycle
              fulldt_s    = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)
              nti = 0
              qoverm = pgroup%sq(is)/pgroup%sm(is)

c             --- Load particles one azimuthal section at a time,
c             --- only if there is a positive number of injected particles.
              do iz=0,nztinj(ij)-1
                zm = ztinjmn(ij) + iz*dztinj(ij)
                ainj = atinject(iz,ij)
                binj = btinject(iz,ij)
                ainjp1 = atinject(iz+1,ij)
                binjp1 = btinject(iz+1,ij)
                ainji = 1./atinject(iz,ij)
                binji = 1./btinject(iz,ij)
                if (ltinj_outward(ij)) then
                  rsign = +1.
                else
                  rsign = -1.
                endif

c               --- Area of section of emitting surface.
                circum = ellipseperimeter(ainj,binj)
                area = dztinj(ij)*circum/nttinj(ij)

                do ith=0,nttinj(ij)-1

c                 --- angle of point in transverse plane
                  aa = ith*2.*pi/nttinj(ij)
                  p1x = ainj*cos(aa)
                  p1y = binj*sin(aa)
c                 rinj = sqrt(p1x**2 + p1y**2)
c                 --- For an ellipse, the larger the axis, the smaller the curvature,
c                 --- so as an adhoc fix for the cylindrical correction to CL,
c                 --- swap x and y to estimate the radius of curvature. This has
c                 --- no effect on a circle.
                  rinj = sqrt((ainj*sin(aa))**2 + (binj*cos(aa))**2)
                  rinji = 1./rinj

c                 --- Find coordinates of the point on the virtual surface
c                 --- along a line perpendicular to the emitting surface.
c                 --- Also calculate the distance of that point from the surface.
c                 --- Note that dx, dy, and dz are used since these relate to
c                 --- the field grid, not the injection grid.
                  p2x = (ainj + rsign*inj_dx)*cos(aa)
                  p2y = (binj + rsign*inj_dx*ainj/binj)*sin(aa)
                  dr = rsign*sqrt((p2x-p1x)**2 + (p2y-p1y)**2)
                  dri = 1./dr

c                 --- Fetch difference between phi at that point and phi on
c                 --- the emitting surface, averaging over z
                  vavez = 0.5*(tinj_phi(ith,iz,ij)+tinj_phi(ith,iz+1,ij))

c                 --- normal velocity
                  vnorm = fulldt_s*qoverm*vavez*dri

c                 --- Axial velocity
c                 --- Calculated using axial E field on virtual surface.
c                 --- This does not take into account any variation in a and btinject
                  vznorm = fulldt_s*qoverm*
     &                     (tinj_phi(ith,iz,ij)-tinj_phi(ith,iz+1,ij))/dztinj(ij)

c                 --- number of particles injected in grid cell
                  if (inject == 1) then
                    rnn = ntinject(jsid+1)*ftinject(jsid+1,ij)/real(nttinj(ij)*nztinj(ij))
                  elseif (inject == 2) then
c                   --- Const for inject = 2.  Denominator includes factor for
c                   --- emitting on concentric cylinders.  See I. Langmuir,
c                   --- K. Blodgett, "Currents Limited by Space Charge Between
c                   --- Concentric Spheres", PhysRev, 1924.
                    const = 4./9.*eps0*(2.*abs(qoverm))**.5*dri**2*
     &                      fulldt_s*area/abs(pgroup%sq(is))/pgroup%sw(is)/
     &                      (1. - 0.8*dr*rinji + 0.66*(dr*rinji)**2)
                    tinjprev(ith,iz,ij) = ftinject(jsid+1,ij)*(
     &                                    inj_param*const*abs(vavez)**1.5 +
     &                                    (1. - inj_param)*tinjprev(ith,iz,ij))
                    rnn = tinjprev(ith,iz,ij)
                    if (ith == 0) tinjprev(nttinj(ij),iz,ij) = tinjprev(ith,iz,ij)
                  elseif (inject == 3) then
c                   --- Find coordinates of the point.
c                   --- This only works with basic 3d and RZ solvers
                    p1x = abs(xtinject(ij) + p1x - xmmin)/dx
                    p1y = abs(ytinject(ij) + p1y - ymmin)/dy
                    p1z = (zm - zmmin)/dz
                    i1x = p1x
                    i1y = p1y
                    i1z = p1z
                    w1x = p1x - i1x
                    w1y = p1y - i1y
                    w1z = p1z - i1z
c                   --- get charge density at injection point
                    rr = rhop(i1x  ,i1y  ,i1z  )*(1.-w1x)*(1.-w1y)*(1.-w1z) +
     &                   rhop(i1x  ,i1y+1,i1z  )*(1.-w1x)*    w1y *(1.-w1z) +
     &                   rhop(i1x+1,i1y  ,i1z  )*    w1x *(1.-w1y)*(1.-w1z) +
     &                   rhop(i1x+1,i1y+1,i1z  )*    w1x *    w1y *(1.-w1z) +
     &                   rhop(i1x  ,i1y  ,i1z+1)*(1.-w1x)*(1.-w1y)*    w1z  +
     &                   rhop(i1x  ,i1y+1,i1z+1)*(1.-w1x)*    w1y *    w1z  +
     &                   rhop(i1x+1,i1y  ,i1z+1)*    w1x *(1.-w1y)*    w1z  +
     &                   rhop(i1x+1,i1y+1,i1z+1)*    w1x *    w1y *    w1z
                    rnn=(eps0*vavez*dri - rr*abs(dr))*area/
     &                   pgroup%sq(is)/pgroup%sw(is)*ftinject(jsid+1,ij)
                  endif
                  rnn = min(rnn,
     &                      abs(jmaxtinj(ij)*area*fulldt_s/pgroup%sq(is)/
     &                      pgroup%sw(is)*ftinject(jsid+1,ij)))
                  nn = int(rnn + wranf())

c                 --- Make sure there is room for more particles.  Estimate
c                 --- total number of particles to be injected from the
c                 --- source.
                  call chckpart(pgroup,is,ntinject(jsid+1) + nn*nttinj(ij)*nztinj(ij),
     &                          0,.true.)
                  call chckpart(pgroup,is,nn,0,.true.)

                  do ip=1,nn
                    ii = pgroup%ins(is) - nti - 1
c                   --- calculate position of new particle
                    wz = wrandom(xrandom,injctcnt,dig1,1,1)
                    pgroup%zp(ii) = zm + wz*dztinj(ij)
                    t = wrandom(xrandom,injctcnt,dig3,1,1)
                    t = aa + (t - .5)*2.*pi/nttinj(ij)
                    pgroup%xp(ii) = xtinject(ij) + (ainj*(1.-wz) + ainjp1*wz)*cos(t)
                    pgroup%yp(ii) = ytinject(ij) + (binj*(1.-wz) + binjp1*wz)*sin(t)
                    vtx = 0.5*vztinject(ij)*emitx_s(jsid+1)/(circum/(2.*pi)) + 
     &                    vthperp_s(jsid+1)
                    vty = 0.5*vztinject(ij)*emity_s(jsid+1)/(circum/(2.*pi)) + 
     &                    vthperp_s(jsid+1)
                    pgroup%uxp(ii) = vnorm*cos(t) +
     &                    vtx*wrandomgauss(vtrandom,injctcnt,dig3,dig4,1,1,.false.)
                    pgroup%uyp(ii) = vnorm*sin(t) +
     &                    vty*wrandomgauss(vtrandom,injctcnt,dig5,dig6,1,1,.false.)
                    pgroup%uzp(ii) = vznorm + vthz_s(jsid+1)*
     &                    wrandomgauss(vzrandom,injctcnt,dig7,dig8,1,1,.false.)
                    if (pgroup%uzp(ii) == 0.) pgroup%uzp(ii) = SMALLPOS
                    if (lrelativ) then
                      pgroup%gaminv(ii) = sqrt(1.-
     &     (pgroup%uxp(ii)**2+pgroup%uyp(ii)**2+pgroup%uzp(ii)**2)*clightsqi)
                      gamma = 1./pgroup%gaminv(ii)
                      pgroup%uxp(ii) = pgroup%uxp(ii)*gamma
                      pgroup%uyp(ii) = pgroup%uyp(ii)*gamma
                      pgroup%uzp(ii) = pgroup%uzp(ii)*gamma
                    else
                      pgroup%gaminv(ii) = 1.
                    endif
c                   --- increment particle counter
                    nti = nti + 1
c                   --- Clear out any old data from pid. This is done here
c                   --- before pid is used anywhere.
                    pgroup%pid(ii,:) = 0.
                    pgroup%pid(ii,injpid) = ninject + ij + rnrev(injctcnt,dig5)
c                   --- Save SSN 
                    if (spid>0) then
                      pgroup%pid(ii,spid) = ssn
                      ssn = ssn+1
                    end if
c                   --- increment random number counter
                    injctcnt = injctcnt + 1
                  enddo
                enddo
              enddo

c             --- Set number of particles that were injected for inject = 2 or 3
              if (inject == 2 .or. inject == 3) then
                ntinject(jsid+1) = ntinject(jsid+1) + nti
              endif

c           --- end loop over species
            enddo

c         --- end loop over transverse injection sources
          enddo

c         --- Call transverse particle scraping routine to force
c         --- removal of particles outside the grid.
          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
            if (.not. pgroup%ldts(is-1)) cycle
            ii = pgroup%ins(is) - ntinject(jsid+1)
            if (ntinject(jsid+1) > 0) then
              call stckxy3d(ntinject(jsid+1),pgroup%xp(ii),xmmax,xmmin,dx,
     &                      pgroup%yp(ii),ymmax,
     &                      ymmin,dy,pgroup%zp(ii),zmmin,dz,pgroup%uxp(ii),
     &                      pgroup%uyp(ii),pgroup%uzp(ii),pgroup%gaminv(ii),
     &                      zgrid,zbeam,l2symtry,l4symtry,absorb,.false.)
              call processlostpart(pgroup,is,clearlostpart,time,zbeam)
            endif
          enddo

          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
            if (.not. pgroup%ldts(is-1)) cycle
c           --- Sum ntinject to get the total number of particles injected.
            npinject = npinject + ntinject(jsid+1)

c           --- make injected particles live particles
            pgroup%ins(is) = pgroup%ins(is) - ntinject(jsid+1)
            pgroup%nps(is) = pgroup%nps(is) + ntinject(jsid+1)
          enddo

c         --- print warning if no particles were injected
          if (npinject == 0) then
            call remark("No particles were injected transversely.")
          endif

c         --- Now, advance the particles off of the emitting surface using
c         --- isochronous leapfrog.  Each particle has it's own time step
c         --- size, uniformly distributed between 0 and dt.
          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
c           --- Skip this if no particles of this species were injected.
            if (.not. pgroup%ldts(is-1)) cycle
            fulldt_s    = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)
            if(ntinject(jsid+1)==0) cycle

c           --- Make sure that there is enough room in the temporary arrays.
            if (ntinject(jsid+1) > npgrp) then
              npgrp = ntinject(jsid+1)
              call gchange("Setpwork3d",0)
            endif

c           --- Make equivalences so names make sense below
            ex => at
            ey => bt
            ez => apt
            bx => bpt
            by => xct
            bz => xpct
            bendres => yct
            bendradi => ypct

c           --- Zero the B field which is accumulated below
            ex = 0.
            ey = 0.
            ez = 0.
            bx = 0.
            by = 0.
            bz = 0.

c           --- Assume that injection will never be done inside a bend.
            bendres = 0.
            bendradi = 1.

c           --- Create equivalences
c           --- In this section of code, this is purely for convenience
c           --- and to make the following code look cleaner.
            nn = ntinject(jsid+1)
            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is)+ntinject(jsid+1)-1
            xx => pgroup%xp(i1:i2)
            yy => pgroup%yp(i1:i2)
            zz => pgroup%zp(i1:i2)
            ux => pgroup%uxp(i1:i2)
            uy => pgroup%uyp(i1:i2)
            uz => pgroup%uzp(i1:i2)
            gg => pgroup%gaminv(i1:i2)
            id => pgroup%pid(i1:i2,injpid)

c           --- get E self-field at initial positions
            if (linj_efromgrid) then
              call fetche3dfrompositions(jsid,pgroup%ndts(is-1),nn,xx,yy,zz,
     &                                   ex,ey,ez)
              call fetchb3dfrompositions(jsid,nn,xx,yy,zz,bx,by,bz)
            endif
            call tinj_sete3d(nn,xx,yy,zz,id,ex,ey,ez)

c           --- get external fields at initial positions
            call exteb3d(nn,xx,yy,zz,uz,gg,0.,fulldt_s*0.5,bx,by,bz,ex,ey,ez,
     &                   pgroup%sm(is),pgroup%sq(is),bendres,bendradi,1.,fulldt_s)

c           --- Fetch the time step size
            tt(1:nn) = fulldt_s*(id - int(id))

c           --- do half velocity advance with E fields
            call epusht3d(nn,ux,uy,uz,ex,ey,ez,pgroup%sq(is),pgroup%sm(is),tt,0.5)

c           --- Advance relativistic Gamma factor
            call gammaadv(nn,gg,ux,uy,uz,gamadv,lrelativ)

c           --- do half velocity advance with B fields
            call bpusht3d(nn,ux,uy,uz,gg,bx,by,bz,pgroup%sq(is),pgroup%sm(is),
     &                    tt,0.5,ibpush)

c           --- do full position advance
            call xpusht3d(nn,xx,yy,zz,ux,uy,uz,gg,tt)

c           --- setrho is now done by the setrho call in padvnc3d.
          enddo

        elseif (itask == 2) then
c Do second part of constant current injection: get new E fields (this is
c done after field solve including injected particles), synchronize velocity
c with position, and then move velocity one half timestep back to match time
c level of rest of particles.
c         --- loop over species
          do is=1,pgroup%ns
            jsid = pgroup%sid(is-1)
            if (.not. pgroup%ldts(is-1)) cycle
            fulldt_s = dt*pgroup%ndts(is-1)*pgroup%dtscale(is)

#ifndef MPIPARALLEL
c           --- Skip this if no particles of this species were injected.
            if(ntinject(jsid+1)==0) cycle
            nn = ntinject(jsid+1)
#else
c           --- Find particles that were freshly injected. These
c           --- are particles which still have noninteger pid(:,injpid).
            allocate(mask(pgroup%nps(is)),stat=allocerror)
            if (allocerror /= 0) then
              print*,"inject3d: allocation error ",allocerror,
     &               ": could not allocate mask to shape ",pgroup%nps(is)
              stop
            endif

            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is)+pgroup%nps(is)-1
            where (pgroup%pid(i1:i2,injpid) > ninject+1)
              mask = ((pgroup%pid(i1:i2,injpid) -
     &             int(pgroup%pid(i1:i2,injpid))))
            elsewhere
              mask = 0.
            endwhere
            nn = COUNT(mask > 0.)
            if (nn == 0) then
              deallocate(mask)
              cycle
            endif
#endif

c           --- Make sure that there is enough room in the temporary arrays.
            if (nn > npgrp) then
              npgrp = nn
              call gchange("Setpwork3d",0)
            endif

c           --- Make equivalences so names make sense below
            ex => at
            ey => bt
            ez => apt
            bx => bpt
            by => xct
            bz => xpct
            bendres => yct
            bendradi => ypct

c           --- Zero the B field which is accumulated below
            ex = 0.
            ey = 0.
            ez = 0.
            bx = 0.
            by = 0.
            bz = 0.

c           --- Assume that injection will never be done inside a bend.
            bendres = 0.
            bendradi = 1.

#ifndef MPIPARALLEL
c           --- Fetch time step size
            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is)+ntinject(jsid+1)-1
            tt(1:nn) = fulldt_s*(pgroup%pid(i1:i2,injpid) -
     &                       int(pgroup%pid(i1:i2,injpid)))
c           --- Create equivalences
            xx => pgroup%xp(i1:i2)
            yy => pgroup%yp(i1:i2)
            zz => pgroup%zp(i1:i2)
            ux => pgroup%uxp(i1:i2)
            uy => pgroup%uyp(i1:i2)
            uz => pgroup%uzp(i1:i2)
            gg => pgroup%gaminv(i1:i2)
            id => pgroup%pid(i1:i2,injpid)
#else
c           --- Create equivalences
            xx => xt(1:nn)
            yy => yt(1:nn)
            zz => zt(1:nn)
            ux => uxt(1:nn)
            uy => uyt(1:nn)
            uz => uzt(1:nn)
            gg => rt(1:nn)
            id => perpscal(1:nn)
c           --- Fetch time step size
            tt(1:nn) = fulldt_s*PACK(mask,mask > 0.)
c           --- Copy particle data into temporary arrays, collecting
c           --- only particles that have been injected this time step.
            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is)+pgroup%nps(is)-1
            xx = PACK(pgroup%xp(i1:i2),mask > 0.)
            yy = PACK(pgroup%yp(i1:i2),mask > 0.)
            zz = PACK(pgroup%zp(i1:i2),mask > 0.)
            ux = PACK(pgroup%uxp(i1:i2),mask > 0.)
            uy = PACK(pgroup%uyp(i1:i2),mask > 0.)
            uz = PACK(pgroup%uzp(i1:i2),mask > 0.)
            gg = PACK(pgroup%gaminv(i1:i2),mask > 0.)
            id = PACK(pgroup%pid(i1:i2,injpid),mask > 0.)
#endif

c           --- calculate new E self-fields
            if (linj_efromgrid) then
              call fetche3dfrompositions(jsid,pgroup%ndts(is-1),nn,xx,yy,zz,
     &                                   ex,ey,ez)
              call fetchb3dfrompositions(jsid,nn,xx,yy,zz,bx,by,bz)
            endif
            call tinj_sete3d(nn,xx,yy,zz,id,ex,ey,ez)

c           --- Get external fields at current positions.
            call exteb3d(nn,xx,yy,zz,uz,gg,-fulldt_s*0.5,0.,bx,by,bz,ex,ey,ez,
     &                pgroup%sm(is),pgroup%sq(is),bendres,bendradi,1.,fulldt_s)

c           --- complete B advance
            call bpusht3d(nn,ux,uy,uz,gg,bx,by,bz,pgroup%sq(is),pgroup%sm(is),
     &                    tt,0.5,ibpush)

c           --- complete the E advance
            call epusht3d(nn,ux,uy,uz,ex,ey,ez,pgroup%sq(is),pgroup%sm(is),
     &                    tt,0.5)

c           --- Advance relativistic Gamma factor
            call gammaadv(nn,gg,ux,uy,uz,gamadv,lrelativ)

c           --- Now, move velocites back one half a step
c           --- first half of a backward B advance
            call bpush3d(nn,ux,uy,uz,gg,bx,by,bz,pgroup%sq(is),pgroup%sm(is),
     &                    -0.5*fulldt_s,ibpush)

c           --- then half of a backward E advance
            call epush3d(nn,ux,uy,uz,ex,ey,ez,pgroup%sq(is),pgroup%sm(is),
     &                    -0.5*fulldt_s)

c           --- Advance relativistic Gamma factor
            call gammaadv(nn,gg,ux,uy,uz,gamadv,lrelativ)

#ifdef MPIPARALLEL
c           --- Copy data back into the particle arrays
            i1 = pgroup%ins(is)
            i2 = pgroup%ins(is)+pgroup%nps(is)-1
            pgroup%xp(i1:i2) = UNPACK(xx,mask>0.,pgroup%xp(i1:i2))
            pgroup%yp(i1:i2) = UNPACK(yy,mask>0.,pgroup%yp(i1:i2))
            pgroup%zp(i1:i2) = UNPACK(zz,mask>0.,pgroup%zp(i1:i2))
            pgroup%uxp(i1:i2) = UNPACK(ux,mask>0.,pgroup%uxp(i1:i2))
            pgroup%uyp(i1:i2) = UNPACK(uy,mask>0.,pgroup%uyp(i1:i2))
            pgroup%uzp(i1:i2) = UNPACK(uz,mask>0.,pgroup%uzp(i1:i2))
            pgroup%pid(i1:i2,injpid) = int(pgroup%pid(i1:i2,injpid))
            deallocate(mask)
#endif

c        --- end of loop over species
         enddo

        endif

      endif
c     --- End of transverse injection

!$OMP MASTER
      if (lw3dtimesubs) timeinject3d = timeinject3d + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine injctint(pgroup)
      use ParticleGroupmodule
      use Subtimers3d
      use Constant
      use Beam_acc
      use InGen
      use InGen3d
      use InPart
      use InMesh3d
      use Particles,Only: npmax,wpid,tpid,rpid
      use InjectVars
      use InjectVars3d
      use Setpwork3d
      use Picglb3d
      use Parallel
      use Obsolete3d,obs_inj_d=>inj_d,obs_inj_f=>inj_f
      type(ParticleGroup):: pgroup

      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: ixmin,ixmax,iymin,iymax
      integer(ISZ):: inx,iny,iz
      integer(ISZ):: js,is,jsid,ip,ij,nn,nnmax,ij1,ij2
      real(kind=8):: total_frac,dd,maxdxdy
      real(kind=8):: maxab1,maxab2,minab1,minab2
      real(kind=8):: circum
      integer(ISZ), external :: oneiftrue, nextpid
      real(kind=8):: ellipseperimeter
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (inject < 1 .or. inject > 3) return

c Set wpid if needed and not set by user
      if((l_inj_rz.or.l_inj_regular) .and. wpid==0) then
        if (wpid == 0) wpid = nextpid()
      endif

c If delay on setting temperature of injected particles, then
c time of particle creation is recorded initially with a negative value
c switched to positive when temperature has been added.
      if(l_inj_delay_temp) l_inj_rec_inittime = .true.

c Set tpid if needed and not set by user
      if(l_inj_rec_inittime .and. tpid==0) then
        if (tpid == 0) tpid = nextpid()
      end if

c Set rpid if needed and not set by user
      if(l_inj_rec_initradius .and. rpid==0) then
        if (rpid == 0) rpid = nextpid()
      end if

c Set pid index for injection information for particles
      if (injpid == 0) injpid = nextpid()

c Make sure that space has been allocated for the arrays. This is needed
c in case the user changes ninject or ns but doesn't explicitly allocate
c the additional needed space.
      call gchange("InjectVars",0)

c Initialize beam size from initial envelope parameters if beam size not set
c One problem with this is that if the user ever wants apinject or bpinject
c to be zero while ap0 or bp0 not be zero.  The user would have to reset
c apinject or bpinject back to zero after the generate.
      do ij=1,ninject
        if (ainject(ij)  == 0.) ainject(ij)  = a0
        if (binject(ij)  == 0.) binject(ij)  = b0
        if (apinject(ij) == 0.) apinject(ij) = ap0
        if (bpinject(ij) == 0.) bpinject(ij) = bp0
      enddo

c     --- Setup grids covering injection surfaces.
      call getparticleextant(1,xmin,xmax,ymin,ymax,xmmin,xmmax,ymmin,ymmax)
      inj_ninj = ninject
      do ij=1,ninject
        ixmin = max(xmin,xinject(ij) - ainject(ij) - 2*dx)*oneiftrue(.not. l_inj_rz)
        ixmax = min(xmax,xinject(ij) + ainject(ij) + 2*dx)
        iymin = max(ymin,yinject(ij) - binject(ij) - 2*dy)*oneiftrue(.not. l_inj_rz)
        iymax = min(ymax,yinject(ij) + binject(ij) + 2*dy)*oneiftrue(.not. l_inj_rz)
        if(inj_dx==0.) then
          inx = nint((ixmax - ixmin)/dx)
        else
          inx = nint((ixmax - ixmin)/inj_dx)
        endif
        if(inj_dy==0.) then
          iny = nint((iymax - iymin)/dy)
        else
          iny = nint((iymax - iymin)/inj_dy)
        endif
        inj_nx = max(inj_nx,inx)
        inj_ny = max(inj_ny,iny)
      enddo
      call gchange("InjectVars3d",0)

      if(inj_dx==0.) inj_dx = dx
      if(inj_dy==0.) inj_dy = dy
      if(inj_dz==0.) inj_dz = dz
      if(l_inj_rz) then
        if(.not. l_inj_area) then
          inj_dx = ainject(1)/(real(nint(ainject(1)/inj_dx))+0.5)
          inj_xwide = 1
          inj_ywide = 1
        end if
      else
        l_inj_area = .true.
      end if
      if(solvergeom==RZgeom) inj_dy = inj_dx
      if(l_inj_rz) inj_ny = 0
      inj_dz0 = inj_dz
      if(inj_nz>1) inj_dz = inj_d(1)*inj_dz/(real(inj_nz-1)*2.**(1./3.)+real(2-inj_nz))**3
      if(inj_nz > 1 .and. ninject > 1) then
        write(0,*) 'inj_nz > 1 not yet supported for ninject > 1'
        write(0,*) 'stop'
        stop
      end if

c     --- Set mins for meshes around each source.
      inj_xmmin = max(xmin, - inj_nx/2*inj_dx)
      inj_ymmin = max(ymin, - inj_ny/2*inj_dy)
      if (l_inj_rz) inj_xmmin = 0.

c     --- Set inject z min and max.
      call updateinj_zmminmax()

c Check if the user has set the obsolete values of inj_d or inj_f.
c If so, then copy values to the new variables in top.
      if (obs_inj_d /= 0.) then
        inj_d = obs_inj_d
        call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        call remark("Notice: The variable inj_d has been moved to the top")
        call remark("        package. Please fix your input files accordingly.")
        call remark("        Your simulation will run ok now, but in future")
        call remark("        versions, it will stop with an error.")
        call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      endif
      if (obs_inj_f /= 0.) then
        inj_f = obs_inj_f
        call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        call remark("Notice: The variable inj_f has been moved to the top")
        call remark("        package. Please fix your input files accordingly.")
        call remark("        Your simulation will run ok now, but in future")
        call remark("        versions, it will stop with an error.")
        call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      endif

c Set value of vzinject for inject==1 if not already set. For other
c types of injection, vzinject is assumed to be zero, since the starting
c velocity is calculated from the local fields.
      do ij=1,ninject
        do is=1,ns
          if (inject == 1 .and. vzinject(ij,is) == 0.) then
            vzinject(ij,is) = vbeam_s(is)
          endif
        enddo
      enddo

c Set the injection length and start, normally the length of one time step
c and at zimin, but can be set by user.
      if (leninjct == 0) leninjct = vbeam*dt

c     --- Calculate number of particles to inject for each species if
c     --- not set by the user.
c     --- npinje_s was being multiplied by sp_fract, but it is also being
c     --- multiplied by finject in inject3d, so the fraction was being
c     --- applied twice. The one here was removed since the finject is
c     --- still needed to control npinje_s per injection source.
c     --- The sp_fract was then also removed from the calculation of sw below.
      do is=1,ns
        if (npinje_s(is) == 0) npinje_s(is) = npinject
        if (rnpinje_s(is) == 0.) rnpinje_s(is) = npinje_s(is)
      enddo

c     --- Print warning if npinject=0
      if (npinject == 0) then
        call remark("ERROR: injctint, npinject is zero, no particles will
     &  be injected.")
      endif

c Each species is stored in seperate parts of the particle arrays, the size of
c which is based off the beam fraction of each species.  Check if there is
c enough space in each of the species blocks.
c If there are no particles yet, then reallocate the arrays to the correct
c size. This is much faster when there are multiple species.
c Also set particle parameters: charge, mass and weight
      total_frac = sum(sp_fract)
      if (maxval(pgroup%nps) == 0) then
        npmax = max(injctspc,npmax)
        np_s = 0
        sp_fract = sp_fract/total_frac
        call alotpart(pgroup)
        sp_fract = sp_fract*total_frac
      endif
      do js=0,pgroup%ns-1
        jsid = pgroup%sid(js)
        if (maxval(pgroup%nps) /= 0) then
          call chckpart(pgroup,js+1,int(injctspc*sp_fract(jsid+1)/total_frac+1),
     &                  0,.true.)
        endif
        if (pgroup%sq(js+1) == 0.) pgroup%sq(js+1) = zion_s(jsid+1)*echarge
        if (pgroup%sm(js+1) == 0.) pgroup%sm(js+1) = aion_s(jsid+1)*amu
        if (pgroup%sw(js+1) == 0.) then
          if (npinje_s(jsid+1) > 0) then
            pgroup%sw(js+1) = (ibeam_s(jsid+1)*dt*pgroup%ndts(js)/
     &                      pgroup%sq(js+1))/npinje_s(jsid+1)
          else if (js > 0 .and.
     &             pgroup%sq(js+1) == pgroup%sq(js) .and.
     &             pgroup%sm(js+1) == pgroup%sm(js) .and.
     &             pgroup%ndts(js) == 2*pgroup%ndts(js-1)) then
c           --- Automatically set the sw for species with larger ndts since in
c           --- general there will be no particles in those groups initially
c           --- (so npinje_s will be zero and the above skipped).
            pgroup%sw(js+1) = pgroup%sw(js)
          endif
        endif
      enddo

c Set up finject array if not set by the user (if sum of finject = 0).
c Give each injection source the same distribution of species.
      total_frac = 0.
      do is=1,ns
        do ij=1,ninject
          total_frac = total_frac + finject(ij,is)
        enddo
      enddo
      if (total_frac == 0.) then
        do is=1,ns
          do ij=1,ninject
            finject(ij,is) = sp_fract(is)
          enddo
        enddo
      endif

c Set the variable INJCTCNT
c This is used for loading the injection particles.  It ensures that the
c injected particles have new random numbers.
      injctcnt = sum(pgroup%nps) + 1 + randoffset

c Set up Setpwork3d arrays
c --- Estimate number of particles injected
      nnmax = 0
      do is=1,ns
        nn = 0
        do ij=1,ninject
          nn = nn + npinje_s(is)*finject(ij,is)
        enddo
        if (nn > nnmax) nnmax = nn
      enddo
      if (nnmax > npgrp) then
        npgrp = nnmax
        call gchange("Setpwork3d",0)
      endif

c Setup transverse injection

c Set the z cell size and number for transverse injection. The user can set nztinj.
      do ij=1,ntinj
        if (nztinj(ij) == 0) nztinj(ij) = nint((ztinjmx(ij) - ztinjmn(ij))/dz)
        if (nztinj(ij) == 0) nztinj(ij) = 1
        dztinj(ij) = (ztinjmx(ij) - ztinjmn(ij))/nztinj(ij)
        nztmax = max(nztmax, nztinj(ij))
      enddo
      call gchange("InjectVars",0)

c Set the default values of atinject and btinject and calculate the number of theta
c points if not already set.
      do ij=1,ntinj
        circum = 0.
        do iz=0,nztmax
          if (atinject(iz,ij) == -1.) atinject(iz,ij) = a0
          if (btinject(iz,ij) == -1.) btinject(iz,ij) = b0
          circum = max(circum,ellipseperimeter(atinject(:,ij),btinject(:,ij)))
        enddo
        if (nttinj(ij) == 0) nttinj(ij) = nint(circum/min(inj_dx,inj_dy))
        nttinjmax = max(nttinjmax,nttinj(ij))
      enddo
      call gchange("InjectVars",0)

c In case there is no axial injection, set the particle weight.
      do is=1,pgroup%ns
        jsid = pgroup%sid(is-1)
        if (pgroup%sw(is) == 0. .and. ntinject(jsid+1) > 0) then
          pgroup%sw(is) = (ibeam_s(jsid+1)*dt*pgroup%ndts(is-1)/
     &                     pgroup%sq(is))/ntinject(jsid+1)
        endif
      enddo

c Set up ftinject array if not set by the user (if sum of ftinject = 0).
c Give each injection source the same distribution of species.
      total_frac = 0.
      do is=1,ns
        do ij=1,ntinj
          total_frac = total_frac + finject(ij,is)
        enddo
      enddo
      if (total_frac == 0.) then
        do is=1,ns
          do ij=1,ntinj
            ftinject(ij,is) = sp_fract(is)
          enddo
        enddo
      endif

c Determine what the value of inj_ns should be. It can either be 1 or
c ns: if only one species will be injected from each source, the it is 1,
c otherwise, ns. This is needed in order that inj_prev and inj_np are set
c correctly. If multiple species are being injected from a single source,
c then inj_prev and inj_np must be saved for each species. In order to
c minimize the bookkeepping, if more than one species is injected from any
c source, than the two arrays are made big enough to save the data for all
c species for each source. If cases come up where this is too much wasted
c space, then the code can later be changed by only making inj_ns as big
c as is needed, i.e. the maximum number of species emitted from any one
c source.  
      inj_ns = 1
      do ij=1,ninject
        nn = 0
        do is=1,ns
          if (finject(ij,is) > 0.) nn = nn + 1
        enddo
        if (nn > 1) inj_ns = ns
      enddo

c Determine what the value if inj_ninj should be. It can either be 1 or
c ninject: if no sources are within 2 grid cells of another, then it is 1,
c otherwise ninject. This is needed since some parameters, inj_prev,
c inj_np, inj_area, and inj_rho, are stored on grids and must be known for
c each emitting source. If the sources are too close, then the data can
c overlap. In that case, the data is saved in seperate planes for each
c source. In order to minimize bookkeepping, if any two source (or more)
c sources are too close, then the arrays are made big enough so that all
c of the data is saved seperately for each source (whether or not they are
c too close to another). If cases comes up where this is too much wasted
c space, then the code can be changed so that inj_ninj is only made as big
c as needed, i.e. the maximum number of sources that are too close to each
c other.  
c The radius of each source is taken to be the maximum of a and b. This
c makes the checks much simpler but overly conservative.
c     inj_ninj = 1
c     maxdxdy = max(inj_dx,inj_dy)
c     do ij1=1,ninject
c       maxab1 = max(ainject(ij1),binject(ij1))
c       minab1 = min(ainjmin(ij1),binjmin(ij1))
c       do ij2=1,ninject
c         if (ij1 == ij2) cycle
c         dd = sqrt((xinject(ij1)-xinject(ij2))**2 +
c    &              (yinject(ij1)-yinject(ij2))**2)
c         maxab2 = max(ainject(ij2),binject(ij2))
c         minab2 = min(ainjmin(ij2),binjmin(ij2))
c         --- Check if the outer edges are near each other. If not
c         --- go to next pair.
c         if (maxab1 + maxab2 < dd - 2.*maxdxdy) cycle
c         --- Check if one is inside the other.
c         if (dd + maxab1 < minab2 - 2.*maxdxdy) cycle
c         if (dd + maxab2 < minab1 - 2.*maxdxdy) cycle
c         --- If it gets through the checks above, these two source are
c         --- too close
c         inj_ninj = ninject
c       enddo
c     enddo
c     --- Also check if sources are at different z locations. If so, then
c     --- save data for each source seperately. Checks if RMS of zinject
c     --- is zero.
c     if (sum(zinject**2) - sum(zinject)**2/ninject /= 0.) then
c       inj_ninj = ninject
c     endif


c Allocate inject arrays and fill with source location information.
c Note that gchange is called instead of gallot since some the of arrays
c in InjectVars3d may have been allocated by the user for specialized
c control of injection.
      call gchange("InjectVars3d",0)
      inj_xwide = max(1,inj_xwide)
      inj_ywide = max(1,inj_ywide)
      if(l_inj_rz) then
        call fill_inj_rz()
      else
        call fill_inj()
      end if

c Check if the voltages were specified.  If so, set so conductor will be
c set up for source.
c      do ij=1,ninject
c        if (vinject(ij) /= 0) lvinject = .true.
c      enddo

c Set vbeamfrm to zero so that the grid does not move off of the source,
c causing a code crash.
      vbeamfrm = 0.

c When using constant current injection, turn off the flag that specifies
c use of the normal E field with the Child-Langmuir profile
      if (inject == 1) linj_enormcl = .false.

!$OMP MASTER
      if (lw3dtimesubs) timeinjctint = timeinjctint + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine updateinj_zmminmax()
      use InjectVars3d,Only: inj_zmmin,inj_zmmax

c Set inject z min and max. This is only really needed for
c the parallel version, but is useful for the serial version to
c prevent particles from being injected off of the grid.
c These are also set in injctint but are reenforced here in case the
c regions have changed.

#ifdef MPIPARALLEL
      use Parallel, Only: my_index,zpslmin,zpslmax
      inj_zmmin = zpslmin(my_index)
      inj_zmmax = zpslmax(my_index)
#else
      use InMesh3d,Only: zmmin,zmmax
      inj_zmmin = zmmin
      inj_zmmax = zmmax
#endif

      return
      end
c=============================================================================
      subroutine fill_inj()
      use Subtimers3d
      use InPart
      use InMesh3d
      use InjectVars
      use InjectVars3d

c     --- Fill the injection grid array.  A padding of 2 grid cells is also
c     --- set.
c     --- The array inj_grid holds the axial location of the injection source
c     --- as a function of the transverse coordinates, in the lab frame.
c     --- The array inj_angl holds the angle due to the curvature of the source
c     --- as a function of the transverse coordinates.


      integer(ISZ):: ij,ix,iy
      real(kind=8):: xx,yy,rrsq,ainj,binj,ai_dx,bi_dy,dxi,dyi
      real(kind=8):: xl,yl,xu,yu,x1,x2,txl,tyl,txu,tyu
      real(kind=8):: ainji,binji
      real(kind=8):: area
      logical(ISZ):: linj_cell
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./inj_dx
      dyi = 1./inj_dy

c Loop over the injection sources.
      do ij=1,ninject

c       --- Provide some simple error checking
        if (ainject(ij) >= abs(rinject(ij)) .or.
     &      binject(ij) >= abs(rinject(ij))) then
          print*,"Error: the radius of curvature of the injection source must"
          print*,"       be greater than its transverse radius"
          stop
        endif

c       --- Calculate the extent of the source plus guard cells.
        ainj    = ainject(ij)
        binj    = binject(ij)
        ainji   = 1./ainj
        binji   = 1./binj
        ai_dx = ainject(ij) + 2.*inj_dx
        bi_dy = binject(ij) + 2.*inj_dy
c       --- Loop over transverse plane.
        do iy = 0,inj_ny
          do ix = 0,inj_nx
c           --- transverse coordinate relative to center of source
            xx = ix*inj_dx + inj_xmmin(ij)
            yy = iy*inj_dy + inj_ymmin(ij)
            
c           --- Check whether the cell is within the injection area.
c           --- This depends on whether the source is elliptical or
c           --- rectangular.
            if(linj_rectangle) then
              linj_cell = (abs(xx) <= ai_dx) .and. (abs(yy) <= bi_dy)
            else
              linj_cell = ((xx*bi_dy)**2 + (yy*ai_dx)**2) <= (ai_dx*bi_dy)**2
            endif

c           --- If the point is within source or guard cells, calculate both
c           --- axial distance and angle of curvature.
            if (linj_cell) then
c             --- Write the expression for inj_grid so as to avoid
c             --- subtraction of large, similar numbers when rinject is large.
c             inj_grid(ix,iy,ij) = (rinject(ij) -
c    &               sqrt(max(0.,rinject(ij)**2 - xx**2 - yy**2)))
              rrsq = xx**2 +  yy**2
              if (rrsq > rinject(ij)**2) rrsq = max(ainject(ij),binject(ij))**2
              inj_grid(ix,iy,ij) = rrsq/
     &             (abs(rinject(ij)) + sqrt(rinject(ij)**2 - rrsq))
              if (rinject(ij) < 0.) inj_grid(ix,iy,ij) = -inj_grid(ix,iy,ij)

              inj_angl(ix,iy,ij) = asin(sqrt(rrsq)/rinject(ij))

c             --- Calculate fraction of grid cell which is within the emitting
c             --- surface.  Used to scale the charge density on the emitting
c             --- surface.
c             --- The calculation is done with the emitting area scaled to
c             --- a unit circle.
c             --- The result is the fraction of the two by two block around
c             --- the grid point which is inside the emitting surface.

c             --- First, put the grid point in the first quadrant and scale it.
c             --- The max(0,) chop off any part that crosses the x or y axis.
              xl = max(0.,(abs(xx) - 0.5*inj_xwide*inj_dx)*ainji)
              yl = max(0.,(abs(yy) - 0.5*inj_ywide*inj_dy)*binji)
              xu =        (abs(xx) + 0.5*inj_xwide*inj_dx)*ainji
              yu =        (abs(yy) + 0.5*inj_ywide*inj_dy)*binji

c             --- If the two by two block is completely outside emitting
c             --- surface, set the area to zero.
              if (xl**2 + yl**2 > 1.) then
                area = 0.

c             --- If two by two block straddles edge of emitting surface, do
c             --- the calculation.
              elseif (xl**2 + yl**2 < 1. .and.
     &                                1. < xu**2 + yu**2) then
c               --- Calculate area of block within first quadrant.
                x1 = max(xl, sqrt(max(0.,1. - yu**2)))
                x2 = min(xu, sqrt(max(0.,1. - yl**2)))
                area = (x1 - xl)*(yu - yl) +
     &                 0.5*x2*sqrt(1. - x2**2) + 0.5*asin(x2) -
     &                 0.5*x1*sqrt(1. - x1**2) - 0.5*asin(x1) -
     &                 (x2 - x1)*yl
c               --- If block extended past x axis, calculate the area of
c               --- that piece.
                if (abs(xx) - inj_dx < 0.) then
                  txu = - (abs(xx) - inj_dx)*ainji
                  txl = 0.
                  x1 = max(txl, sqrt(max(0.,1. - yu**2)))
                  x2 = min(txu, sqrt(max(0.,1. - yl**2)))
                  area = area + (x1 - txl)*(yu - yl) +
     &                          0.5*x2*sqrt(1. - x2**2) + 0.5*asin(x2) -
     &                          0.5*x1*sqrt(1. - x1**2) - 0.5*asin(x1) -
     &                          (x2 - x1)*yl
                endif
c               --- If block extended past y axis, calculate the area of
c               --- that piece.
                if (abs(yy) - inj_dy < 0.) then
                  tyu = - (abs(yy) - inj_dy)*binji
                  tyl = 0.
                  x1 = max(xl, sqrt(max(0.,1. - tyu**2)))
                  x2 = min(xu, sqrt(max(0.,1. - tyl**2)))
                  area = area + (x1 - xl)*(tyu - tyl) +
     &                          0.5*x2*sqrt(1. - x2**2) + 0.5*asin(x2) -
     &                          0.5*x1*sqrt(1. - x1**2) - 0.5*asin(x1)
                endif

c               --- Normalize the area so that the total area of the block
c               --- would be 1 (since area was calculated above on a
c               --- unit circle).
                area = area*dxi*dyi*ainj*binj/(inj_xwide*inj_ywide)

c               --- If the block is completely within the emitting surface,
c               --- set the area to 1.
              else
                area = 1.
              endif

c             --- Take into account the tent function of the charge density,
c             --- i.e. convert the integral of the area to an integral
c             --- of the tent function over the part of the two by two block
c             --- within the emitting surface.  This is an approximate
c             --- calculation.  The conversion was obtained by completely
c             --- integrating over one dimension and partially intregrating
c             --- over the other.
              if (area <= 0.5) then
                area = 2.*area**2
              else
                area = 4.*area - 2.*area**2 -1.
              endif

c             --- Now assign to inj_area. This is done here to simplify the
c             --- coding above, primarily because of the third index here.
c             --- Note that the min works, since inj_ninj is either 1 or
c             --- ninject: if 1, then the min is always 1, if ninject, then
c             --- the min is always ij.
              inj_area(ix,iy,ij) = area

            endif
          enddo
        enddo
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timefill_inj = timefill_inj + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine fill_inj_rz()
      use Subtimers3d
      use InPart
      use InMesh3d
      use InjectVars
      use InjectVars3d

c     --- Fill the injection grid array.  A padding of 2 grid cells is also
c     --- set.
c     --- The array inj_grid holds the axial location of the injection source
c     --- as a function of the transverse coordinates, in the lab frame.
c     --- The array inj_angl holds the angle due to the curvature of the source
c     --- as a function of the transverse coordinates.


      integer(ISZ):: ij,ix
      real(kind=8):: xx,ainj,ai_dx,dxi
      real(kind=8):: xl,xu
      real(kind=8):: ainji
      real(kind=8):: area
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./inj_dx

c Loop over the injection sources.
      do ij=1,ninject
c       --- Calculate the extent of the source plus guard cells.
        ainj    = ainject(ij)
        ainji   = 1./ainj
        ai_dx = ainject(ij) + 2.*inj_dx
c       --- Loop over transverse plane.
        do ix = 0,inj_nx
c           --- transverse coordinate relative to center of source
            xx = ix*inj_dx + inj_xmmin(ij)
c           --- If the point is within source or guard cells, calculate both
c           --- axial distance and angle of curvature.
            if (xx - ai_dx <= 0.) then
c             --- Write the expression for inj_grid so as to avoid
c             --- subtraction of large, similar numbers when rinject is large.
c             inj_grid(ix,0,ij) = zinject(ij) + (rinject(ij) -
c    &               sqrt(max(0.,rinject(ij)**2 - xx**2)))
              inj_grid(ix,0,ij) = xx**2 /
     &         (abs(rinject(ij)) + sqrt(max(0.,rinject(ij)**2 - xx**2)))
              if (rinject(ij) < 0.) inj_grid(ix,0,ij) = -inj_grid(ix,0,ij)
              if (xx**2 > rinject(ij)**2) inj_grid(ix,0,ij) = rinject(ij)

              inj_angl(ix,0,ij) = asin(max(-1.,min(1.,abs(xx)/rinject(ij))))

c             --- Calculate fraction of grid cell which is within the emitting
c             --- surface.  Used to scale the charge density on the emitting
c             --- surface.
c             --- The calculation is done with the emitting area scaled to
c             --- a unit circle.
c             --- The result is the fraction of the two by two block around
c             --- the grid point which is inside the emitting surface.

c             --- First, put the grid point in the first quadrant and scale it.
c             --- The max(0,) chop off any part that crosses the x or y axis.
              xl = max(0.,(abs(xx) - 0.5*inj_xwide*inj_dx)*ainji)
              xu =        (abs(xx) + 0.5*inj_xwide*inj_dx)*ainji

c             --- If the two by two block is completely outside emitting
c             --- surface, set the area to zero.
              if (xl >= 1.) then
                area = 0.

c             --- If two by two block straddles edge of emitting surface, do
c             --- the calculation.
              elseif (xl < 1. .and. 1. < xu) then
                if(l_inj_area) then
                  area = (0.5*(xl*ainj+ainject(ij))*(ainject(ij)-xl*ainj))
     &                 / (xx*inj_dx*inj_xwide)
                else
                  area = 0.
                end if

c               --- If the block is completely within the emitting surface,
c               --- set the area to 1.
              else
                area = 1.
              endif

c             --- Now assign to inj_area. This is done here to simplify the
c             --- coding above, primarily because of the third index here.
c             --- Note that the min works, since inj_ninj is either 1 or
c             --- ninject: if 1, then the min is always 1, if ninject, then
c             --- the min is always ij.
              inj_area(ix,0,min(ij,inj_ninj)) = area

            endif
        enddo
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timefill_inj = timefill_inj + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine inj_transform(np,x,y,z,ni,ijp,tsign,ishift)
      use Subtimers3d
      use InjectVars
      use InjectVars3d
      integer(ISZ):: np,ni
      real(kind=8):: x(np),y(np),z(np)
      integer(ISZ):: ijp(ni),tsign
      integer(ISZ):: ishift

c Transform coordinates from frame of source to frame of lab, or
c vice versa.
c When tsign = +1, transform from source frame to lab frame
c      tsign = -1, transform from lab frame to source frame
c ishift specifies what kind of coordinate shift is done.
c  ishift = 1, data is shifted by grid position
c  ishift = 2, data is shifted by velocity angle
c  for any other value, the shift is skipped.

      real(kind=8):: theta,phi,ct,st,cp,sp
      real(kind=8):: xp,yp,zp
      integer(ISZ):: i,ij
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- If transforming to coordinates of source, subtract source center
      if (tsign < 0 .and. ishift > 0) then
        if (ishift == 1) then
          do i=1,np
            ij = ijp(min(i,ni))
            if (0 < ij .and. ij <= ninject  ) then
              x(i) = x(i) - xinject(ij)
              y(i) = y(i) - yinject(ij)
              z(i) = z(i) - zinject(ij)
            endif
          enddo
        else if (ishift == 2) then
          print*,"Error: invalid ishift value - if you get this error,"
          print*,"       report it to DPGrote@lbl.gov"
          stop
c         --- This code is never used and so is commented out. Note also that
c         --- vzinject now has a second dimension, ns.
c         do i=1,np
c           ij = ijp(min(i,ni))
c           if (0 < ij .and. ij <= ninject  ) then
c             x(i) = x(i) - xpinject(ij)*vzinject(ij)
c             y(i) = y(i) - ypinject(ij)*vzinject(ij)
c             z(i) = z(i) - vzinject(ij)
c           endif
c         enddo
        endif
      endif

      if (maxval(abs(thetainject)) > 0. .or. maxval(abs(phiinject)) > 0.) then
        if (maxval(ijp) == minval(ijp)) then
          ij = ijp(1)
c         --- Convert xpinject and ypinject into the theta and phi angles.
c         --- In the transformation from the source to the lab frame,
c         --- for points along the axis in the source frame, the slopes
c         --- in the lab frame must be preserved. So, for two points
c         --- (0,0,-z) and (0,0,z) in the source frame, the slopes are given by
c         --- xpinject = (x2 - x1)/(z2 - z1) and
c         --- ypinject = (y2 - y1)/(z2 - z1)
c         --- where (x1,y1,z1) and (x2,y2,z2) are the two points rotated into
c         --- the lab frame. Plugging in the transformation equations below
c         --- produces and inverting produces the following result.
c         theta = atan(xpinject(ij))
c         phi = atan(ypinject(ij)*cos(theta))
          theta = thetainject(ij)
          phi = phiinject(ij)
          ct = cos(theta)
          st = sin(theta)
          cp = cos(phi)
          sp = sin(phi)
        endif
        do i=1,np
          ij = ijp(min(i,ni))
          if (0 < ij .and. ij <= ninject  ) then
            if (ni > 1) then
c             --- See comments above.
              theta = thetainject(ij)
              phi = phiinject(ij)
              ct = cos(theta)
              st = sin(theta)
              cp = cos(phi)
              sp = sin(phi)
            endif
            if (tsign == -1) then
c             --- Transform from lab frame to source frame
c             --- This transformation is a rotation about the y axis by theta,
c             --- followed by a rotation about the new x axis by phi.
              xp = +x(i)*ct              - z(i)*st
              yp = -x(i)*st*sp + y(i)*cp - z(i)*ct*sp
              zp = +x(i)*st*cp + y(i)*sp + z(i)*ct*cp
            else
c             --- Transform from source frame to lab frame
c             --- This transformation is a rotation about the x axis by phi,
c             --- followed by a rotation about the new y axis by theta.
              xp = +x(i)*ct - y(i)*st*sp + z(i)*st*cp
              yp =          + y(i)*cp    + z(i)*sp
              zp = -x(i)*st - y(i)*ct*sp + z(i)*ct*cp
            endif
            x(i) = xp
            y(i) = yp
            z(i) = zp
          endif
        enddo
      endif

c     --- If transforming from coordinates of source, add source center
      if (tsign > 0 .and. ishift > 0) then
        if (ishift == 1) then
          do i=1,np
            ij = ijp(min(i,ni))
            if (0 < ij .and. ij <= ninject  ) then
              x(i) = x(i) + xinject(ij)
              y(i) = y(i) + yinject(ij)
              z(i) = z(i) + zinject(ij)
            endif
          enddo
        else if (ishift == 2) then
          print*,"Error: invalid ishift value - if you get this error,"
          print*,"       report it to DPGrote@lbl.gov"
          stop
c         --- This code is never used and so is commented out. Note also that
c         --- vzinject now has a second dimension, ns.
c         do i=1,np
c           ij = ijp(min(i,ni))
c           if (0 < ij .and. ij <= ninject  ) then
c             x(i) = x(i) + xpinject(ij)*vzinject(ij)
c             y(i) = y(i) + ypinject(ij)*vzinject(ij)
c             z(i) = z(i) + vzinject(ij)
c           endif
c         enddo
        endif
      endif

!$OMP MASTER
      if (lw3dtimesubs) timeinj_transform = timeinj_transform + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine getinj_phi()
      use Subtimers3d
      use InjectVars
      use InjectVars3d

c Calculate potential drop from emitting surface at distance
c of dz from the surface.  This is only done for grid cells
c within two grid cells of the elliptical emitting surface, and
c within the axial extent of the grid. This is done over the
c full axial extent since points which are outside of the
c injection region maybe needed for the interlopation below
c to get zp.

      integer(ISZ):: ij,ix,iy,i2x,i2y,i2z
      real(kind=8):: dxi,dyi,inj_dz_tmp
      real(kind=8):: xxsq,yysq,r1sq,r2sq
      real(kind=8):: ainjsqi,binjsqi
      real(kind=8):: xm,ym,rrsq,aa,az
      real(kind=8),allocatable:: xx(:),yy(:),zz(:),pp(:)
      integer(ISZ),allocatable:: xi(:),yi(:),in(:)
      integer(ISZ):: nn,ii,i1,i2
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (ninject == 0) return
      if (inject < 1 .or. inject > 3) return

      allocate(xx((1+inj_nx)*(1+inj_ny)*ninject))
      allocate(yy((1+inj_nx)*(1+inj_ny)*ninject))
      allocate(zz((1+inj_nx)*(1+inj_ny)*ninject))
      allocate(xi((1+inj_nx)*(1+inj_ny)*ninject))
      allocate(yi((1+inj_nx)*(1+inj_ny)*ninject))
      allocate(in((1+inj_nx)*(1+inj_ny)*ninject))
      allocate(pp((1+inj_nx)*(1+inj_ny)*ninject))

      call updateinj_zmminmax()

      dxi = 1./inj_dx
      dyi = 1./inj_dy

      if(inj_nz>1) then
        inj_dz_tmp = inj_dz
        inj_dz = inj_dz0
      end if

c     --- loop over injection sources
c     --- i1 and i2 keep track of the data as new data from each injection
c     --- source is added. i1 is the starting location for the each set of
c     --- data and i2 is the location of the last data point.
      i1 = 1
      do ij=1,ninject
        i2 = i1 - 1

c       --- Set some temporaries.
        ainjsqi = 1./(ainject(ij) + 2.*inj_dx)**2
        binjsqi = 1./(binject(ij) + 2.*inj_dy)**2

        do iy = 0,inj_ny
          do ix = 0,inj_nx
            xm = inj_xmmin(ij) + ix*inj_dx
            ym = inj_ymmin(ij) + iy*inj_dy
            xxsq = xm**2
            yysq = ym**2
            if(linj_rectangle) then
              r1sq = xxsq*ainjsqi
              r2sq = yysq*binjsqi
            else
              r1sq = xxsq*ainjsqi + yysq*binjsqi
              r2sq = 0. 
            endif
            if (r1sq < 1.0 .and. r2sq < 1.0) then

              i2 = i2 + 1

c             --- Save coordinates relative to injection grid.
              xi(i2) = ix
              yi(i2) = iy

c             --- angle of point in transverse plane
              aa = atan2(ym,dvnz(xm))
              az = inj_angl(ix,iy,ij)

c             --- Find coordinates of the point a distance dz in front
c             --- of the source along a line perpendicular to the
c             --- emitting surface.
              xx(i2) = xm - inj_dz*cos(aa)*sin(az)*inj_d(ij)
              yy(i2) = ym - inj_dz*sin(aa)*sin(az)*inj_d(ij)
              zz(i2) = inj_grid(ix,iy,ij) + cos(az)*inj_d(ij)*inj_dz
              in(i2) = ij

            endif
          enddo
        enddo

        call inj_transform(i2-i1+1,xx(i1),yy(i1),zz(i1),1,ij,1,1)

c       --- Select only those data points which are within the grid.
c       --- This also removes particles not within the domain of this
c       --- process.
        nn = i1
        do ii=i1,i2
          if (inj_zmmin <= zz(ii) .and. zz(ii) < inj_zmmax) then
            if (nn < ii) then
              xi(nn) = xi(ii)
              yi(nn) = yi(ii)
              xx(nn) = xx(ii)
              yy(nn) = yy(ii)
              zz(nn) = zz(ii)
              in(nn) = in(ii)
            endif
            nn = nn + 1
          endif
        enddo
        i2 = nn - 1

        i1 = i2 + 1

      enddo
      nn = i1 - 1

      pp = 0.
      call fetchphi(nn,xx,yy,zz,pp)

c     --- Calculate the difference between phi at that point and phi on
c     --- the emitting surface.
      inj_phi = 0.
      do ii=1,nn
        inj_phi(xi(ii),yi(ii),in(ii)) = vinject(in(ii)) - pp(ii)
      enddo
        
      if(inj_nz>1) then
        inj_dz = inj_dz_tmp
        call getinj_phi_mr(nn,xx,yy,zz,xi,yi)
      end if

#ifdef MPIPARALLEL
c     --- The result calculated on each processor is gathered on all
      call parallelnonzerorealarray(inj_phi,(1+inj_nx)*(1+inj_ny)*ninject)
#endif

c     --- Smooth out inj_phi for this emitting surface.
      if (inj_nsmooth > 0) then
        do ij=1,ninject
          call inj_smoother(inj_nx,inj_ny,inj_phi(:,:,ij),inj_dx,inj_dy,
     &                      inj_xmmin(ij),inj_ymmin(ij),
     &                      xinject(ij),yinject(ij),ainject(ij),binject(ij),
     &                      inj_nsmooth)
        enddo
      endif

c     --- Calculate the transverse fields
c     --- Note that sign of inj_ex and inj_ey is not what you think it might
c     --- be since inj_phi is actually vinject - phi, with an extra minus sign.
      if (linj_eperp) then
        do ij = 1,inj_ninj
          do iy = 1,inj_ny-1
            do ix = 1,inj_nx-1
              inj_ex(ix,iy,ij)=(inj_phi(ix+1,iy  ,ij)-inj_phi(ix-1,iy  ,ij))*dxi*0.5
              inj_ey(ix,iy,ij)=(inj_phi(ix  ,iy+1,ij)-inj_phi(ix  ,iy-1,ij))*dyi*0.5
            enddo
          enddo
          do ix = 1,inj_nx-1
            inj_ex(ix,0,ij) = (inj_phi(ix+1,0,ij) - inj_phi(ix-1,0,ij))*dxi*0.5
          enddo
          do iy = 1,inj_ny-1
            inj_ey(0,iy,ij) = (inj_phi(0,iy+1,ij) - inj_phi(0,iy-1,ij))*dyi*0.5
          enddo
        enddo
      endif

      deallocate(xx,yy,zz,pp)
      deallocate(xi,yi,in)

!$OMP MASTER
      if (lw3dtimesubs) timegetinj_phi = timegetinj_phi + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine inj_smoother(inj_nx,inj_ny,inj_phi,inj_dx,inj_dy,
     &                        inj_xmmin,inj_ymmin,
     &                        xinj,yinj,ainj,binj,inj_nsmooth)
      integer(ISZ):: inj_nx,inj_ny
      real(kind=8):: inj_phi(0:inj_nx,0:inj_ny)
      real(kind=8):: inj_dx,inj_dy,inj_xmmin,inj_ymmin,xinj,yinj,ainj,binj
      integer(ISZ):: inj_nsmooth

c This routine smoothes the normal electric field in front of the emitting
c surface, hopefully removing defects  at places where the emitting surface
c crosses grid lines.
c It holds the field at the center point fixed. For points outside of the
c emitter, a linear extrapolation is done using the results of a least
c square fit of the field as a function of radius for the points within
c a grid cell of the edge. The points outside are also held fixed.

      integer(ISZ):: ix,iy,ix0,iy0,is
      integer(ISZ):: ixm1,iym1,ixp1,iyp1
      real(kind=8):: rrsq1,rrsq2
      real(kind=8):: isum,fsum,xsum,xsqsum,fxsum,c1,c2
      real(kind=8):: s(0:inj_nx,0:inj_ny)
      real(kind=8):: ainjsqi1,binjsqi1
      real(kind=8):: ainjsqi2,binjsqi2

c     --- Calculate the parameters for the least square fit.
      isum = 0.
      fsum = 0.
      xsum = 0.
      xsqsum = 0.
      fxsum = 0.
      ainjsqi1 = 1./ainj**2
      binjsqi1 = 1./binj**2
      ainjsqi2 = 1./(ainj - inj_dx)**2
      binjsqi2 = 1./(binj - inj_dx)**2
      do iy = 0,inj_ny
        do ix = 0,inj_nx
          rrsq1 = (ix*inj_dx + inj_xmmin)**2*ainjsqi1 +
     &            (iy*inj_dy + inj_ymmin)**2*binjsqi1
          rrsq2 = (ix*inj_dx + inj_xmmin)**2*ainjsqi2 +
     &            (iy*inj_dy + inj_ymmin)**2*binjsqi2
          if (rrsq1 < 1. .and. rrsq2 > 1.) then
            isum = isum + 1.
            fsum = fsum + inj_phi(ix,iy)
            xsum = xsum + sqrt(rrsq1)
            xsqsum = xsqsum + rrsq1
            fxsum = fxsum + sqrt(rrsq1)*inj_phi(ix,iy)
          endif
        enddo
      enddo
      c1 = (fsum*xsqsum/xsum - fxsum)/(isum*xsqsum/xsum - xsum)
      c2 = (fsum - c1*isum)/xsum

c     --- Set values for points just outside of the emitter edge, using
c     --- the equation from the fit above.
      ainjsqi2 = 1./(ainj + 2.*inj_dx)**2
      binjsqi2 = 1./(binj + 2.*inj_dx)**2
      do iy = 0,inj_ny
        do ix = 0,inj_nx
          rrsq1 = (ix*inj_dx + inj_xmmin)**2*ainjsqi1 +
     &            (iy*inj_dy + inj_ymmin)**2*binjsqi1
          rrsq2 = (ix*inj_dx + inj_xmmin)**2*ainjsqi2 +
     &            (iy*inj_dy + inj_ymmin)**2*binjsqi2
          if (rrsq1 > 1. .and. rrsq2 < 1.) then
            inj_phi(ix,iy) = c1 + c2*sqrt(rrsq1)
          endif
        enddo
      enddo

c     --- Do the smoothing.
      ix0 = -inj_xmmin/inj_dx
      iy0 = -inj_ymmin/inj_dy
      do is = 1,inj_nsmooth
        s = inj_phi
        do iy = 0,inj_ny
          do ix = 0,inj_nx
            rrsq1 = (ix*inj_dx + inj_xmmin)**2*ainjsqi1 +
     &              (iy*inj_dy + inj_ymmin)**2*binjsqi1
            if (rrsq1 < 1. .and. (ix /= ix0 .or. iy /= iy0)) then
              ixm1 = ix - 1
              ixp1 = ix + 1
              iym1 = iy - 1
              iyp1 = iy + 1
              if (ix == 0) ixm1 = 1
              if (ix == inj_nx) ixp1 = inj_nx-1
              if (iy == 0) iym1 = 1
              if (iy == inj_ny) iyp1 = inj_ny-1
              inj_phi(ix,iy) = 0.0625*(s(ixm1,iym1) + s(ixm1,iyp1)  +
     &                                 s(ixp1,iym1) + s(ixp1,iyp1)) +
     &                         0.1250*(s(ixm1,iy  ) + s(ixp1,iy  )  +
     &                                 s(ix  ,iym1) + s(ix  ,iyp1)) +
     &                         0.2500*(s(ix  ,iy  ))
            endif
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine inj_sete(pgroup,ipmin,np,ex,ey,ez)
      use ParticleGroupmodule
      use InjectVars
      type(ParticleGroup):: pgroup
      integer(ISZ):: ipmin,np
      real(kind=8):: ex(np),ey(np),ez(np)

      if (inject <= 0) return

      call inj_sete3d(np,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                pgroup%pid(ipmin,injpid),ex,ey,ez)

      return
      end
c=============================================================================
      subroutine inj_sete3d(npart,xp,yp,zp,pid,ex,ey,ez)
      use Subtimers3d
      use InGen3d
      use InPart
      use InjectVars
      use InjectVars3d
      integer(ISZ):: npart
      real(kind=8):: xp(npart),yp(npart),zp(npart),pid(npart)
      real(kind=8):: ex(npart),ey(npart),ez(npart)

c Calculate the electric field from the grid for particles near the emitting
c surface.  The normal electric field is calculated from the potential drop
c across a length of 'dz' along a line normal to the emitting surface.  The
c field components are then obtained from the normal field. Optionally,
c the tangential fields near the emitting surface can also be included.
c Eventually this routine could be expanded to calculate the E field near
c any conductor surface.

      real(kind=8),allocatable:: xx(:),yy(:),zz(:)
      real(kind=8),allocatable:: tex(:),tey(:),tez(:)
      integer(ISZ),allocatable:: ijp(:)
      integer(ISZ):: ip,ix,iy,iz,ij
      real(kind=8):: rinj,rinji,rr2,wr,dd,ca,zii
      real(kind=8):: aa,az,atx,aty,px,py,pz,wx,wy,wz,en,etx,ety
      real(kind=8):: dxi,dyi,dzi,inj_di
      real(kind=8):: sx,sy
      real(kind=8):: fourthirds,onethird
      real(kind=8):: substarttime,wtime
      integer(ISZ):: spreadx,spready
      if (lw3dtimesubs) substarttime = wtime()

c     --- The user can force the code to skip the special calculation of the 
c     --- E fields for particles behind the virtual surface.
      if (linj_efromgrid) return

      allocate(xx(npart),yy(npart),zz(npart))
      allocate(tex(npart),tey(npart),tez(npart))
      allocate(ijp(npart))

      if(inj_nz>1) call inj_sete_3darray(npart,xp,yp,zp,pid,inj_dz,ex,ey,ez)

      if(solvergeom==XYZgeom .or. solvergeom==RZgeom) then
        if(.not.l_inj_rz) then
          spreadx = 1
          spready = 1
        else
          spreadx = 1
          spready = 0
        end if
      elseif(solvergeom==XZgeom) then
        spreadx = 1
        spready = 0
      elseif(solvergeom==Zgeom) then
        spreadx = 0
        spready = 0
      end if
      if (l_inj_rz) spready = 0
      tex = 0.
      tey = 0.
      tez = 0.

      dxi = 1./inj_dx
      dyi = 1./inj_dy
      dzi = 1./inj_dz
      fourthirds = 4./3.
      onethird = 1./3.

c     --- Transform the particles into the frame of the injection source,
c     --- copying particle data to temporary arrays.
      xx = xp
      yy = yp
      zz = zp
      ijp = int(pid)
      call inj_transform(npart,xx,yy,zz,npart,ijp,-1,1)

c     --- Loop over particles
      do ip=1,npart

c       --- Get number of injection source of particle.  The particle id
c       --- has the number of the injection source the particle was
c       --- emitted from. The int is needed since it is doing double
c       --- duty and has other information stored in the fractional part.
        ij = ijp(ip)
        if (ij <= 0 .or. ij > ninject) cycle
        zii = 1./(inj_dz*inj_d(ij))
        inj_di = 1./inj_d(ij)

c       --- set temporaries
        rinj = rinject(ij)
        rinji = 1./rinj

c       --- Calculate distance of particle from the emitting surface
        rr2 = xx(ip)**2 + yy(ip)**2
c       if (abs(zz(ip)) < abs(rinj)) then
        if (-abs(rinj) < zz(ip) .and. zz(ip) < abs(rinj)) then
c         --- The expression below for zz<rinj is mathematically identical
c         --- to the commented out expression below.  That form is
c         --- used so that for large values of rinj, i.e. a flat emitting
c         --- surface, the correct value of the distance is calculated,
c         --- namely zz.
c         dd = rinj - sqrt(rr2 + (rinj - zz(ip))**2)
          dd = (2.*zz(ip) - (zz(ip)**2 + rr2)*rinji)/
     &         (1. + sqrt(rr2*rinji**2 + (1. - zz(ip)*rinji)**2))
        else
c         --- When zz>rinj, a seperate equation is needed.
          dd = abs(rinj) + sqrt(rr2 + (zz(ip) - rinj)**2)
        endif
        dd = dd*dzi

c       --- Only calculate E-field if particle close to emitting surface.
        if (0.0 <= inj_d(ij)*dd .and. abs(dd) <= abs(inj_d(ij))) then
c         --- The check if a particle is within the radial area of the
c         --- source is probably not needed. It has been removed since there
c         --- are some cases which this check is incorrect - for example
c         --- in a source with a diverging beam. The particles spread out
c         --- beyond the edge of the source and would not pass this check.
c    &    ((xx(ip)*binject(ij))**2 + (yy(ip)*ainject(ij))**2 <=
c    &    (ainject(ij)*binject(ij))**2) .and.
c    &    ((xx(ip)*binjmin(ij))**2 + (yy(ip)*ainjmin(ij))**2 >=
c    &    (ainjmin(ij)*binjmin(ij))**2) ) then

c         --- angle of point in transverse plane
          if(l_inj_rz) then
            aa = 0.
          else
            aa = atan2(yy(ip),dvnz(xx(ip)))
          end if

c         --- angle relative to z axis
          az = atan2(sqrt(rr2)*sign(1.,rinj),(rinj - zz(ip))*sign(1.,rinj))

c         --- Find coordinates of the point on the secondary surface in front
c         --- of the emitting surface along a line perpendicular to the
c         --- emitting surface.
          wr = 1./(1. - inj_dz*dd*rinji)
          px = abs(xx(ip)*wr - inj_xmmin(ij))*dxi
          py = abs(yy(ip)*wr - inj_ymmin(ij))*dyi
          ix = spreadx*int(px)
          iy = spready*int(py)
          wx = spreadx*(px - ix)
          wy = spready*(py - iy)

          if (spreadx*px < 0. .or. px*spreadx > inj_nx .or.
     &        spready*py < 0. .or. py*spready > inj_ny) cycle

c         --- Calculate the normal electric field from the potential drop in
c         --- front of the emitting surface. The normal field given by the
c         --- Child-Langmuir solution (for planar surfaces) is used.
          en = zii*(inj_phi(ix        ,iy        ,ij)*(1.-wx)*(1.-wy) +
     &              inj_phi(ix+spreadx,iy        ,ij)*    wx *(1.-wy) +
     &              inj_phi(ix        ,iy+spready,ij)*(1.-wx)*    wy  +
     &              inj_phi(ix+spreadx,iy+spready,ij)*    wx *    wy   )
          if (linj_enormcl) then
c           --- Scale the normal E field to match the Child-Langmuir solution
            en = en*fourthirds*(dd*inj_di)**onethird
          endif

c         --- Set the particle's electric field based off of the normal
c         --- electric fields.
          if(l_inj_rz) aa = atan2(yy(ip),dvnz(xx(ip)))
          tex(ip) = -en*sin(az)*cos(aa)
          tey(ip) = -en*sin(az)*sin(aa)
          tez(ip) =  en*cos(az)

c         --- Zero out the original data
          ex(ip) = 0.
          ey(ip) = 0.
          ez(ip) = 0.

          if (linj_eperp .and. .not.l_inj_rz) then
c           --- Add in the tangential fields if requested.
c           --- Fetch the transverse field components. The field is linearly
c           --- interpolated between the value of the emitting surface
c           --- (i.e. zero) and the value on the secondary surface.
            atx = atan2(xx(ip)*sign(1.,rinj),(rinj - zz(ip))*sign(1.,rinj))
            aty = atan2(yy(ip)*sign(1.,rinj),(rinj - zz(ip))*sign(1.,rinj))
            sx = 1.
            sy = 1.
            if (xx(ip) < inj_xmmin(ij)) sx = -1.
            if (yy(ip) < inj_ymmin(ij)) sy = -1.
            etx = (dd*inj_di)*sx*(inj_ex(ix        ,iy        ,ij)*(1.-wx)*(1.-wy) +
     &                            inj_ex(ix+spreadx,iy        ,ij)*    wx *(1.-wy) +
     &                            inj_ex(ix        ,iy+spready,ij)*(1.-wx)*    wy  +
     &                            inj_ex(ix+spreadx,iy+spready,ij)*    wx *    wy   )
            ety = (dd*inj_di)*sy*(inj_ey(ix        ,iy        ,ij)*(1.-wx)*(1.-wy) +
     &                            inj_ey(ix+spreadx,iy        ,ij)*    wx *(1.-wy) +
     &                            inj_ey(ix        ,iy+spready,ij)*(1.-wx)*    wy  +
     &                            inj_ey(ix+spreadx,iy+spready,ij)*    wx *    wy   )

            tex(ip) = tex(ip) + etx*cos(atx)
            tey(ip) = tey(ip) + ety*cos(aty)
            tez(ip) = tez(ip) + etx*sin(atx) + ety*sin(aty)
          endif

        endif
      enddo

c     --- Transform E field from injection source frame to lab frame
      call inj_transform(npart,tex,tey,tez,npart,ijp,1,0)

c     --- Copy the temporary arrays into the original. Note that for
c     --- particles near the source, the original was zero-ed out.
c     --- For particles not near the source, the temp arrays are zero.
      ex = ex + tex
      ey = ey + tey
      ez = ez + tez

!$OMP MASTER
      if (lw3dtimesubs) timeinj_sete3d = timeinj_sete3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine inj_addtemp3d(pgroup,npart,ipmin,dz)
      use ParticleGroupmodule
      use Subtimers3d
      use InGen3d
      use InPart
      use InPart3d
      use InjectVars
      use InjectVars3d
      use Particles,Only: tpid
      type(ParticleGroup):: pgroup
      integer(ISZ):: npart,ipmin
      real(kind=8):: dz

c Add temperature to injected particles when they cross virtual surface
c at distance inj_dtemp from emitter. The time of creation array pid(:,tpid)
c is used to flag particles that not have had temperature added yet.

      real(kind=8),allocatable:: xx(:),yy(:),zz(:)
      real(kind=8),allocatable:: vx(:),vy(:),vz(:)
      integer(ISZ),allocatable:: ijp(:)
      integer(ISZ):: ip,ij
      real(kind=8):: rinj,rinji,rr2,dd,ca,zii
      real(kind=8):: aa,az,atx,aty,tx,ty,tz
      real(kind=8):: dzi,inj_di
      real(kind=8):: rnorm
c     real(kind=8):: wrandomgauss
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (inject <= 0) return
      if (.not. l_inj_delay_temp) return

      allocate(xx(npart),yy(npart),zz(npart))
      allocate(vx(npart),vy(npart),vz(npart))
      allocate(ijp(npart))

      dzi = 1./dz

c     --- Transform the particles into the frame of the injection source
      xx = pgroup%xp(ipmin:ipmin+npart-1)
      yy = pgroup%yp(ipmin:ipmin+npart-1)
      zz = pgroup%zp(ipmin:ipmin+npart-1)
      ijp = int(pgroup%pid(ipmin:ipmin+npart-1,injpid))
      call inj_transform(npart,xx,yy,zz,npart,ijp,-1,1)

c     --- Loop over particles
      do ip=1,npart

c       --- if pid(pi,tpid)>0., temperature has already been added
        if(pgroup%pid(ip+ipmin-1,tpid)>0.) cycle

c       --- Get number of injection source of particle.  The particle id
c       --- has the number of the injection source the particle was
c       --- emitted from. The int is needed since it is doing double
c       --- duty and has other information stored in the fractional part.
        ij = ijp(ip)
        if (ij <= 0 .or. ij > ninject) cycle
        zii = 1./(dz*inj_d(ij))
        inj_di = 1./inj_d(ij)

c       --- set temporaries
        rinj = rinject(ij)
        rinji = 1./rinj

c       --- Calculate distance of particle from the emitting surface
        rr2 = xx(ip)**2 + yy(ip)**2
        if (abs(zz(ip)) < abs(rinj)) then
c         --- The expression below for zz<rinj is mathematically identical
c         --- to the commented out expression below.  That form is
c         --- used so that for large values of rinj, i.e. a flat emitting
c         --- surface, the correct value of the distance is calculated,
c         --- namely zz.
c         dd = rinj - sqrt(rr2 + (rinj - zz(ip))**2)
          dd = (2.*zz(ip) - (zz(ip)**2 + rr2)*rinji)/
     &         (1. + sqrt(rr2*rinji**2 + (1. - zz(ip)*rinji)**2))
        else
c         --- When zz>rinj, a seperate equation is needed.
          dd = rinj + sqrt(rr2 + (zz(ip) - rinj)**2)
        endif
        dd = dd*dzi

c       --- Only add temperature if particle close to emitting surface.
        if (.not. (0.0 <= inj_dtemp(ij)*dd .and. abs(dd) <= abs(inj_dtemp(ij)))) then
c         --- The check if a particle is within the radial area of the
c         --- source is probably not needed. It has been removed since there
c         --- are some cases which this check is incorrect - for example
c         --- in a source with a diverging beam. The particles spread out
c         --- beyond the edge of the source and would not pass this check.
c    &    ((xx(ip)*binject(ij))**2 + (yy(ip)*ainject(ij))**2 <=
c    &    (ainject(ij)*binject(ij))**2) .and.
c    &    ((xx(ip)*binjmin(ij))**2 + (yy(ip)*ainjmin(ij))**2 >=
c    &    (ainjmin(ij)*binjmin(ij))**2) ) then

c         --- angle of point in transverse plane
          aa = atan2(yy(ip),dvnz(xx(ip)))

c         --- angle relative to z axis
          az = asin(sqrt(rr2)*rinji)
c         az = atan2(sqrt(rr2),(rinj - zz(ip))*sign(1.,rinj)) XXX

          atx = asin(xx(ip)/sqrt(rinj**2 - yy(ip)**2))
          aty = asin(yy(ip)/sqrt(rinj**2 - xx(ip)**2))

c         --- Add the particle's temperature
c         --- These should be used but the index ii to pass in is not easy
c         --- to obtain.
c         wrandomgauss(vtrandom,ii,dig3,dig4,1,1,.false.)
c         wrandomgauss(vtrandom,ii,dig5,dig6,1,1,.false.)
c         wrandomgauss(vzrandom,ii,dig7,dig8,1,1,.false.)
          tx = vthperp_s(ij)*rnorm()
          ty = vthperp_s(ij)*rnorm()
          tz = vthz_s(ij)*rnorm()
          if(l_inj_addtempz_abs) tz=abs(tz)
          vx(ip) = - tz*sin(az)*cos(aa) + tx*cos(atx)
          vy(ip) = - tz*sin(az)*sin(aa) + ty*cos(aty)
          vz(ip) = + tz*cos(az)         + tx*sin(atx) + ty*sin(aty)

c         --- reverse sign of pid(,tpid) indicating that temperature has been added
          pgroup%pid(ip+ipmin-1,tpid) = -pgroup%pid(ip+ipmin-1,tpid)

        endif
      enddo

c     --- Transform thermal velocities to lab frame
      call inj_transform(npart,vx,vy,vz,npart,ijp,1,0)
      pgroup%uxp(ipmin:ipmin+npart-1) = pgroup%uxp(ipmin:ipmin+npart-1) + vx
      pgroup%uyp(ipmin:ipmin+npart-1) = pgroup%uyp(ipmin:ipmin+npart-1) + vy
      pgroup%uzp(ipmin:ipmin+npart-1) = pgroup%uzp(ipmin:ipmin+npart-1) + vz

      deallocate(xx,yy,zz)
      deallocate(vx,vy,vz)
      deallocate(ijp)

!$OMP MASTER
      if (lw3dtimesubs) timeinj_addtemp3d = timeinj_addtemp3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine loadrho_inject(pgroup)
      use ParticleGroupmodule
      use InjectVars
      use InjectVars3d,Only: inj_nz
      use InGen3d

      type(ParticleGroup):: pgroup
      if (inject == 0) return

      if(inj_nz>1) call inj_setrhomr(pgroup)

      return
      end
c=============================================================================
      subroutine inj_setrho3d(pgroup,dz,l2symtry,l4symtry)
      use ParticleGroupmodule
      use Subtimers3d
      use InPart
      use InjectVars
      use InjectVars3d
      type(ParticleGroup):: pgroup
      real(kind=8):: dz
      logical(ISZ):: l2symtry,l4symtry

c Calculate the charge density on the emitting surface.  The particles are
c mapped onto the surface assuming a spherical emitter.  The charge density
c scale to account for symmetries and to account for the fraction of the
c grid cell within emitting surface.
c The scaling is done so that the charge density accurately represents
c the charge divided by the volume filled by that charge.  This only
c affects the edge of the emitting surface where only part of a grid cell
c will have charge in it, but the charge is divided by the area of all of
c the grid cell.

      integer(ISZ),allocatable:: ijp(:)
      integer(ISZ):: is,ip,ix,iy,ij,ij1,i1,i2
      real(kind=8):: rinj,rinji,rr2,dd
      real(kind=8):: aa,az,px,py,wx1,wy1,wx0,wy0,ww
      real(kind=8):: dxi,dyi,dzi
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      allocate(ijp(pgroup%npmax))

      dxi = 1./inj_dx
      dyi = 1./inj_dy
      dzi = 1./dz

c     --- zero out the array
      call zeroarry(inj_rho(0,0,1),(1+inj_nx)*(1+inj_ny)*inj_ninj)

c     --- Loop over species and particles
      do is=1,ns
        i1 = pgroup%ins(is)
        i2 = pgroup%ins(is)+pgroup%nps(is)-1

c       --- Transform the particles into the frame of the injection source
        if (pgroup%nps(is) > 0) then
          ijp(i1:i2) = int(pgroup%pid(i1:i2,injpid))
          call inj_transform(pgroup%nps(is),
     &                    pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                    pgroup%nps(is),ijp(i1:i2),-1,1)
        endif

        ww = pgroup%sq(is)*pgroup%sw(is)*dxi*dyi*dzi
        do ip=i1,i2
          if (pgroup%uzp(ip) == 0.) cycle

c         --- Get number of injection source of particle.  The particle id
c         --- has the number of the injection source the particle was
c         --- emitted from. The int is needed since it is doing double
c         --- duty and has other information stored in the fractional part.
          ij = int(pgroup%pid(ip,injpid))

c         --- Skip particles that did not come from an injection source.
c         --- This is not necessarily the correct thing to do, but at least
c         --- for now, it avoids memory problems and seg faults.
c         --- Ultimately, a particle should be able to check if its near
c         --- any conductor and do special coding for the field near it.
          if (ij <= 0 .or. ij > ninject) cycle

          ix = int(abs(pgroup%xp(ip) - inj_xmmin(ij))*dxi)
          iy = int(abs(pgroup%yp(ip) - inj_ymmin(ij))*dyi)

c         --- set temporaries
          rinji = 1./rinject(ij)
          rinj = rinject(ij)

c         --- Calculate distance of particle from the emitting surface
          rr2 = pgroup%xp(ip)**2 + pgroup%yp(ip)**2
          if (abs(pgroup%zp(ip)) < abs(rinj)) then
c           --- The expression below for zp(ip)<rinj is mathematically identical
c           --- to the commented out expression below.  That form is
c           --- used so that for large values of rinj, i.e. a flat emitting
c           --- surface, the correct value of the distance is calculated,
c           --- namely zp(ip).
c           dd = rinj - sqrt(rr2 + (rinj - zp(ip))**2)
            dd = (2.*pgroup%zp(ip) - (pgroup%zp(ip)**2 + rr2)*rinji)/
     &           (1. + sqrt(rr2*rinji**2 + (1. - pgroup%zp(ip)*rinji)**2))
            if (rinj < 0.) dd = -dd
          else
c           --- When zp(ip)>rinj, a seperate equation is needed.
            dd = abs(rinj) + sqrt(rr2 + (pgroup%zp(ip) - rinj)**2)
          endif
          dd = dd*dzi

c         --- Only deposit charge of particles close to emitting surface.
          if (abs(dd) < 1. .and. pgroup%uzp(ip) /= 0.) then

c           --- angle of point in transverse plane
            aa = atan2(pgroup%yp(ip),dvnz(pgroup%xp(ip)))

c           --- angle relative to z axis
            az = atan2(sqrt(rr2),rinj-pgroup%zp(ip))

c           --- Map the particle position onto the emitting surface.
            px = abs(rinj*sin(az)*cos(aa) - inj_xmmin(ij))*dxi
            py = abs(rinj*sin(az)*sin(aa) - inj_ymmin(ij))*dyi
            ix = int(px)
            iy = int(py)
            wx1 = px - ix
            wy1 = py - iy
            wx0 = 1. - wx1
            wy0 = 1. - wy1

c           --- Deposit the particle's charge onto the emitting surface.
            ij1 = min(ij,inj_ninj)
            inj_rho(ix  ,iy  ,ij1)=inj_rho(ix  ,iy  ,ij1) + ww*wx0*wy0*(1.-dd)
            inj_rho(ix+1,iy  ,ij1)=inj_rho(ix+1,iy  ,ij1) + ww*wx1*wy0*(1.-dd)
            inj_rho(ix  ,iy+1,ij1)=inj_rho(ix  ,iy+1,ij1) + ww*wx0*wy1*(1.-dd)
            inj_rho(ix+1,iy+1,ij1)=inj_rho(ix+1,iy+1,ij1) + ww*wx1*wy1*(1.-dd)

          endif
        enddo

c       --- Transform the particles back into the lab frame
        if (pgroup%nps(is) > 0) then
          call inj_transform(pgroup%nps(is),
     &                    pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                    pgroup%nps(is),ijp(i1:i2),1,1)
        endif

      enddo

c     --- Loop over unique injection sources
      do ij = 1,inj_ninj

c       --- Scale the charge density by the inverse of the fraction of the
c       --- contributing area which is within the emitting surface.
c       --- The factor of two is needed since particles are contributing to
c       --- rho only on the right side of the emitting surface.
        where (inj_area(:,:,ij) > 0.)
          inj_rho(:,:,ij) = 2.*inj_rho(:,:,ij)/inj_area(:,:,ij)
        end where

      enddo

      deallocate(ijp)

!$OMP MASTER
      if (lw3dtimesubs) timeinj_setrho3d = timeinj_setrho3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine inj_setrho3d_z(pgroup,dz,nz)
      use ParticleGroupmodule
      use Subtimers3d
      use InPart
      use InjectVars
      use InjectVars3d
      use Constant
      type(ParticleGroup):: pgroup
      integer(ISZ):: nz
      real(kind=8):: dz

c Calculate the charge density on the emitting surface.  The particles are
c mapped onto the surface assuming a spherical emitter.  The charge density
c scale to account for symmetries and to account for the fraction of the
c grid cell within emitting surface.
c The scaling is done so that the charge density accurately represents
c the charge divided by the volume filled by that charge.  This only
c affects the edge of the emitting surface where only part of a grid cell
c will have charge in it, but the charge is divided by the area of all of
c the grid cell.

      integer(ISZ),allocatable:: ijp(:)
      integer(ISZ):: is,ip,ij,j,ij1,i1,i2
      real(kind=8):: dd
      real(kind=8):: ww
      real(kind=8):: dzi
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      allocate(ijp(pgroup%npmax))

      dzi = 1./dz

c     --- zero out the array
      call zeroarry(inj_rho(0,0,1),inj_ninj)

c     --- Loop over species and particles
      do is=1,ns
        i1 = pgroup%ins(is)
        i2 = pgroup%ins(is)+pgroup%nps(is)-1

c       --- Transform the particles into the frame of the injection source
c       --- This is only really here to be consistent with the other routines.
c       --- All that this should do is add zinject to zp.
        if (pgroup%nps(is) > 0) then
          ijp(i1:i2) = int(pgroup%pid(i1:i2,injpid))
          call inj_transform(pgroup%nps(is),
     &                    pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                    pgroup%nps(is),ijp(i1:i2),-1,1)
        endif

        ww = 2.*pgroup%sq(is)*pgroup%sw(is)*dzi
        do ip=i1,i2
          if (pgroup%uzp(ip) == 0.) cycle

c         --- Get number of injection source of particle.  The particle id
c         --- has the number of the injection source the particle was
c         --- emitted from. The int is needed since it is doing double
c         --- duty and has other information stored in the fractional part.
          ij = int(pgroup%pid(ip,injpid))

c         --- Skip particles that did not come from an injection source.
c         --- This is not necessarily the correct thing to do, but at least
c         --- for now, it avoids memory problems and seg faults.
c         --- Ultimately, a particle should be able to check if its near
c         --- any conductor and do special coding for the field near it.
          if (ij <= 0 .or. ij > ninject) cycle

c         --- Calculate distance of particle from the emitting surface
          dd = pgroup%zp(ip)*dzi

c         --- Only deposit charge of particles close to emitting surface.
          if (dd < 1.) then

c           --- Deposit the particle's charge onto the emitting surface.
            ij1 = min(ij,inj_ninj)
            inj_rho(0,0,ij1)=inj_rho(0,0,ij1)+ww*(1.-dd)

          endif
        enddo

c       --- Transform the particles back into the lab frame
        if (pgroup%nps(is) > 0) then
          call inj_transform(pgroup%nps(is),
     &                    pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                    pgroup%nps(is),ijp(i1:i2),1,1)
        endif

      enddo

      deallocate(ijp)

!$OMP MASTER
      if (lw3dtimesubs) timeinj_setrho3d = timeinj_setrho3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine getinj_phi_3d()
      use Subtimers3d
      use InjectVars
      use InjectVars3d

c Calculate potential from emitting surface at distance
c of dz from the surface.  This is only done for grid cells
c within two grid cells of the elliptical emitting surface, and
c within the axial extent of the grid. This is done over the
c full axial extent since points which are outside of the
c injection region maybe needed for the interlopation below
c to get zp.

      integer(ISZ):: ij,ix,iy,iz,i2x,i2y,i2z
      real(kind=8):: dxi,dyi,inj_dz_tmp
      real(kind=8):: xxsq,yysq,r1sq,r2sq
      real(kind=8):: ainjsqi,binjsqi
      real(kind=8):: xm,ym,rrsq,aa,az
      real(kind=8),allocatable:: xx(:),yy(:),zz(:),pp(:)
      integer(ISZ),allocatable:: xi(:),yi(:),zi(:),in(:)
      real(kind=8):: xg(0:inj_nx),yg(0:inj_ny),zg(0:inj_nz)
      real(kind=8):: zinj_grid,onethird,dz0,dz1,dzmin,zdist
      integer(ISZ):: nn,ii,i1,i2

      allocate(xx((1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject))
      allocate(yy((1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject))
      allocate(zz((1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject))
      allocate(xi((1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject))
      allocate(yi((1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject))
      allocate(zi((1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject))
      allocate(in((1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject))
      allocate(pp((1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject))

      call updateinj_zmminmax()

      dxi = 1./inj_dx
      dyi = 1./inj_dy

c     --- loop over injection sources
c     --- i1 and i2 keep track of the data as new data from each injection
c     --- source is added. i1 is the starting location for the each set of
c     --- data and i2 is the location of the last data point.
      i1 = 1
      do ij=1,ninject
        i2 = i1 - 1

c       --- Set some temporaries.
        ainjsqi = 1./(ainject(ij) + 2.*inj_dx)**2
        binjsqi = 1./(binject(ij) + 2.*inj_dy)**2

        onethird = 1./3.
        dzmin = inj_d(ij)*inj_dz0/(real(inj_nz-1)*2.**onethird+2.-real(inj_nz))**3

        zg = 0.
        zg(0) = inj_d(ij)*inj_dz0
        do iz = 1, inj_nz-1
          if(iz==1) then
            dz0 = dzmin
            dz1 = dzmin
            zinj_grid = dz0
          else
            dz0 = dz1
            zinj_grid = zinj_grid+dz0
            dz1 = inj_d(ij)*inj_dz0*(1.-(real(inj_nz-iz-1)/real(inj_nz-1))*(1.-dzmin**onethird))**3
     &          - zinj_grid
          end if
          zg(iz) = zg(iz-1)-dz0
        end do

        do ix = 0,inj_nx
          xg(ix) = inj_xmmin(ij) + ix*inj_dx
        end do
        do iy = 0,inj_ny
          yg(iy) = inj_ymmin(ij) + iy*inj_dy
        end do

        do iz = 0, inj_nz
          do iy = 0,inj_ny
            ym = yg(iy)
            do ix = 0,inj_nx
              xm = xg(ix)

                i2 = i2 + 1

c               --- Save coordinates relative to injection grid.
                xi(i2) = ix
                yi(i2) = iy
                zi(i2) = iz

c               --- angle of point in transverse plane
                aa = atan2(ym,dvnz(xm))
                az = inj_angl(ix,iy,ij)

                zdist = inj_dz0*inj_d(ij)-zg(iz)
c               --- Find coordinates of the point a distance dz in front
c               --- of the source along a line perpendicular to the
c               --- emitting surface.
                xx(i2) = xm - cos(aa)*sin(az)*zdist
                yy(i2) = ym - sin(aa)*sin(az)*zdist
                zz(i2) = inj_grid(ix,iy,ij) + cos(az)*zdist
                in(i2) = ij

c              endif
            enddo
          enddo
        enddo

        call inj_transform(i2-i1+1,xx(i1),yy(i1),zz(i1),1,ij,1,1)

c       --- Select only those data points which are within the grid.
c       --- This also removes particles not within the domain of this
c       --- process.
        nn = i1
        do ii=i1,i2
          if (inj_zmmin <= zz(ii) .and. zz(ii) < inj_zmmax .and. nn < ii) then
            xi(nn) = xi(ii)
            yi(nn) = yi(ii)
            zi(nn) = zi(ii)
            xx(nn) = xx(ii)
            yy(nn) = yy(ii)
            zz(nn) = zz(ii)
            in(nn) = in(ii)
            nn = nn + 1
          endif
        enddo
        i2 = nn - 1

        i1 = i2 + 1

      enddo
      nn = i1 - 1

      pp = 0.
      call fetchphi(nn,xx,yy,zz,pp)

c     --- Calculate the difference between phi at that point and phi on
c     --- the emitting surface.
      inj_phi_3d = 0.
      do ii=1,nn
        inj_phi_3d(xi(ii),yi(ii),zi(ii),in(ii)) = pp(ii)
      enddo
        
      deallocate(xx,yy,zz,xi,yi,zi,in,pp)

#ifdef MPIPARALLEL
c     --- The result calculated on each processor is gathered on all
      call parallelnonzerorealarray(inj_phi_3d,(1+inj_nx)*(1+inj_ny)*(1+inj_nz)*ninject)
#endif

      return
      end
c=============================================================================
      subroutine inj_setrhomr(pgroup)
      use ParticleGroupmodule
      use Subtimers3d
      use Constant,Only: pi
      use InPart,Only: ns
      use Particles,Only: wpid
      use InGen3d
      use InjectVars,Only: ninject,inj_d,injpid,rinject
      use InjectVars3d
      type(ParticleGroup):: pgroup

c Calculate the charge density in the region of the injector mesh refinement.

      integer(ISZ), allocatable:: ijp(:)
      integer(ISZ):: i1,i2
      real(kind=8), allocatable, dimension(:) :: dz_local
      integer(ISZ):: sx,sy
      real(kind=8):: dxi,dyi,dzmin,dz0,dz1,zinj_grid,onethird,rinj,rinji
      real(kind=8):: px,py,wx0,wx1,wy0,wy1,wz0,wz1,rr2,dd,aa,az
      real(kind=8):: deltaz(ninject)
      integer(ISZ):: ix,iy,iz,is,ij,ij1,ip
      real(kind=8):: wws,ww
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      allocate(ijp(pgroup%npmax))
      allocate(dz_local(0:inj_nz))

      if(solvergeom==XYZgeom) then
        sx = 1
        sy = 1
      elseif(solvergeom==RZgeom) then
        sx = 1
        if(l_inj_rz) then
          sy = 0
        else
          sy = 1
        end if
      elseif(solvergeom==XZgeom) then
        sx = 1
        sy = 0
      elseif(solvergeom==Zgeom) then
        sx = 0
        sy = 0
      endif

      onethird = 1./3.

      dzmin = 1./(real(inj_nz-1)*2.**onethird+2.-real(inj_nz))**3
      deltaz = inj_d*inj_dz0

      dxi = 1./inj_dx
      dyi = 1./inj_dy

      dz_local(0) = dzmin
      do iz = 1, inj_nz-1
        if(iz==1) then
          dz0 = dzmin
          dz1 = dzmin
          zinj_grid = dz0
        else
          dz0 = dz1
          zinj_grid = zinj_grid+dz0
          dz1 = (1.-(real(inj_nz-iz-1)/real(inj_nz-1))*(1.-dzmin**onethird))**3
     &        - zinj_grid
        end if
        dz_local(iz) = dz1
      end do
      dz_local(inj_nz) = dz_local(inj_nz-1)

c     --- zero out the array
      inj_q = 0.

c     --- Loop over species and particles
      if(.not. l_inj_rz) then
        do is=1,pgroup%ns
          i1 = pgroup%ins(is)
          i2 = pgroup%ins(is)+pgroup%nps(is)-1
c         --- Transform the particles into the frame of the injection source
          if (pgroup%nps(is) > 0) then
            ijp(i1:i2) = int(pgroup%pid(i1:i2,injpid))
            call inj_transform(pgroup%nps(is),
     &                   pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                   pgroup%nps(is),ijp(i1:i2),-1,1)
          endif
          wws = pgroup%sq(is)*pgroup%sw(is)*dxi*dyi
          do ip=i1,i2
            if (pgroup%uzp(ip) == 0.) cycle
            if(wpid>0) then
              ww=wws*pgroup%pid(ip,wpid)
            else
              ww=wws
            end if

c           --- Get number of injection source of particle.  The particle id
c           --- has the number of the injection source the particle was
c           --- emitted from. The int is needed since it is doing double
c           --- duty and has other information stored in the fractional part.
            ij = int(pgroup%pid(ip,injpid))

c           --- Skip particles that did not come from an injection source.
c           --- This is not necessarily the correct thing to do, but at least
c           --- for now, it avoids memory problems and seg faults.
c           --- Ultimately, a particle should be able to check if its near
c           --- any conductor and do special coding for the field near it.
            if (ij <= 0 .or. ij > ninject) cycle

            ix = int(abs(pgroup%xp(ip) - inj_xmmin(ij))*dxi)
            iy = int(abs(pgroup%yp(ip) - inj_ymmin(ij))*dyi)

c           --- set temporaries
            rinji = 1./rinject(ij)
            rinj = rinject(ij)

c           --- Calculate distance of particle from the emitting surface
            rr2 = pgroup%xp(ip)**2 + pgroup%yp(ip)**2
            if (abs(pgroup%zp(ip)) < abs(rinj)) then
c             --- The expression below for zp(ip)<rinj is mathematically identical
c             --- to the commented out expression below.  That form is
c             --- used so that for large values of rinj, i.e. a flat emitting
c             --- surface, the correct value of the distance is calculated,
c             --- namely zp(ip).
c             dd = rinj - sqrt(rr2 + (rinj - zp(ip))**2)
              dd = (2.*pgroup%zp(ip) - (pgroup%zp(ip)**2 + rr2)*rinji)/
     &             (1. + sqrt(rr2*rinji**2 + (1. - pgroup%zp(ip)*rinji)**2))
              if (rinj < 0.) dd = -dd
            else
c             --- When zp(ip)>rinj, a seperate equation is needed.
              dd = abs(rinj) + sqrt(rr2 + (pgroup%zp(ip) - rinj)**2)
            endif
            dd = dd/deltaz(ij)

c           --- Only deposit charge of particles in emitting area
            if (abs(dd) < 1.) then

c             --- angle of point in transverse plane
              aa = atan2(pgroup%yp(ip),dvnz(pgroup%xp(ip)))

c             --- angle relative to z axis
              az = atan2(sqrt(rr2),rinj-pgroup%zp(ip))

c             --- Map the particle position onto the emitting surface.
              px = abs(rinj*sin(az)*cos(aa) - inj_xmmin(ij))*dxi
              py = abs(rinj*sin(az)*sin(aa) - inj_ymmin(ij))*dyi
              ix = int(px)
              iy = int(py)
              wx1 = px - ix
              wy1 = py - iy
              wx0 = 1. - wx1
              wy0 = 1. - wy1
              if(dd>dzmin) then
                iz  = 1+int(real(inj_nz-1)*(1.-(1.-dd**onethird)/(1.-dzmin**onethird)))
                wz1 = (dd-(1.-((real(inj_nz-iz))/real(inj_nz-1))*(1.-dzmin**onethird))**3)/dz_local(iz)
              else
c               prevent deposition in cell closer to emitter (gives better result on 1D Lampel-Tiefenback test)
                IF(l_inj_no_rho_on_emit) cycle
                iz  = 0
                wz1 = dd/dzmin
              end if
              wz0 = 1.-wz1

c             --- Deposit the particle's charge onto the emitting region.
              ij1 = min(ij,inj_ninj)
              inj_q(ix   ,iy   ,iz  , ij1) = inj_q(ix   ,iy   ,iz  ,ij1) + ww*wx0*wy0*wz0
              inj_q(ix+sx,iy   ,iz  , ij1) = inj_q(ix+sx,iy   ,iz  ,ij1) + ww*wx1*wy0*wz0
              inj_q(ix   ,iy+sy,iz  , ij1) = inj_q(ix   ,iy+sy,iz  ,ij1) + ww*wx0*wy1*wz0
              inj_q(ix+sx,iy+sy,iz  , ij1) = inj_q(ix+sx,iy+sy,iz  ,ij1) + ww*wx1*wy1*wz0
              inj_q(ix   ,iy   ,iz+1, ij1) = inj_q(ix   ,iy   ,iz+1,ij1) + ww*wx0*wy0*wz1
              inj_q(ix+sx,iy   ,iz+1, ij1) = inj_q(ix+sx,iy   ,iz+1,ij1) + ww*wx1*wy0*wz1
              inj_q(ix   ,iy+sy,iz+1, ij1) = inj_q(ix   ,iy+sy,iz+1,ij1) + ww*wx0*wy1*wz1
              inj_q(ix+sx,iy+sy,iz+1, ij1) = inj_q(ix+sx,iy+sy,iz+1,ij1) + ww*wx1*wy1*wz1

            endif
          enddo

c         --- Transform the particles back into the lab frame
          if (pgroup%nps(is) > 0) then
            call inj_transform(pgroup%nps(is),
     &                    pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                    pgroup%nps(is),ijp(i1:i2),1,1)
          endif
        enddo

        if (solvergeom/=Zgeom) then
c         --- Loop over unique injection sources
         do ij1 = 1,inj_ninj

c         --- Scale the charge density by the inverse of the fraction of the/
c         --- contributing area which is within the emitting surface.
            do iy=0,inj_ny
              do ix=0,inj_nx
                if (inj_area(ix,iy,ij1) > 0.)
     &            inj_q(ix,iy,:,ij1) = inj_q(ix,iy,:,ij1)/inj_area(ix,iy,ij1)
              enddo
            enddo
          enddo
        end if

      else !l_inj_rz=.true.
        do is=1,pgroup%ns
          i1 = pgroup%ins(is)
          i2 = pgroup%ins(is)+pgroup%nps(is)-1
c         --- Transform the particles into the frame of the injection source
          if (pgroup%nps(is) > 0) then
            ijp(i1:i2) = int(pgroup%pid(i1:i2,injpid))
            call inj_transform(pgroup%nps(is),
     &                    pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                    pgroup%nps(is),ijp(i1:i2),-1,1)
          endif
          wws = pgroup%sq(is)*pgroup%sw(is)
          do ip=i1,i2
            if (pgroup%uzp(ip) == 0.) cycle
            if(wpid>0) then
              ww=wws*pgroup%pid(ip,wpid)
            else
              ww=wws
            end if
            rr2 = pgroup%xp(ip)**2 + pgroup%yp(ip)**2
            ix = int(sqrt(rr2)*dxi)

c           --- Get number of injection source of particle.  The particle id
c           --- has the number of the injection source the particle was
c           --- emitted from. The int is needed since it is doing double
c           --- duty and has other information stored in the fractional part.
            ij = int(pgroup%pid(ip,injpid))

c           --- Skip particles that did not come from an injection source.
c           --- This is not necessarily the correct thing to do, but at least
c           --- for now, it avoids memory problems and seg faults.
c           --- Ultimately, a particle should be able to check if its near
c           --- any conductor and do special coding for the field near it.
            if (ij <= 0 .or. ij > ninject) cycle

c           --- set temporaries
            rinji = 1./rinject(ij)
            rinj = rinject(ij)

c           --- Calculate distance of particle from the emitting surface
            if (abs(pgroup%zp(ip)) < abs(rinj)) then
c             --- The expression below for zp(ip)<rinj is mathematically identical
c             --- to the commented out expression below.  That form is
c             --- used so that for large values of rinj, i.e. a flat emitting
c             --- surface, the correct value of the distance is calculated,
c             --- namely zp(ip).
c             dd = rinj - sqrt(rr2 + (rinj - zp(ip))**2)
              dd = (2.*pgroup%zp(ip) - (pgroup%zp(ip)**2 + rr2)*rinji)/
     &             (1. + sqrt(rr2*rinji**2 + (1. - pgroup%zp(ip)*rinji)**2))
            else
c             --- When zp(ip)>rinj, a seperate equation is needed.
              dd = abs(rinj) + sqrt(rr2 + (pgroup%zp(ip) - rinj)**2)
            endif
            dd = dd/deltaz(ij)

c           --- Only deposit charge of particles in emitting area
            if (abs(dd) < 1.) then

c             --- angle relative to z axis
              az = atan2(sqrt(rr2),rinj-pgroup%zp(ip))

c             --- Map the particle position onto the emitting surface.
              px = abs(rinj*sin(az))*dxi
              ix = int(px)
              wx1 = px - ix
              wx0 = 1. - wx1
              if(dd>dzmin) then
                iz  = 1+int(real(inj_nz-1)*(1.-(1.-dd**onethird)/(1.-dzmin**onethird)))
                wz1 = (dd-(1.-((real(inj_nz-iz))/real(inj_nz-1))*(1.-dzmin**onethird))**3)/dz_local(iz)
              else
c               prevent deposition in cell closer to emitter (gives better result on 1D Lampel-Tiefenback test)
                IF(l_inj_no_rho_on_emit) cycle
                iz  = 0
                wz1 = dd/dzmin
              end if
              wz0 = 1.-wz1

c             --- Deposit the particle's charge onto the emitting region.
              ij1 = min(ij,inj_ninj)
              inj_q(ix   ,0   ,iz  , ij1) = inj_q(ix   ,0   ,iz  ,ij1) + ww*wx0*wz0
              inj_q(ix+sx,0   ,iz  , ij1) = inj_q(ix+sx,0   ,iz  ,ij1) + ww*wx1*wz0
              inj_q(ix   ,0   ,iz+1, ij1) = inj_q(ix   ,0   ,iz+1,ij1) + ww*wx0*wz1
              inj_q(ix+sx,0   ,iz+1, ij1) = inj_q(ix+sx,0   ,iz+1,ij1) + ww*wx1*wz1

            endif
          enddo

c         --- Transform the particles back into the lab frame
          if (pgroup%nps(is) > 0) then
            call inj_transform(pgroup%nps(is),
     &                      pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                      pgroup%nps(is),ijp(i1:i2),1,1)
          endif
        enddo

c       --- Loop over unique injection sources
        do ij1 = 1,inj_ninj
c         --- Scale the charge density by the inverse of the fraction of the
c         --- contributing area which is within the emitting surface.
          do ix=0,inj_nx
            if (inj_area(ix,0,ij1) > 0.)
     &        inj_q(ix,0,:,ij1) = inj_q(ix,0,:,ij1) / inj_area(ix,0,ij1)
            if(ix==0) then
              inj_q(ix,0,:,ij1) = inj_q(ix,0,:,ij1) / (pi*0.25*inj_dx**2)
            else
              inj_q(ix,0,:,ij1) = inj_q(ix,0,:,ij1) / (2.*pi*ix*inj_dx**2)
            end if
          enddo
        enddo

      end if

c     --- assign charge deposited on emitter surface (slice 0) to slice 1
      IF(.not. l_inj_no_rho_on_emit) then
        inj_q(:,:,1,:) = inj_q(:,:,1,:) + inj_q(:,:,0,:)
        inj_q(:,:,0,:) = 0.
      end if

!$OMP MASTER
      if (lw3dtimesubs) timeinj_setrhomr = timeinj_setrhomr + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine getinj_phi_mr(nn,xxp,yyp,zzp,xi,yi)
      use Subtimers3d
      use InPart
      use InGen
      use InGen3d
      use InjectVars
      use InjectVars3d
      use InMesh3d
      use Constant
      use Fields3d
      use Fields3dParticles
      use Efields3d
      use Picglb
      use Picglb3d
      integer(ISZ), intent(in) :: nn
      real(kind=8), dimension(nn), intent(in) :: xxp, yyp, zzp
      integer(ISZ), dimension(nn), intent(in) :: xi, yi

c Calculate the phi in the emitting region by doing 1-d tridiagonal solves
c of Poisson's equation along lines normal to the emitting surface.

      integer(ISZ):: i,is,ip,ix,iy,iz,ij,ij1,nrhs,irhs,kd,ldab,ldb,n,info
      real(kind=8):: rinj,rinji,rr2,dd,deltaz(ninject)
      real(kind=8):: aa,az,px,py,wx1,wy1,wz1,wx0,wy0,wz0,ww,wws,dz0,dz1,zinj_grid
      real(kind=8):: dxi,dyi,dzmin,ainjsqi,binjsqi,rrsq,onethird,ext,eyt,ezt,atx,aty
      real(kind=8), allocatable, dimension(:,:) :: rhs, ab
      real(kind=8), allocatable, dimension(:) :: ex_tmp, ey_tmp, ez_tmp
      real(kind=8), allocatable, dimension(:) :: dz_local
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      allocate(ex_tmp(nn),ey_tmp(nn),ez_tmp(nn))
      allocate(dz_local(0:inj_nz))

      onethird = 1./3.

      dzmin = 1./(real(inj_nz-1)*2.**onethird+2.-real(inj_nz))**3
      deltaz = inj_d*inj_dz0

      dxi = 1./inj_dx
      dyi = 1./inj_dy

      dz_local(0) = dzmin
      do iz = 1, inj_nz-1
        if(iz==1) then
          dz0 = dzmin
          dz1 = dzmin
          zinj_grid = dz0
        else
          dz0 = dz1
          zinj_grid = zinj_grid+dz0
          dz1 = (1.-(real(inj_nz-iz-1)/real(inj_nz-1))*(1.-dzmin**onethird))**3
     &        - zinj_grid
        end if
        dz_local(iz) = dz1
      end do
      dz_local(inj_nz) = dz_local(inj_nz-1)

c     --- initialize inj_phi_3d
      call getinj_phi_3d()

      allocate(rhs(1:inj_nz-1,nn),ab(2,1:inj_nz-1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"getinj_phi_mr: allocation error ",allocerror,
     &         ": could not allocate rhs and ab to shape ",inj_nz,nn
        stop
      endif

      do ij = 1, inj_ninj
        nrhs = nn
c       --- set temporaries
        rinji = 1./rinject(ij)
        rinj = rinject(ij)
        do iz = 1, inj_nz-1
          if(iz==1) then
            dz0 = dzmin*deltaz(ij)
            dz1 = dzmin*deltaz(ij)
            zinj_grid = dz0
          else
            dz0 = dz1
            zinj_grid = zinj_grid+dz0
            dz1 = deltaz(ij)*(1.-(real(inj_nz-iz-1)/real(inj_nz-1))*(1.-dzmin**onethird))**3
     &          - zinj_grid
          end if
          ab(1,iz) = -1./dz0
          ab(2,iz) = 1./dz0+1./dz1
          do irhs = 1, nrhs
            ix = xi(irhs)
            iy = yi(irhs)
            if (l_inj_use_rho_with_mr) then
              rhs(iz,irhs) = inj_q(ix,iy,iz,ij)/eps0
            else
              rhs(iz,irhs) = 0.
            endif
          end do
          if(iz==1) then
            do irhs = 1, nrhs
              ix = xi(irhs)
              iy = yi(irhs)
              inj_phi_3d(ix,iy,iz-1,ij) = vinject(ij)
              rhs(iz,irhs) = rhs(iz,irhs) + inj_phi_3d(ix,iy,iz-1,ij)/dz0
            end do
          end if
          if(iz==inj_nz-1) then
            do irhs = 1, nrhs
              ix = xi(irhs)
              iy = yi(irhs)
              inj_phi_3d(ix,iy,iz+1,ij) = vinject(ij)-inj_phi(ix,iy,ij)
              rhs(iz,irhs) = rhs(iz,irhs) + inj_phi_3d(ix,iy,iz+1,ij)/dz1
            end do
          end if
        end do
        kd = 1
        ldab = 1+kd
        ldb = inj_nz-1
        n = inj_nz-1
#ifdef CYGWIN
        call dpbsv_('u',n,kd,nrhs,ab,ldab,rhs,ldb,info)
#else
        call dpbsv('u',n,kd,nrhs,ab,ldab,rhs,ldb,info)
#endif
        if(info/=0) then
          write(0,*) 'ERROR in subroutine getinj_phi_mr at exit of call to LAPACK routine dpbsv'
          write(0,*) 'Info /= 0, info = ',info
          write(0,*) 'Stop'
          stop
        end if
        do iz = 1, inj_nz-1
          do irhs = 1, nrhs
            ix = xi(irhs)
            iy = yi(irhs)
            inj_phi_3d(ix,iy,iz,ij) = rhs(iz,irhs)
          end do
        end do

        call inj_transform(nrhs,xxp,yyp,zzp,1,ij,-1,1)
        do iz = 1, inj_nz-1
          if(iz==1) then
            dz0 = dzmin*deltaz(ij)
            dz1 = dzmin*deltaz(ij)
            zinj_grid = dz0
          else
            dz0 = dz1
            zinj_grid = zinj_grid+dz0
            dz1 = deltaz(ij)*(1.-(real(inj_nz-iz-1)/real(inj_nz-1))*(1.-dzmin**onethird))**3
     &          - zinj_grid
          end if
          do irhs = 1, nrhs
            ix = xi(irhs)
            iy = yi(irhs)

c           --- computes components of electric field in source curved coordinates
            ezt = (inj_phi_3d(ix,iy,iz-1,ij)-inj_phi_3d(ix,iy,iz+1,ij))/(dz0+dz1)
            if(linj_eperp .and. ix>0 .and. ix<inj_nx) then
              ext = 0.5*(inj_phi_3d(ix-1,iy,iz,ij)-inj_phi_3d(ix+1,iy,iz,ij))*dxi
            else
              ext = 0.
            end if
            if(linj_eperp .and. iy>0 .and. iy<inj_ny) then
              eyt = 0.5*(inj_phi_3d(ix,iy-1,iz,ij)-inj_phi_3d(ix,iy+1,iz,ij))*dyi
            else
              eyt = 0.
            end if

c           --- angle of point in transverse plane
            aa = atan2(yyp(irhs),dvnz(xxp(irhs)))

c           --- angle relative to z axis
            rr2 = xxp(irhs)**2 + yyp(irhs)**2
            az = asin(sqrt(rr2)*rinji)
c         az = atan2(sqrt(rr2),(rinj - zz(ip))*sign(1.,rinj)) XXX

            atx = asin(xxp(irhs)/sqrt(rinj**2 - yyp(irhs)**2))
            aty = asin(yyp(irhs)/sqrt(rinj**2 - xxp(irhs)**2))

c           --- transform field vector from source coordinates to lab coordinates
            inj_ex_3d(ix,iy,iz,ij) = - ezt*sin(az)*cos(aa) + ext*cos(atx)
            inj_ey_3d(ix,iy,iz,ij) = - ezt*sin(az)*sin(aa) + eyt*cos(aty)
            inj_ez_3d(ix,iy,iz,ij) = ezt*cos(az) + ext*sin(atx) + eyt*sin(aty)
          end do
          inj_ez_3d(:,:,0,ij) = inj_ez_3d(:,:,1,ij)
        end do
        call inj_transform(nrhs,xxp,yyp,zzp,1,ij,1,1)
        call inj_transform(nrhs,inj_ex_3d,inj_ey_3d,inj_ez_3d,1,ij,1,0)
c        inj_phi_3d = inj_phi_tmp
c           --- Obtain the self-field from the electrostatic potential
        call fetche3dfrompositions(0,1,nn,xxp,yyp,zzp,ex_tmp,ey_tmp,ez_tmp)

        dd = inj_d(ij)*dzmin
        if(dd>dzmin) then
          iz  = 1+int(real(inj_nz-1)*(1.-(1.-dd**onethird)/(1.-dzmin**onethird)))
          wz1 = (dd-(1.-((real(inj_nz-iz))/real(inj_nz-1))*(1.-dzmin**onethird))**3)/dz_local(iz)
c          iz = int(real(inj_nz) - log(1./dd)/log(2.))
c          wz1 = dd*2.**(inj_nz-iz)-1.
        else
          iz = 0
          wz1 = dd/dzmin
        end if
        wz0 = 1.-wz1
        do irhs = 1, nrhs
          ix = xi(irhs)
          iy = yi(irhs)
          inj_ex_3d(ix,iy,inj_nz,ij) = ex_tmp(irhs)
          inj_ey_3d(ix,iy,inj_nz,ij) = ey_tmp(irhs)
          inj_ez_3d(ix,iy,inj_nz,ij) = ez_tmp(irhs)
          inj_phi(ix,iy,ij) = vinject(ij) - (inj_phi_3d(ix,iy,iz,ij)*wz0+inj_phi_3d(ix,iy,iz+1,ij)*wz1)
        end do
      end do

      deallocate(rhs,ab)
      deallocate(ex_tmp,ey_tmp,ez_tmp)
      deallocate(dz_local)

!$OMP MASTER
      if (lw3dtimesubs) timeinj_setrho3d = timeinj_setrho3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine inj_sete_3darray(npart,xp,yp,zp,pid,dz,ex,ey,ez)
      use Subtimers3d
      use InGen3d
      use InPart
      use InjectVars
      use InjectVars3d
      integer(ISZ):: npart
      real(kind=8):: xp(npart),yp(npart),zp(npart),pid(npart)
      real(kind=8):: ex(npart),ey(npart),ez(npart)
      real(kind=8):: dz

c Calculate the electric field from the grid for particles near the emitting
c surface.  The normal electric field is calculated from the potential drop
c across a length of 'dz' along a line normal to the emitting surface.  The
c field components are then obtained from the normal field. Optionally,
c the tangential fields near the emitting surface can also be included.
c Eventually this routine could be expanded to calculate the E field near
c any conductor surface.

      real(kind=8),allocatable:: xx(:),yy(:),zz(:)
      real(kind=8),allocatable:: tex(:),tey(:),tez(:)
      integer(ISZ),allocatable:: ijp(:)
      integer(ISZ):: ip,ix,iy,iz,ij
      real(kind=8):: rinj,rinji,rr2,wr,dd,ca,dzmin,deltaz(ninject)
      real(kind=8):: aa,az,atx,aty,px,py,pz,wx,wy,wz,en
      real(kind=8):: dxi,dyi
      real(kind=8):: sx,sy,s0,s1,s2,s3,s4,s5,s6,s7
      real(kind=8):: fourthirds,onethird
      real(kind=8):: substarttime,wtime
      integer(ISZ):: spreadx,spready
      real(kind=8):: dz0, dz1, zinj_grid
      real(kind=8), dimension(0:inj_nz) :: dz_local
      if (lw3dtimesubs) substarttime = wtime()

      allocate(xx(npart),yy(npart),zz(npart))
      allocate(tex(npart),tey(npart),tez(npart))
      allocate(ijp(npart))

      if(solvergeom==XYZgeom) then
        spreadx = 1
        spready = 1
      elseif(solvergeom==RZgeom) then
        spreadx = 1
        if(l_inj_rz) then
          spready = 0
        else
          spready = 1
        end if
      elseif(solvergeom==XZgeom) then
        spreadx = 1
        spready = 0
      elseif(solvergeom==Zgeom) then
        spreadx = 0
        spready = 0
      end if
      tex = 0.
      tey = 0.
      tez = 0.

      dxi = 1./inj_dx
      dyi = 1./inj_dy
      fourthirds = 4./3.
      onethird = 1./3.

      dzmin = 1./(real(inj_nz-1)*2.**onethird+real(2-inj_nz))**3
c      dzmin = 1./(2.**(inj_nz-1))
      deltaz(:) = inj_d(1)*inj_dz0

      dz_local(0) = dzmin
      do iz = 1, inj_nz-1
        if(iz==1) then
          dz0 = dzmin
          dz1 = dzmin
          zinj_grid = dz0
        else
          dz0 = dz1
c         dz1 = 2.*dz1
          zinj_grid = zinj_grid+dz0
          dz1 = (1.-(real(inj_nz-iz-1)/real(inj_nz-1))*(1.-dzmin**onethird))**3
     &        - zinj_grid
        end if
        dz_local(iz) = dz1
      end do
      dz_local(inj_nz) = dz_local(inj_nz-1)

c     --- Transform the particles into the frame of the injection source
      xx = xp
      yy = yp
      zz = zp
      ijp = int(pid)
      call inj_transform(npart,xx,yy,zz,npart,ijp,-1,1)

c     --- Loop over particles
      do ip=1,npart

c       --- Get number of injection source of particle.  The particle id
c       --- has the number of the injection source the particle was
c       --- emitted from. The int is needed since it is doing double
c       --- duty and has other information stored in the fractional part.
        ij = ijp(ip)
        if (ij <= 0 .or. ij > ninject) cycle

c       --- set temporaries
        rinj = rinject(ij)
        rinji = 1./rinj

c       --- Calculate distance of particle from the emitting surface
        rr2 = xx(ip)**2 + yy(ip)**2
        if (abs(zz(ip)) < abs(rinj)) then
c         --- The expression below for zz<rinj is mathematically identical
c         --- to the commented out expression below.  That form is
c         --- used so that for large values of rinj, i.e. a flat emitting
c         --- surface, the correct value of the distance is calculated,
c         --- namely zz.
c         dd = rinj - sqrt(rr2 + (rinj - zz(ip))**2)
          dd = (2.*zz(ip) - (zz(ip)**2 + rr2)*rinji)/
     &         (1. + sqrt(rr2*rinji**2 + (1. - zz(ip)*rinji)**2))
        else
c         --- When zz>rinj, a seperate equation is needed.
          dd = abs(rinj) + sqrt(rr2 + (zz(ip) - rinj)**2)
        endif
        dd = dd/deltaz(ij)

c       --- Only calculate E-field if particle close to emitting surface.
c        if (0.0 <= inj_d(ij)*dd .and. abs(dd) <= 0.5) then
        if (0.0 <= inj_d(ij)*dd .and. abs(dd) <= 1.) then
c         --- The check if a particle is within the radial area of the
c         --- source is probably not needed. It has been removed since there
c         --- are some cases which this check is incorrect - for example
c         --- in a source with a diverging beam. The particles spread out
c         --- beyond the edge of the source and would not pass this check.
c    &    ((xx(ip)*binject(ij))**2 + (yy(ip)*ainject(ij))**2 <=
c    &    (ainject(ij)*binject(ij))**2) .and.
c    &    ((xx(ip)*binjmin(ij))**2 + (yy(ip)*ainjmin(ij))**2 >=
c    &    (ainjmin(ij)*binjmin(ij))**2) ) then

c         --- angle of point in transverse plane
          aa = atan2(yy(ip),dvnz(xx(ip)))

c         --- angle relative to z axis
          az = asin(sqrt(rr2)*rinji)
c         az = atan2(sqrt(rr2),(rinj - zz(ip))*sign(1.,rinj)) XXX

c         --- Find coordinates of the point on the secondary surface in front
c         --- of the emitting surface along a line perpendicular to the
c         --- emitting surface.
          wr = 1./(1. - deltaz(ij)*dd*rinji)
          if(l_inj_rz) then
            px = abs(rinj*sin(az) - inj_xmmin(ij))*dxi
            px = abs(xx(ip)*wr - inj_xmmin(ij))*dxi
            py = 0.
          else
            px = abs(rinj*sin(az)*cos(aa) - inj_xmmin(ij))*dxi
            py = abs(rinj*sin(az)*sin(aa) - inj_ymmin(ij))*dyi
            px = abs(xx(ip)*wr - inj_xmmin(ij))*dxi
            py = abs(yy(ip)*wr - inj_ymmin(ij))*dyi
          end if

          if (spreadx*px < 0. .or. (px+1)*spreadx > inj_nx .or.
     &        spready*py < 0. .or. (py+1)*spready > inj_ny) cycle

          ix = spreadx*int(px)
          iy = spready*int(py)
          wx = spreadx*(px - ix)
          wy = spready*(py - iy)
          if(dd>dzmin) then
            iz = 1+int(real(inj_nz-1)*(1.-(1.-dd**onethird)/(1.-dzmin**onethird)))
            wz = (dd-(1.-((real(inj_nz)-real(iz))/real(inj_nz-1))*(1.-dzmin**onethird))**3)/dz_local(iz)
c            iz = int(real(inj_nz) - log(1./dd)/log(2.))
c            wz = dd*2.**(inj_nz-iz)-1.
          else
            iz = 0
            wz = dd/dzmin
          end if


          ex(ip) = 0.
          ey(ip) = 0.
          ez(ip) = 0.

          s0 = (1.-wx)*(1.-wy)*(1.-wz)
          s1 =     wx *(1.-wy)*(1.-wz)
          s2 = (1.-wx)*    wy *(1.-wz)
          s3 =     wx *    wy *(1.-wz)
          s4 = (1.-wx)*(1.-wy)*    wz
          s5 =     wx *(1.-wy)*    wz
          s6 = (1.-wx)*    wy *    wz
          s7 =     wx *    wy *    wz 

          if(l_inj_rz) then
            en = inj_ex_3d(ix        ,iy        ,iz  ,ij)*s0 +
     &           inj_ex_3d(ix+spreadx,iy        ,iz  ,ij)*s1 +
     &           inj_ex_3d(ix        ,iy+spready,iz  ,ij)*s2 +
     &           inj_ex_3d(ix+spreadx,iy+spready,iz  ,ij)*s3 +
     &           inj_ex_3d(ix        ,iy        ,iz+1,ij)*s4 +
     &           inj_ex_3d(ix+spreadx,iy        ,iz+1,ij)*s5 +
     &           inj_ex_3d(ix        ,iy+spready,iz+1,ij)*s6 +
     &           inj_ex_3d(ix+spreadx,iy+spready,iz+1,ij)*s7
            tex(ip) = en*cos(aa)
            tey(ip) = en*sin(aa)
          else
            sx = 1.
            sy = 1.
            if (xx(ip) < inj_xmmin(ij)) sx = -1.
            if (yy(ip) < inj_ymmin(ij)) sy = -1.

            tex(ip) = sx*(inj_ex_3d(ix        ,iy        ,iz  ,ij)*s0 +
     &                    inj_ex_3d(ix+spreadx,iy        ,iz  ,ij)*s1 +
     &                    inj_ex_3d(ix        ,iy+spready,iz  ,ij)*s2 +
     &                    inj_ex_3d(ix+spreadx,iy+spready,iz  ,ij)*s3 +
     &                    inj_ex_3d(ix        ,iy        ,iz+1,ij)*s4 +
     &                    inj_ex_3d(ix+spreadx,iy        ,iz+1,ij)*s5 +
     &                    inj_ex_3d(ix        ,iy+spready,iz+1,ij)*s6 +
     &                    inj_ex_3d(ix+spreadx,iy+spready,iz+1,ij)*s7 )

            tey(ip) = sy*(inj_ey_3d(ix        ,iy        ,iz  ,ij)*s0 +
     &                    inj_ey_3d(ix+spreadx,iy        ,iz  ,ij)*s1 +
     &                    inj_ey_3d(ix        ,iy+spready,iz  ,ij)*s2 +
     &                    inj_ey_3d(ix+spreadx,iy+spready,iz  ,ij)*s3 +
     &                    inj_ey_3d(ix        ,iy        ,iz+1,ij)*s4 +
     &                    inj_ey_3d(ix+spreadx,iy        ,iz+1,ij)*s5 +
     &                    inj_ey_3d(ix        ,iy+spready,iz+1,ij)*s6 +
     &                    inj_ey_3d(ix+spreadx,iy+spready,iz+1,ij)*s7 )
          end if

          tez(ip) =     inj_ez_3d(ix        ,iy        ,iz  ,ij)*s0 +
     &                  inj_ez_3d(ix+spreadx,iy        ,iz  ,ij)*s1 +
     &                  inj_ez_3d(ix        ,iy+spready,iz  ,ij)*s2 +
     &                  inj_ez_3d(ix+spreadx,iy+spready,iz  ,ij)*s3 +
     &                  inj_ez_3d(ix        ,iy        ,iz+1,ij)*s4 +
     &                  inj_ez_3d(ix+spreadx,iy        ,iz+1,ij)*s5 +
     &                  inj_ez_3d(ix        ,iy+spready,iz+1,ij)*s6 +
     &                  inj_ez_3d(ix+spreadx,iy+spready,iz+1,ij)*s7

        endif
      enddo

      call inj_transform(npart,tex,tey,tez,npart,ijp,1,0)
      ex = ex + tex
      ey = ey + tey
      ez = ez + tez

      deallocate(xx,yy,zz)
      deallocate(tex,tey,tez)
      deallocate(ijp)

!$OMP MASTER
      if (lw3dtimesubs) timeinj_sete3d = timeinj_sete3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c===========================================================================
      subroutine gettinj_phi()
      use Subtimers3d
      use Constant
      use InMesh3d
      use InjectVars
      use InjectVars3d

c Calculate potential drop from transverse emitting surface at distance
c of dx from the surface.

      integer(ISZ):: ij,iz,ith
      real(kind=8):: aa,rsign,aoverb
      real(kind=8):: xx((1+nttinjmax)*(1+nztmax)*ntinj)
      real(kind=8):: yy((1+nttinjmax)*(1+nztmax)*ntinj)
      real(kind=8):: zz((1+nttinjmax)*(1+nztmax)*ntinj)
      integer(ISZ):: ti((1+nttinjmax)*(1+nztmax)*ntinj)
      integer(ISZ):: zi((1+nttinjmax)*(1+nztmax)*ntinj)
      integer(ISZ):: in((1+nttinjmax)*(1+nztmax)*ntinj)
      real(kind=8):: pp((1+nttinjmax)*(1+nztmax)*ntinj)
      integer(ISZ):: nn,ii,i1,i2
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (ntinj == 0) return

      if (inject < 1 .or. inject > 3) return

      call updateinj_zmminmax()

c     --- loop over injection sources
c     --- i1 and i2 keep track of the data as new data from each injection
c     --- source is added. i1 is the starting location for the each set of
c     --- data and i2 is the location of the last data point.
      i1 = 1
      do ij=1,ntinj
        i2 = i1 - 1

        do iz=0,nztmax
          do ith=0,nttinj(ij)

            i2 = i2 + 1

            aa = ith*2.*pi/nttinj(ij)
            if (ltinj_outward(ij)) then
              rsign = +1.
            else
              rsign = -1.
            endif
            aoverb = atinject(iz,ij)/btinject(iz,ij)
            xx(i2) = xtinject(ij) + (atinject(iz,ij) + rsign*inj_dx)*cos(aa)
            yy(i2) = ytinject(ij) + (btinject(iz,ij) + rsign*inj_dx*aoverb)*sin(aa)
            zz(i2) = ztinjmn(ij) + iz*dztinj(ij)

c           --- Save coordinates relative to injection grid.
            zi(i2) = iz
            ti(i2) = ith
            in(i2) = ij

          enddo
        enddo

c       --- Select only those data points which are within the grid.
c       --- This also removes particles not within the domain of this
c       --- process.
        nn = i1
        do ii=i1,i2
          if (inj_zmmin <= zz(ii) .and. zz(ii) <= inj_zmmax) then
            if (nn < ii) then
              zi(nn) = zi(ii)
              ti(nn) = ti(ii)
              xx(nn) = xx(ii)
              yy(nn) = yy(ii)
              zz(nn) = zz(ii)
              in(nn) = in(ii)
            endif
            nn = nn + 1
          endif
        enddo
        i2 = nn - 1

        i1 = i2 + 1

      enddo
      nn = i1 - 1

      pp = 0.
      call fetchphi(nn,xx,yy,zz,pp)

c     --- Calculate the difference between phi at that point and phi on
c     --- the emitting surface.
      tinj_phi = 0.
      do ii=1,nn
        tinj_phi(ti(ii),zi(ii),in(ii)) = vtinject(in(ii)) - pp(ii)
      enddo
        
#ifdef MPIPARALLEL
c     --- The result calculated on each processor is gathered on all
      call parallelnonzerorealarray(inj_phi,(1+inj_nx)*(1+inj_ny)*ninject)
#endif


!$OMP MASTER
      if (lw3dtimesubs) timegettinj_phi = timegettinj_phi + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine tinj_sete3d(npart,xp,yp,zp,pid,ex,ey,ez)
      use Subtimers3d
      use Constant
      use InGen3d
      use InMesh3d
      use Picglb3d
      use InPart
      use InjectVars
      use InjectVars3d
      integer(ISZ):: npart
      real(kind=8):: xp(npart),yp(npart),zp(npart),pid(npart)
      real(kind=8):: ex(npart),ey(npart),ez(npart)

c Calculate the electric field from the grid for particles near the transverse
c emitting surfaces.  The normal electric field is calculated from the potential drop
c across a length of 'dx' along a line normal to the emitting surface.  The
c field components are then obtained from the normal field.
c Eventually this routine could be expanded to calculate the E field near
c any conductor surface.

      real(kind=8):: xx,yy,zz
      integer(ISZ):: ijp
      integer(ISZ):: ip,iz,ith,ij
      real(kind=8):: ainj,binj,ainjp1,binjp1,aoverb
      real(kind=8):: p1x,p1y,p2x,p2y
      real(kind=8):: dr,rr,dd,rsign,dzi,rinj
      real(kind=8):: aa,wt,wz,en
      real(kind=8):: sx,sy
      real(kind=8):: fourthirds,onethird
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- The user can force the code to skip the special calculation of the 
c     --- E fields for particles behind the virtual surface.
      if (linj_efromgrid) return

      fourthirds = 4./3.
      onethird = 1./3.

c     --- Loop over particles
      do ip=1,npart

c       --- Get number of injection source of particle.  The particle id
c       --- has the number of the injection source the particle was
c       --- emitted from. The int is needed since it is doing double
c       --- duty and has other information stored in the fractional part.
        ij = int(pid(ip)) - ninject
        if (ij <= 0 .or. ij > ntinj) cycle

        xx = xp(ip) - xtinject(ij)
        yy = yp(ip) - ytinject(ij)
        zz = zp(ip) - ztinjmn(ij)

        iz = zz/dztinj(ij)
        if (iz < 0. .or. iz >= nztinj(ij)) cycle
        wz = zz/dztinj(ij) - iz
        
        ainj = atinject(iz,ij)
        binj = btinject(iz,ij)
        ainjp1 = atinject(iz+1,ij)
        binjp1 = btinject(iz+1,ij)
        aoverb = ainj/binj

        if(.not.l_inj_rz) then
          aa = atan2(yy*aoverb,xx)
        else
          aa = 0.
        endif
        if (aa < 0.) aa = aa + 2*pi
        ith = aa/(2.*pi)*nttinj(ij)

        if (ltinj_outward(ij)) then
          rsign = +1.
        else
          rsign = -1.
        endif

c       --- Calculate distance of particle from the emitting surface
        rr = sqrt(xx**2 + (yy*aoverb)**2)
        p1x = (ainj*(1.-wz) + ainjp1*wz)*cos(aa)
        p1y = (binj*(1.-wz) + binjp1*wz)*sin(aa)
        p2x = p1x + rsign*inj_dx*cos(aa)
        p2y = p1y + rsign*inj_dx*aoverb*sin(aa)
        dr = sqrt((p2x - p1x)**2 + (p2y - p1y)**2)
        dd = (rr - rinj)/dr

c       --- Only calculate E-field if particle close to emitting surface.
        if (0.0 <= dd .and. abs(dd) <= 1.) then

c         --- Find coordinates of the point on the secondary surface in front
c         --- of the emitting surface along a line perpendicular to the
c         --- emitting surface.
          wt = aa/(2.*pi)*nttinj(ij) - ith

c         --- Calculate the normal electric field from the potential drop in
c         --- front of the emitting surface. The normal field given by the
c         --- Child-Langmuir solution (for planar surfaces) is used.
          en = (tinj_phi(ith  ,iz  ,ij)*(1.-wt)*(1.-wz) +
     &          tinj_phi(ith+1,iz  ,ij)*    wt *(1.-wz) +
     &          tinj_phi(ith  ,iz+1,ij)*(1.-wt)*    wz  +
     &          tinj_phi(ith+1,iz+1,ij)*    wt *    wz   )/dr
          if (linj_enormcl) then
c           --- Scale the normal E field to match the Child-Langmuir solution
            en = en*fourthirds*dd**onethird
          endif

c         --- Set the particle's electric field based off of the normal
c         --- electric fields.
          if(l_inj_rz) aa = atan2(yy*aoverb,dvnz(xx))
          ex(ip) = -en*cos(aa)
          ey(ip) = -en*sin(aa)
          if (linj_eperp) then
            ez(ip) = dd*((tinj_phi(ith  ,iz  ,ij)*(1.-wt)*(1.-wz) +
     &                    tinj_phi(ith+1,iz  ,ij)*    wt *(1.-wz) +
     &                    tinj_phi(ith  ,iz+1,ij)*(1.-wt)*    wz  +
     &                    tinj_phi(ith+1,iz+1,ij)*    wt *    wz) -
     &                   (tinj_phi(ith  ,iz+1,ij)*(1.-wt)*(1.-wz) +
     &                    tinj_phi(ith+1,iz+1,ij)*    wt *(1.-wz) +
     &                    tinj_phi(ith  ,iz+2,ij)*(1.-wt)*    wz  +
     &                    tinj_phi(ith+1,iz+2,ij)*    wt *    wz))/dztinj(ij)
          else
            ez(ip) =  0.
          endif

        endif
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timetinj_sete3d = timetinj_sete3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      real(kind=8) function ellipseperimeter(a,b)
      use Constant
      real(kind=8):: a,b
c Use several terms of Gauss-Kummer series to get good accuracy.
      real(kind=8):: h,p

      h = (a-b)**2/(a+b)**2
      p = pi*(a+b)*(1 + h/4 + h**2/64 + h**3/256 + 25*h**4/16384 + 49*h**5/65536)
      ellipseperimeter = p

      return
      end
c=============================================================================
