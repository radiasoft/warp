#include "top.h"
module em2d_depos
contains
subroutine depose_jxjy_esirkepov_linear_serial(j,np,xp,yp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,dt,dx,dy,l_particles_weight)
   implicit none
   real(kind=8), dimension(-1:,-1:,1:), intent(in out) :: j
   integer(ISZ) :: np
   real(kind=8), dimension(np) :: xp,yp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,xmin,ymin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dtsdx,dtsdy,sd(18),xint,yint,wx(1:4,1:5),wy(1:5,1:4)
   real(kind=8) :: xold,yold,xmid,ymid,x,y,wq,wqx,wqy,tmp,vx,vy,dts2dx,dts2dy,s1x,s2x,s1y,s2y,invsurf,invdtdx,invdtdy
   real(kind=8), DIMENSION(6) :: sx, sy, sx0, sy0, dsx, dsy
   integer(ISZ) :: iixp0,ijxp0,iixp,ijxp,ip,dix,diy,idx,idy

      dxi = 1./dx
      dyi = 1./dy
      dtsdx = dt*dxi
      dtsdy = dt*dyi
      dts2dx = 0.5*dtsdx
      dts2dy = 0.5*dtsdy
      invsurf = 1./(dx*dy)
      invdtdx = 1./(dt*dx)
      invdtdy = 1./(dt*dy)

      do ip=1,np
      
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        
        xold=x-dtsdx*vx
        yold=y-dtsdy*vy

        if (l_particles_weight) then
          wq=q*w(ip)
        else
          wq=q
        end if
        wqx = wq*invdtdy
        wqy = wq*invdtdx

!       computation of current at x n+1/2 v n+1/2

        iixp0=int(xold)
        ijxp0=int(yold)

        xint=xold-iixp0
        yint=yold-ijxp0

        sx(1)=0.;sx(2)=0.;sx(5)=0.;sx(6)=0.
        sy(1)=0.;sy(2)=0.;sy(5)=0.;sy(6)=0.

!        sx0(1,iv) = 0.
        sx0(2) = 0.
        sx0(3) = 1.-xint
        sx0(4) = xint
        sx0(5) = 0.
!        sx0(6,iv) = 0.

!        sy0(1,iv) = 0.
        sy0(2) = 0.
        sy0(3) = 1.-yint
        sy0(4) = yint
        sy0(5) = 0.
!        sy0(6,iv) = 0.

        iixp=int(x)
        ijxp=int(y)
        xint = x-iixp
        yint = y-ijxp

        dix = iixp-iixp0
        diy = ijxp-ijxp0

        sx(2+dix) = 0.
        sx(3+dix) = 1.-xint
        sx(4+dix) = xint
        sx(5+dix) = 0.       

        sy(2+diy) = 0.
        sy(3+diy) = 1.-yint
        sy(4+diy) = yint
        sy(5+diy) = 0.       
       

        idx = MIN(0,dix)
        idy = MIN(0,diy)

        iixp0 = iixp0+idx
        ijxp0 = ijxp0+idy

        dsx(2)=sx(2)-sx0(2);dsx(3)=sx(3)-sx0(3)
        dsx(4)=sx(4)-sx0(4);dsx(5)=sx(5)-sx0(5)
        dsy(2)=sy(2)-sy0(2);dsy(3)=sy(3)-sy0(3)
        dsy(4)=sy(4)-sy0(4);dsy(5)=sy(5)-sy0(5)

        tmp = (sy0(3+idy)+0.5*dsy(3+idy))*wqx; wx(2,1) = dsx(3+idx)*tmp
                                               wx(3,1) = dsx(4+idx)*tmp
                                               wx(4,1) = dsx(5+idx)*tmp
        tmp = (sy0(4+idy)+0.5*dsy(4+idy))*wqx; wx(2,2) = dsx(3+idx)*tmp
                                               wx(3,2) = dsx(4+idx)*tmp
                                               wx(4,2) = dsx(5+idx)*tmp
        tmp = (sy0(5+idy)+0.5*dsy(5+idy))*wqx; wx(2,4) = dsx(3+idx)*tmp
                                               wx(3,4) = dsx(4+idx)*tmp
                                               wx(4,4) = dsx(5+idx)*tmp

        tmp = (sx0(3+idx)+0.5*dsx(3+idx))*wqy; wy(2,1) = dsy(3+idy)*tmp
                                               wy(2,2) = dsy(4+idy)*tmp
                                               wy(2,4) = dsy(5+idy)*tmp
        tmp = (sx0(4+idx)+0.5*dsx(4+idx))*wqy; wy(3,1) = dsy(3+idy)*tmp
                                               wy(3,2) = dsy(4+idy)*tmp
                                               wy(3,4) = dsy(5+idy)*tmp
        tmp = (sx0(5+idx)+0.5*dsx(5+idx))*wqy; wy(4,1) = dsy(3+idy)*tmp
                                               wy(4,2) = dsy(4+idy)*tmp
                                               wy(4,4) = dsy(5+idy)*tmp

        sd(1) = wx(2,1)
        sd(2) = wx(3,1)+sd(1)
        sd(3) = wx(4,1)+sd(2)

        sd(4) = wx(2,2)
        sd(5) = wx(3,2)+sd(4)
        sd(6) = wx(4,2)+sd(5)

        sd(7) = wx(2,4)
        sd(8) = wx(3,4)+sd(7)
        sd(9) = wx(4,4)+sd(8)

        sd(10) = wy(2,1)
        sd(13) = wy(2,2)+sd(10)
        sd(16) = wy(2,4)+sd(13)

        sd(11) = wy(3,1)
        sd(14) = wy(3,2)+sd(11)
        sd(17) = wy(3,4)+sd(14)

        sd(12) = wy(4,1)
        sd(15) = wy(4,2)+sd(12)
        sd(18) = wy(4,4)+sd(15)

        j(iixp0,  ijxp0  ,1)=j(iixp0  ,ijxp0  ,1)-sd(1)
        j(iixp0+1,ijxp0  ,1)=j(iixp0+1,ijxp0  ,1)-sd(2)
        j(iixp0+2,ijxp0  ,1)=j(iixp0+2,ijxp0  ,1)-sd(3)
        j(iixp0,  ijxp0+1,1)=j(iixp0  ,ijxp0+1,1)-sd(4)
        j(iixp0+1,ijxp0+1,1)=j(iixp0+1,ijxp0+1,1)-sd(5)
        j(iixp0+2,ijxp0+1,1)=j(iixp0+2,ijxp0+1,1)-sd(6)
        j(iixp0,  ijxp0+2,1)=j(iixp0  ,ijxp0+2,1)-sd(7)
        j(iixp0+1,ijxp0+2,1)=j(iixp0+1,ijxp0+2,1)-sd(8)
        j(iixp0+2,ijxp0+2,1)=j(iixp0+2,ijxp0+2,1)-sd(9)
            
        j(iixp0,  ijxp0  ,2)=j(iixp0  ,ijxp0  ,2)-sd(10)
        j(iixp0+1,ijxp0  ,2)=j(iixp0+1,ijxp0  ,2)-sd(11)
        j(iixp0+2,ijxp0  ,2)=j(iixp0+2,ijxp0  ,2)-sd(12)
        j(iixp0,  ijxp0+1,2)=j(iixp0  ,ijxp0+1,2)-sd(13)
        j(iixp0+1,ijxp0+1,2)=j(iixp0+1,ijxp0+1,2)-sd(14)
        j(iixp0+2,ijxp0+1,2)=j(iixp0+2,ijxp0+1,2)-sd(15)
        j(iixp0,  ijxp0+2,2)=j(iixp0  ,ijxp0+2,2)-sd(16)
        j(iixp0+1,ijxp0+2,2)=j(iixp0+1,ijxp0+2,2)-sd(17)
        j(iixp0+2,ijxp0+2,2)=j(iixp0+2,ijxp0+2,2)-sd(18)
      
        ! Esirkepov deposition of Jx and Jy is over; now starts linear deposition of Jz
        xmid=x-dts2dx*vx
        ymid=y-dts2dy*vy

        ! there is a shift of 0.5 since Ez/Jz are aligned with Bz, at the center of the cell
        x = x-0.5
        y = y-0.5

        wq = wq*uzp(ip)*gaminv(ip)*invsurf
      
        iixp=int(x)
        ijxp=int(y)

        xint = x-iixp
        yint = y-ijxp

        s1x = 1.-xint
        s2x = xint

        s1y = 1.-yint
        s2y = yint

        j(iixp  ,ijxp  ,3)=j(iixp  ,ijxp  ,3)+s1x*s1y*wq
        j(iixp+1,ijxp  ,3)=j(iixp+1,ijxp  ,3)+s2x*s1y*wq
        j(iixp  ,ijxp+1,3)=j(iixp  ,ijxp+1,3)+s1x*s2y*wq
        j(iixp+1,ijxp+1,3)=j(iixp+1,ijxp+1,3)+s2x*s2y*wq
        
      
    end do

  return
end subroutine depose_jxjy_esirkepov_linear_serial

 subroutine geteb2d_linear_serial(np,xp,yp,ex,ey,ez,bx,by,bz,xmin,ymin,dx,dy,exg,eyg,ezg,bxg,byg,bzg)
   
      integer(ISZ) :: np
      real(kind=8), dimension(np) :: xp,yp,ex,ey,ez,bx,by,bz
      real(kind=8), dimension(-1:,-1:) :: exg,eyg,ezg,bxg,byg,bzg 
      real(kind=8) :: xmin,ymin,dx,dy
      integer(ISZ) :: ip, iixp, ijxp
      real(kind=8) :: dxi, dyi, x, y, xint, yint, s1x, s2x, s1y, s2y, w1, w2, w3, w4

      dxi = 1./dx
      dyi = 1./dy

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi

        iixp=int(x)
        ijxp=int(y)

        xint=x-iixp
        yint=y-ijxp

        s1x=xint
        s2x=1.-xint
        s1y=yint
        s2y=1.-yint

        w1 = s1x*s1y
        w2 = s2x*s1y 
        w3 = s1x*s2y
        w4 = s2x*s2y 
          
        ex(ip) = ex(ip)+(w1*exg(iixp+1,ijxp+1)+w2*exg(iixp,ijxp+1)+w3*exg(iixp+1,ijxp)+w4*exg(iixp,ijxp))
        ey(ip) = ey(ip)+(w1*eyg(iixp+1,ijxp+1)+w2*eyg(iixp,ijxp+1)+w3*eyg(iixp+1,ijxp)+w4*eyg(iixp,ijxp))
        ez(ip) = ez(ip)+(w1*ezg(iixp+1,ijxp+1)+w2*ezg(iixp,ijxp+1)+w3*ezg(iixp+1,ijxp)+w4*ezg(iixp,ijxp))

        bx(ip) = bx(ip)+(w1*bxg(iixp+1,ijxp+1)+w2*bxg(iixp,ijxp+1)+w3*bxg(iixp+1,ijxp)+w4*bxg(iixp,ijxp))
        by(ip) = by(ip)+(w1*byg(iixp+1,ijxp+1)+w2*byg(iixp,ijxp+1)+w3*byg(iixp+1,ijxp)+w4*byg(iixp,ijxp))
        bz(ip) = bz(ip)+(w1*bzg(iixp+1,ijxp+1)+w2*bzg(iixp,ijxp+1)+w3*bzg(iixp+1,ijxp)+w4*bzg(iixp,ijxp))

     end do

   return
 end subroutine geteb2d_linear_serial
end module em2d_depos

subroutine depose_current_em2d(np,xp,yp,uxp,uyp,uzp,gaminv,w,q,dt,l_particles_weight)
   use EM2D_FIELDobjects
   use em2d_depos
   implicit none
   integer(ISZ) :: np
   real(kind=8), dimension(np) :: xp,yp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: dt,q
   logical(ISZ) :: l_particles_weight

   integer(ISZ) :: ip, np_inpatch
   logical(ISZ) :: l_inpatch(np)
   
   TYPE(EM2D_FIELDtype), POINTER :: f

   f => field
   if (l_onegrid) then
     call depose_jxjy_esirkepov_linear_serial(f%J,np,xp,yp,uxp,uyp,uzp,gaminv, &
                                              w,q,f%xmin,f%ymin,dt,f%dx,f%dy,l_particles_weight)
   else
     np_inpatch = 0
     do ip=1,np
       IF(xp(ip)>xminpatch_scatter .and. xp(ip)<xmaxpatch_scatter .and. &
          yp(ip)>yminpatch_scatter .and. yp(ip)<ymaxpatch_scatter) then
         l_inpatch(ip) = .true.
         np_inpatch = np_inpatch+1
       else
         l_inpatch(ip) = .false.
       END if
     end do
     if (np-np_inpatch>0) then
       call depose_jxjy_esirkepov_linear_serial(f%J,np-np_inpatch, &
                                                pack(xp,.not. l_inpatch), &
                                                pack(yp,.not. l_inpatch), &
                                                pack(uxp,.not. l_inpatch), &
                                                pack(uyp,.not. l_inpatch), &
                                                pack(uzp,.not. l_inpatch), &
                                                pack(gaminv,.not. l_inpatch), &
                                                w,q,f%xmin,f%ymin,dt,f%dx,f%dy,l_particles_weight)
     end if
     if (np_inpatch>0) then
       f => fpatchfine
       call depose_jxjy_esirkepov_linear_serial(f%J,np_inpatch, &
                                                pack(xp,l_inpatch), &
                                                pack(yp,l_inpatch), &
                                                pack(uxp,l_inpatch), &
                                                pack(uyp,l_inpatch), &
                                                pack(uzp,l_inpatch), &
                                                pack(gaminv,l_inpatch), &
                                                w,q,f%xmin,f%ymin,dt,f%dx,f%dy,l_particles_weight)
     end if     
   endif
end subroutine depose_current_em2d

subroutine geteb_em2d(np,xp,yp,ex,ey,ez,bx,by,bz)
   use EM2D_FIELDobjects
   use em2d_depos
   implicit none
   
   integer(ISZ) :: np
   real(kind=8), dimension(np) :: xp,yp,ex,ey,ez,bx,by,bz

   integer(ISZ) :: ip, ipt, np_inpatch, np_outpatch
   logical(ISZ) :: l_inpatch(np)
   real(kind=8), allocatable, dimension(:) :: ext,eyt,ezt,bxt,byt,bzt
   real(kind=8) :: d1,d2,d3,d4,d,wtz(np)

   TYPE(EM2D_FIELDtype), POINTER :: f

   f => field
   if (l_onegrid) then
     call geteb2d_linear_serial(np,xp,yp,ex,ey,ez,bx,by,bz,f%xmin,f%ymin,f%dx,f%dy, &
          f%ex,f%ey,f%ez,f%bx,f%by,f%bz)
   else
     np_inpatch = 0
     do ip=1,np
       IF(xp(ip)>xminpatch_gather .and. xp(ip)<xmaxpatch_gather .and. &
          yp(ip)>yminpatch_gather .and. yp(ip)<ymaxpatch_gather) then
         l_inpatch(ip) = .true.
         np_inpatch = np_inpatch+1
       else
         l_inpatch(ip) = .false.
       END if
     end do
     np_outpatch = np-np_inpatch
     if (np_outpatch>0) then
       allocate(ext(np_outpatch),eyt(np_outpatch),ezt(np_outpatch), &
                bxt(np_outpatch),byt(np_outpatch),bzt(np_outpatch))
       ext=0.; eyt=0.; ezt=0.; 
       bxt=0.; byt=0.; bzt=0.; 
       call geteb2d_linear_serial(np,pack(xp,.not. l_inpatch), &
                                     pack(yp,.not. l_inpatch), &
                                     ext,eyt,ezt,bxt,byt,bzt, &
                                     f%xmin,f%ymin,f%dx,f%dy, &
                                     f%ex,f%ey,f%ez,f%bx,f%by,f%bz)
       ipt = 1
       do ip=1,np
         if (.not. l_inpatch(ip)) then
           ex(ip)=ex(ip)+ext(ipt)
           ey(ip)=ey(ip)+eyt(ipt)
           ez(ip)=ez(ip)+ezt(ipt)
           bx(ip)=bx(ip)+bxt(ipt)
           by(ip)=by(ip)+byt(ipt)
           bz(ip)=bz(ip)+bzt(ipt)
         end if
         ipt=ipt+1
       enddo
       deallocate(ext,eyt,ezt,bxt,byt,bzt)
     end if
     if (np_inpatch>0) then
       f => fpatchfine
       allocate(ext(np_inpatch),eyt(np_inpatch),ezt(np_inpatch), &
                bxt(np_inpatch),byt(np_inpatch),bzt(np_inpatch))
       ext=0.; eyt=0.; ezt=0.; 
       bxt=0.; byt=0.; bzt=0.; 
       call geteb2d_linear_serial(np,pack(xp,l_inpatch), &
                                     pack(yp,l_inpatch), &
                                     ext,eyt,ezt,bxt,byt,bzt, &
                                     f%xmin,f%ymin,f%dx,f%dy, &
                                     exfsum,eyfsum,ezfsum,bxfsum,byfsum,bzfsum)
       ipt = 1
       do ip=1,np
         if (l_inpatch(ip)) then
           ex(ip)=ex(ip)+ext(ipt)
           ey(ip)=ey(ip)+eyt(ipt)
           ez(ip)=ez(ip)+ezt(ipt)
           bx(ip)=bx(ip)+bxt(ipt)
           by(ip)=by(ip)+byt(ipt)
           bz(ip)=bz(ip)+bzt(ipt)
         end if
         ipt=ipt+1
       enddo
       deallocate(ext,eyt,ezt,bxt,byt,bzt)
     end if
   endif
end subroutine geteb_em2d

subroutine smooth2d_lindman(q,nx,ny)
 implicit none

 integer(ISZ) :: nx,ny,ns,i1,i2,j1,j2,is,i,j,ntemp

 real(kind=8), dimension(0:nx+3,0:ny+2) :: q
 real(kind=8), dimension(5) :: cs,ds,dc
 real(kind=8), dimension(:), ALLOCATABLE :: temp

 data cs /4*.25,-1.25/,ds/4*.5,3.5/,ns/5/
 data dc /4*2.,-2.8/

     ntemp = 2*max(nx,ny)+4
     ALLOCATE(temp(0:ntemp))


      i1=0
      i2=nx+2
      j1=0
      j2=ny+2

      temp=0.

!     x smoothing

      do 110 is=1,ns
      do  i=2,nx-1,2
!cdir nodep
      do  j=1,ny+1
      temp(j+j1)=q(i-1,j)+dc(is)*q(i,j)+q(i+1,j)
      q(i-1,j)=cs(is)*temp(j+j2)
      temp(j+j2)=q(i,j)+dc(is)*q(i+1,j)+q(i+2,j)
      q(i,j)=cs(is)*temp(j+j1)

      enddo
      enddo

      do  j=1,ny+1
      q(nx,j)=cs(is)*temp(j+j2)
      q(1,j)=0.
      q(nx+1,j)=0.
      enddo

 110  continue

!     y smoothing
!     -----------
      do 160 is=1,ns

      do j=2,ny-1,2

!cdir nodep
      do  i=1,nx+1
      temp(i+i1)=q(i,j-1)+dc(is)*q(i,j)+q(i,j+1)
      q(i,j-1)=cs(is)*temp(i+i2)
      temp(i+i2)=q(i,j)+dc(is)*q(i,j+1)+q(i,j+2)
      q(i,j)=cs(is)*temp(i+i1)
      enddo
      enddo

      do  i=1,nx+1
      q(i,ny)=cs(is)*temp(i+i2)
      q(i,1)=0.
      q(i,ny+1)=0.
      enddo

 160  continue
      DEALLOCATE(temp)

      return
      end subroutine smooth2d_lindman

subroutine em2d_smoothdensity()
   use EM2D_FIELDobjects
   implicit none
   
   TYPE(EM2D_FIELDtype), POINTER :: f
   real(kind=8), dimension(:,:), pointer :: jaux
   integer(ISZ) :: i,ngrids

   if (l_onegrid) then
     ngrids = 1
   else
     ngrids = 2
   end if
   do i=1, ngrids
     if(i==1) then
       f => field
     else
       f => fpatchfine
     endif
    jaux => f%J(:,:,1)
    call smooth2d_lindman(jaux,f%nx,f%ny)
    jaux => f%J(:,:,2)
    call smooth2d_lindman(jaux,f%nx,f%ny)
    jaux => f%J(:,:,3)
    call smooth2d_lindman(jaux,f%nx,f%ny)
  end do
  
  return
end subroutine em2d_smoothdensity

subroutine em2d_step()
      use InGen
      use Constant
      use GlobalVars
      use Picglb
      use Particles
      use Beam_acc
      use DKInterptmp
      use EM2D_FIELDobjects
      use em2d_depos
      implicit None
      
!     --- Create local pointers to the arrays in pgroup.
      real(kind=8),pointer:: xp(:),yp(:),zp(:),uxp(:),uyp(:),uzp(:)
      real(kind=8),pointer:: gaminv(:),pid(:,:)
      real(kind=8),pointer:: sm(:),sq(:),sw(:),dtscale(:)
      integer(ISZ),pointer:: ins(:),nps(:)
      real(kind=8) :: wtmp(nparpgrp)

      integer(ISZ) :: is, ipmin, ip

!  Set storage for field arrays if not adequate already
      if (npfield < nparpgrp) then
         npfield = nparpgrp
         call gchange("DKInterptmp",0)
      endif

!     --- Create local pointers to the arrays in pgroup.
      xp => pgroup%xp
      yp => pgroup%yp
      zp => pgroup%zp
      uxp => pgroup%uxp
      uyp => pgroup%uyp
      uzp => pgroup%uzp
      gaminv => pgroup%gaminv
      if (pgroup%npid > 0) pid => pgroup%pid

      sm => pgroup%sm
      sq => pgroup%sq
      sw => pgroup%sw
      ins => pgroup%ins
      nps => pgroup%nps
      dtscale => pgroup%dtscale

      wtmp = 0.

      field%J = 0.      

      ! put fields back on staggered grid
      do is=1,pgroup%ns
         do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)

            ex=0.; ey=0.; ez=0.; bx=0.; by=0.; bz=0.

            call geteb_em2d(ip,xp(ipmin),yp(ipmin),ex,ey,ez,bx,by,bz)

            call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin), &
                          bx, by, bz, sq(is), sm(is), 0.5*dt, ibpush)
            call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin), &
                          ex, ey, ez, sq(is), sm(is), 0.5*dt)
!              --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin), &
                          gamadv,lrelativ)
            call xpush3d(ip,xp(ipmin),yp(ipmin),zp(ipmin), &
                         uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),dt)
            
         end do
      end do
      
      call particleboundaries3d(pgroup)
      
      do is=1,pgroup%ns
         do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)

            ! we assume that all particles have same weight
            call depose_current_em2d(ip,xp(ipmin),yp(ipmin), &
                                     uxp(ipmin),uyp(ipmin),uzp(ipmin), &
                                     gaminv(ipmin),wtmp,sq(is)*sw(is),dt, &
                                     .false.)
                               
         end do
      end do
      if(l_smoothdensity) call em2d_smoothdensity()

      call grimax(field) 
      
      call push_em_b(field,0.5*dt)
      call push_em_e(field,dt,clight,mu0)
      call push_em_b(field,0.5*dt)

!     put fields values at nodes
      call griuni(field) 

      do is=1,pgroup%ns
         do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)

            ex=0.; ey=0.; ez=0.; bx=0.; by=0.; bz=0.

            call geteb_em2d(ip,xp(ipmin),yp(ipmin),ex,ey,ez,bx,by,bz)

            call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin), &
                          ex, ey, ez, sq(is), sm(is), 0.5*dt)
!              --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin), &
                          gamadv,lrelativ)
            call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin), &
                          bx, by, bz, sq(is), sm(is), 0.5*dt, ibpush)
         end do
      end do

      it=it+1
      time=time+dt
      
  return
end subroutine em2d_step





