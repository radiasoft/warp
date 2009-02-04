#include "top.h"
module em2d_depos
use EM2D_FIELDtypemodule
contains
subroutine depose_jxjy_esirkepov_linear_serial(j,np,xp,yp,xpold,ypold,uzp,gaminv,w,q,xmin,ymin,dt,dx,dy,nx,ny,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny
   real(kind=8), dimension(-1:nx+2,-1:ny+1,3), intent(in out) :: j
   real(kind=8), dimension(np) :: xp,yp,xpold,ypold,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,xmin,ymin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dtsdx,dtsdy,sd(18),xint,yint,wx(1:4,1:5),wy(1:5,1:4)
   real(kind=8) :: xold,yold,xmid,ymid,x,y,wq,wqx,wqy,tmp,vx,vy,vz,dts2dx,dts2dy,s1x,s2x,s1y,s2y,invsurf,invdtdx,invdtdy
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
      
      dsx = 0.
      dsy = 0.
      
      do ip=1,np
      
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        
        vx = (xp(ip)-xpold(ip))/dt
        vy = (yp(ip)-ypold(ip))/dt
        vz = gaminv(ip)*uzp(ip)
        
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

        iixp0=floor(xold)
        ijxp0=floor(yold)

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

        iixp=floor(x)
        ijxp=floor(y)
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
!        write(0,*) dix,idx,diy,idy,dsx(2+idx),dsy(2+idy)
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

        wq = wq*vz*invsurf
      
        iixp=floor(x)
        ijxp=floor(y)

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

subroutine depose_jxjy_esirkepov_n_serial(cj,np,xp,yp,xpold,ypold,uzp,gaminv,w,q,xmin,ymin, &
                                                 dt,dx,dy,nx,ny,nox,noy,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nox,noy
   real(kind=8), dimension(-1:nx+1,-1:ny+1,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,xpold,ypold,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,xmin,ymin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dtsdx,dtsdy,xint,yint
   real(kind=8),dimension(-nox:nox+1,-noy:noy+1) :: wx,wy,wz,sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, oxint, oyint, onesixth, twothird
   real(kind=8), DIMENSION(-nox:nox+1) :: sx, sy, sz, sx0, sy0, sz0, dsx, dsy, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k

      sx0=0.;sy0=0.;sz0=0.
      sdz=0.
      onesixth = 1./6.
      twothird = 2./3.
      
      dxi = 1./dx
      dyi = 1./dy
      dtsdx = dt*dxi
      dtsdy = dt*dyi
      dts2dx = 0.5*dtsdx
      dts2dy = 0.5*dtsdy
      invvol = 1./(dx*dy)
      invdtdx = 1./(dt*dy)
      invdtdy = 1./(dt*dx)
      invdtdz = 1./(dt*dx*dy)

      do ip=1,np
      
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        
        xold = (xpold(ip)-xmin)*dxi
        yold = (ypold(ip)-ymin)*dyi
        
        vx = (xp(ip)-xpold(ip))/dt
        vy = (yp(ip)-ypold(ip))/dt
        vz = uzp(ip)*gaminv(ip)

        if (l_particles_weight) then
          wq=q*w(ip)
        else
          wq=q*w(1)
        end if
        wqx = wq*invdtdx
        wqy = wq*invdtdy
        wqz = wq*invdtdz

!       computation of current at x(n+1/2),v(n+1/2)

        iixp0=floor(x)
        ijxp0=floor(y)

        xint=x-iixp0
        yint=y-ijxp0

        if (nox==1) then
          sx0(0) = 1.-xint
          sx0(1) = xint
        else if (nox==2) then
          oxint = 1.-xint
          sx0(-1) = onesixth*oxint**3
          sx0(0)  = twothird-xint**2*(1.-0.5*xint)
          sx0(1)  = twothird-oxint**2*(1.-0.5*oxint)
          sx0(2)  = onesixth*xint**3
        end if
        
        if (noy==1) then
          sy0(0) = 1.-yint
          sy0(1) = yint
        else if (noy==2) then
          oyint = 1.-yint
          sy0(-1) = onesixth*oyint**3
          sy0(0)  = twothird-yint**2*(1.-0.5*yint)
          sy0(1)  = twothird-oyint**2*(1.-0.5*oyint)
          sy0(2)  = onesixth*yint**3
        end if

        iixp=floor(xold)
        ijxp=floor(yold)
        xint = xold-iixp
        yint = yold-ijxp

        dix = iixp-iixp0
        diy = ijxp-ijxp0

        sx=0.;sy=0.;sz=0.

        if (nox==1) then
          sx(0+dix) = 1.-xint
          sx(1+dix) = xint
        else if (nox==2) then
          oxint = 1.-xint
          sx(-1+dix) = onesixth*oxint**3
          sx(0+dix)  = twothird-xint**2*(1.-0.5*xint)
          sx(1+dix)  = twothird-oxint**2*(1.-0.5*oxint)
          sx(2+dix)  = onesixth*xint**3
        end if
        
        if (noy==1) then
          sy(0+diy) = 1.-yint
          sy(1+diy) = yint
        else if (noy==2) then
          oyint = 1.-yint
          sy(-1+diy) = onesixth*oyint**3
          sy(0+diy)  = twothird-yint**2*(1.-0.5*yint)
          sy(1+diy)  = twothird-oyint**2*(1.-0.5*oyint)
          sy(2+diy)  = onesixth*yint**3
        end if

        dsx = sx - sx0
        dsy = sy - sy0

        do k=-nox, nox+1
          do j=-nox, noy+1
              wx(i,j) = wqx*dsx(i)*(sy0(j)+0.5*dsy(j))
              wy(i,j) = wqy*dsy(j)*(sx0(i)+0.5*dsx(i))
          end do
        end do

        do i = -nox, nox
          sdx(i,:)  = wx(i,:)
          if (i>-1) sdx(i,:)=sdx(i,:)+sdx(i-1,:)
          cj(iixp0+i,ijxp0-1:ijxp0+2,1) = cj(iixp0+i,ijxp0-1:ijxp0+2,1)+sdx(i,:)
        end do        

        do j = -noy, noy
          sdy(:,j)  = wy(:,j)
          if (j>-1) sdy(:,j)=sdy(:,j)+sdy(:,j-1)
          cj(iixp0-1:iixp0+2,ijxp0+j,2) = cj(iixp0-1:iixp0+2,ijxp0+j,2)+sdy(:,j)
        end do        

     
        ! Esirkepov deposition of Jx and Jy is over; now starts linear deposition of Jz
        xmid=x-dts2dx*vx
        ymid=y-dts2dy*vy

        ! there is a shift of 0.5 since Ez/Jz are aligned with Bz, at the center of the cell
        x = x-0.5
        y = y-0.5

        wq = wq*vz*invvol
      
        iixp=floor(x)
        ijxp=floor(y)

        xint = x-iixp
        yint = y-ijxp


        if (nox==1) then
          sx0(0) = 1.-xint
          sx0(1) = xint
        else if (nox==2) then
          oxint = 1.-xint
          sx0(-1) = onesixth*oxint**3
          sx0(0)  = twothird-xint**2*(1.-0.5*xint)
          sx0(1)  = twothird-oxint**2*(1.-0.5*oxint)
          sx0(2)  = onesixth*xint**3
        end if
        
        if (noy==1) then
          sy0(0) = 1.-yint
          sy0(1) = yint
        else if (noy==2) then
          oyint = 1.-yint
          sy0(-1) = onesixth*oyint**3
          sy0(0)  = twothird-yint**2*(1.-0.5*yint)
          sy0(1)  = twothird-oyint**2*(1.-0.5*oyint)
          sy0(2)  = onesixth*yint**3
        end if

        do j = -noy+1,noy
          do i = -nox+1,nox
            cj(iixp+i  ,ijxp+j  ,3)=cj(iixp+i  ,ijxp+j  ,3)+sx0(j)*sy0(k)*wq
          end do
        end do
        
    end do

  return
end subroutine depose_jxjy_esirkepov_n_serial



subroutine depose_rho_esirkepov_linear_serial(rho,np,xp,yp,w,q,xmin,ymin,dx,dy,nx,ny,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny
   real(kind=8), dimension(-1:nx+2,-1:ny+1), intent(in out) :: rho
   real(kind=8), dimension(np) :: xp,yp,w
   real(kind=8) :: q,dt,dx,dy,xmin,ymin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,xint,yint
   real(kind=8) :: x,y,wq,invsurf,s1x,s2x,s1y,s2y
   integer(ISZ) :: iixp,ijxp,ip,dix,diy

      dxi = 1./dx
      dyi = 1./dy
      invsurf = dxi*dyi

      do ip=1,np
      
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        
        if (l_particles_weight) then
          wq=q*w(ip)*invsurf
        else
          wq=q*invsurf
        end if
      
        iixp=floor(x)
        ijxp=floor(y)

        xint = x-iixp
        yint = y-ijxp

        s1x = 1.-xint
        s2x = xint

        s1y = 1.-yint
        s2y = yint

        rho(iixp  ,ijxp  )=rho(iixp  ,ijxp  )+s1x*s1y*wq
        rho(iixp+1,ijxp  )=rho(iixp+1,ijxp  )+s2x*s1y*wq
        rho(iixp  ,ijxp+1)=rho(iixp  ,ijxp+1)+s1x*s2y*wq
        rho(iixp+1,ijxp+1)=rho(iixp+1,ijxp+1)+s2x*s2y*wq
      
    end do

  return
end subroutine depose_rho_esirkepov_linear_serial
subroutine geteb2d_linear_serial(np,xp,yp,ex,ey,ez,bx,by,bz,xmin,ymin,dx,dy,nx,ny,exg,eyg,ezg,bxg,byg,bzg)
   
      integer(ISZ) :: np,nx,ny
      real(kind=8), dimension(np) :: xp,yp,ex,ey,ez,bx,by,bz
      real(kind=8), dimension(-1:nx+2,-1:ny+1) :: exg,eyg,ezg,bxg,byg,bzg 
      real(kind=8) :: xmin,ymin,dx,dy
      integer(ISZ) :: ip, iixp, ijxp
      real(kind=8) :: dxi, dyi, x, y, xint, yint, s1x, s2x, s1y, s2y, w1, w2, w3, w4

      dxi = 1./dx
      dyi = 1./dy

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi

        iixp=floor(x)
        ijxp=floor(y)

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
 
 subroutine getf2d_linear_serial(np,xp,yp,fx,fy,fz,xmin,ymin,dx,dy,nx,ny,fxg,fyg,fzg)
   
      integer(ISZ) :: np,nx,ny
      real(kind=8), dimension(np) :: xp,yp,fx,fy,fz
      real(kind=8), dimension(-1:nx+2,-1:ny+1) :: fxg,fyg,fzg
      real(kind=8) :: xmin,ymin,dx,dy
      integer(ISZ) :: ip, iixp, ijxp
      real(kind=8) :: dxi, dyi, x, y, xint, yint, s1x, s2x, s1y, s2y, w1, w2, w3, w4

      dxi = 1./dx
      dyi = 1./dy

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi

        iixp=floor(x)
        ijxp=floor(y)

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
          
        fx(ip) = fx(ip)+(w1*fxg(iixp+1,ijxp+1)+w2*fxg(iixp,ijxp+1)+w3*fxg(iixp+1,ijxp)+w4*fxg(iixp,ijxp))
        fy(ip) = fy(ip)+(w1*fyg(iixp+1,ijxp+1)+w2*fyg(iixp,ijxp+1)+w3*fyg(iixp+1,ijxp)+w4*fyg(iixp,ijxp))
        fz(ip) = fz(ip)+(w1*fzg(iixp+1,ijxp+1)+w2*fzg(iixp,ijxp+1)+w3*fzg(iixp+1,ijxp)+w4*fzg(iixp,ijxp))

     end do

   return
 end subroutine getf2d_linear_serial
end module em2d_depos

subroutine depose_current_em2d(np,xp,yp,uxp,uyp,uzp,gaminv,w,q,dt,l_particles_weight,field,fpatchfine)
   use EM2D_FIELDtypemodule
   use EM2D_FIELDobjects
   use em2d_depos
   implicit none
   integer(ISZ) :: np
   real(kind=8), dimension(np) :: xp,yp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: dt,q
   logical(ISZ) :: l_particles_weight
   TYPE(EM2D_FIELDtype), target :: field,fpatchfine

   integer(ISZ) :: ip, np_inpatch
   logical(ISZ) :: l_inpatch(np)
   
   TYPE(EM2D_FIELDtype), POINTER :: f,ff

   f => field
   ff => fpatchfine

   if (l_onegrid) then
     call depose_jxjy_esirkepov_linear_serial(f%J,np,xp,yp,uxp,uyp,uzp,gaminv, &
                                              w,q,f%xmin,f%ymin,dt,f%dx,f%dy, &
                                              f%nx,f%ny,l_particles_weight)
   else
     np_inpatch = 0
     do ip=1,np
       IF(xp(ip)>ff%xminpatch_scatter .and. xp(ip)<ff%xmaxpatch_scatter .and. &
          yp(ip)>ff%yminpatch_scatter .and. yp(ip)<ff%ymaxpatch_scatter) then
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
                                                pack(w,.not. l_inpatch), &
                                                q,f%xmin,f%ymin,dt, &
                                                f%dx,f%dy,f%nx,f%ny,l_particles_weight)
     end if
     if (np_inpatch>0) then
       call depose_jxjy_esirkepov_linear_serial(ff%J,np_inpatch, &
                                                pack(xp,l_inpatch), &
                                                pack(yp,l_inpatch), &
                                                pack(uxp,l_inpatch), &
                                                pack(uyp,l_inpatch), &
                                                pack(uzp,l_inpatch), &
                                                pack(gaminv,l_inpatch), &
                                                pack(w,l_inpatch), &
                                                q,ff%xmin,ff%ymin,dt, &
                                                ff%dx,ff%dy,ff%nx,ff%ny,l_particles_weight)
     end if     
   endif
end subroutine depose_current_em2d

subroutine geteb_em2d(np,xp,yp,ex,ey,ez,bx,by,bz,field,fpatchfine)
   use EM2D_FIELDtypemodule
   use EM2D_FIELDobjects
   use em2d_depos
   implicit none
   
   integer(ISZ) :: np
   real(kind=8), dimension(np) :: xp,yp,ex,ey,ez,bx,by,bz
   TYPE(EM2D_FIELDtype), target :: field,fpatchfine

   integer(ISZ) :: ip, ipt, np_inpatch, np_outpatch
   logical(ISZ) :: l_inpatch(np)
   real(kind=8), allocatable, dimension(:) :: ext,eyt,ezt,bxt,byt,bzt
   real(kind=8) :: d1,d2,d3,d4,d,wtz(np)

   TYPE(EM2D_FIELDtype), POINTER :: f, ff

   f => field
   ff => fpatchfine
   
   if (l_onegrid) then
     call geteb2d_linear_serial(np,xp,yp,ex,ey,ez,bx,by,bz,f%xmin,f%ymin,f%dx,f%dy,f%nx,f%ny, &
          f%ex,f%ey,f%ez,f%bx,f%by,f%bz)
   else
     np_inpatch = 0
     do ip=1,np
       IF(xp(ip)>ff%xminpatch_gather .and. xp(ip)<ff%xmaxpatch_gather .and. &
          yp(ip)>ff%yminpatch_gather .and. yp(ip)<ff%ymaxpatch_gather) then
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
                                     f%xmin,f%ymin,f%dx,f%dy,f%nx,f%ny, &
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
       allocate(ext(np_inpatch),eyt(np_inpatch),ezt(np_inpatch), &
                bxt(np_inpatch),byt(np_inpatch),bzt(np_inpatch))
       ext=0.; eyt=0.; ezt=0.; 
       bxt=0.; byt=0.; bzt=0.; 
       call geteb2d_linear_serial(np,pack(xp,l_inpatch), &
                                     pack(yp,l_inpatch), &
                                     ext,eyt,ezt,bxt,byt,bzt, &
                                     ff%xmin,ff%ymin,ff%dx,ff%dy,ff%nx,ff%ny, &
                                     ff%exfsum,ff%eyfsum,ff%ezfsum,ff%bxfsum,ff%byfsum,ff%bzfsum)
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

subroutine getf_em2d(np,xp,yp,fx,fy,fz,field,fpatchfine,which)
   use EM2D_FIELDobjects
   use em2d_depos
   implicit none
   
   integer(ISZ) :: np,which
   real(kind=8), dimension(np) :: xp,yp,fx,fy,fz
   TYPE(EM2D_FIELDtype), target :: field,fpatchfine

   integer(ISZ) :: ip, ipt, np_inpatch, np_outpatch
   logical(ISZ) :: l_inpatch(np)
   real(kind=8), allocatable, dimension(:) :: fxt,fyt,fzt
   real(kind=8), pointer, dimension(:,:) :: fxg,fyg,fzg
   real(kind=8) :: d1,d2,d3,d4,d,wtz(np)

   TYPE(EM2D_FIELDtype), POINTER :: f, ff
   
   f => field
   ff => fpatchfine
   
   if(which==1) then
     fxg => f%Ex
     fyg => f%Ey
     fzg => f%Ez
   else
     fxg => f%Bx
     fyg => f%By
     fzg => f%Bz
   end if
   if (l_onegrid) then
     call getf2d_linear_serial(np,xp,yp,fx,fy,fz,f%xmin,f%ymin,f%dx,f%dy,f%nx,f%ny, &
          fxg,fyg,fzg)
   else
     np_inpatch = 0
     do ip=1,np
       IF(xp(ip)>ff%xminpatch_gather .and. xp(ip)<ff%xmaxpatch_gather .and. &
          yp(ip)>ff%yminpatch_gather .and. yp(ip)<ff%ymaxpatch_gather) then
         l_inpatch(ip) = .true.
         np_inpatch = np_inpatch+1
       else
         l_inpatch(ip) = .false.
       END if
     end do
     np_outpatch = np-np_inpatch
     if (np_outpatch>0) then
       allocate(fxt(np_outpatch),fyt(np_outpatch),fzt(np_outpatch))
       fxt=0.; fyt=0.; fzt=0.; 
       call getf2d_linear_serial(np_outpatch,pack(xp,.not. l_inpatch), &
                                     pack(yp,.not. l_inpatch), &
                                     fxt,fyt,fzt, &
                                     f%xmin,f%ymin,f%dx,f%dy,f%nx,f%ny, &
                                     fxg,fyg,fzg)
       ipt = 1
       do ip=1,np
         if (.not. l_inpatch(ip)) then
           fx(ip)=fx(ip)+fxt(ipt)
           fy(ip)=fy(ip)+fyt(ipt)
           fz(ip)=fz(ip)+fzt(ipt)
           ipt=ipt+1
         end if
       enddo
       deallocate(fxt,fyt,fzt)
     end if
     if (np_inpatch>0) then
       if(which==1) then
         fxg => ff%Exfsum
         fyg => ff%Eyfsum
         fzg => ff%Ezfsum
       else
         fxg => ff%Bxfsum
         fyg => ff%Byfsum
         fzg => ff%Bzfsum
       end if
       allocate(fxt(np_inpatch),fyt(np_inpatch),fzt(np_inpatch))
       fxt=0.; fyt=0.; fzt=0.; 
       call getf2d_linear_serial(np_inpatch,pack(xp,l_inpatch), &
                                     pack(yp,l_inpatch), &
                                     fxt,fyt,fzt, &
                                     ff%xmin,ff%ymin,ff%dx,ff%dy,ff%nx,ff%ny, &
                                     fxg,fyg,fzg)
       ipt = 1
       do ip=1,np
         if (l_inpatch(ip)) then
           fx(ip)=fx(ip)+fxt(ipt)
           fy(ip)=fy(ip)+fyt(ipt)
           fz(ip)=fz(ip)+fzt(ipt)
           ipt=ipt+1
         end if
       enddo
       deallocate(fxt,fyt,fzt)
     end if
   endif
end subroutine getf_em2d

subroutine em2d_geteb2d_linear_serial(np,xp,yp,ex,ey,ez,bx,by,bz,xmin,ymin,dx,dy,nx,ny,exg,eyg,ezg,bxg,byg,bzg)
  use em2d_depos,Only: geteb2d_linear_serial
   
  integer(ISZ) :: np,nx,ny
  real(kind=8), dimension(np) :: xp,yp,ex,ey,ez,bx,by,bz
  real(kind=8), dimension(-1:nx+2,-1:ny+1) :: exg,eyg,ezg,bxg,byg,bzg 
  real(kind=8) :: xmin,ymin,dx,dy

  call geteb2d_linear_serial(np,xp,yp,ex,ey,ez,bx,by,bz,xmin,ymin,dx,dy,nx,ny,exg,eyg,ezg,bxg,byg,bzg)

  return
end subroutine em2d_geteb2d_linear_serial

subroutine em2d_getf2d_linear_serial(np,xp,yp,fx,fy,fz,xmin,ymin,dx,dy,nx,ny,fxg,fyg,fzg)
  use em2d_depos,Only: getf2d_linear_serial
   
  integer(ISZ) :: np,nx,ny
  real(kind=8), dimension(np) :: xp,yp,fx,fy,fz
  real(kind=8), dimension(-1:nx+2,-1:ny+1) :: fxg,fyg,fzg
  real(kind=8) :: xmin,ymin,dx,dy

  call getf2d_linear_serial(np,xp,yp,fx,fy,fz,xmin,ymin,dx,dy,nx,ny,fxg,fyg,fzg)

  return
end subroutine em2d_getf2d_linear_serial

subroutine em2d_depose_jxjy_esirkepov_linear_serial(j,np,xp,yp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,dt,dx,dy,nx,ny,l_particles_weight)
  use em2d_depos,Only: depose_jxjy_esirkepov_linear_serial
  integer(ISZ) :: np,nx,ny
  real(kind=8), dimension(-1:nx+2,-1:ny+1,3), intent(in out) :: j
  real(kind=8), dimension(np) :: xp,yp,uxp,uyp,uzp,gaminv,w
  real(kind=8) :: q,dt,dx,dy,xmin,ymin
  logical(ISZ) :: l_particles_weight

  call depose_jxjy_esirkepov_linear_serial(j,np,xp,yp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,dt,dx,dy,nx,ny,l_particles_weight)

  return
end subroutine em2d_depose_jxjy_esirkepov_linear_serial

subroutine em2d_depose_rho_esirkepov_linear_serial(rho,np,xp,yp,w,q,xmin,ymin,dx,dy,nx,ny,l_particles_weight)
  use em2d_depos,Only: depose_rho_esirkepov_linear_serial
  integer(ISZ) :: np,nx,ny
  real(kind=8), dimension(-1:nx+2,-1:ny+1), intent(in out) :: rho
  real(kind=8), dimension(np) :: xp,yp,w
  real(kind=8) :: q,dx,dy,xmin,ymin
  logical(ISZ) :: l_particles_weight

  call depose_rho_esirkepov_linear_serial(rho,np,xp,yp,w,q,xmin,ymin,dx,dy,nx,ny,l_particles_weight)

  return
end subroutine em2d_depose_rho_esirkepov_linear_serial

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

subroutine smooth2d_121(q,nx,ny)
 implicit none

 integer(ISZ) :: nx,ny,i,j

 real(kind=8), dimension(0:nx+3,0:ny+2) :: q
 real(kind=8) :: temp0, temp1

!     x smoothing

      do  j=0,ny+2
        temp0 = q(0,j)
        do  i=1,nx+2
          temp1 = q(i,j)
          q(i,j) = 0.5*q(i,j)+0.25*(temp0+q(i+1,j))
          temp0 = temp1
        end do
      end do

!     y smoothing

      do  i=0,nx+3
        temp0 = q(i,0)
        do  j=1,ny+1
          temp1 = q(i,j)
          q(i,j) = 0.5*q(i,j)+0.25*(temp0+q(i,j+1))
          temp0 = temp1
        end do
      end do

      return
      end subroutine smooth2d_121

!subroutine em2d_smoothdensity()
!   use EM2D_FIELDobjects
!   implicit none
!   
!   TYPE(EM2D_FIELDtype), POINTER :: f
!   real(kind=8), dimension(:,:), pointer :: jaux
!   integer(ISZ) :: i,ngrids
!
!   if (l_onegrid) then
!     ngrids = 1
!   else
!     ngrids = 2
!   end if
!   do i=1, ngrids
!     if(i==1) then
!       f => field
!     else
!       f => fpatchfine
!     endif
!    jaux => f%J(:,:,1)
!    call smooth2d_lindman(jaux,f%nx,f%ny)
!    jaux => f%J(:,:,2)
!    call smooth2d_lindman(jaux,f%nx,f%ny)
!    jaux => f%J(:,:,3)
!    call smooth2d_lindman(jaux,f%nx,f%ny)
!  end do
  
!  return
!end subroutine em2d_smoothdensity

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
      real(kind=8),pointer:: ex(:),ey(:),ez(:),bx(:),by(:),bz(:)
      real(kind=8),pointer:: gaminv(:),pid(:,:)
      real(kind=8),pointer:: sm(:),sq(:),sw(:),dtscale(:)
      integer(ISZ),pointer:: ins(:),nps(:)
      real(kind=8) :: wtmp(nparpgrp)

      integer(ISZ) :: is, ipmin, ip

!     --- Create local pointers to the arrays in pgroup.
      xp => pgroup%xp
      yp => pgroup%yp
      zp => pgroup%zp
      uxp => pgroup%uxp
      uyp => pgroup%uyp
      uzp => pgroup%uzp
      gaminv => pgroup%gaminv
      ex => pgroup%ex
      ey => pgroup%ey
      ez => pgroup%ez
      bx => pgroup%bx
      by => pgroup%by
      bz => pgroup%bz
      if (pgroup%npid > 0) pid => pgroup%pid

      sm => pgroup%sm
      sq => pgroup%sq
      sw => pgroup%sw
      ins => pgroup%ins
      nps => pgroup%nps
      dtscale => pgroup%dtscale

      wtmp = 0.

!      field%J = 0.      

      ! put fields back on staggered grid
      do is=1,pgroup%ns
         do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)

            ex(ipmin:ipmin+ip-1) = 0.
            ey(ipmin:ipmin+ip-1) = 0.
            ez(ipmin:ipmin+ip-1) = 0.
            bx(ipmin:ipmin+ip-1) = 0.
            by(ipmin:ipmin+ip-1) = 0.
            bz(ipmin:ipmin+ip-1) = 0.

!           call geteb_em2d(ip,xp(ipmin),yp(ipmin),ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin))

            call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin), &
                          bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is), 0.5*dt, ibpush)
            call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin), &
                          ex(ipmin), ey(ipmin), ez(ipmin), sq(is), sm(is), 0.5*dt)
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
!            call depose_current_em2d(ip,xp(ipmin),yp(ipmin), &
!                                     uxp(ipmin),uyp(ipmin),uzp(ipmin), &
!                                     gaminv(ipmin),wtmp,sq(is)*sw(is),dt, &
!                                     .false.)
                               
         end do
      end do
!      if(l_smoothdensity) call em2d_smoothdensity()

!      call grimax(field) 
      
!      call push_em_b(field,0.5*dt)
!      call push_em_e(field,dt)
!      call push_em_b(field,0.5*dt)
!      call move_window_fields()

!     put fields values at nodes
!      call griuni(field) 

      do is=1,pgroup%ns
         do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)

            ex=0.; ey=0.; ez=0.; bx=0.; by=0.; bz=0.

!            call geteb_em2d(ip,xp(ipmin),yp(ipmin),ex,ey,ez,bx,by,bz)

            call epush3d (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin), &
                          ex(ipmin), ey(ipmin), ez(ipmin), sq(is), sm(is), 0.5*dt)
!              --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin), &
                          gamadv,lrelativ)
            call bpush3d (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin), &
                          bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is), 0.5*dt, ibpush)
         end do
      end do

      it=it+1
      time=time+dt
      
  return
end subroutine em2d_step

