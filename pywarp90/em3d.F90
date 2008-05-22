#include "top.h"
subroutine depose_jxjyjz_esirkepov_linear_serial(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz
   real(kind=8), dimension(-1:nx+1,-1:ny+1,-1:nz+1,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(-1:2,-1:2,-1:2) :: wx,wy,wz,sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz
   real(kind=8), DIMENSION(-1:2) :: sx, sy, sz, sx0, sy0, sz0, dsx, dsy, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k

      sx0=0.;sy0=0.;sz0=0.
      sdz=0.
      
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      dtsdx = dt*dxi
      dtsdy = dt*dyi
      dtsdz = dt*dzi
      dts2dx = 0.5*dtsdx
      dts2dy = 0.5*dtsdy
      dts2dz = 0.5*dtsdz
      invvol = 1./(dx*dy*dz)
      invdtdx = 1./(dt*dy*dz)
      invdtdy = 1./(dt*dx*dz)
      invdtdz = 1./(dt*dx*dy)

      do ip=1,np
      
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi
        
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)
        
        xold=x-dtsdx*vx
        yold=y-dtsdy*vy
        zold=z-dtsdz*vz

        if (l_particles_weight) then
          wq=q*w(ip)
        else
          wq=q
        end if
        wqx = wq*invdtdx
        wqy = wq*invdtdy
        wqz = wq*invdtdz

!       computation of current at x(n+1/2),v(n+1/2)

        iixp0=floor(x)
        ijxp0=floor(y)
        ikxp0=floor(z)

        xint=x-iixp0
        yint=y-ijxp0
        zint=z-ikxp0

        sx0(0) = 1.-xint
        sx0(1) = xint
        sy0(0) = 1.-yint
        sy0(1) = yint
        sz0(0) = 1.-zint
        sz0(1) = zint

        iixp=floor(xold)
        ijxp=floor(yold)
        ikxp=floor(zold)
        xint = xold-iixp
        yint = yold-ijxp
        zint = zold-ikxp

!        if (iixp/=iixp0 .and. ijxp/=ijxp0) write(0,*) '**************************** WARNING 1:',iixp,iixp0,ijxp,ijxp0
!        if (iixp/=iixp0 .and. ikxp/=ikxp0) write(0,*) '**************************** WARNING 2:',iixp,iixp0,ikxp,ikxp0
!        if (ikxp/=ikxp0 .and. ijxp/=ijxp0) write(0,*) '**************************** WARNING 3:',ikxp,ikxp0,ijxp,ijxp0
!        if (iixp/=iixp0) write(0,*) '****** WARNING x:',ip,iixp,iixp0
!        if (ijxp/=ijxp0) write(0,*) '****** WARNING y:',ip,ijxp,ijxp0
!        if (ikxp/=ikxp0) write(0,*) '****** WARNING z:',ip,ikxp,ikxp0

        dix = iixp-iixp0
        diy = ijxp-ijxp0
        diz = ikxp-ikxp0

        sx=0.;sy=0.;sz=0.
        sx(0+dix) = 1.-xint
        sx(1+dix) = xint
        sy(0+diy) = 1.-yint
        sy(1+diy) = yint
        sz(0+diz) = 1.-zint
        sz(1+diz) = zint

!        idx = MIN(0,dix)
!        idy = MIN(0,diy)
!        idz = MIN(0,diz)

!        iixp0 = iixp0-idx
!        ijxp0 = ijxp0-idy
!        ikxp0 = ikxp0-idz

        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0

        do k=-1, 2
          do j=-1, 2
            do i=-1, 2
              wx(i,j,k) = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
              wy(i,j,k) = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
              wz(i,j,k) = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
            end do
          end do
        end do

        do i = -1, 1
          sdx(i,:,:)  = wx(i,:,:)
          if (i>-1) sdx(i,:,:)=sdx(i,:,:)+sdx(i-1,:,:)
          cj(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1) = cj(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1)+sdx(i,:,:)
        end do        

        do j = -1, 1
          sdy(:,j,:)  = wy(:,j,:)
          if (j>-1) sdy(:,j,:)=sdy(:,j,:)+sdy(:,j-1,:)
          cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2) = cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2)+sdy(:,j,:)
        end do        

        do k = -1, 1
          sdz(:,:,k)  = wz(:,:,k)
          if (k>-1) sdz(:,:,k)=sdz(:,:,k)+sdz(:,:,k-1)
          cj(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+K,3) = cj(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+k,3)+sdz(:,:,k)
        end do        

    end do

  return
end subroutine depose_jxjyjz_esirkepov_linear_serial

subroutine depose_jxjyjz_esirkepov_linear_serial_old(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                     dt,dx,dy,dz,nx,ny,nz,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz
   real(kind=8), dimension(-1:nx+1,-1:ny+1,-1:nz+1,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(-1:2,-1:2,-1:2) :: wx,wy,wz,sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz
   real(kind=8), DIMENSION(-1:2) :: sx, sy, sz, sx0, sy0, sz0, dsx, dsy, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k

      sx0=0.;sy0=0.;sz0=0.
      
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      dtsdx = dt*dxi
      dtsdy = dt*dyi
      dtsdz = dt*dzi
      dts2dx = 0.5*dtsdx
      dts2dy = 0.5*dtsdy
      dts2dz = 0.5*dtsdz
      invvol = 1./(dx*dy*dz)
      invdtdx = 1./(dt*dy*dz)
      invdtdy = 1./(dt*dx*dz)
      invdtdz = 1./(dt*dx*dy)

      write(0,*) (minval(xp)-xmin)/dx,(minval(yp)-ymin)/dy,(minval(zp)-zmin)/dz
      write(0,*) (minval(xp-uxp*gaminv*dt)-xmin)/dx,(minval(yp-uyp*gaminv*dt)-ymin)/dy,(minval(zp-uzp*gaminv*dt)-zmin)/dz

      do ip=1,np
      
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi
        
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)
        
        xold=x-dtsdx*vx
        yold=y-dtsdy*vy
        zold=z-dtsdz*vz

        if (l_particles_weight) then
          wq=q*w(ip)
        else
          wq=q
        end if
        wqx = -wq*invdtdx
        wqy = -wq*invdtdy
        wqz = -wq*invdtdz

!       computation of current at x(n+1/2),v(n+1/2)

        iixp0=floor(xold)
        ijxp0=floor(yold)
        ikxp0=floor(zold)

        xint=xold-iixp0
        yint=yold-ijxp0
        zint=zold-ikxp0

        sx0(0) = 1.-xint
        sx0(1) = xint
        sy0(0) = 1.-yint
        sy0(1) = yint
        sz0(0) = 1.-zint
        sz0(1) = zint

        iixp=floor(x)
        ijxp=floor(y)
        ikxp=floor(z)
        xint = x-iixp
        yint = y-ijxp
        zint = z-ikxp

        dix = iixp-iixp0
        diy = ijxp-ijxp0
        diz = ikxp-ikxp0

        sx=0.;sy=0.;sz=0.
        sx(0+dix) = 1.-xint
        sx(1+dix) = xint
        sy(0+diy) = 1.-yint
        sy(1+diy) = yint
        sz(0+diz) = 1.-zint
        sz(1+diz) = zint

!        idx = MIN(0,dix)
!        idy = MIN(0,diy)
!        idz = MIN(0,diz)

!        iixp0 = iixp0-idx
!        ijxp0 = ijxp0-idy
!        ikxp0 = ikxp0-idz

        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0

        do k=-1, 2
          do j=-1, 2
            do i=-1, 2
              wx(i,j,k) = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
              wy(i,j,k) = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
              wz(i,j,k) = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
            end do
          end do
        end do

        do i = -1, 1
          sdx(i,:,:)  = wx(i,:,:)
          if (i>-1) sdx(i,:,:)=sdx(i,:,:)+sdx(i-1,:,:)
          cj(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1) = cj(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1)+sdx(i,:,:)
        end do        

        do j = -1, 1
          sdy(:,j,:)  = wy(:,j,:)
          if (j>-1) sdy(:,j,:)=sdy(:,j,:)+sdy(:,j-1,:)
          cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2) = cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2)+sdy(:,j,:)
        end do        

        do k = -1, 1
          sdz(:,:,k)  = wz(:,:,k)
          if (k>-1) sdz(:,:,k)=sdz(:,:,k)+sdz(:,:,k-1)
          cj(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+K,3) = cj(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+k,3)+sdz(:,:,k)
        end do        

    end do

  return
end subroutine depose_jxjyjz_esirkepov_linear_serial_old 

subroutine depose_rho_linear_serial(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz
   real(kind=8), dimension(-1:nx+1,-1:ny+1,-1:nz+1), intent(in out) :: rho
   real(kind=8), dimension(np) :: xp,yp,zp,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dzi,xint,yint,zint
   real(kind=8) :: x,y,z,wq,invvol,s1x,s2x,s1y,s2y,s1z,s2z
   integer(ISZ) :: j,k,l,ip,dix,diy,diz
   
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      invvol = dxi*dyi*dzi

      do ip=1,np
      
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi
        
        if (l_particles_weight) then
          wq=q*w(ip)*invvol
        else
          wq=q*invvol
        end if
      
        j=floor(x)
        k=floor(y)
        l=floor(z)

        xint = x-j
        yint = y-k
        zint = z-l

        s1x = 1.-xint
        s2x = xint

        s1y = 1.-yint
        s2y = yint

        s1z = 1.-zint
        s2z = zint

        rho(j  ,k  ,l  )=rho(j  ,k  ,l  )+s1x*s1y*s1z*wq
        rho(j+1,k  ,l  )=rho(j+1,k  ,l  )+s2x*s1y*s1z*wq
        rho(j  ,k+1,l  )=rho(j  ,k+1,l  )+s1x*s2y*s1z*wq
        rho(j+1,k+1,l  )=rho(j+1,k+1,l  )+s2x*s2y*s1z*wq
        rho(j  ,k  ,l+1)=rho(j  ,k  ,l+1)+s1x*s1y*s2z*wq
        rho(j+1,k  ,l+1)=rho(j+1,k  ,l+1)+s2x*s1y*s2z*wq
        rho(j  ,k+1,l+1)=rho(j  ,k+1,l+1)+s1x*s2y*s2z*wq
        rho(j+1,k+1,l+1)=rho(j+1,k+1,l+1)+s2x*s2y*s2z*wq
      
    end do

  return
end subroutine depose_rho_linear_serial

 subroutine getf3d_linear(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,exg,eyg,ezg)
   
      integer(ISZ) :: np,nx,ny,nz
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      real(kind=8), dimension(-1:nx+1,-1:ny+1,-1:nz+1) :: exg,eyg,ezg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      integer(ISZ) :: ip, j, k, l
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, s1x, s2x, s1y, s2y, s1z, s2z
      real(kind=8) :: w1, w2, w3, w4, w5, w6, w7, w8

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        j=floor(x)
        k=floor(y)
        l=floor(z)

        xint=x-j
        yint=y-k
        zint=z-l

        s1x=xint
        s2x=1.-xint
        s1y=yint
        s2y=1.-yint
        s1z=zint
        s2z=1.-zint

        w1 = s1x*s1y*s1z
        w2 = s2x*s1y*s1z
        w3 = s1x*s2y*s1z
        w4 = s2x*s2y*s1z
        w5 = s1x*s1y*s2z
        w6 = s2x*s1y*s2z
        w7 = s1x*s2y*s2z
        w8 = s2x*s2y*s2z
          
        ex(ip) = ex(ip)+w1*exg(j+1,k+1,l+1)+w2*exg(j,k+1,l+1)+w3*exg(j+1,k,l+1)+w4*exg(j,k,l+1) &
                       +w5*exg(j+1,k+1,l  )+w6*exg(j,k+1,l  )+w7*exg(j+1,k,l  )+w8*exg(j,k,l  )
        ey(ip) = ey(ip)+w1*eyg(j+1,k+1,l+1)+w2*eyg(j,k+1,l+1)+w3*eyg(j+1,k,l+1)+w4*eyg(j,k,l+1) &
                       +w5*eyg(j+1,k+1,l  )+w6*eyg(j,k+1,l  )+w7*eyg(j+1,k,l  )+w8*eyg(j,k,l  )
        ez(ip) = ez(ip)+w1*ezg(j+1,k+1,l+1)+w2*ezg(j,k+1,l+1)+w3*ezg(j+1,k,l+1)+w4*ezg(j,k,l+1) &
                       +w5*ezg(j+1,k+1,l  )+w6*ezg(j,k+1,l  )+w7*ezg(j+1,k,l  )+w8*ezg(j,k,l  )
     end do

   return
 end subroutine getf3d_linear
