#include "top.h"
subroutine depose_jxjyjz_esirkepov_linear_serial(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(-1:2,-1:2,-1:2) :: wx,wy,wz,sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz
   real(kind=8), DIMENSION(-1:2) :: sx, sy, sz, sx0, sy0, sz0, dsx, dsy, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k, &
                   ixmin, ixmax, iymin, iymax, izmin, izmax

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
          wq=q*w(1)
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

        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0

        ixmin = min(0,dix)
        ixmax = max(0,dix)
        iymin = min(0,diy)
        iymax = max(0,diy)
        izmin = min(0,diz)
        izmax = max(0,diz)

        do k=izmin, izmax+1
          do j=iymin, iymax+1
            do i=ixmin, ixmax+1
              wx(i,j,k) = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
              wy(i,j,k) = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
              wz(i,j,k) = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
            end do
          end do
        end do

        do i = ixmin, ixmax
          sdx(i,:,:)  = wx(i,:,:)
          if (i>ixmin) sdx(i,:,:)=sdx(i,:,:)+sdx(i-1,:,:)
          cj(iixp0+i,ijxp0+iymin:ijxp0+iymax+1,ikxp0+izmin:ikxp0+izmax+1,1) = &
          cj(iixp0+i,ijxp0+iymin:ijxp0+iymax+1,ikxp0+izmin:ikxp0+izmax+1,1) + sdx(i,iymin:iymax+1,izmin:izmax+1)
        end do        
        
        do j = iymin, iymax
          sdy(:,j,:)  = wy(:,j,:)
          if (j>iymin) sdy(:,j,:)=sdy(:,j,:)+sdy(:,j-1,:)
          cj(iixp0+ixmin:iixp0+ixmax+1,ijxp0+j,ikxp0+izmin:ikxp0+izmax+1,2) = &
          cj(iixp0+ixmin:iixp0+ixmax+1,ijxp0+j,ikxp0+izmin:ikxp0+izmax+1,2) + sdy(ixmin:ixmax+1,j,izmin:izmax+1)
        end do        

        do k = izmin, izmax
          sdz(:,:,k)  = wz(:,:,k)
          if (k>izmin) sdz(:,:,k)=sdz(:,:,k)+sdz(:,:,k-1)
          cj(iixp0+ixmin:iixp0+ixmax+1,ijxp0+iymin:ijxp0+iymax+1,ikxp0+K,3) = &
          cj(iixp0+ixmin:iixp0+ixmax+1,ijxp0+iymin:ijxp0+iymax+1,ikxp0+k,3) + sdz(ixmin:ixmax+1,iymin:iymax+1,k)
        end do        

    end do

  return
end subroutine depose_jxjyjz_esirkepov_linear_serial

subroutine depose_jxjyjz_esirkepov_n(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                                 nox,noy,noz,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(-int(nox/2)-1:int((nox+1)/2)+1, &
                          -int(noy/2)-1:int((noy+1)/2)+1, &
                          -int(noz/2)-1:int((noz+1)/2)+1) :: wx,wy,wz,sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8), DIMENSION(-int(nox/2)-1:int((nox+1)/2)+1) :: sx, sx0, dsx
   real(kind=8), DIMENSION(-int(noy/2)-1:int((noy+1)/2)+1) :: sy, sy0, dsy
   real(kind=8), DIMENSION(-int(noz/2)-1:int((noz+1)/2)+1) :: sz, sz0, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k, &
                   ixmin, ixmax, iymin, iymax, izmin, izmax

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
          wq=q*w(1)
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

        if (nox==1) then
          sx0( 0) = 1.-xint
          sx0( 1) = xint
        elseif (nox==2) then
          xint=xint-0.5
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
        end if

        if (noy==1) then
          sy0( 0) = 1.-yint
          sy0( 1) = yint
        elseif (noy==2) then
          yint=yint-0.5
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
        elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy0(-1) = onesixth*oyintsq*oyint
          sy0( 0) = twothird-yintsq*(1.-yint/2)
          sy0( 1) = twothird-oyintsq*(1.-oyint/2)
          sy0( 2) = onesixth*yintsq*yint
        end if

        if (noz==1) then
          sz0( 0) = 1.-zint
          sz0( 1) = zint
        elseif (noz==2) then
          zint=zint-0.5
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
        end if

        iixp=floor(xold)
        ijxp=floor(yold)
        ikxp=floor(zold)
        xint = xold-iixp
        yint = yold-ijxp
        zint = zold-ikxp

        dix = iixp-iixp0
        diy = ijxp-ijxp0
        diz = ikxp-ikxp0

        sx=0.;sy=0.;sz=0.

        if (nox==1) then
          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
        elseif (nox==2) then
          xint=xint-0.5
          xintsq = xint*xint
          sx(-1+dix) = 0.5*(0.5-xint)**2
          sx( 0+dix) = 0.75-xintsq
          sx( 1+dix) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1+dix) = onesixth*oxintsq*oxint
          sx( 0+dix) = twothird-xintsq*(1.-xint/2)
          sx( 1+dix) = twothird-oxintsq*(1.-oxint/2)
          sx( 2+dix) = onesixth*xintsq*xint
        end if

        if (noy==1) then
          sy( 0+diy) = 1.-yint
          sy( 1+diy) = yint
        elseif (noy==2) then
          yint=yint-0.5
          yintsq = yint*yint
          sy(-1+diy) = 0.5*(0.5-yint)**2
          sy( 0+diy) = 0.75-yintsq
          sy( 1+diy) = 0.5*(0.5+yint)**2
        elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1+diy) = onesixth*oyintsq*oyint
          sy( 0+diy) = twothird-yintsq*(1.-yint/2)
          sy( 1+diy) = twothird-oyintsq*(1.-oyint/2)
          sy( 2+diy) = onesixth*yintsq*yint
        end if

        if (noz==1) then
          sz( 0+diz) = 1.-zint
          sz( 1+diz) = zint
        elseif (noz==2) then
          zint=zint-0.5
          zintsq = zint*zint
          sz(-1+diz) = 0.5*(0.5-zint)**2
          sz( 0+diz) = 0.75-zintsq
          sz( 1+diz) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1+diz) = onesixth*ozintsq*ozint
          sz( 0+diz) = twothird-zintsq*(1.-zint/2)
          sz( 1+diz) = twothird-ozintsq*(1.-ozint/2)
          sz( 2+diz) = onesixth*zintsq*zint
        end if

        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0
        
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox-1)/2)
        iymin = min(0,diy)-int(noy/2)
        iymax = max(0,diy)+int((noy-1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz-1)/2)
!        ixmin = min(0,dix)-(nox-1)
!        ixmax = max(0,dix)+(nox-1)
!        iymin = min(0,diy)-(noy-1)
!        iymax = max(0,diy)+(noy-1)
!        izmin = min(0,diz)-(noz-1)
!        izmax = max(0,diz)+(noz-1)

        do k=izmin, izmax+1
          do j=iymin, iymax+1
            do i=ixmin, ixmax+1
              wx(i,j,k) = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
              wy(i,j,k) = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
              wz(i,j,k) = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
            end do
          end do
        end do

        do i = ixmin, ixmax
          sdx(i,:,:)  = wx(i,:,:)
          if (i>ixmin) sdx(i,:,:)=sdx(i,:,:)+sdx(i-1,:,:)
          cj(iixp0+i,ijxp0+iymin:ijxp0+iymax+1,ikxp0+izmin:ikxp0+izmax+1,1) = &
          cj(iixp0+i,ijxp0+iymin:ijxp0+iymax+1,ikxp0+izmin:ikxp0+izmax+1,1) + sdx(i,iymin:iymax+1,izmin:izmax+1)
        end do        
        
        do j = iymin, iymax
          sdy(:,j,:)  = wy(:,j,:)
          if (j>iymin) sdy(:,j,:)=sdy(:,j,:)+sdy(:,j-1,:)
          cj(iixp0+ixmin:iixp0+ixmax+1,ijxp0+j,ikxp0+izmin:ikxp0+izmax+1,2) = &
          cj(iixp0+ixmin:iixp0+ixmax+1,ijxp0+j,ikxp0+izmin:ikxp0+izmax+1,2) + sdy(ixmin:ixmax+1,j,izmin:izmax+1)
        end do        

        do k = izmin, izmax
          sdz(:,:,k)  = wz(:,:,k)
          if (k>izmin) sdz(:,:,k)=sdz(:,:,k)+sdz(:,:,k-1)
          cj(iixp0+ixmin:iixp0+ixmax+1,ijxp0+iymin:ijxp0+iymax+1,ikxp0+K,3) = &
          cj(iixp0+ixmin:iixp0+ixmax+1,ijxp0+iymin:ijxp0+iymax+1,ikxp0+k,3) + sdz(ixmin:ixmax+1,iymin:iymax+1,k)
        end do        

    end do

  return
end subroutine depose_jxjyjz_esirkepov_n

subroutine depose_jxjyjz_pxpypz_esirkepov_linear_serial(cj,mp,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,m,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,l_particles_weight,l_relativ)
   ! mp is the the kinetic energy density
    use Constant
  implicit none
   integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj,mp
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,m,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight, l_relativ

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint,wp,vxsq,vysq,vzsq,vsq
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

        vxsq = vx**2
        vysq = vy**2
        vzsq = vz**2

        if (l_particles_weight) then
          wq=q*w(ip)
          if (l_relativ) then
            vsq = vxsq+vysq+vzsq
            wp=m*w(ip)*(1./gaminv(ip)-1.)*clight**2/vsq
          else
            wp=0.5*m*w(ip)
          end if
        else
          wq=q*w(1)
          if (l_relativ) then
            vsq = vxsq+vysq+vzsq
            wp=m*w(1)*(1./gaminv(ip)-1.)*clight**2/vsq
          else
            wp=0.5*m*w(1)
          end if
        end if

        wqx = invdtdx
        wqy = invdtdy
        wqz = invdtdz

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
          cj(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1) = cj(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1)+wq*sdx(i,:,:)
          mp(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1) = mp(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1)+wp*vx*sdx(i,:,:)
        end do        

        do j = -1, 1
          sdy(:,j,:)  = wy(:,j,:)
          if (j>-1) sdy(:,j,:)=sdy(:,j,:)+sdy(:,j-1,:)
          cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2) = cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2)+wq*sdy(:,j,:)
          mp(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2) = mp(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2)+wp*vy*sdy(:,j,:)
        end do        

        do k = -1, 1
          sdz(:,:,k)  = wz(:,:,k)
          if (k>-1) sdz(:,:,k)=sdz(:,:,k)+sdz(:,:,k-1)
          cj(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+K,3) = cj(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+k,3)+wq*sdz(:,:,k)
          mp(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+K,3) = mp(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+k,3)+wp*vz*sdz(:,:,k)
        end do        

    end do

  return
end subroutine depose_jxjyjz_pxpypz_esirkepov_linear_serial

subroutine deposcor_jxjyjz_pxpypz_esirkepov_linear_serial(cj,mp,bx,by,bz, &
                                                          np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,m,xmin,ymin,zmin, &
                                                          dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,l_particles_weight,l_relativ)
   ! mp is the the kinetic energy density
    use Constant
  implicit none
   integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj,mp
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: bx,by,bz
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,m,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight, l_relativ

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint,wp,vxsq,vysq,vzsq,vsq
   real(kind=8),dimension(-1:2,-1:2,-1:2) :: wx,wy,wz,sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz,vol,dtcoef,total_field_density, &
                   thispart_field_density,total_kinetic_energy,thispart_kinetic_energy
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
      vol = dx*dy*dz
      invvol = 1./vol
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

        vxsq = vx**2
        vysq = vy**2
        vzsq = vz**2

        if (l_particles_weight) then
          wq=q*w(ip)
          if (l_relativ) then
            vsq = vxsq+vysq+vzsq
            wp=m*w(ip)*(1./gaminv(ip)-1.)*clight**2/vsq
          else
            wp=0.5*m*w(ip)
          end if
        else
          wq=q*w(1)
          if (l_relativ) then
            vsq = vxsq+vysq+vzsq
            wp=m*w(1)*(1./gaminv(ip)-1.)*clight**2/vsq
          else
            wp=0.5*m*w(1)
          end if
        end if

        wqx = invdtdx
        wqy = invdtdy
        wqz = invdtdz

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
        
        ! first, estimate dteff
!    print 'sum',sum(self.field.Mp[...,-1]*w3d.dz), \
!          sum(emass*(1./where(top.pgroup.gaminv==0.,1.,top.pgroup.gaminv)-1.)*top.pgroup.sw[0]*clight**2), \
!          sum(self.field.J[...,-1]**2)*(top.dt/eps0)**2*eps0*w3d.dz/2
        dtcoef = 1.
        
        do k=-1, 2
         do j=-1, 2
          do i = -1, 1
           sdx(i,:,:)  = wx(i,:,:)
           if (i>-1) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
           total_field_density = (0.5*(cj(iixp0+i,ijxp0+j,ikxp0+k,1)*dt)**2/eps0)*vol
           thispart_field_density = (0.5*(wq*sdx(i,j,k)*dt)**2/eps0)*vol
           total_kinetic_energy = mp(iixp0+i,ijxp0+j,ikxp0+k,1)
           thispart_kinetic_energy = wp*vx*sdx(i,j,k)
           if (total_field_density>total_kinetic_energy) then
           end if
          end do        
         end do        
        end do        

        do j = -1, 1
          sdy(:,j,:)  = wy(:,j,:)
          if (j>-1) sdy(:,j,:)=sdy(:,j,:)+sdy(:,j-1,:)
          cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2) = cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2)+wq*sdy(:,j,:)
          mp(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2) = mp(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2)+wp*vy*sdy(:,j,:)
        end do        
        
        do k=-1, 1
         do j=-1, 2
          do i = -1, 2
           sdz(i,j,k)  = wz(i,j,k)
           if (k>-1) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)

           total_field_density     = (0.5*(cj(iixp0+i,ijxp0+j,ikxp0+k,3)*dt)**2/eps0)*vol

           thispart_field_density  = (0.5*(wq*sdz(i,j,k)*dt)**2/eps0)*vol

           total_kinetic_energy    = mp(iixp0+i,ijxp0+j,ikxp0+k,3)

           thispart_kinetic_energy = wp*vx*sdz(i,j,k)

!           if ( (total_field_density>total_kinetic_energy) .and. (sign(1.,wq*sdz(i,j,k))==sign(1.,cj(iixp0+i,ijxp0+j,ikxp0+k,3))) ) then
!             dtcoef = 
!           end if
          end do        
         end do        
        end do        

        
        ! second, fix new positions and current

        do i = -1, 1
          sdx(i,:,:)  = wx(i,:,:)
          if (i>-1) sdx(i,:,:)=sdx(i,:,:)+sdx(i-1,:,:)
          cj(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1) = cj(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1)+wq*sdx(i,:,:)
          mp(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1) = mp(iixp0+i,ijxp0-1:ijxp0+2,ikxp0-1:ikxp0+2,1)+wp*vx*sdx(i,:,:)
        end do        

        do j = -1, 1
          sdy(:,j,:)  = wy(:,j,:)
          if (j>-1) sdy(:,j,:)=sdy(:,j,:)+sdy(:,j-1,:)
          cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2) = cj(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2)+wq*sdy(:,j,:)
          mp(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2) = mp(iixp0-1:iixp0+2,ijxp0+j,ikxp0-1:ikxp0+2,2)+wp*vy*sdy(:,j,:)
        end do        

        do k = -1, 1
          sdz(:,:,k)  = wz(:,:,k)
          if (k>-1) sdz(:,:,k)=sdz(:,:,k)+sdz(:,:,k-1)
          cj(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+K,3) = cj(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+k,3)+wq*sdz(:,:,k)
          mp(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+K,3) = mp(iixp0-1:iixp0+2,ijxp0-1:ijxp0+2,ikxp0+k,3)+wp*vz*sdz(:,:,k)
        end do        

    end do

  return
end subroutine deposcor_jxjyjz_pxpypz_esirkepov_linear_serial

subroutine depose_rho_linear_serial(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: rho
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

subroutine depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: rho
   real(kind=8), dimension(np) :: xp,yp,zp,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight

   real(kind=8) :: dxi,dyi,dzi,xint,yint,zint, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq
   real(kind=8) :: x,y,z,wq,invvol
   real(kind=8) :: sx(-int(nox/2):int((nox+1)/2)), &
                   sy(-int(noy/2):int((noy+1)/2)), &
                   sz(-int(noz/2):int((noz+1)/2))
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   integer(ISZ) :: j,k,l,ip,jj,kk,ll,ixmin, ixmax, iymin, iymax, izmin, izmax
   
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      invvol = dxi*dyi*dzi

      ixmin = -int(nox/2)
      ixmax = int((nox+1)/2)
      iymin = -int(noy/2)
      iymax = int((noy+1)/2)
      izmin = -int(noz/2)
      izmax = int((noz+1)/2)

      do ip=1,np
        
        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi
        
        j=floor(x)
        k=floor(y)
        l=floor(z)

        xint = x-j
        yint = y-k
        zint = z-l

        if (l_particles_weight) then
          wq=q*w(ip)*invvol
        else
          wq=q*invvol
        end if
      
        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
          xint=xint-0.5
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end if

        if (noy==1) then
          sy( 0) = 1.-yint
          sy( 1) = yint
        elseif (noy==2) then
          yint=yint-0.5
          yintsq = yint*yint
          sy(-1) = 0.5*(0.5-yint)**2
          sy( 0) = 0.75-yintsq
          sy( 1) = 0.5*(0.5+yint)**2
        elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1) = onesixth*oyintsq*oyint
          sy( 0) = twothird-yintsq*(1.-yint/2)
          sy( 1) = twothird-oyintsq*(1.-oyint/2)
          sy( 2) = onesixth*yintsq*yint
        end if

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
          zint=zint-0.5
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end if

         do ll = izmin, izmax
          do kk = iymin, iymax
            do jj = ixmin, ixmax
              rho(j+jj,k+kk,l+ll)=rho(j+jj,k+kk,l+ll)+sx(jj)*sy(kk)*sz(ll)*wq
            end do
          end do
        end do

    end do

  return
end subroutine depose_rho_n

 subroutine getf3d_linear(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,exg,eyg,ezg)
   
      integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
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

 subroutine gete3d_linear_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz, &
                                            nxguard,nyguard,nzguard,exg,eyg,ezg)
   
      integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      integer(ISZ) :: ip, j, k, l
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, s1x, s2x, s1y, s2y, s1z, s2z

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

        ex(ip) = ex(ip) + s1y*s1z*exg(j,k+1,l+1) &
                        + s2y*s1z*exg(j,k  ,l+1) &
                        + s1y*s2z*exg(j,k+1,l  ) &
                        + s2y*s2z*exg(j,k  ,l  )
                       
        ey(ip) = ey(ip) + s1x*s1z*eyg(j+1,k,l+1) &
                        + s2x*s1z*eyg(j  ,k,l+1) &
                        + s1x*s2z*eyg(j+1,k,l  ) &
                        + s2x*s2z*eyg(j  ,k,l  )
                       
        ez(ip) = ez(ip) + s1x*s1y*ezg(j+1,k+1,l) &
                        + s2x*s1y*ezg(j  ,k+1,l) &
                        + s1x*s2y*ezg(j+1,k  ,l) &
                        + s2x*s2y*ezg(j  ,k  ,l)
                       
     end do

   return
 end subroutine gete3d_linear_energy_conserving

  subroutine gete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                       nox,noy,noz,exg,eyg,ezg)
   
      integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      integer(ISZ) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
                      ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
                      xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
      real(kind=8), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(kind=8), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
      real(kind=8), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(kind=8), DIMENSION(-int((nox)/2):int((nox-1)/2)) :: sx0
      real(kind=8), DIMENSION(-int((noy)/2):int((noy-1)/2)) :: sy0
      real(kind=8), DIMENSION(-int((noz)/2):int((noz-1)/2)) :: sz0
      real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox-1)/2)
      iymin = -int(noy/2)
      iymax =  int((noy-1)/2)
      izmin = -int(noz/2)
      izmax =  int((noz-1)/2)

      ixmin0 = -int((nox)/2)
      ixmax0 =  int((nox-1)/2)
      iymin0 = -int((noy)/2)
      iymax0 =  int((noy-1)/2)
      izmin0 = -int((noz)/2)
      izmax0 =  int((noz-1)/2)

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

        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
          xint=xint-0.5
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end if

        if (noy==1) then
          sy( 0) = 1.-yint
          sy( 1) = yint
        elseif (noy==2) then
          yint=yint-0.5
          yintsq = yint*yint
          sy(-1) = 0.5*(0.5-yint)**2
          sy( 0) = 0.75-yintsq
          sy( 1) = 0.5*(0.5+yint)**2
        elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1) = onesixth*oyintsq*oyint
          sy( 0) = twothird-yintsq*(1.-yint/2)
          sy( 1) = twothird-oyintsq*(1.-oyint/2)
          sy( 2) = onesixth*yintsq*yint
        end if

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
          zint=zint-0.5
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end if

        if (nox==1) then
          sx0( 0) = 1.
        elseif (nox==2) then
          xint=xint+0.5
          sx0(-1) = 1.-xint
          sx0( 0) = xint
        elseif (nox==3) then
          xint=xint-0.5
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
        end if

        if (noy==1) then
          sy0( 0) = 1.
        elseif (noy==2) then
          yint=yint+0.5
          sy0(-1) = 1.-yint
          sy0( 0) = yint
        elseif (noy==3) then
          yint=yint-0.5
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
        end if

        if (noz==1) then
          sz0( 0) = 1.
        elseif (noz==2) then
          zint=zint+0.5
          sz0(-1) = 1.-zint
          sz0( 0) = zint
        elseif (noz==3) then
          zint=zint-0.5
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
        end if

        do ll = izmin, izmax+1
          do kk = iymin, iymax+1
            do jj = ixmin0, ixmax0
              ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j+jj,k+kk,l+ll)
            end do
          end do
        end do

        do ll = izmin, izmax+1
          do kk = iymin0, iymax0
            do jj = ixmin, ixmax+1
              ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k+kk,l+ll)
            end do
          end do
        end do

        do ll = izmin0, izmax0
          do kk = iymin, iymax+1
            do jj = ixmin, ixmax+1
              ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj,k+kk,l+ll)
            end do
          end do
        end do
                     
     end do

   return
 end subroutine gete3d_n_energy_conserving

 subroutine gete3d_linear_energy_conserving2(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz, &
                                             nxguard,nyguard,nzguard,exg,eyg,ezg)
   
      integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      integer(ISZ) :: ip, j, k, l
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, s1x, s2x, s1y, s2y, s1z, s2z

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi-0.5

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

        ex(ip) = ex(ip) + s1y*s1z*exg(j,k+1,l+1) &
                        + s2y*s1z*exg(j,k  ,l+1) &
                        + s1y*s2z*exg(j,k+1,l  ) &
                        + s2y*s2z*exg(j,k  ,l  )
                       
        ey(ip) = ey(ip) + s1x*s1z*eyg(j+1,k,l+1) &
                        + s2x*s1z*eyg(j  ,k,l+1) &
                        + s1x*s2z*eyg(j+1,k,l  ) &
                        + s2x*s2z*eyg(j  ,k,l  )
                       
        ez(ip) = ez(ip) + s1x*s1y*s2z*ezg(j+1,k+1,l) &
                        + s2x*s1y*s2z*ezg(j  ,k+1,l) &
                        + s1x*s2y*s2z*ezg(j+1,k  ,l) &
                        + s2x*s2y*s2z*ezg(j  ,k  ,l) &
                        + s1x*s1y*s1z*ezg(j+1,k+1,l+1) &
                        + s2x*s1y*s1z*ezg(j  ,k+1,l+1) &
                        + s1x*s2y*s1z*ezg(j+1,k  ,l+1) &
                        + s2x*s2y*s1z*ezg(j  ,k  ,l+1)
                       
     end do

   return
 end subroutine gete3d_linear_energy_conserving2

 subroutine getb3d_linear_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz, &
                                            nxguard,nyguard,nzguard,bxg,byg,bzg)
   
      integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,bx,by,bz
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      integer(ISZ) :: ip, j, k, l
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, s1x, s2x, s1y, s2y, s1z, s2z

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

        bx(ip) = bx(ip) + s1x*bxg(j  ,k,l) &
                        + s2x*bxg(j+1,k,l) 
                       
        by(ip) = by(ip) + s1y*byg(j,k  ,l) &
                        + s2y*byg(j,k+1,l) 
                       
        bz(ip) = bz(ip) + s1z*bzg(j,k,l  ) &
                        + s2z*bzg(j,k,l+1) 
                       
     end do

   return
 end subroutine getb3d_linear_energy_conserving

subroutine getb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                       nox,noy,noz,bxg,byg,bzg)
   
      integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,bx,by,bz
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      integer(ISZ) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
                      ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
                      xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq
      real(kind=8), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(kind=8), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
      real(kind=8), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(kind=8), DIMENSION(-int((nox-1)/2):int(nox/2)) :: sx0
      real(kind=8), DIMENSION(-int((noy-1)/2):int(noy/2)) :: sy0
      real(kind=8), DIMENSION(-int((noz-1)/2):int(noz/2)) :: sz0
      real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox-1)/2)
      iymin = -int(noy/2)
      iymax =  int((noy-1)/2)
      izmin = -int(noz/2)
      izmax =  int((noz-1)/2)

      ixmin0 = -int((nox)/2)
      ixmax0 =  int((nox-1)/2)
      iymin0 = -int((noy)/2)
      iymax0 =  int((noy-1)/2)
      izmin0 = -int((noz)/2)
      izmax0 =  int((noz-1)/2)

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

        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
          xint=xint-0.5
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
        elseif (nox==3) then
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end if

        if (noy==1) then
          sy( 0) = 1.-yint
          sy( 1) = yint
        elseif (noy==2) then
          yint=yint-0.5
          yintsq = yint*yint
          sy(-1) = 0.5*(0.5-yint)**2
          sy( 0) = 0.75-yintsq
          sy( 1) = 0.5*(0.5+yint)**2
        elseif (noy==3) then
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1) = onesixth*oyintsq*oyint
          sy( 0) = twothird-yintsq*(1.-yint/2)
          sy( 1) = twothird-oyintsq*(1.-oyint/2)
          sy( 2) = onesixth*yintsq*yint
        end if

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
          zint=zint-0.5
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
        elseif (noz==3) then
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end if

        if (nox==1) then
          sx0( 0) = 1.
        elseif (nox==2) then
          xint=xint+0.5
          sx0(-1) = 1.-xint
          sx0( 0) = xint
        elseif (nox==3) then
          xint=xint-0.5
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
        end if

        if (noy==1) then
          sy0( 0) = 1.
        elseif (noy==2) then
          yint=yint+0.5
          sy0(-1) = 1.-yint
          sy0( 0) = yint
        elseif (noy==3) then
          yint=yint-0.5
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
        end if

        if (noz==1) then
          sz0( 0) = 1.
        elseif (noz==2) then
          zint=zint+0.5
          sz0(-1) = 1.-zint
          sz0( 0) = zint
        elseif (noz==3) then
          zint=zint-0.5
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
        end if

        do ll = izmin0, izmax0
          do kk = iymin0, iymax0
            do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k+kk,l+ll)
            end do
          end do
        end do

        do ll = izmin0, izmax0
          do kk = iymin, iymax+1
            do jj = ixmin0, ixmax0
              by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j+jj,k+kk,l+ll)
            end do
          end do
        end do

        do ll = izmin, izmax+1
          do kk = iymin0, iymax0
            do jj = ixmin0, ixmax0
              bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j+jj,k+kk,l+ll)
            end do
          end do
        end do
                     
     end do

   return
 end subroutine getb3d_n_energy_conserving

