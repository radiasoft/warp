#include "top.h"

subroutine depose_jxjy_esirkepov_linear_serial_2d(j,np,xp,yp,xpold,ypold,uzp,gaminv,w,q,xmin,ymin,dt,dx,dy,nx,ny,l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny
   real(kind=8), dimension(-1:nx+1,-1:ny+1,3), intent(in out) :: j
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

        xmid = 0.5*(x+xold)
        ymid = 0.5*(y+yold)

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

        sx=0.;sy=0.
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
            
        j(iixp0,  ijxp0  ,3)=j(iixp0  ,ijxp0  ,3)-sd(10)
        j(iixp0+1,ijxp0  ,3)=j(iixp0+1,ijxp0  ,3)-sd(11)
        j(iixp0+2,ijxp0  ,3)=j(iixp0+2,ijxp0  ,3)-sd(12)
        j(iixp0,  ijxp0+1,3)=j(iixp0  ,ijxp0+1,3)-sd(13)
        j(iixp0+1,ijxp0+1,3)=j(iixp0+1,ijxp0+1,3)-sd(14)
        j(iixp0+2,ijxp0+1,3)=j(iixp0+2,ijxp0+1,3)-sd(15)
        j(iixp0,  ijxp0+2,3)=j(iixp0  ,ijxp0+2,3)-sd(16)
        j(iixp0+1,ijxp0+2,3)=j(iixp0+1,ijxp0+2,3)-sd(17)
        j(iixp0+2,ijxp0+2,3)=j(iixp0+2,ijxp0+2,3)-sd(18)
      
        ! Esirkepov deposition of Jx and Jz is over; now starts linear deposition of Jy
        wq = wq*vz*invsurf
      
        iixp=floor(xmid)
        ijxp=floor(ymid)

        xint = xmid-iixp
        yint = ymid-ijxp

        s1x = 1.-xint
        s2x = xint

        s1y = 1.-yint
        s2y = yint

        j(iixp  ,ijxp  ,2)=j(iixp  ,ijxp  ,2)+s1x*s1y*wq
        j(iixp+1,ijxp  ,2)=j(iixp+1,ijxp  ,2)+s2x*s1y*wq
        j(iixp  ,ijxp+1,2)=j(iixp  ,ijxp+1,2)+s1x*s2y*wq
        j(iixp+1,ijxp+1,2)=j(iixp+1,ijxp+1,2)+s2x*s2y*wq
        
      
    end do
    write(0,*) '*** sum j',sum(j(:,:,1)),sum(j(:,:,2)),sum(j(:,:,3))
  return
end subroutine depose_jxjy_esirkepov_linear_serial_2d

subroutine depose_jxjyjz_esirkepov_n_2d(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry,l_2drz)
   implicit none
   integer(ISZ) :: np,nx,nz,nox,noz,nxguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dz,xmin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry,l_2drz

   real(kind=8) :: dxi,dzi,dtsdx,dtsdz,xint,yint,zint
   real(kind=8),dimension(-int(nox/2)-1:int((nox+1)/2)+1, &
                          -int(noz/2)-1:int((noz+1)/2)+1) :: sdx,sdz
   real(kind=8) :: xold,yold,zold,rold,xmid,zmid,x,y,z,r,wq,wqx,wqz,tmp,vx,vy,vz,dts2dx,dts2dz, &
                   s1x,s2x,s1z,s2z,invvol,invdtdx,invdtdz, &
                   oxint,ozint,xintsq,zintsq,oxintsq,ozintsq, &
                   dtsdx0,dtsdz0,dts2dx0,dts2dz0
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8), DIMENSION(-int(nox/2)-1:int((nox+1)/2)+1) :: sx, sx0, dsx
   real(kind=8), DIMENSION(-int(noz/2)-1:int((noz+1)/2)+1) :: sz, sz0, dsz
   integer(ISZ) :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc, &
                   ixmin, ixmax, izmin, izmax, icell, ncells

    sx0=0.;sz0=0.
    sdz=0.
      
      dxi = 1./dx
      dzi = 1./dz
      invvol = 1./(dx*dz)
      dtsdx0 = dt*dxi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dz0 = 0.5*dtsdz0
      invdtdx = 1./(dt*dz)
      invdtdz = 1./(dt*dx)

      do ip=1,np
      
        x = xp(ip)
        if (l_2drz) then
          y = yp(ip)
          x=sqrt(x*x+y*y)
        end if
        x=x*dxi
        z = zp(ip)*dzi
          
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)

        if (l_2drz) then
          xold = xp(ip)-dt*vx
          yold = yp(ip)-dt*vy
          rold = sqrt(xold*xold+yold*yold)
          xold=rold*dxi
          vx = (x-xold)/dtsdx0
        else
          xold=x-dtsdx0*vx
        end if
        zold=z-dtsdz0*vz
 
        if (l4symtry) then
          x=abs(x)
          xold=abs(xold)
          vx = (x-xold)/dtsdx0
        end if

        x = x-xmin*dxi
        z = z-zmin*dzi
        xold = xold-xmin*dxi
        zold = zold-zmin*dzi
        
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1+max( int(abs(x-xold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        z=zold
        
        do icell = 1,ncells

        xold = x
        zold = z
        
        x = x+dtsdx*vx
        z = z+dtsdz*vz
        
        if (l_particles_weight) then
          wq=q*w(ip)
        else
          wq=q*w(1)
        end if
        wqx = wq*invdtdx
        wqz = wq*invdtdz
!       computation of current at x(n+1/2),v(n+1/2)

        iixp0=floor(x)
        ikxp0=floor(z)

        xint=x-iixp0
        zint=z-ikxp0

        select case(nox)
         case(0)
          sx0( 0) = 1.
         case(1)
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
        end select        

        select case(noz)
         case(0)
          sz0( 0) = 1.
         case(1)
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
        end select        

        iixp=floor(xold)
        ikxp=floor(zold)
        xint = xold-iixp
        zint = zold-ikxp

        dix = iixp-iixp0
        diz = ikxp-ikxp0

        sx=0.;sz=0.

        select case(nox)
         case(0)
          sx( 0+dix) = 1.
         case(1)
          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx(-1+dix) = 0.5*(0.5-xint)**2
          sx( 0+dix) = 0.75-xintsq
          sx( 1+dix) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1+dix) = onesixth*oxintsq*oxint
          sx( 0+dix) = twothird-xintsq*(1.-xint/2)
          sx( 1+dix) = twothird-oxintsq*(1.-oxint/2)
          sx( 2+dix) = onesixth*xintsq*xint
        end select        

        select case(noz)
         case(0)
          sz( 0+diz) = 1.
         case(1)
          sz( 0+diz) = 1.-zint
          sz( 1+diz) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz(-1+diz) = 0.5*(0.5-zint)**2
          sz( 0+diz) = 0.75-zintsq
          sz( 1+diz) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1+diz) = onesixth*ozintsq*ozint
          sz( 0+diz) = twothird-zintsq*(1.-zint/2)
          sz( 1+diz) = twothird-ozintsq*(1.-ozint/2)
          sz( 2+diz) = onesixth*zintsq*zint
        end select        

        dsx = sx - sx0
        dsz = sz - sz0
        
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox+1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz+1)/2)


        do k=izmin, izmax
            do i=ixmin, ixmax
              ic = iixp0+i
              kc = ikxp0+k
              if(i<ixmax) then
                sdx(i,k)  = wqx*dsx(i)*( sz0(k) + 0.5*dsz(k))
                if (i>ixmin) sdx(i,k)=sdx(i,k)+sdx(i-1,k)
                cj(ic,kc,1) = cj(ic,kc,1) + sdx(i,k)
              end if
              cj(ic,kc,2) = cj(ic,kc,2) + wq*vy*invvol*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))/ncells
              if(k<izmax) then
                sdz(i,k)  = wqz*dsz(k)*(sx0(i)+0.5*dsx(i)) 
                if (k>izmin) sdz(i,k)=sdz(i,k)+sdz(i,k-1)
                cj(ic,kc,3) = cj(ic,kc,3) + sdz(i,k)
              end if
          end do        
        end do        

      end do

    end do

  return
end subroutine depose_jxjyjz_esirkepov_n_2d

subroutine depose_jxjyjz_villasenor_n_2d(cj,np,xp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
                                                 dt,dx,dz,nx,nz,nxguard,nzguard, &
                                                 nox,noz,l_particles_weight,l4symtry)
   implicit none
   integer(ISZ) :: np,nx,nz,nox,noz,nxguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dz,xmin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry

   real(kind=8) :: dxi,dzi,dtsdx,dtsdz,xint,yint,zint
   real(kind=8),dimension(-int(nox/2)-1:int((nox+1)/2)+1, &
                          -int(noz/2)-1:int((noz+1)/2)+1) :: sdx,sdz
   real(kind=8) :: xold,zold,xmid,zmid,x,z,wq,wqx,wqz,tmp,vx,vy,vz,deltx,deltz
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8)::xnode(-nxguard:nx+nxguard),znode(-nzguard:nz+nzguard),length, &
                 x1,x2,z1,z2,qq(np),xx1,zz1,xx2,zz2,invvol,deltsx,deltsz,dtsdx0,dtsdz0, &
                 lengood,wr1,wr2,wz1,wz2,curx,curz,xmoy,zmoy,lengthx,lengthz
   integer(ISZ) :: iixp0,ikxp0,iixp,ikxp,ip,dix,diz,idx,idz,i,k,ic,kc, interx, interz, jn2,ln2,jni,lni,&
                   ixmin, ixmax, izmin, izmax, icell, ncells,jn,ln,jn1,ln1,j,jtot,nloop
   logical(ISZ) :: doit
   
      dxi = 1./dx
      dzi = 1./dz
      invvol = 1./(dx*dz)
      dtsdx0 = dt*dxi
      dtsdz0 = dt*dzi

      deltsx = 1. / 1.e10  
      deltsz = 1. / 1.e10  
      
      do jn = -nxguard, nx+nxguard
        xnode(jn) = xmin/dx + jn!*dx
      end do

      do ln = -nzguard, nz+nzguard
        znode(ln) = zmin/dz + ln!*dz
      end do

      do ip=1,np
      
        x = (xp(ip)-xmin)*dxi
        z = (zp(ip)-zmin)*dzi
        
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)
        
        xold=x-dtsdx0*vx
        zold=z-dtsdz0*vz
 
        if (l4symtry) then
          x=abs(x)
          xold=abs(xold)
          vx = (x-xold)/dtsdx0
        end if
       
        if (l_particles_weight) then
          wq=q*w(ip)*invvol
        else
          wq=q*w(1)*invvol
        end if

        jn1 = floor(xold)
        ln1 = floor(zold)
        jn2 = floor(x)
        ln2 = floor(z)
        jni = jn2
        lni = ln2
        
        x1 = xold
        z1 = zold
        x2 = x
        z2 = z

        doit = .true.

        nloop = 0
        do while(doit)
        
          deltx = x2-x1
          deltz = z2-z1

          interx = 0
          interz = 0
          
          if (jn1 .ne. jn2) then
            if (jn2.gt.jn1) then
              lengthx = abs ( (xnode (jn2) - xold ) / deltx)
              interx = 1
            else  
              lengthx = abs ( (xnode (jn1) - xold ) / deltx)
              interx = -1
            endif
          end if
             
          if (ln1 .ne. ln2) then
            if (ln2.gt.ln1) then
              lengthz = abs ( (znode (ln2) - zold ) / deltz)
              interz = 1
            else  
              lengthz = abs ( (znode (ln1) - zold ) / deltz)
              interz = -1
            endif
          end if

          if (interx/=0 .and. interz/=0) then
            if (lengthx<=lengthz) then
              interz = 0
            else 
              interx = 0
            end if
          end if

          if (interx/=0 .or. interz/=0) then
            if (interx/=0) then
              length = lengthx
            else 
              length = lengthz
            end if
            x2 = x1+length*deltx
            z2 = z1+length*deltz
            jn2=jn1+interx
            ln2=ln1+interz
          end if
          
        jn = jn1
        ln = ln1

        wr1 = (x1-xnode(jn))
        wr2 = (x2-xnode(jn))
!        wr1 = (x1**2-xnode(jn)**2)/((xnode(jn)+xnode(jn+1)))
!        wr2 = (x2**2-xnode(jn)**2)/((xnode(jn)+xnode(jn+1)))

        wz1 = (z1-znode(ln))
        wz2 = (z2-znode(ln))

        curx = (wr2-wr1) * dx * wq/dt
        curz = (wz2-wz1) * dz * wq/dt

        xmoy = 0.5*(wr1+wr2)
        zmoy = 0.5*(wz1+wz2)

        write(0,*) jn,ln,curx,curz,xmoy,zmoy,interx,interz,length

        cj (jn,  ln,  1) = cj (jn  ,ln  ,1) + curx * (1. - zmoy)
        cj (jn,  ln+1,1) = cj (jn  ,ln+1,1) + curx * zmoy
        cj (jn,  ln,  3) = cj (jn  ,ln  ,3) + curz * (1. - xmoy)
        cj (jn+1,ln,  3) = cj (jn+1,ln  ,3) + curz * xmoy

        if (interx/=0 .or. interz/=0) then
          x1 = x2
          z1 = z2
          x2 = x
          z2 = z
          jn1 = jn2
          ln1 = ln2
          jn2 = jni
          ln2 = lni
        endif
        nloop = nloop+1
        doit = nloop<3 .and. (interx/=0 .or. interz/=0) 
        write(0,*) 'doit',doit,nloop,interx,interz
      enddo  
    enddo  

  return
end subroutine depose_jxjyjz_villasenor_n_2d

subroutine depose_jxjyjz_esirkepov_linear_serialold(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
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
end subroutine depose_jxjyjz_esirkepov_linear_serialold

subroutine depose_jxjyjz_esirkepov_linear_serialnew(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                                 l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(-1:2, &
                          -1:2, &
                          -1:2) :: sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                   dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8), DIMENSION(-1:2) :: sx, sx0, dsx
   real(kind=8), DIMENSION(-1:2) :: sy, sy0, dsy
   real(kind=8), DIMENSION(-1:2) :: sz, sz0, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                   ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ixmin2, ixmax2, iymin2, iymax2, izmin2, izmax2 

    sx0=0.;sy0=0.;sz0=0.
    sdz=0.
      
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      dtsdx0 = dt*dxi
      dtsdy0 = dt*dyi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dy0 = 0.5*dtsdy0
      dts2dz0 = 0.5*dtsdz0
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
        
        xold=x-dtsdx0*vx
        yold=y-dtsdy0*vy
        zold=z-dtsdz0*vz
 
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdy = dtsdy0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dy = dts2dy0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        y=yold
        z=zold
        
        do icell = 1,ncells

        xold = x
        yold = y
        zold = z
        
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz
        
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

          sx0( 0) = 1.-xint
          sx0( 1) = xint
          sy0( 0) = 1.-yint
          sy0( 1) = yint
          sz0( 0) = 1.-zint
          sz0( 1) = zint

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

          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
          sy( 0+diy) = 1.-yint
          sy( 1+diy) = yint
          sz( 0+diz) = 1.-zint
          sz( 1+diz) = zint

        dsx = sx - sx0

        dsy = sy - sy0

        dsz = sz - sz0
        
        ixmin = min(0,dix)
        ixmax = max(0,dix)+1
        iymin = min(0,diy)
        iymax = max(0,diy)+1
        izmin = min(0,diz)
        izmax = max(0,diz)+1

        do k=izmin, izmax
          i=ixmin
          do j=iymin, iymax
            sdx(i,j,k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
          end do        
          do i=ixmin+1, ixmax-1
            do j=iymin, iymax
               sdx(i,j,k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
               sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
            end do        
          end do        
        end do        

        do k=izmin, izmax
          j=iymin
          do i=ixmin, ixmax
            sdy(i,j,k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
          end do
          do j=iymin+1, iymax-1
            do i=ixmin, ixmax
              sdy(i,j,k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
              sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
            end do        
          end do        
        end do        


        do j=iymin, iymax
          k=izmin
          do i=ixmin, ixmax
            sdz(i,j,k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
          end do
          do k=izmin+1, izmax-1
            do i=ixmin, ixmax
              sdz(i,j,k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
              sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
            end do        
          end do        
        end do        
        
        ixmin2 = iixp0+ixmin
        ixmax2 = iixp0+ixmax
        iymin2 = ijxp0+iymin
        iymax2 = ijxp0+iymax
        izmin2 = ikxp0+izmin
        izmax2 = ikxp0+izmax
        cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,1) = cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,1) &
                                                        + sdx(ixmin:ixmax,iymin:iymax,izmin:izmax)
        cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,2) = cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,2) &
                                                        + sdy(ixmin:ixmax,iymin:iymax,izmin:izmax)
        cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,3) = cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,3) &
                                                        + sdz(ixmin:ixmax,iymin:iymax,izmin:izmax)

      end do
 
    end do

  return
end subroutine depose_jxjyjz_esirkepov_linear_serialnew

subroutine depose_jxjyjz_esirkepov_linear_serial(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                                 l_particles_weight)
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(-1:2, &
                          -1:2, &
                          -1:2) :: sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                   dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8), DIMENSION(-1:2) :: sx, sx0, dsx
   real(kind=8), DIMENSION(-1:2) :: sy, sy0, dsy
   real(kind=8), DIMENSION(-1:2) :: sz, sz0, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                   ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells 

    sx0=0.;sy0=0.;sz0=0.
    sdz=0.
      
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      dtsdx0 = dt*dxi
      dtsdy0 = dt*dyi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dy0 = 0.5*dtsdy0
      dts2dz0 = 0.5*dtsdz0
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
        
        xold=x-dtsdx0*vx
        yold=y-dtsdy0*vy
        zold=z-dtsdz0*vz
 
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdy = dtsdy0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dy = dts2dy0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        y=yold
        z=zold
        
        do icell = 1,ncells

        xold = x
        yold = y
        zold = z
        
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz
        
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

          sx0( 0) = 1.-xint
          sx0( 1) = xint
          sy0( 0) = 1.-yint
          sy0( 1) = yint
          sz0( 0) = 1.-zint
          sz0( 1) = zint

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

          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
          sy( 0+diy) = 1.-yint
          sy( 1+diy) = yint
          sz( 0+diz) = 1.-zint
          sz( 1+diz) = zint

        dsx = sx - sx0

        dsy = sy - sy0

        dsz = sz - sz0
        
        ixmin = min(0,dix)
        ixmax = max(0,dix)+1
        iymin = min(0,diy)
        iymax = max(0,diy)+1
        izmin = min(0,diz)
        izmax = max(0,diz)+1

        do k=izmin, izmax
          do j=iymin, iymax
            do i=ixmin, ixmax
              ic = iixp0+i
              jc = ijxp0+j
              kc = ikxp0+k
              if(i<ixmax) then
                sdx(i,j,k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
                if (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                cj(ic,jc,kc,1) = cj(ic,jc,kc,1) + sdx(i,j,k)
              end if
              if(j<iymax) then
!        write(0,*) ip,np,i,j,k,wqy,dsy(j),yold,y,dtsdy,vy, ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
                sdy(i,j,k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
                if (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                cj(ic,jc,kc,2) = cj(ic,jc,kc,2) + sdy(i,j,k)
              end if
              if(k<izmax) then
                sdz(i,j,k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
                if (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                cj(ic,jc,kc,3) = cj(ic,jc,kc,3) + sdz(i,j,k)
              end if
            end do        
          end do        
        end do        

      end do
 
    end do

  return
end subroutine depose_jxjyjz_esirkepov_linear_serial

subroutine depose_jxjyjz_esirkepov_nnew(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                                 nox,noy,noz,l_particles_weight,l4symtry)
! although it vectorizes better, this version is slower than the old one
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(-int(nox/2)-1:int((nox+1)/2)+1, &
                          -int(noy/2)-1:int((noy+1)/2)+1, &
                          -int(noz/2)-1:int((noz+1)/2)+1) :: sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                   dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8), DIMENSION(-int(nox/2)-1:int((nox+1)/2)+1) :: sx, sx0, dsx
   real(kind=8), DIMENSION(-int(noy/2)-1:int((noy+1)/2)+1) :: sy, sy0, dsy
   real(kind=8), DIMENSION(-int(noz/2)-1:int((noz+1)/2)+1) :: sz, sz0, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                   ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells, ixmin2, ixmax2, iymin2, iymax2, izmin2, izmax2 

    sx0=0.;sy0=0.;sz0=0.
    sdz=0.
      
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      dtsdx0 = dt*dxi
      dtsdy0 = dt*dyi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dy0 = 0.5*dtsdy0
      dts2dz0 = 0.5*dtsdz0
      invvol = 1./(dx*dy*dz)
      invdtdx = 1./(dt*dy*dz)
      invdtdy = 1./(dt*dx*dz)
      invdtdz = 1./(dt*dx*dy)

      do ip=1,np
      
        x = xp(ip)*dxi
        y = yp(ip)*dyi
        z = zp(ip)*dzi
        
        vx = uxp(ip)*gaminv(ip)
        vy = uyp(ip)*gaminv(ip)
        vz = uzp(ip)*gaminv(ip)
        
        xold=x-dtsdx0*vx
        yold=y-dtsdy0*vy
        zold=z-dtsdz0*vz
 
        if (l4symtry) then
          x=abs(x)
          y=abs(y)
          xold=abs(xold)
          yold=abs(yold)
          vx = (x-xold)/dtsdx0
          vy = (y-yold)/dtsdy0
        end if
        
        x = x - xmin*dxi
        y = y - ymin*dyi
        z = z - zmin*dzi
        
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdy = dtsdy0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dy = dts2dy0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        y=yold
        z=zold
        
        do icell = 1,ncells

        xold = x
        yold = y
        zold = z
        
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz
        
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

        select case(nox)
         case(0)
          sx0( 0) = 1.
         case(1)
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy0( 0) = 1.
         case(1)
          sy0( 0) = 1.-yint
          sy0( 1) = yint
         case(2)
!          yint=yint-0.5
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy0(-1) = onesixth*oyintsq*oyint
          sy0( 0) = twothird-yintsq*(1.-yint/2)
          sy0( 1) = twothird-oyintsq*(1.-oyint/2)
          sy0( 2) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz0( 0) = 1.
         case(1)
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
        end select        

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

        select case(nox)
         case(0)
          sx( 0+dix) = 1.
         case(1)
          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx(-1+dix) = 0.5*(0.5-xint)**2
          sx( 0+dix) = 0.75-xintsq
          sx( 1+dix) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1+dix) = onesixth*oxintsq*oxint
          sx( 0+dix) = twothird-xintsq*(1.-xint/2)
          sx( 1+dix) = twothird-oxintsq*(1.-oxint/2)
          sx( 2+dix) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy( 0+diy) = 1.
         case(1)
          sy( 0+diy) = 1.-yint
          sy( 1+diy) = yint
         case(2)
!          yint=yint-0.5
          yintsq = yint*yint
          sy(-1+diy) = 0.5*(0.5-yint)**2
          sy( 0+diy) = 0.75-yintsq
          sy( 1+diy) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1+diy) = onesixth*oyintsq*oyint
          sy( 0+diy) = twothird-yintsq*(1.-yint/2)
          sy( 1+diy) = twothird-oyintsq*(1.-oyint/2)
          sy( 2+diy) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz( 0+diz) = 1.
         case(1)
          sz( 0+diz) = 1.-zint
          sz( 1+diz) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz(-1+diz) = 0.5*(0.5-zint)**2
          sz( 0+diz) = 0.75-zintsq
          sz( 1+diz) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1+diz) = onesixth*ozintsq*ozint
          sz( 0+diz) = twothird-zintsq*(1.-zint/2)
          sz( 1+diz) = twothird-ozintsq*(1.-ozint/2)
          sz( 2+diz) = onesixth*zintsq*zint
        end select        

        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0
        
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox+1)/2)
        iymin = min(0,diy)-int(noy/2)
        iymax = max(0,diy)+int((noy+1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz+1)/2)

         do k=izmin, izmax
          i=ixmin
          do j=iymin, iymax
            sdx(i,j,k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
          end do        
          do i=ixmin+1, ixmax-1
            do j=iymin, iymax
               sdx(i,j,k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
               sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
            end do        
          end do        
        end do        

        do k=izmin, izmax
          j=iymin
          do i=ixmin, ixmax
            sdy(i,j,k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
          end do
          do j=iymin+1, iymax-1
            do i=ixmin, ixmax
              sdy(i,j,k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
              sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
            end do        
          end do        
        end do        


        do j=iymin, iymax
          k=izmin
          do i=ixmin, ixmax
            sdz(i,j,k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
          end do
          do k=izmin+1, izmax-1
            do i=ixmin, ixmax
              sdz(i,j,k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
              sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
            end do        
          end do        
        end do        
        
        ixmin2 = iixp0+ixmin
        ixmax2 = iixp0+ixmax
        iymin2 = ijxp0+iymin
        iymax2 = ijxp0+iymax
        izmin2 = ikxp0+izmin
        izmax2 = ikxp0+izmax
        cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,1) = cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,1) &
                                                        + sdx(ixmin:ixmax,iymin:iymax,izmin:izmax)
        cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,2) = cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,2) &
                                                        + sdy(ixmin:ixmax,iymin:iymax,izmin:izmax)
        cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,3) = cj(ixmin2:ixmax2,iymin2:iymax2,izmin2:izmax2,3) &
                                                        + sdz(ixmin:ixmax,iymin:iymax,izmin:izmax)


      end do
 
    end do

  return
end subroutine depose_jxjyjz_esirkepov_nnew

subroutine depose_jxjyjz_esirkepov_n(cj,np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin, &
                                                 dt,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                                 nox,noy,noz,l_particles_weight,l4symtry)
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,yp,zp,uxp,uyp,uzp,gaminv,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry

   real(kind=8) :: dxi,dyi,dzi,dtsdx,dtsdy,dtsdz,xint,yint,zint
   real(kind=8),dimension(-int(nox/2)-1:int((nox+1)/2)+1, &
                          -int(noy/2)-1:int((noy+1)/2)+1, &
                          -int(noz/2)-1:int((noz+1)/2)+1) :: sdx,sdy,sdz
   real(kind=8) :: xold,yold,zold,xmid,ymid,zmid,x,y,z,wq,wqx,wqy,wqz,tmp,vx,vy,vz,dts2dx,dts2dy,dts2dz, &
                   s1x,s2x,s1y,s2y,s1z,s2z,invvol,invdtdx,invdtdy,invdtdz, &
                   oxint,oyint,ozint,xintsq,yintsq,zintsq,oxintsq,oyintsq,ozintsq, &
                   dtsdx0,dtsdy0,dtsdz0,dts2dx0,dts2dy0,dts2dz0
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   real(kind=8), DIMENSION(-int(nox/2)-1:int((nox+1)/2)+1) :: sx, sx0, dsx
   real(kind=8), DIMENSION(-int(noy/2)-1:int((noy+1)/2)+1) :: sy, sy0, dsy
   real(kind=8), DIMENSION(-int(noz/2)-1:int((noz+1)/2)+1) :: sz, sz0, dsz
   integer(ISZ) :: iixp0,ijxp0,ikxp0,iixp,ijxp,ikxp,ip,dix,diy,diz,idx,idy,idz,i,j,k,ic,jc,kc, &
                   ixmin, ixmax, iymin, iymax, izmin, izmax, icell, ncells

    sx0=0.;sy0=0.;sz0=0.
    sdz=0.
      
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      dtsdx0 = dt*dxi
      dtsdy0 = dt*dyi
      dtsdz0 = dt*dzi
      dts2dx0 = 0.5*dtsdx0
      dts2dy0 = 0.5*dtsdy0
      dts2dz0 = 0.5*dtsdz0
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
        
        xold=x-dtsdx0*vx
        yold=y-dtsdy0*vy
        zold=z-dtsdz0*vz
 
        if (l4symtry) then
          x=abs(x)
          y=abs(y)
          xold=abs(xold)
          yold=abs(yold)
          vx = (x-xold)/dtsdx0
          vy = (y-yold)/dtsdy0
        end if
        
        ! computes maximum number of cells traversed by particle in a given dimension
        ncells = 1+max( int(abs(x-xold)), int(abs(y-yold)), int(abs(z-zold)))
        
        dtsdx = dtsdx0/ncells
        dtsdy = dtsdy0/ncells
        dtsdz = dtsdz0/ncells
        dts2dx = dts2dx0/ncells
        dts2dy = dts2dy0/ncells
        dts2dz = dts2dz0/ncells
        
        x=xold
        y=yold
        z=zold
        
        do icell = 1,ncells

        xold = x
        yold = y
        zold = z
        
        x = x+dtsdx*vx
        y = y+dtsdy*vy
        z = z+dtsdz*vz
        
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

        select case(nox)
         case(0)
          sx0( 0) = 1.
         case(1)
          sx0( 0) = 1.-xint
          sx0( 1) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx0(-1) = onesixth*oxintsq*oxint
          sx0( 0) = twothird-xintsq*(1.-xint/2)
          sx0( 1) = twothird-oxintsq*(1.-oxint/2)
          sx0( 2) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy0( 0) = 1.
         case(1)
          sy0( 0) = 1.-yint
          sy0( 1) = yint
         case(2)
!          yint=yint-0.5
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy0(-1) = onesixth*oyintsq*oyint
          sy0( 0) = twothird-yintsq*(1.-yint/2)
          sy0( 1) = twothird-oyintsq*(1.-oyint/2)
          sy0( 2) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz0( 0) = 1.
         case(1)
          sz0( 0) = 1.-zint
          sz0( 1) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz0(-1) = onesixth*ozintsq*ozint
          sz0( 0) = twothird-zintsq*(1.-zint/2)
          sz0( 1) = twothird-ozintsq*(1.-ozint/2)
          sz0( 2) = onesixth*zintsq*zint
        end select        

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

        select case(nox)
         case(0)
          sx( 0+dix) = 1.
         case(1)
          sx( 0+dix) = 1.-xint
          sx( 1+dix) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx(-1+dix) = 0.5*(0.5-xint)**2
          sx( 0+dix) = 0.75-xintsq
          sx( 1+dix) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1+dix) = onesixth*oxintsq*oxint
          sx( 0+dix) = twothird-xintsq*(1.-xint/2)
          sx( 1+dix) = twothird-oxintsq*(1.-oxint/2)
          sx( 2+dix) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy( 0+diy) = 1.
         case(1)
          sy( 0+diy) = 1.-yint
          sy( 1+diy) = yint
         case(2)
!          yint=yint-0.5
          yintsq = yint*yint
          sy(-1+diy) = 0.5*(0.5-yint)**2
          sy( 0+diy) = 0.75-yintsq
          sy( 1+diy) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1+diy) = onesixth*oyintsq*oyint
          sy( 0+diy) = twothird-yintsq*(1.-yint/2)
          sy( 1+diy) = twothird-oyintsq*(1.-oyint/2)
          sy( 2+diy) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz( 0+diz) = 1.
         case(1)
          sz( 0+diz) = 1.-zint
          sz( 1+diz) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz(-1+diz) = 0.5*(0.5-zint)**2
          sz( 0+diz) = 0.75-zintsq
          sz( 1+diz) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1+diz) = onesixth*ozintsq*ozint
          sz( 0+diz) = twothird-zintsq*(1.-zint/2)
          sz( 1+diz) = twothird-ozintsq*(1.-ozint/2)
          sz( 2+diz) = onesixth*zintsq*zint
        end select        

        dsx = sx - sx0
        dsy = sy - sy0
        dsz = sz - sz0
        
        ixmin = min(0,dix)-int(nox/2)
        ixmax = max(0,dix)+int((nox+1)/2)
        iymin = min(0,diy)-int(noy/2)
        iymax = max(0,diy)+int((noy+1)/2)
        izmin = min(0,diz)-int(noz/2)
        izmax = max(0,diz)+int((noz+1)/2)

        do k=izmin, izmax
          do j=iymin, iymax
            do i=ixmin, ixmax
              ic = iixp0+i
              jc = ijxp0+j
              kc = ikxp0+k
              if(i<ixmax) then
                sdx(i,j,k)  = wqx*dsx(i)*( (sy0(j)+0.5*dsy(j))*sz0(k) + (0.5*sy0(j)+1./3.*dsy(j))*dsz(k))
                if (i>ixmin) sdx(i,j,k)=sdx(i,j,k)+sdx(i-1,j,k)
                cj(ic,jc,kc,1) = cj(ic,jc,kc,1) + sdx(i,j,k)
              end if
              if(j<iymax) then
!        write(0,*) ip,np,i,j,k,wqy,dsy(j),yold,y,dtsdy,vy, ( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
                sdy(i,j,k)  = wqy*dsy(j)*( (sz0(k)+0.5*dsz(k))*sx0(i) + (0.5*sz0(k)+1./3.*dsz(k))*dsx(i))
                if (j>iymin) sdy(i,j,k)=sdy(i,j,k)+sdy(i,j-1,k)
                cj(ic,jc,kc,2) = cj(ic,jc,kc,2) + sdy(i,j,k)
              end if
              if(k<izmax) then
                sdz(i,j,k)  = wqz*dsz(k)*( (sx0(i)+0.5*dsx(i))*sy0(j) + (0.5*sx0(i)+1./3.*dsx(i))*dsy(j))
                if (k>izmin) sdz(i,j,k)=sdz(i,j,k)+sdz(i,j,k-1)
                cj(ic,jc,kc,3) = cj(ic,jc,kc,3) + sdz(i,j,k)
              end if
            end do        
          end do        
        end do        

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

subroutine depose_rho_n_2dxz(rho,np,xp,yp,zp,w,q,xmin,zmin,dx,dz,nx,nz,nxguard,nzguard,nox,noz, &
                        l_particles_weight,l4symtry,l_2drz)
   implicit none
   integer(ISZ) :: np,nx,nz,nox,noz,nxguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,0:0,-nzguard:nz+nzguard), intent(in out) :: rho
   real(kind=8), dimension(np) :: xp,yp,zp,w
   real(kind=8) :: q,dt,dx,dz,xmin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry,l_2drz

   real(kind=8) :: dxi,dzi,xint,zint, &
                   oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
   real(kind=8) :: x,z,r,wq,invvol
   real(kind=8) :: sx(-int(nox/2):int((nox+1)/2)), &
                   sz(-int(noz/2):int((noz+1)/2))
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   integer(ISZ) :: j,l,ip,jj,ll,ixmin, ixmax, izmin, izmax
   
      dxi = 1./dx
      dzi = 1./dz
      invvol = dxi*dzi

      ixmin = -int(nox/2)
      ixmax = int((nox+1)/2)
      izmin = -int(noz/2)
      izmax = int((noz+1)/2)

      do ip=1,np
        
        if (l_2drz) then
          r = sqrt(xp(ip)*xp(ip)+yp(ip)*yp(ip))
          x = (r-xmin)*dxi
          z = (zp(ip)-zmin)*dzi
        else
          x = (xp(ip)-xmin)*dxi
          z = (zp(ip)-zmin)*dzi
        end if
        
        if (l4symtry) then
          x=abs(x)
        end if
        
        j=floor(x)
        l=floor(z)

        xint = x-j
        zint = z-l

        if (l_particles_weight) then
          wq=q*w(ip)*invvol
        else
          wq=q*invvol
        end if
      
        select case(nox)
         case(0)
          sx( 0) = 1.
         case(1)
          sx( 0) = 1.-xint
          sx( 1) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end select        

        select case(noz)
         case(0)
          sz( 0) = 1.
         case(1)
          sz( 0) = 1.-zint
          sz( 1) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end select        

         do ll = izmin, izmax
            do jj = ixmin, ixmax
              rho(j+jj,0,l+ll)=rho(j+jj,0,l+ll)+sx(jj)*sz(ll)*wq
            end do
        end do

    end do

  return
end subroutine depose_rho_n_2dxz

subroutine depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz, &
                        l_particles_weight,l4symtry)
   implicit none
   integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in out) :: rho
   real(kind=8), dimension(np) :: xp,yp,zp,w
   real(kind=8) :: q,dt,dx,dy,dz,xmin,ymin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry

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
        
        if (l4symtry) then
          x=abs(x)
          y=abs(y)
        end if
        
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
      
        select case(nox)
         case(0)
          sx( 0) = 1.
         case(1)
          sx( 0) = 1.-xint
          sx( 1) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy( 0) = 1.
         case(1)
          sy( 0) = 1.-yint
          sy( 1) = yint
         case(2)
!          yint=yint-0.5
          yintsq = yint*yint
          sy(-1) = 0.5*(0.5-yint)**2
          sy( 0) = 0.75-yintsq
          sy( 1) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1) = onesixth*oyintsq*oyint
          sy( 0) = twothird-yintsq*(1.-yint/2)
          sy( 1) = twothird-oyintsq*(1.-oyint/2)
          sy( 2) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz( 0) = 1.
         case(1)
          sz( 0) = 1.-zint
          sz( 1) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end select        

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

subroutine depose_j_n_2dxz(cj,np,xp,zp,ux,uy,uz,gaminv,w,q,xmin,zmin,dt,dx,dz,nx,nz,nxguard,nzguard,nox,noz, &
                        l_particles_weight,l4symtry)
   implicit none
   integer(ISZ) :: np,nx,nz,nox,noz,nxguard,nzguard
   real(kind=8), dimension(-nxguard:nx+nxguard,0:0,-nzguard:nz+nzguard,3), intent(in out) :: cj
   real(kind=8), dimension(np) :: xp,zp,w,ux,uy,uz,gaminv
   real(kind=8) :: q,dt,dx,dz,xmin,zmin
   logical(ISZ) :: l_particles_weight,l4symtry

   real(kind=8) :: dxi,dzi,xint,zint, &
                   oxint,ozint,xintsq,zintsq,oxintsq,ozintsq
   real(kind=8) :: x,z,wq,invvol,vx,vy,vz
   real(kind=8) :: sx(-int(nox/2):int((nox+1)/2)), &
                   sz(-int(noz/2):int((noz+1)/2))
   real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.
   integer(ISZ) :: j,l,ip,jj,ll,ixmin, ixmax, izmin, izmax
   
      dxi = 1./dx
      dzi = 1./dz
      invvol = dxi*dzi

      ixmin = -int(nox/2)
      ixmax = int((nox+1)/2)
      izmin = -int(noz/2)
      izmax = int((noz+1)/2)

      do ip=1,np
      
        vx = ux(ip)*gaminv(ip)
        vy = uy(ip)*gaminv(ip)
        vz = uz(ip)*gaminv(ip)
        
        x = (xp(ip)-0.5*vx*dt-xmin)*dxi
        z = (zp(ip)-0.5*vz*dt-zmin)*dzi
        
        if (l4symtry) then
          x=abs(x)
        end if
        
        j=floor(x)
        l=floor(z)

        xint = x-j
        zint = z-l

        if (l_particles_weight) then
          wq=q*w(ip)*invvol
        else
          wq=q*invvol
        end if
      
        select case(nox)
         case(0)
          sx( 0) = 1.
         case(1)
          sx( 0) = 1.-xint
          sx( 1) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end select        

        select case(noz)
         case(0)
          sz( 0) = 1.
         case(1)
          sz( 0) = 1.-zint
          sz( 1) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end select        

         do ll = izmin, izmax
            do jj = ixmin, ixmax
              cj(j+jj  ,0,l+ll  ,1) = cj(j+jj  ,0,l+ll  ,1) + sx(jj)*sz(ll)*wq*vx*0.5
              cj(j+jj-1,0,l+ll  ,1) = cj(j+jj-1,0,l+ll  ,1) + sx(jj)*sz(ll)*wq*vx*0.5
              cj(j+jj  ,0,l+ll  ,2) = cj(j+jj  ,0,l+ll  ,2) + sx(jj)*sz(ll)*wq*vy
              cj(j+jj  ,0,l+ll  ,3) = cj(j+jj  ,0,l+ll  ,3) + sx(jj)*sz(ll)*wq*vz*0.5
              cj(j+jj  ,0,l+ll-1,3) = cj(j+jj  ,0,l+ll-1,3) + sx(jj)*sz(ll)*wq*vz*0.5
            end do
        end do

    end do

  return
end subroutine depose_j_n_2dxz

 subroutine getf3d_linear(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,exg,eyg,ezg)
   
 implicit none
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

 subroutine getf3d_n(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz, &
                     nxguard,nyguard,nzguard,nox,noy,noz,exg,eyg,ezg,l4symtry)
   
 implicit none
      integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard,nox,noy,noz
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      logical(ISZ) :: l4symtry
      integer(ISZ) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
                      ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint
      real(kind=8) :: xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
      real(kind=8), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(kind=8), DIMENSION(-int(noy/2):int((noy+1)/2)) :: sy
      real(kind=8), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox+1)/2)
      iymin = -int(noy/2)
      iymax =  int((noy+1)/2)
      izmin = -int(noz/2)
      izmax =  int((noz+1)/2)

      signx = 1.
      signy = 1.

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
          if (y<0.) then
            y = -y
            signy = -1.
          else
            signy = 1.
          end if
        end if

        j=floor(x)
        k=floor(y)
        l=floor(z)

        xint=x-j
        yint=y-k
        zint=z-l

        select case(nox)
         case(0)
          sx( 0) = 1.
         case(1)
          sx( 0) = 1.-xint
          sx( 1) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end select        

        select case(noy)
         case(0)
          sy( 0) = 1.
         case(1)
          sy( 0) = 1.-yint
          sy( 1) = yint
         case(2)
!          yint=yint-0.5
          yintsq = yint*yint
          sy(-1) = 0.5*(0.5-yint)**2
          sy( 0) = 0.75-yintsq
          sy( 1) = 0.5*(0.5+yint)**2
         case(3)
          oyint = 1.-yint
          yintsq = yint*yint
          oyintsq = oyint*oyint
          sy(-1) = onesixth*oyintsq*oyint
          sy( 0) = twothird-yintsq*(1.-yint/2)
          sy( 1) = twothird-oyintsq*(1.-oyint/2)
          sy( 2) = onesixth*yintsq*yint
        end select        

        select case(noz)
         case(0)
          sz( 0) = 1.
         case(1)
          sz( 0) = 1.-zint
          sz( 1) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end select        

        do ll = izmin, izmax
          do kk = iymin, iymax
            do jj = ixmin, ixmax
              ex(ip) = ex(ip) + sx(jj)*sy(kk)*sz(ll)*exg(j+jj,k+kk,l+ll)*signx
            end do
          end do
        end do

        do ll = izmin, izmax
          do kk = iymin, iymax
            do jj = ixmin, ixmax
              ey(ip) = ey(ip) + sx(jj)*sy(kk)*sz(ll)*eyg(j+jj,k+kk,l+ll)*signy
            end do
          end do
        end do

        do ll = izmin, izmax
          do kk = iymin, iymax
            do jj = ixmin, ixmax
              ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz(ll)*ezg(j+jj,k+kk,l+ll)
            end do
          end do
        end do

     end do

   return
 end subroutine getf3d_n

subroutine getf2dxz_n(np,xp,yp,zp,ex,ey,ez,xmin,zmin,dx,dz,nx,ny,nz, &
                     nxguard,nyguard,nzguard,nox,noz,exg,eyg,ezg,l4symtry,l_2drz)
   
 implicit none
      integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard,nox,noz
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(kind=8) :: xmin,zmin,dx,dz
      logical(ISZ) :: l4symtry,l_2drz
      integer(ISZ) :: ip, j, l, ixmin, ixmax, izmin, izmax, &
                      ixmin0, ixmax0, izmin0, izmax0, jj, ll
      real(kind=8) :: dxi, dzi, x, y, z, xint, zint, r, costheta, sintheta
      real(kind=8) :: xintsq,oxint,zintsq,ozint,oxintsq,ozintsq,signx
      real(kind=8), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(kind=8), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox+1)/2)
      izmin = -int(noz/2)
      izmax =  int((noz+1)/2)

      signx = 1.
      
      do ip=1,np

        if (l_2drz) then
          x = xp(ip)
          y = yp(ip)
          r=sqrt(x*x+y*y)
          if (r*dxi>1.e-20) then
            costheta=x/r
            sintheta=y/r
          else  
            costheta=1.
            sintheta=0.
          end if
          x = (r-xmin)*dxi
        else
          x = (xp(ip)-xmin)*dxi
        end if

        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
        end if

        j=floor(x)
        l=floor(z)

        xint=x-j
        zint=z-l

        select case(nox)
         case(0)
          sx( 0) = 1.
         case(1)
          sx( 0) = 1.-xint
          sx( 1) = xint
         case(2)
!          xint=xint-0.5
          xintsq = xint*xint
          sx(-1) = 0.5*(0.5-xint)**2
          sx( 0) = 0.75-xintsq
          sx( 1) = 0.5*(0.5+xint)**2
         case(3)
          oxint = 1.-xint
          xintsq = xint*xint
          oxintsq = oxint*oxint
          sx(-1) = onesixth*oxintsq*oxint
          sx( 0) = twothird-xintsq*(1.-xint/2)
          sx( 1) = twothird-oxintsq*(1.-oxint/2)
          sx( 2) = onesixth*xintsq*xint
        end select        

        select case(noz)
         case(0)
          sz( 0) = 1.
         case(1)
          sz( 0) = 1.-zint
          sz( 1) = zint
         case(2)
!          zint=zint-0.5
          zintsq = zint*zint
          sz(-1) = 0.5*(0.5-zint)**2
          sz( 0) = 0.75-zintsq
          sz( 1) = 0.5*(0.5+zint)**2
         case(3)
          ozint = 1.-zint
          zintsq = zint*zint
          ozintsq = ozint*ozint
          sz(-1) = onesixth*ozintsq*ozint
          sz( 0) = twothird-zintsq*(1.-zint/2)
          sz( 1) = twothird-ozintsq*(1.-ozint/2)
          sz( 2) = onesixth*zintsq*zint
        end select        

        if (l_2drz) then

          do ll = izmin, izmax
            do jj = ixmin, ixmax
              ex(ip) = ex(ip) + sx(jj)*sz(ll)*(exg(j+jj,0,l+ll)*costheta-eyg(j+jj,0,l+ll)*sintheta)
            end do
          end do

          do ll = izmin, izmax
            do jj = ixmin, ixmax
              ey(ip) = ey(ip) + sx(jj)*sz(ll)*(exg(j+jj,0,l+ll)*sintheta+eyg(j+jj,0,l+ll)*costheta)
            end do
          end do
        
        else

          do ll = izmin, izmax
            do jj = ixmin, ixmax
              ex(ip) = ex(ip) + sx(jj)*sz(ll)*exg(j+jj,0,l+ll)*signx
            end do
          end do

          do ll = izmin, izmax
            do jj = ixmin, ixmax
              ey(ip) = ey(ip) + sx(jj)*sz(ll)*eyg(j+jj,0,l+ll)
            end do
          end do

        end if

        do ll = izmin, izmax
          do jj = ixmin, ixmax
            ez(ip) = ez(ip) + sx(jj)*sz(ll)*ezg(j+jj,0,l+ll)
          end do
        end do

     end do

   return
 end subroutine getf2dxz_n

 subroutine gete3d_linear_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz, &
                                            nxguard,nyguard,nzguard,exg,eyg,ezg)
   
 implicit none
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

 subroutine geteb3d_linear_energy_conserving(np,xp,yp,zp,ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz, &
                                            nxguard,nyguard,nzguard,exg,eyg,ezg,bxg,byg,bzg)
   
 implicit none
      integer(ISZ) :: np,nx,ny,nz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez,bx,by,bz
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg,bxg,byg,bzg
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

        bx(ip) = bx(ip) + s1x*bxg(j  ,k,l) &
                        + s2x*bxg(j+1,k,l) 
                       
        by(ip) = by(ip) + s1y*byg(j,k  ,l) &
                        + s2y*byg(j,k+1,l) 
                       
        bz(ip) = bz(ip) + s1z*bzg(j,k,l  ) &
                        + s2z*bzg(j,k,l+1) 
                       
     end do

   return
 end subroutine geteb3d_linear_energy_conserving

  subroutine gete2dxz_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,zmin,dx,dz,nx,nz,nxguard,nzguard, &
                                       nox,noz,exg,eyg,ezg,l4symtry,l_2drz)
   
   implicit none
     integer(ISZ) :: np,nx,nz,nox,noz,nxguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      logical(ISZ) :: l4symtry,l_2drz
      real(kind=8), dimension(-nxguard:nx+nxguard,1,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(kind=8) :: xmin,zmin,dx,dz,costheta,sintheta
      integer(ISZ) :: ip, j, l, ixmin, ixmax, izmin, izmax, &
                      ixmin0, ixmax0, izmin0, izmax0, jj, ll
      real(kind=8) :: dxi, dzi, x, y, z, r, xint, zint, &
                      xintsq,oxint,zintsq,ozint,oxintsq,ozintsq,signx
      real(kind=8), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(kind=8), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(kind=8), DIMENSION(-int((nox)/2):int((nox-1)/2)) :: sx0
      real(kind=8), DIMENSION(-int((noz)/2):int((noz-1)/2)) :: sz0
      real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox-1)/2)
      izmin = -int(noz/2)
      izmax =  int((noz-1)/2)

      ixmin0 = -int((nox)/2)
      ixmax0 =  int((nox-1)/2)
      izmin0 = -int((noz)/2)
      izmax0 =  int((noz-1)/2)

      signx = 1.

      do ip=1,np

        if (l_2drz) then
          x = xp(ip)
          y = yp(ip)
          r=sqrt(x*x+y*y)
          if (r*dxi>1.e-20) then
            costheta=x/r
            sintheta=y/r
          else  
            costheta=1.
            sintheta=0.
          end if
          x = (r-xmin)*dxi
        else
          x = (xp(ip)-xmin)*dxi
        end if

        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
        end if
        
        j=floor(x)
        l=floor(z)

        xint=x-j
        zint=z-l

        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
!          xint=xint-0.5
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

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
!          zint=zint-0.5
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

        xint=xint-0.5
        if (nox==1) then
          sx0( 0) = 1.
        elseif (nox==2) then
!          xint=xint+0.5
          sx0(-1) = 1.-xint
          sx0( 0) = xint
        elseif (nox==3) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
        end if

        zint=zint-0.5
        if (noz==1) then
          sz0( 0) = 1.
        elseif (noz==2) then
!          zint=zint+0.5
          sz0(-1) = 1.-zint
          sz0( 0) = zint
        elseif (noz==3) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
        end if

        if (l_2drz) then
       
!          write(0,*) 'field gathering needs to be done for fstype=4 in EM-RZ'
!          stop
          do ll = izmin, izmax+1
            do jj = ixmin0, ixmax0
              ex(ip) = ex(ip) + sz(ll)*sx0(jj)*(exg(j+jj,1,l+ll)*costheta-eyg(j+jj,1,l+ll)*sintheta)
            end do
          end do

          do ll = izmin, izmax+1
            do jj = ixmin, ixmax+1
              ey(ip) = ey(ip) + sz(ll)*sx(jj)*(exg(j+jj,1,l+ll)*sintheta+eyg(j+jj,1,l+ll)*costheta)
            end do
          end do

        else

          do ll = izmin, izmax+1
            do jj = ixmin0, ixmax0
              ex(ip) = ex(ip) + sx0(jj)*sz(ll)*exg(j+jj,1,l+ll)*signx
            end do
          end do

          do ll = izmin, izmax+1
            do jj = ixmin, ixmax+1
              ey(ip) = ey(ip) + sx(jj)*sz(ll)*eyg(j+jj,1,l+ll)
            end do
          end do

        end if

          do ll = izmin0, izmax0
            do jj = ixmin, ixmax+1
              ez(ip) = ez(ip) + sx(jj)*sz0(ll)*ezg(j+jj,1,l+ll)
            end do
          end do
                     
     end do

   return
 end subroutine gete2dxz_n_energy_conserving

  subroutine gete3d_n_energy_conserving(np,xp,yp,zp,ex,ey,ez,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                       nox,noy,noz,exg,eyg,ezg,l4symtry)
   
   implicit none
     integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,ex,ey,ez
      logical(ISZ) :: l4symtry
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      integer(ISZ) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
                      ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
                      xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
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

      signx = 1.
      signy = 1.

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
          if (y<0.) then
            y = -y
            signy = -1.
          else
            signy = 1.
          end if
        end if
        
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
!          xint=xint-0.5
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
!          yint=yint-0.5
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
!          zint=zint-0.5
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

        xint=xint-0.5
        if (nox==1) then
          sx0( 0) = 1.
        elseif (nox==2) then
!          xint=xint+0.5
          sx0(-1) = 1.-xint
          sx0( 0) = xint
        elseif (nox==3) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
        end if

        yint=yint-0.5
        if (noy==1) then
          sy0( 0) = 1.
        elseif (noy==2) then
!          yint=yint+0.5
          sy0(-1) = 1.-yint
          sy0( 0) = yint
        elseif (noy==3) then
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
        end if

        zint=zint-0.5
        if (noz==1) then
          sz0( 0) = 1.
        elseif (noz==2) then
!          zint=zint+0.5
          sz0(-1) = 1.-zint
          sz0( 0) = zint
        elseif (noz==3) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
        end if

        do ll = izmin, izmax+1
          do kk = iymin, iymax+1
            do jj = ixmin0, ixmax0
              ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j+jj,k+kk,l+ll)*signx
            end do
          end do
        end do

        do ll = izmin, izmax+1
          do kk = iymin0, iymax0
            do jj = ixmin, ixmax+1
              ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj,k+kk,l+ll)*signy
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


 subroutine getb3d_linear_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz, &
                                            nxguard,nyguard,nzguard,bxg,byg,bzg)
   
 implicit none
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

subroutine getb2dxz_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,zmin,dx,dz,nx,nz,nxguard,nzguard, &
                                       nox,noz,bxg,byg,bzg,l4symtry,l_2drz)
   
      implicit none
      integer(ISZ) :: np,nx,nz,nox,noz,nxguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,bx,by,bz
      logical(ISZ) :: l4symtry,l_2drz
      real(kind=8), dimension(-nxguard:nx+nxguard,1,-nzguard:nz+nzguard) :: bxg,byg,bzg
      real(kind=8) :: xmin,zmin,dx,dz
      integer(ISZ) :: ip, j, l, ixmin, ixmax, izmin, izmax, &
                      ixmin0, ixmax0, izmin0, izmax0, jj, ll
      real(kind=8) :: dxi, dzi, x, y, z, xint, zint, &
                      xintsq,oxint,zintsq,ozint,oxintsq,ozintsq,signx, &
                      r, costheta, sintheta
      real(kind=8), DIMENSION(-int(nox/2):int((nox+1)/2)) :: sx
      real(kind=8), DIMENSION(-int(noz/2):int((noz+1)/2)) :: sz
      real(kind=8), DIMENSION(-int((nox)/2):int((nox-1)/2)) :: sx0
      real(kind=8), DIMENSION(-int((noz)/2):int((noz-1)/2)) :: sz0
      real(kind=8), parameter :: onesixth=1./6.,twothird=2./3.

      dxi = 1./dx
      dzi = 1./dz

      ixmin = -int(nox/2)
      ixmax =  int((nox-1)/2)
      izmin = -int(noz/2)
      izmax =  int((noz-1)/2)

      ixmin0 = -int((nox)/2)
      ixmax0 =  int((nox-1)/2)
      izmin0 = -int((noz)/2)
      izmax0 =  int((noz-1)/2)

      signx = 1.

      do ip=1,np

        if (l_2drz) then
          x = xp(ip)
          y = yp(ip)
          r=sqrt(x*x+y*y)
          if (r*dxi>1.e-20) then
            costheta=x/r
            sintheta=y/r
          else  
            costheta=1.
            sintheta=0.
          end if
          x = (r-xmin)*dxi
        else
          x = (xp(ip)-xmin)*dxi
        end if

        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
        end if
        
        j=floor(x)
        l=floor(z)

        xint=x-j
        zint=z-l

        if (nox==1) then
          sx( 0) = 1.-xint
          sx( 1) = xint
        elseif (nox==2) then
!          xint=xint-0.5
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

        if (noz==1) then
          sz( 0) = 1.-zint
          sz( 1) = zint
        elseif (noz==2) then
!          zint=zint-0.5
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

        xint=xint-0.5
        if (nox==1) then
          sx0( 0) = 1.
        elseif (nox==2) then
!          xint=xint+0.5
          sx0(-1) = 1.-xint
          sx0( 0) = xint
        elseif (nox==3) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
        end if

        zint=zint-0.5
        if (noz==1) then
          sz0( 0) = 1.
        elseif (noz==2) then
!          zint=zint+0.5
          sz0(-1) = 1.-zint
          sz0( 0) = zint
        elseif (noz==3) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
        end if

        if (l_2drz) then

          do ll = izmin0, izmax0
            do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sz0(ll)*sx(jj)*(bxg(j+jj,1,l+ll)*costheta-byg(j+jj,1,l+ll)*sintheta)
            end do
          end do

          do ll = izmin0, izmax0
            do jj = ixmin0, ixmax0
              by(ip) = by(ip) + sz0(ll)*sx0(jj)*(bxg(j+jj,1,l+ll)*sintheta+byg(j+jj,1,l+ll)*costheta)
            end do
          end do

        else

          do ll = izmin0, izmax0
            do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sx(jj)*sz0(ll)*bxg(j+jj,1,l+ll)*signx
            end do
          end do

          do ll = izmin0, izmax0
            do jj = ixmin0, ixmax0
              by(ip) = by(ip) + sx0(jj)*sz0(ll)*byg(j+jj,1,l+ll)
            end do
          end do

        end if

        do ll = izmin, izmax+1
            do jj = ixmin0, ixmax0
              bz(ip) = bz(ip) + sx0(jj)*sz(ll)*bzg(j+jj,1,l+ll)
            end do
        end do
                
     end do

   return
 end subroutine getb2dxz_n_energy_conserving

subroutine getb3d_n_energy_conserving(np,xp,yp,zp,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
                                       nox,noy,noz,bxg,byg,bzg,l4symtry)
   
      implicit none
      integer(ISZ) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
      real(kind=8), dimension(np) :: xp,yp,zp,bx,by,bz
      logical(ISZ) :: l4symtry
      real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
      real(kind=8) :: xmin,ymin,zmin,dx,dy,dz
      integer(ISZ) :: ip, j, k, l, ixmin, ixmax, iymin, iymax, izmin, izmax, &
                      ixmin0, ixmax0, iymin0, iymax0, izmin0, izmax0, jj, kk, ll
      real(kind=8) :: dxi, dyi, dzi, x, y, z, xint, yint, zint, &
                      xintsq,oxint,yintsq,oyint,zintsq,ozint,oxintsq,oyintsq,ozintsq,signx,signy
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

      signx = 1.
      signy = 1.

      do ip=1,np

        x = (xp(ip)-xmin)*dxi
        y = (yp(ip)-ymin)*dyi
        z = (zp(ip)-zmin)*dzi

        if (l4symtry) then
          if (x<0.) then
            x = -x
            signx = -1.
          else
            signx = 1.
          end if
          if (y<0.) then
            y = -y
            signy = -1.
          else
            signy = 1.
          end if
        end if
        
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
!          xint=xint-0.5
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
!          yint=yint-0.5
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
!          zint=zint-0.5
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

        xint=xint-0.5
        if (nox==1) then
          sx0( 0) = 1.
        elseif (nox==2) then
!          xint=xint+0.5
          sx0(-1) = 1.-xint
          sx0( 0) = xint
        elseif (nox==3) then
          xintsq = xint*xint
          sx0(-1) = 0.5*(0.5-xint)**2
          sx0( 0) = 0.75-xintsq
          sx0( 1) = 0.5*(0.5+xint)**2
        end if

        yint=yint-0.5
        if (noy==1) then
          sy0( 0) = 1.
        elseif (noy==2) then
!          yint=yint+0.5
          sy0(-1) = 1.-yint
          sy0( 0) = yint
        elseif (noy==3) then
          yintsq = yint*yint
          sy0(-1) = 0.5*(0.5-yint)**2
          sy0( 0) = 0.75-yintsq
          sy0( 1) = 0.5*(0.5+yint)**2
        end if

        zint=zint-0.5
        if (noz==1) then
          sz0( 0) = 1.
        elseif (noz==2) then
!          zint=zint+0.5
          sz0(-1) = 1.-zint
          sz0( 0) = zint
        elseif (noz==3) then
          zintsq = zint*zint
          sz0(-1) = 0.5*(0.5-zint)**2
          sz0( 0) = 0.75-zintsq
          sz0( 1) = 0.5*(0.5+zint)**2
        end if

        do ll = izmin0, izmax0
          do kk = iymin0, iymax0
            do jj = ixmin, ixmax+1
              bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj,k+kk,l+ll)*signx
            end do
          end do
        end do

        do ll = izmin0, izmax0
          do kk = iymin, iymax+1
            do jj = ixmin0, ixmax0
              by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j+jj,k+kk,l+ll)*signy
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

 subroutine project_jxjyjz(jfine,jcoarse,jcoarse_mother,nxf,nyf,nzf,nxc,nyc,nzc, &
                           nxguard,nyguard,nzguard,rapx,rapy,rapz,ixc,iyc,izc,l_2dxz,icycle,novercycle)
 ! Projection of J from one fine grid onto a coarse grid
 implicit none
 logical(ISZ) :: l_2dxz
 integer(ISZ) :: nxf,nyf,nzf,nxc,nyc,nzc,ixc,iyc,izc,rapx,rapy,rapz,nxguard,nyguard,nzguard,icycle,novercycle
 real(kind=8), DIMENSION(-nxguard:nxf+nxguard,-nyguard:nyf+nyguard,-nzguard:nzf+nzguard,3) :: jfine
 real(kind=8), DIMENSION(-nxguard:nxf/rapx+nxguard,-nyguard:nyf/rapy+nyguard,-nzguard:nzf/rapz+nzguard,3) :: jcoarse
 real(kind=8), DIMENSION(-nxguard:nxc+nxguard,-nyguard:nyc+nyguard,-nzguard:nzc+nzguard,3) :: jcoarse_mother

 INTEGER :: j, k, l, jg, kg, lg, ixmin, ixmax, iymin, iymax, izmin, izmax
 real(kind=8) :: wx, wy, wz, owx, owy, owz, invrapvol, irapx, irapy, irapz

   irapx = 1./rapx
   irapy = 1./rapy
   irapz = 1./rapz
   invrapvol = irapx*irapy*irapz/novercycle

   if(icycle==0) jcoarse(:,:,:,:) = 0.
   
   ixmin = -nxguard
   ixmax = nxf+nxguard
   iymin = -nyguard
   iymax = nyf+nyguard
   izmin = -nzguard
   izmax = nzf+nzguard

   if (.not.l_2dxz) then

   do l = izmin, izmax
      lg = floor(l*irapz)  
      wz = REAL(MOD(l+nzguard*rapz,rapz))*irapz
      owz= 1.-wz
      do k = iymin, iymax
         kg = floor(k*irapy)
         wy = REAL(MOD(k+nyguard*rapy,rapy))*irapy
         owy= 1.-wy
         do j = ixmin, ixmax-1
            jg = floor(j*irapx)
            jcoarse(jg,kg  ,lg  ,1) = jcoarse(jg,kg  ,lg  ,1) + owy*owz*jfine(j,k,l,1)*invrapvol
            if (kg<iymax) &
            jcoarse(jg,kg+1,lg  ,1) = jcoarse(jg,kg+1,lg  ,1) +  wy*owz*jfine(j,k,l,1)*invrapvol
            if (lg<izmax) &
            jcoarse(jg,kg  ,lg+1,1) = jcoarse(jg,kg  ,lg+1,1) + owy* wz*jfine(j,k,l,1)*invrapvol
            if (kg<iymax .and. lg<izmax) &
            jcoarse(jg,kg+1,lg+1,1) = jcoarse(jg,kg+1,lg+1,1) +  wy* wz*jfine(j,k,l,1)*invrapvol
         end do
      end do
   end do

   do l = izmin, izmax
      lg = floor(l*irapz)  
      wz = REAL(MOD(l+nzguard*rapz,rapz))*irapz
      owz= 1.-wz
      do k = iymin, iymax-1
         kg = floor(k*irapy)
         do j = ixmin, ixmax
            jg = floor(j*irapx)
            wx = REAL(MOD(j+nxguard*rapx,rapx))*irapx
            owx= 1.-wx
            jcoarse(jg  ,kg,lg  ,2) = jcoarse(jg  ,kg,lg  ,2) + owx*owz*jfine(j,k,l,2)*invrapvol
            if (jg<ixmax) &
            jcoarse(jg+1,kg,lg  ,2) = jcoarse(jg+1,kg,lg  ,2) +  wx*owz*jfine(j,k,l,2)*invrapvol
            if (lg<izmax) &
            jcoarse(jg  ,kg,lg+1,2) = jcoarse(jg  ,kg,lg+1,2) + owx* wz*jfine(j,k,l,2)*invrapvol
            if (jg<ixmax .and. lg<izmax) &
            jcoarse(jg+1,kg,lg+1,2) = jcoarse(jg+1,kg,lg+1,2) +  wx* wz*jfine(j,k,l,2)*invrapvol
         end do
      end do
   end do

   do l = izmin, izmax-1
      lg = floor(l*irapz)  
      do k = iymin, iymax
         kg = floor(k*irapy)
         wy = REAL(MOD(k+nyguard*rapy,rapy))*irapy
         owy= 1.-wy
         do j = ixmin, ixmax
            jg = floor(j*irapx)
            wx = REAL(MOD(j+nxguard*rapx,rapx))*irapx
            owx= 1.-wx
            jcoarse(jg  ,kg  ,lg,3) = jcoarse(jg  ,kg  ,lg,3) + owy*owx*jfine(j,k,l,3)*invrapvol
            if (kg<iymax) &
            jcoarse(jg  ,kg+1,lg,3) = jcoarse(jg  ,kg+1,lg,3) +  wy*owx*jfine(j,k,l,3)*invrapvol
            if (jg<ixmax) &
            jcoarse(jg+1,kg  ,lg,3) = jcoarse(jg+1,kg  ,lg,3) + owy* wx*jfine(j,k,l,3)*invrapvol
            if (jg<ixmax .and. kg<iymax) &
            jcoarse(jg+1,kg+1,lg,3) = jcoarse(jg+1,kg+1,lg,3) +  wy* wx*jfine(j,k,l,3)*invrapvol
         end do
      end do
   end do

   else
   
   k=0
   kg=0
   do l = izmin, izmax
      lg = floor(l*irapz)  
      wz = REAL(MOD(l+nzguard*rapz,rapz))*irapz
      owz= 1.-wz
         do j = ixmin, ixmax-1
            jg = floor(j*irapx)
            jcoarse(jg,kg  ,lg  ,1) = jcoarse(jg,kg  ,lg  ,1) + owz*jfine(j,k,l,1)*invrapvol
            if (lg<izmax) &
            jcoarse(jg,kg  ,lg+1,1) = jcoarse(jg,kg  ,lg+1,1) +  wz*jfine(j,k,l,1)*invrapvol
         end do
   end do

   do l = izmin, izmax
      lg = floor(l*irapz)  
      wz = REAL(MOD(l+nzguard*rapz,rapz))*irapz
      owz= 1.-wz
         do j = ixmin, ixmax
            jg = floor(j*irapx)
            wx = REAL(MOD(j+nxguard*rapx,rapx))*irapx
            owx= 1.-wx
            jcoarse(jg  ,kg,lg  ,2) = jcoarse(jg  ,kg,lg  ,2) + owx*owz*jfine(j,k,l,2)*invrapvol
            if (jg<ixmax) &
            jcoarse(jg+1,kg,lg  ,2) = jcoarse(jg+1,kg,lg  ,2) +  wx*owz*jfine(j,k,l,2)*invrapvol
            if (lg<izmax) &
            jcoarse(jg  ,kg,lg+1,2) = jcoarse(jg  ,kg,lg+1,2) + owx* wz*jfine(j,k,l,2)*invrapvol
            if (jg<ixmax .and. lg<izmax) &
            jcoarse(jg+1,kg,lg+1,2) = jcoarse(jg+1,kg,lg+1,2) +  wx* wz*jfine(j,k,l,2)*invrapvol
         end do
   end do

   do l = izmin, izmax-1
      lg = floor(l*irapz)  
         do j = ixmin, ixmax
            jg = floor(j*irapx)
            wx = REAL(MOD(j+nxguard*rapx,rapx))*irapx
            owx= 1.-wx
            jcoarse(jg  ,kg  ,lg,3) = jcoarse(jg  ,kg  ,lg,3) + owx*jfine(j,k,l,3)*invrapvol
            if (jg<ixmax) &
            jcoarse(jg+1,kg  ,lg,3) = jcoarse(jg+1,kg  ,lg,3) +  wx*jfine(j,k,l,3)*invrapvol
         end do
   end do

   endif
   
   jcoarse_mother(ixc-nxguard:ixc+nxf/rapx+nxguard,iyc-nyguard:iyc+nyf/rapy+nyguard,izc-nzguard:izc+nzf/rapz+nzguard,:) = &
   jcoarse_mother(ixc-nxguard:ixc+nxf/rapx+nxguard,iyc-nyguard:iyc+nyf/rapy+nyguard,izc-nzguard:izc+nzf/rapz+nzguard,:) + &
   jcoarse(:,:,:,:)

   return
 end subroutine project_jxjyjz

 subroutine project_rho(rhofine,rhocoarse,rhocoarse_mother,nxf,nyf,nzf,nxc,nyc,nzc,nxguard,nyguard,nzguard, &
                        rapx,rapy,rapz,ixc,iyc,izc,l_2dxz)
 ! Projection of J from one fine grid onto a coarse grid
 implicit none
 logical(ISZ) :: l_2dxz
 integer(ISZ) :: nxf,nyf,nzf,nxc,nyc,nzc,ixc,iyc,izc,rapx,rapy,rapz,nxguard,nyguard,nzguard
 real(kind=8), DIMENSION(-nxguard:nxf+nxguard,-nyguard:nyf+nyguard,-nzguard:nzf+nzguard) :: rhofine
 real(kind=8), DIMENSION(-nxguard:nxf/rapx+nxguard,-nyguard:nyf/rapy+nyguard,-nzguard:nzf/rapz+nzguard) :: rhocoarse
 real(kind=8), DIMENSION(-nxguard:nxc+nxguard,-nyguard:nyc+nyguard,-nzguard:nzc+nzguard) :: rhocoarse_mother

 INTEGER :: j, k, l, jg, kg, lg, ixmin, ixmax, iymin, iymax, izmin, izmax
 real(kind=8) :: wx, wy, wz, owx, owy, owz, invrapvol, irapx, irapy, irapz

   irapx = 1./rapx
   irapy = 1./rapy
   irapz = 1./rapz
   invrapvol = irapx*irapy*irapz
   
   ixmin = -nxguard
   ixmax = nxf+nxguard
   iymin = -nyguard
   iymax = nyf+nyguard
   izmin = -nzguard
   izmax = nzf+nzguard

   rhocoarse(:,:,:) = 0.

   if (.not.l_2dxz) then

   do l = -nzguard, nzf+nzguard
      lg = floor(l*irapz)  
      wz = REAL(MOD(l+nzguard*rapz,rapz))*irapz
      owz= 1.-wz
      do k = -nyguard, nyf+nyguard
         kg = floor(k*irapy)
         wy = REAL(MOD(k+nyguard*rapy,rapy))*irapy
         owy= 1.-wy
         do j = -nxguard, nxf+nxguard
            jg = floor(j*irapx)
            wx = REAL(MOD(j+nxguard*rapx,rapx))*irapx
            owx= 1.-wx
            rhocoarse(jg,kg  ,lg  ) = rhocoarse(jg,kg  ,lg  ) + owx*owy*owz*rhofine(j,k,l)*invrapvol
            if (kg<iymax) &
            rhocoarse(jg,kg+1,lg  ) = rhocoarse(jg,kg+1,lg  ) + owx* wy*owz*rhofine(j,k,l)*invrapvol
            if (lg<izmax) &
            rhocoarse(jg,kg  ,lg+1) = rhocoarse(jg,kg  ,lg+1) + owx*owy* wz*rhofine(j,k,l)*invrapvol
            if (lg<izmax .and. kg<iymax) &
            rhocoarse(jg,kg+1,lg+1) = rhocoarse(jg,kg+1,lg+1) + owx* wy* wz*rhofine(j,k,l)*invrapvol
            if (jg<ixmax) &
            rhocoarse(jg+1,kg  ,lg  ) = rhocoarse(jg+1,kg  ,lg  ) + wx*owy*owz*rhofine(j,k,l)*invrapvol
            if (jg<ixmax .and. kg<iymax) &
            rhocoarse(jg+1,kg+1,lg  ) = rhocoarse(jg+1,kg+1,lg  ) + wx* wy*owz*rhofine(j,k,l)*invrapvol
            if (jg<ixmax .and. lg<izmax) &
            rhocoarse(jg+1,kg  ,lg+1) = rhocoarse(jg+1,kg  ,lg+1) + wx*owy* wz*rhofine(j,k,l)*invrapvol
            if (jg<ixmax .and. kg<iymax .and. lg<izmax) &
            rhocoarse(jg+1,kg+1,lg+1) = rhocoarse(jg+1,kg+1,lg+1) + wx* wy* wz*rhofine(j,k,l)*invrapvol
         end do
      end do
   end do
   
   else
   
   k=0
   kg=0
   do l = -nzguard, nzf+nzguard
      lg = floor(l*irapz)  
      wz = REAL(MOD(l+nzguard*rapz,rapz))*irapz
      owz= 1.-wz
         do j = -nxguard, nxf+nxguard
            jg = floor(j*irapx)
            wx = REAL(MOD(j+nxguard*rapx,rapx))*irapx
            owx= 1.-wx
            rhocoarse(jg,kg  ,lg  ) = rhocoarse(jg,kg  ,lg  )     + owx*owz*rhofine(j,k,l)*invrapvol
            if (lg<izmax) &
            rhocoarse(jg,kg  ,lg+1) = rhocoarse(jg,kg  ,lg+1)     + owx* wz*rhofine(j,k,l)*invrapvol
            if (jg<ixmax) &
            rhocoarse(jg+1,kg  ,lg  ) = rhocoarse(jg+1,kg  ,lg  ) + wx *owz*rhofine(j,k,l)*invrapvol
            if (jg<ixmax .and. lg<izmax) &
            rhocoarse(jg+1,kg  ,lg+1) = rhocoarse(jg+1,kg  ,lg+1) + wx * wz*rhofine(j,k,l)*invrapvol
         end do
   end do

   endif
   
   rhocoarse_mother(ixc-nxguard:ixc+nxf/rapx+nxguard,iyc-nyguard:iyc+nyf/rapy+nyguard,izc-nzguard:izc+nzf/rapz+nzguard) = &
   rhocoarse_mother(ixc-nxguard:ixc+nxf/rapx+nxguard,iyc-nyguard:iyc+nyf/rapy+nyguard,izc-nzguard:izc+nzf/rapz+nzguard) + &
   rhocoarse(:,:,:)

   return
 end subroutine project_rho

subroutine apply_dmask(rho,jc,dmaskx,dmasky,dmaskz,bounds,nguarddepos,ntrans,nx,ny,nz,nxguard,nyguard,nzguard,l_pushf,l_2dxz)
 ! Projection of J from one fine grid onto a coarse grid
 use EM3D_FIELDobjects, only : otherproc
 implicit none
 logical(ISZ) :: l_2dxz, l_pushf
 integer(ISZ) :: nx,ny,nz,nxguard,nyguard,nzguard,bounds(10),nguarddepos(3),ntrans(3)
 real(kind=8), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: rho
 real(kind=8), DIMENSION(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: jc
 real(kind=8), DIMENSION(-nxguard:nx+nxguard) :: dmaskx
 real(kind=8), DIMENSION(-nyguard:ny+nyguard) :: dmasky
 real(kind=8), DIMENSION(-nzguard:nz+nzguard) :: dmaskz

 INTEGER :: j, k, l, i

if (.true.) then
 do j = 0, nx
   i = -ntrans(1)!/2
   if (j<(nguarddepos(1)-i) .and. bounds(1)/=otherproc) then
     if (j>(nguarddepos(1)-ntrans(1)-i)) then
       dmaskx(j) = real(j-nguarddepos(1)+ntrans(1)+i)/ntrans(1)
     end if
   else if (j>(nx-nguarddepos(1)+i) .and. bounds(2)/=otherproc) then
     if (j<nx-nguarddepos(1)+i+ntrans(1)) then
       dmaskx(j) = real(nx-nguarddepos(1)+i+ntrans(1)-j)/ntrans(1)
     end if
   else
     dmaskx(j) = 1.
   end if
 end do 
else
 do j = 0, nx
   if (j<nguarddepos(1) .and. bounds(1)/=otherproc) then
     if (j>nguarddepos(1)-ntrans(1)) then
       dmaskx(j) = real(j-nguarddepos(1)+ntrans(1))/ntrans(1)
     end if
   else if (j>(nx-nguarddepos(1)) .and. bounds(2)/=otherproc) then
     if (j<nx-nguarddepos(1)+ntrans(1)) then
       dmaskx(j) = real(nx-nguarddepos(1)+ntrans(1)-j)/ntrans(1)
     end if
   else
     dmaskx(j) = 1.
   end if
 end do 
endif

 if (.not.l_2dxz) then
 do k = 0, ny
   if (k<nguarddepos(2) .and. bounds(3)/=otherproc) then
     if (k>nguarddepos(2)-ntrans(2)) then
       dmasky(k) = real(k-nguarddepos(2)+ntrans(2))/ntrans(2)
     end if
   else if (k>(ny-nguarddepos(2)) .and. bounds(4)/=otherproc) then
     if (k<ny-nguarddepos(2)+ntrans(2)) then
       dmasky(k) = real(ny-nguarddepos(2)+ntrans(2)-j)/ntrans(2)
     end if
   else
     dmasky(k) = 1.
   end if
 end do 
 endif

 do l = 0, nz
   if (l<nguarddepos(3) .and. bounds(5)/=otherproc) then
     if (l>nguarddepos(3)-ntrans(3)) then
       dmaskz(l) = real(l-nguarddepos(3)+ntrans(3))/ntrans(3)
     end if
   else if (l>(nz-nguarddepos(3)) .and. bounds(6)/=otherproc) then
     if (l<nz-nguarddepos(3)+ntrans(3)) then
       dmaskz(l) = real(nz-nguarddepos(3)+ntrans(3)-l)/ntrans(3)
     end if
   else
     dmaskz(l) = 1.
   end if
 end do 

!dmaskx=1.
dmaskz=1.

 if (.not.l_2dxz) then

   do l = 0, nz
      do k = 0, ny
         do j = 0, nx
           jc(j,k,l,1) = jc(j,k,l,1) * 0.5*(dmaskx(j)+dmaskx(j+1)) * dmasky(k) * dmaskz(l)
       end do
      end do
   end do
   do l = 0, nz
      do k = 0, ny
         do j = 0, nx
           jc(j,k,l,2) = jc(j,k,l,2) * 0.5*(dmasky(k)+dmasky(k+1)) * dmaskx(j) * dmaskz(l)
       end do
      end do
   end do
   do l = 0, nz
      do k = 0, ny
         do j = 0, nx
           jc(j,k,l,3) = jc(j,k,l,3) * 0.5*(dmaskz(l)+dmaskz(l+1)) * dmaskx(j) * dmasky(k)
       end do
      end do
   end do
   if (l_pushf) then
     do l = 0, nz
        do k = 0, ny
           do j = 0, nx
             rho(j,k,l) = rho(j,k,l) * dmaskx(j) * dmasky(k) * dmaskz(l)
         end do
       end do
     end do
   end if

 else
   k = 0

   do l = 0, nz
         do j = 0, nx
           jc(j,k,l,1) = jc(j,k,l,1) * 0.5*(dmaskx(j)+dmaskx(j+1)) * dmaskz(l)
       end do
   end do
   do l = 0, nz
         do j = 0, nx
           jc(j,k,l,2) = jc(j,k,l,2) * dmaskx(j) * dmaskz(l)
       end do
   end do
   do l = 0, nz
         do j = 0, nx
           jc(j,k,l,3) = jc(j,k,l,3) * 0.5*(dmaskz(l)+dmaskz(l+1)) * dmaskx(j) 
       end do
   end do
   if (l_pushf) then
     do l = 0, nz
           do j = 0, nx
             rho(j,k,l) = rho(j,k,l) * dmaskx(j) * dmaskz(l)
         end do
     end do
   end if
 end if

 return
end subroutine apply_dmask

subroutine setebp(emfield,icycle,novercycle)

 use EM3D_YEEFIELDtypemodule
 implicit none

 TYPE(EM3D_YEEFIELDtype) :: emfield
 real(kind=8) :: w,ow
 integer(ISZ) :: icycle,novercycle
   
 if (novercycle==1) then
   emfield%exp = emfield%ex
   emfield%eyp = emfield%ey
   emfield%ezp = emfield%ez
   emfield%bxp = emfield%bx
   emfield%byp = emfield%by
   emfield%bzp = emfield%bz 
 else
   if (icycle==0) then
     emfield%expnext = emfield%ex
     emfield%eypnext = emfield%ey
     emfield%ezpnext = emfield%ez
     emfield%bxpnext = emfield%bx
     emfield%bypnext = emfield%by
     emfield%bzpnext = emfield%bz 
   end if
   w = 1./(novercycle-icycle)
   ow = 1.-w
   emfield%exp = ow*emfield%exp + w*emfield%expnext
   emfield%eyp = ow*emfield%eyp + w*emfield%eypnext
   emfield%ezp = ow*emfield%ezp + w*emfield%ezpnext
   emfield%bxp = ow*emfield%bxp + w*emfield%bxpnext
   emfield%byp = ow*emfield%byp + w*emfield%bypnext
   emfield%bzp = ow*emfield%bzp + w*emfield%bzpnext
 end if

 return
end subroutine setebp

subroutine addsubstractfields(child,child_coarse,parent,lc,ref,l_2dxz)
! Add own field and field from parent, substracting field from core_coarse, and 
! putting the result in Exp, Eyp, Ezp, Bxp, Byp and Bzp.

 use EM3D_BLOCKtypemodule
 use EM3D_YEEFIELDtypemodule
 implicit none

 TYPE(EM3D_BLOCKtype) :: child, child_coarse, parent
 integer(ISZ) :: lc(3),ref(3) ! lower bounds of child grid in parent grid, refinement
 logical(ISZ) :: l_2dxz

 TYPE(EM3D_YEEFIELDtype), pointer :: cf, cc, p
 INTEGER :: j, k, l, jg, kg, lg, jgp, kgp, lgp, rapx, rapy, rapz, nxf, nyf, nzf
 real(kind=8) :: wx, wy, wz, owx, owy, owz, irapx, irapy, irapz

   cf => child%core%yf
   cc => child_coarse%core%yf
   p  => parent%core%yf

   rapx = ref(1)
   rapy = ref(2)
   rapz = ref(3)
   irapx = 1./rapx
   irapy = 1./rapy
   irapz = 1./rapz
   
   nxf = cf%nx
   nyf = cf%ny
   nzf = cf%nz

   if (.not.l_2dxz) then
   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      wz = REAL(MOD(l,rapz))*irapz
      owz= 1.-wz
      do k = 0, nyf
         kg = k*irapy
         kgp = kg+lc(2)
         wy = REAL(MOD(k,rapy))*irapy
         owy= 1.-wy
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            cf%exp(j,k,l) = cf%exp(j,k,l) &
                          - owy * owz * cc%exp(jg   ,kg   ,lg   ) &
                          -  wy * owz * cc%exp(jg   ,kg+1 ,lg   ) &
                          - owy *  wz * cc%exp(jg   ,kg   ,lg+1 ) &
                          -  wy *  wz * cc%exp(jg   ,kg+1 ,lg+1 ) &
                          + owy * owz * p%exp(jgp  ,kgp  ,lgp  ) &
                          +  wy * owz * p%exp(jgp  ,kgp+1,lgp  ) &
                          + owy *  wz * p%exp(jgp  ,kgp  ,lgp+1) &
                          +  wy *  wz * p%exp(jgp  ,kgp+1,lgp+1) 
         end do
      end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      wz = REAL(MOD(l,rapz))*irapz
      owz= 1.-wz
      do k = 0, nyf
         kg = k*irapy
         kgp = kg+lc(2)
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            wx = REAL(MOD(j,rapx))*irapx
            owx= 1.-wx
            cf%eyp(j,k,l) = cf%eyp(j,k,l) &
                          - owx * owz * cc%eyp(jg   ,kg   ,lg   ) &
                          -  wx * owz * cc%eyp(jg+1 ,kg   ,lg   ) &
                          - owx *  wz * cc%eyp(jg   ,kg   ,lg+1 ) &
                          -  wx *  wz * cc%eyp(jg+1 ,kg   ,lg+1 ) &
                          + owx * owz * p%eyp(jgp  ,kgp  ,lgp  ) &
                          +  wx * owz * p%eyp(jgp+1,kgp  ,lgp  ) &
                          + owx *  wz * p%eyp(jgp  ,kgp  ,lgp+1) &
                          +  wx *  wz * p%eyp(jgp+1,kgp  ,lgp+1) 
         end do
      end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      do k = 0, nyf
         kg = k*irapy
         kgp = kg+lc(2)
         wy = REAL(MOD(k,rapy))*irapy
         owy= 1.-wy
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            wx = REAL(MOD(j,rapx))*irapx
            owx= 1.-wx
            cf%ezp(j,k,l) = cf%ezp(j,k,l) &
                          - owx * owy * cc%ezp(jg   ,kg   ,lg   ) &
                          -  wx * owy * cc%ezp(jg+1 ,kg   ,lg   ) &
                          - owx *  wy * cc%ezp(jg   ,kg+1 ,lg   ) &
                          -  wx *  wy * cc%ezp(jg+1 ,kg+1 ,lg   ) &
                          + owx * owy * p%ezp(jgp  ,kgp  ,lgp  ) &
                          +  wx * owy * p%ezp(jgp+1,kgp  ,lgp  ) &
                          + owx *  wy * p%ezp(jgp  ,kgp+1,lgp  ) &
                          +  wx *  wy * p%ezp(jgp+1,kgp+1,lgp  ) 
         end do
      end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      do k = 0, nyf
         kg = k*irapy
         kgp = kg+lc(2)
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            wx = REAL(MOD(j,rapx))*irapx
            owx= 1.-wx
            cf%bxp(j,k,l) = cf%bxp(j,k,l) &
                          - owx * cc%bxp(jg   ,kg   ,lg   ) &
                          -  wx * cc%bxp(jg+1 ,kg   ,lg   ) &
                          + owx * p%bxp(jgp  ,kgp  ,lgp  ) &
                          +  wx * p%bxp(jgp+1,kgp  ,lgp  ) 
         end do
      end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      do k = 0, nyf
         kg = k*irapy
         kgp = kg+lc(2)
         wy = REAL(MOD(k,rapy))*irapy
         owy= 1.-wy
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            cf%byp(j,k,l) = cf%byp(j,k,l) &
                          - owy * cc%byp(jg   ,kg   ,lg   ) &
                          -  wy * cc%byp(jg   ,kg+1 ,lg   ) &
                          + owy * p%byp(jgp  ,kgp  ,lgp  ) &
                          +  wy * p%byp(jgp  ,kgp+1,lgp  ) 
         end do
      end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      wz = REAL(MOD(l,rapz))*irapz
      owz= 1.-wz
      do k = 0, nyf
         kg = k*irapy
         kgp = kg+lc(2)
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            cf%bzp(j,k,l) = cf%bzp(j,k,l) &
                          - owz * cc%bzp(jg   ,kg   ,lg   ) &
                          -  wz * cc%bzp(jg   ,kg   ,lg+1 ) &
                          + owz * p%bzp(jgp  ,kgp  ,lgp  ) &
                          +  wz * p%bzp(jgp  ,kgp  ,lgp+1) 
         end do
      end do
   end do

   else

   k=0
   kg=0
   kgp=0
   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      wz = REAL(MOD(l,rapz))*irapz
      owz= 1.-wz
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            cf%exp(j,k,l) = cf%exp(j,k,l) &
                          - owz * cc%exp(jg   ,kg   ,lg   ) &
                          -  wz * cc%exp(jg   ,kg   ,lg+1 ) &
                          + owz * p%exp(jgp  ,kgp  ,lgp  ) &
                          +  wz * p%exp(jgp  ,kgp  ,lgp+1) 
         end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      wz = REAL(MOD(l,rapz))*irapz
      owz= 1.-wz
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            wx = REAL(MOD(j,rapx))*irapx
            owx= 1.-wx
            cf%eyp(j,k,l) = cf%eyp(j,k,l) &
                          - owx * owz * cc%eyp(jg   ,kg   ,lg   ) &
                          -  wx * owz * cc%eyp(jg+1 ,kg   ,lg   ) &
                          - owx *  wz * cc%eyp(jg   ,kg   ,lg+1 ) &
                          -  wx *  wz * cc%eyp(jg+1 ,kg   ,lg+1 ) &
                          + owx * owz * p%eyp(jgp  ,kgp  ,lgp  ) &
                          +  wx * owz * p%eyp(jgp+1,kgp  ,lgp  ) &
                          + owx *  wz * p%eyp(jgp  ,kgp  ,lgp+1) &
                          +  wx *  wz * p%eyp(jgp+1,kgp  ,lgp+1) 
         end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            wx = REAL(MOD(j,rapx))*irapx
            owx= 1.-wx
            cf%ezp(j,k,l) = cf%ezp(j,k,l) &
                          - owx * cc%ezp(jg   ,kg   ,lg   ) &
                          -  wx * cc%ezp(jg+1 ,kg   ,lg   ) &
                          + owx * p%ezp(jgp  ,kgp  ,lgp  ) &
                          +  wx * p%ezp(jgp+1,kgp  ,lgp  ) 
         end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            wx = REAL(MOD(j,rapx))*irapx
            owx= 1.-wx
            cf%bxp(j,k,l) = cf%bxp(j,k,l) &
                          - owx * cc%bxp(jg   ,kg   ,lg   ) &
                          -  wx * cc%bxp(jg+1 ,kg   ,lg   ) &
                          + owx * p%bxp(jgp  ,kgp  ,lgp  ) &
                          +  wx * p%bxp(jgp+1,kgp  ,lgp  ) 
         end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            cf%byp(j,k,l) = cf%byp(j,k,l) &
                          - cc%byp(jg   ,kg   ,lg   ) &
                          + p%byp(jgp  ,kgp  ,lgp  ) 
         end do
   end do

   do l = 0, nzf
      lg = l*irapz
      lgp = lg+lc(3)
      wz = REAL(MOD(l,rapz))*irapz
      owz= 1.-wz
         do j = 0, nxf-1
            jg = j*irapx
            jgp = jg+lc(1)
            cf%bzp(j,k,l) = cf%bzp(j,k,l) &
                          - owz * cc%bzp(jg   ,kg   ,lg   ) &
                          -  wz * cc%bzp(jg   ,kg   ,lg+1 ) &
                          + owz * p%bzp(jgp  ,kgp  ,lgp  ) &
                          +  wz * p%bzp(jgp  ,kgp  ,lgp+1) 
         end do
   end do
   endif

   return
 end subroutine addsubstractfields


subroutine addsubstractfields_nodal(child,child_coarse,parent,lc,ref,l_2dxz)
! Add own field and field from parent, substracting field from core_coarse, and 
! putting the result in Exp, Eyp, Ezp, Bxp, Byp and Bzp.

 use EM3D_BLOCKtypemodule
 use EM3D_YEEFIELDtypemodule
 implicit none

 TYPE(EM3D_BLOCKtype) :: child, child_coarse, parent
 integer(ISZ) :: lc(3),ref(3) ! lower bounds of child grid in parent grid, refinement
 logical(ISZ) :: l_2dxz

 TYPE(EM3D_YEEFIELDtype), pointer :: cf, cc, p
 INTEGER :: j, k, l, jg, kg, lg, jgp, kgp, lgp, rapx, rapy, rapz, nxf, nyf, nzf, incx, incy, incz
 real(kind=8) :: wx, wy, wz, owx, owy, owz, irapx, irapy, irapz

   cf => child%core%yf
   cc => child_coarse%core%yf
   p  => parent%core%yf

   rapx = ref(1)
   rapy = ref(2)
   rapz = ref(3)
   irapx = 1./rapx
   irapy = 1./rapy
   irapz = 1./rapz
   
   nxf = cf%nx
   nyf = cf%ny
   nzf = cf%nz

   if (.not.l_2dxz) then

   do l = -cf%nzguard, cf%nz+cf%nzguard
      lg = floor(l*irapz)  
      lgp = lg+lc(3)
      wz = REAL(MOD(l+cf%nzguard*rapz,rapz))*irapz
      owz= 1.-wz
      if (lg<cc%nz+cc%nzguard) then
        incz = 1
      else
        incz = 0
      end if
      do k = -cf%nyguard, cf%ny+cf%nyguard
         kg = floor(k*irapy)  
         wy = REAL(MOD(k+cf%nyguard*rapy,rapy))*irapy
         kgp = kg+lc(2)
         owy= 1.-wy
         if (kg<cc%ny+cc%nyguard) then
           incy = 1
         else
           incy = 0
         end if
         do j = -cf%nxguard, cf%nx+cf%nxguard
            jg = floor(j*irapx)  
            jgp = jg+lc(1)
            wx = REAL(MOD(j+cf%nxguard*rapx,rapx))*irapx
            owx= 1.-wx
            if (jg<cc%nx+cc%nxguard) then
              incx = 1
            else
              incx = 0
            end if
            cf%exp(j,k,l) = cf%exp(j,k,l) &
                          - owx * owy * owz * cc%exp(jg   ,kg   ,lg   ) &
                          - owx *  wy * owz * cc%exp(jg   ,kg+incy ,lg   ) &
                          - owx * owy *  wz * cc%exp(jg   ,kg   ,lg+incz ) &
                          - owx *  wy *  wz * cc%exp(jg   ,kg+incy ,lg+incz ) &
                          -  wx * owy * owz * cc%exp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wy * owz * cc%exp(jg+incx ,kg+incy ,lg   ) &
                          -  wx * owy *  wz * cc%exp(jg+incx ,kg   ,lg+incz ) &
                          -  wx *  wy *  wz * cc%exp(jg+incx ,kg+incy ,lg+incz ) &
                          + owx * owy * owz * p%exp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wy * owz * p%exp(jgp  ,kgp+incy,lgp  ) &
                          + owx * owy *  wz * p%exp(jgp  ,kgp  ,lgp+incz) &
                          + owx *  wy *  wz * p%exp(jgp  ,kgp+incy,lgp+incz) &
                          +  wx * owy * owz * p%exp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wy * owz * p%exp(jgp+incx,kgp+incy,lgp  ) &
                          +  wx * owy *  wz * p%exp(jgp+incx,kgp  ,lgp+incz) &
                          +  wx *  wy *  wz * p%exp(jgp+incx,kgp+incy,lgp+incz) 
            cf%eyp(j,k,l) = cf%eyp(j,k,l) &
                          - owx * owy * owz * cc%eyp(jg   ,kg   ,lg   ) &
                          - owx *  wy * owz * cc%eyp(jg   ,kg+incy ,lg   ) &
                          - owx * owy *  wz * cc%eyp(jg   ,kg   ,lg+incz ) &
                          - owx *  wy *  wz * cc%eyp(jg   ,kg+incy ,lg+incz ) &
                          -  wx * owy * owz * cc%eyp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wy * owz * cc%eyp(jg+incx ,kg+incy ,lg   ) &
                          -  wx * owy *  wz * cc%eyp(jg+incx ,kg   ,lg+incz ) &
                          -  wx *  wy *  wz * cc%eyp(jg+incx ,kg+incy ,lg+incz ) &
                          + owx * owy * owz * p%eyp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wy * owz * p%eyp(jgp  ,kgp+incy,lgp  ) &
                          + owx * owy *  wz * p%eyp(jgp  ,kgp  ,lgp+incz) &
                          + owx *  wy *  wz * p%eyp(jgp  ,kgp+incy,lgp+incz) &
                          +  wx * owy * owz * p%eyp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wy * owz * p%eyp(jgp+incx,kgp+incy,lgp  ) &
                          +  wx * owy *  wz * p%eyp(jgp+incx,kgp  ,lgp+incz) &
                          +  wx *  wy *  wz * p%eyp(jgp+incx,kgp+incy,lgp+incz) 
            cf%ezp(j,k,l) = cf%ezp(j,k,l) &
                          - owx * owy * owz * cc%ezp(jg   ,kg   ,lg   ) &
                          - owx *  wy * owz * cc%ezp(jg   ,kg+incy ,lg   ) &
                          - owx * owy *  wz * cc%ezp(jg   ,kg   ,lg+incz ) &
                          - owx *  wy *  wz * cc%ezp(jg   ,kg+incy ,lg+incz ) &
                          -  wx * owy * owz * cc%ezp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wy * owz * cc%ezp(jg+incx ,kg+incy ,lg   ) &
                          -  wx * owy *  wz * cc%ezp(jg+incx ,kg   ,lg+incz ) &
                          -  wx *  wy *  wz * cc%ezp(jg+incx ,kg+incy ,lg+incz ) &
                          + owx * owy * owz * p%ezp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wy * owz * p%ezp(jgp  ,kgp+incy,lgp  ) &
                          + owx * owy *  wz * p%ezp(jgp  ,kgp  ,lgp+incz) &
                          + owx *  wy *  wz * p%ezp(jgp  ,kgp+incy,lgp+incz) &
                          +  wx * owy * owz * p%ezp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wy * owz * p%ezp(jgp+incx,kgp+incy,lgp  ) &
                          +  wx * owy *  wz * p%ezp(jgp+incx,kgp  ,lgp+incz) &
                          +  wx *  wy *  wz * p%ezp(jgp+incx,kgp+incy,lgp+incz) 
            cf%bxp(j,k,l) = cf%bxp(j,k,l) &
                          - owx * owy * owz * cc%bxp(jg   ,kg   ,lg   ) &
                          - owx *  wy * owz * cc%bxp(jg   ,kg+incy ,lg   ) &
                          - owx * owy *  wz * cc%bxp(jg   ,kg   ,lg+incz ) &
                          - owx *  wy *  wz * cc%bxp(jg   ,kg+incy ,lg+incz ) &
                          -  wx * owy * owz * cc%bxp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wy * owz * cc%bxp(jg+incx ,kg+incy ,lg   ) &
                          -  wx * owy *  wz * cc%bxp(jg+incx ,kg   ,lg+incz ) &
                          -  wx *  wy *  wz * cc%bxp(jg+incx ,kg+incy ,lg+incz ) &
                          + owx * owy * owz * p%bxp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wy * owz * p%bxp(jgp  ,kgp+incy,lgp  ) &
                          + owx * owy *  wz * p%bxp(jgp  ,kgp  ,lgp+incz) &
                          + owx *  wy *  wz * p%bxp(jgp  ,kgp+incy,lgp+incz) &
                          +  wx * owy * owz * p%bxp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wy * owz * p%bxp(jgp+incx,kgp+incy,lgp  ) &
                          +  wx * owy *  wz * p%bxp(jgp+incx,kgp  ,lgp+incz) &
                          +  wx *  wy *  wz * p%bxp(jgp+incx,kgp+incy,lgp+incz) 
            cf%byp(j,k,l) = cf%byp(j,k,l) &
                          - owx * owy * owz * cc%byp(jg   ,kg   ,lg   ) &
                          - owx *  wy * owz * cc%byp(jg   ,kg+incy ,lg   ) &
                          - owx * owy *  wz * cc%byp(jg   ,kg   ,lg+incz ) &
                          - owx *  wy *  wz * cc%byp(jg   ,kg+incy ,lg+incz ) &
                          -  wx * owy * owz * cc%byp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wy * owz * cc%byp(jg+incx ,kg+incy ,lg   ) &
                          -  wx * owy *  wz * cc%byp(jg+incx ,kg   ,lg+incz ) &
                          -  wx *  wy *  wz * cc%byp(jg+incx ,kg+incy ,lg+incz ) &
                          + owx * owy * owz * p%byp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wy * owz * p%byp(jgp  ,kgp+incy,lgp  ) &
                          + owx * owy *  wz * p%byp(jgp  ,kgp  ,lgp+incz) &
                          + owx *  wy *  wz * p%byp(jgp  ,kgp+incy,lgp+incz) &
                          +  wx * owy * owz * p%byp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wy * owz * p%byp(jgp+incx,kgp+incy,lgp  ) &
                          +  wx * owy *  wz * p%byp(jgp+incx,kgp  ,lgp+incz) &
                          +  wx *  wy *  wz * p%byp(jgp+incx,kgp+incy,lgp+incz) 
            cf%bzp(j,k,l) = cf%bzp(j,k,l) &
                          - owx * owy * owz * cc%bzp(jg   ,kg   ,lg   ) &
                          - owx *  wy * owz * cc%bzp(jg   ,kg+incy ,lg   ) &
                          - owx * owy *  wz * cc%bzp(jg   ,kg   ,lg+incz ) &
                          - owx *  wy *  wz * cc%bzp(jg   ,kg+incy ,lg+incz ) &
                          -  wx * owy * owz * cc%bzp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wy * owz * cc%bzp(jg+incx ,kg+incy ,lg   ) &
                          -  wx * owy *  wz * cc%bzp(jg+incx ,kg   ,lg+incz ) &
                          -  wx *  wy *  wz * cc%bzp(jg+incx ,kg+incy ,lg+incz ) &
                          + owx * owy * owz * p%bzp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wy * owz * p%bzp(jgp  ,kgp+incy,lgp  ) &
                          + owx * owy *  wz * p%bzp(jgp  ,kgp  ,lgp+incz) &
                          + owx *  wy *  wz * p%bzp(jgp  ,kgp+incy,lgp+incz) &
                          +  wx * owy * owz * p%bzp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wy * owz * p%bzp(jgp+incx,kgp+incy,lgp  ) &
                          +  wx * owy *  wz * p%bzp(jgp+incx,kgp  ,lgp+incz) &
                          +  wx *  wy *  wz * p%bzp(jgp+incx,kgp+incy,lgp+incz) 
         end do
      end do
   end do

   else
   k=0
   kg=0
   kgp=0
   do l = -cf%nzguard, cf%nz+cf%nzguard
      lg = floor(l*irapz)  
      lgp = lg+lc(3)
      wz = REAL(MOD(l+cf%nzguard*rapz,rapz))*irapz
      owz= 1.-wz
      if (lg<cc%nz+cc%nzguard) then
        incz = 1
      else
        incz = 0
      end if
         do j = -cf%nxguard, cf%nx+cf%nxguard
            jg = floor(j*irapx)  
            jgp = jg+lc(1)
            wx = REAL(MOD(j+cf%nxguard*rapx,rapx))*irapx
            owx= 1.-wx
            if (jg<cc%nx+cc%nxguard) then
              incx = 1
            else
              incx = 0
            end if
            cf%exp(j,k,l) = cf%exp(j,k,l) &
                          - owx * owz * cc%exp(jg   ,kg   ,lg   ) &
                          - owx *  wz * cc%exp(jg   ,kg   ,lg+incz ) &
                          -  wx * owz * cc%exp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wz * cc%exp(jg+incx ,kg   ,lg+incz ) &
                          + owx * owz * p%exp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wz * p%exp(jgp  ,kgp  ,lgp+incz) &
                          +  wx * owz * p%exp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wz * p%exp(jgp+incx,kgp  ,lgp+incz) 
            cf%eyp(j,k,l) = cf%eyp(j,k,l) &
                          - owx * owz * cc%eyp(jg   ,kg   ,lg   ) &
                          - owx *  wz * cc%eyp(jg   ,kg   ,lg+incz ) &
                          -  wx * owz * cc%eyp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wz * cc%eyp(jg+incx ,kg   ,lg+incz ) &
                          + owx * owz * p%eyp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wz * p%eyp(jgp  ,kgp  ,lgp+incz) &
                          +  wx * owz * p%eyp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wz * p%eyp(jgp+incx,kgp  ,lgp+incz) 
            cf%ezp(j,k,l) = cf%ezp(j,k,l) &
                          - owx * owz * cc%ezp(jg   ,kg   ,lg   ) &
                          - owx *  wz * cc%ezp(jg   ,kg   ,lg+incz ) &
                          -  wx * owz * cc%ezp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wz * cc%ezp(jg+incx ,kg   ,lg+incz ) &
                          + owx * owz * p%ezp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wz * p%ezp(jgp  ,kgp  ,lgp+incz) &
                          +  wx * owz * p%ezp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wz * p%ezp(jgp+incx,kgp  ,lgp+incz) 
            cf%bxp(j,k,l) = cf%bxp(j,k,l) &
                          - owx * owz * cc%bxp(jg   ,kg   ,lg   ) &
                          - owx *  wz * cc%bxp(jg   ,kg   ,lg+incz ) &
                          -  wx * owz * cc%bxp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wz * cc%bxp(jg+incx ,kg   ,lg+incz ) &
                          + owx * owz * p%bxp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wz * p%bxp(jgp  ,kgp  ,lgp+incz) &
                          +  wx * owz * p%bxp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wz * p%bxp(jgp+incx,kgp  ,lgp+incz) 
            cf%byp(j,k,l) = cf%byp(j,k,l) &
                          - owx * owz * cc%byp(jg   ,kg   ,lg   ) &
                          - owx *  wz * cc%byp(jg   ,kg   ,lg+incz ) &
                          -  wx * owz * cc%byp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wz * cc%byp(jg+incx ,kg   ,lg+incz ) &
                          + owx * owz * p%byp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wz * p%byp(jgp  ,kgp  ,lgp+incz) &
                          +  wx * owz * p%byp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wz * p%byp(jgp+incx,kgp  ,lgp+incz) 
            cf%bzp(j,k,l) = cf%bzp(j,k,l) &
                          - owx * owz * cc%bzp(jg   ,kg   ,lg   ) &
                          - owx *  wz * cc%bzp(jg   ,kg   ,lg+incz ) &
                          -  wx * owz * cc%bzp(jg+incx ,kg   ,lg   ) &
                          -  wx *  wz * cc%bzp(jg+incx ,kg   ,lg+incz ) &
                          + owx * owz * p%bzp(jgp  ,kgp  ,lgp  ) &
                          + owx *  wz * p%bzp(jgp  ,kgp  ,lgp+incz) &
                          +  wx * owz * p%bzp(jgp+incx,kgp  ,lgp  ) &
                          +  wx *  wz * p%bzp(jgp+incx,kgp  ,lgp+incz) 
         end do
   end do

   endif
    
   return
 end subroutine addsubstractfields_nodal

