from warp import *
from grid_1d import *
from grid_2d import *
emittance_version = "$Id: emittance.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

def emit1(threshold=0.05,js=0):
  (gg,ggmesh) = gather_1d(top.xp[top.ins[js]-1:top.ins[js]+top.nps[js]-1])
  gg = gg/max(gg)
  dd = scatter_1d(gg,ggmesh,top.xp[top.ins[js]-1:top.ins[js]+top.nps[js]-1])

  xx = compress(greater(dd,threshold),
                top.xp[top.ins[js]-1:top.ins[js]+top.nps[js]-1])
  xp = compress(greater(dd,threshold),
                top.uxp[top.ins[js]-1:top.ins[js]+top.nps[js]-1]/
                top.uzp[top.ins[js]-1:top.ins[js]+top.nps[js]-1])

  (gg,ggmesh) = gather_1d(xp)
  gg = gg/max(gg)
  dd = scatter_1d(gg,ggmesh,xp)

  xx = compress(greater(dd,threshold),xx)
  xp = compress(greater(dd,threshold),xp)

  txe = ((ave(xx**2) - ave(xx)**2)*(ave(xp**2) - ave(xp)**2) -
         (ave(xx*xp) - ave(xx)*ave(xp))**2)
  epsx = 4.*sqrt(txe)

  (gg,ggmesh) = gather_1d(top.yp[top.ins[js]-1:top.ins[js]+top.nps[js]-1])
  gg = gg/max(gg)
  dd = scatter_1d(gg,ggmesh,top.yp[top.ins[js]-1:top.ins[js]+top.nps[js]-1])

  yy = compress(greater(dd,threshold),
                top.yp[top.ins[js]-1:top.ins[js]+top.nps[js]-1])
  yp = compress(greater(dd,threshold),
                top.uyp[top.ins[js]-1:top.ins[js]+top.nps[js]-1]/
                top.uzp[top.ins[js]-1:top.ins[js]+top.nps[js]-1])

  (gg,ggmesh) = gather_1d(yp)
  gg = gg/max(gg)
  dd = scatter_1d(gg,ggmesh,yp)

  yy = compress(greater(dd,threshold),yy)
  yp = compress(greater(dd,threshold),yp)

  tye = ((ave(yy**2) - ave(yy)**2)*(ave(yp**2) - ave(yp)**2) -
         (ave(yy*yp) - ave(yy)*ave(yp))**2)
  epsy = 4.*sqrt(tye)

  return (epsx,epsy)

def emit2(threshold=0.05,js=0,iw=0,ngridx=20,ngridy=20):
  """
Calculates the emittance with thesholding. Particles in regions where
the density is less than threshold are not included. The thresholding is
done in a way which mimics experiment - the particles are binned onto a
grid which has the slope subtracted (to minimize empty space).
Input:
  - threshold=0.05 is threshold value
  - js=0 species of particles to include
  - iw=0 z-window to select particles
  - ngridx=20 number of grid points position is binned into
  - ngridy=20 number of grid points velocity is binned into
Output:
  - a tuple containing epsx and epsy
  """
  # --- Select out live particles within the z window.
  i1 = top.ins[js]-1
  i2 = top.ins[js]+top.nps[js]-1
  zw0 = top.zwindows[0,iw]+top.zbeam
  zw1 = top.zwindows[1,iw]+top.zbeam
  ii = compress(logical_and(greater(top.uzp[i1:i2],0.),
                 logical_and(less(zw0,top.zp[i1:i2]),less(top.zp[i1:i2],zw1))),
                arange(i1,i2))
  xx = take(top.xp,ii)
  yy = take(top.yp,ii)
  zz = take(top.zp,ii)
  vx = take(top.uxp,ii)
  vy = take(top.uyp,ii)
  vz = take(top.uzp,ii)

  # --- Set grid ranges
  wmin = min(xx)
  wmax = max(xx)
  slope = (top.xxpbar[iw]-top.xbar[iw]*top.xpbar[iw])/top.xrms[iw]**2
  vxms = vx/vz - xx*slope
  hmin = min(vxms)
  hmax = max(vxms)

  # --- Bin up the data onto a 2-D grid
  xxpmesh = fzeros((ngridx+1,ngridy+1),'d')
  getpsgrd(top.nps[js],xx,vx,ngridx,ngridy,xxpmesh,wmin,wmax,hmin,hmax,
           top.zwindows[0,iw]+top.zbeam,top.zwindows[1,iw]+top.zbeam,
           zz,vz,slope)
  xxpmesh = xxpmesh/maxnd(xxpmesh)
  # --- Scatter the resulting normalized density back to the particle
  # --- locations.
  dd = scatter_2d(xxpmesh,wmin,hmin,wmax,hmax,xx,vxms)

  # --- Select particles that are in the regions where the threshold is met.
  xx = compress(greater(dd,threshold),xx)
  xp = compress(greater(dd,threshold),vx/vz)

  # --- Calculate emittance from the remaining particles.
  txe = ((ave(xx**2) - ave(xx)**2)*(ave(xp**2) - ave(xp)**2) -
         (ave(xx*xp) - ave(xx)*ave(xp))**2)
  epsx = 4.*sqrt(txe)

  # --- Repeat for y
  wmin = min(yy)
  wmax = max(yy)
  slope = (top.yypbar[iw]-top.ybar[iw]*top.ypbar[iw])/top.yrms[iw]**2
  vyms = vy/vz - yy*slope
  hmin = min(vyms)
  hmax = max(vyms)
  yypmesh = fzeros((ngridx+1,ngridy+1),'d')
  getpsgrd(top.nps[js],yy,vy,ngridx,ngridy,yypmesh,wmin,wmax,hmin,hmax,
           top.zwindows[0,iw]+top.zbeam,top.zwindows[1,iw]+top.zbeam,
           zz,vz,slope)
  yypmesh = yypmesh/maxnd(yypmesh)
  dd = scatter_2d(yypmesh,wmin,hmin,wmax,hmax,yy,vyms)
  yy = compress(greater(dd,threshold),yy)
  yp = compress(greater(dd,threshold),vy/vz)
  tye = ((ave(yy**2) - ave(yy)**2)*(ave(yp**2) - ave(yp)**2) -
         (ave(yy*yp) - ave(yy)*ave(yp))**2)
  epsy = 4.*sqrt(tye)

  # --- return the resulting emittances.
  return (epsx,epsy)



def emitn2(threshold=0.05,js=0,iw=0,ngridx=20,ngridy=20):
  """
Calculates the normalized emittance with thesholding. Particles in
regions where the density is less than threshold are not included. The
thresholding is done in a way which mimics experiment - the particles
are binned onto a grid which has the slope subtracted (to minimize empty
space).
Input:
  - threshold=0.05 is threshold value
  - js=0 species of particles to include
  - iw=0 z-window to select particles
  - ngridx=20 number of grid points position is binned into
  - ngridy=20 number of grid points velocity is binned into
Output:
  - a tuple containing epsnx and epsny
  """
  # --- Select out live particles within the z window.
  i1 = top.ins[js]-1
  i2 = top.ins[js]+top.nps[js]-1
  zw0 = top.zwindows[0,iw]+top.zbeam
  zw1 = top.zwindows[1,iw]+top.zbeam
  ii = compress(logical_and(greater(top.uzp[i1:i2],0.),
                 logical_and(less(zw0,top.zp[i1:i2]),less(top.zp[i1:i2],zw1))),
                arange(i1,i2))
  xx = take(top.xp,ii)
  yy = take(top.yp,ii)
  zz = take(top.zp,ii)
  vx = take(top.uxp,ii)
  vy = take(top.uyp,ii)
  vz = ones(len(ii),'d')

  # --- Set grid ranges
  wmin = min(xx)
  wmax = max(xx)
  slope = ((top.xxpbar[iw]-top.xbar[iw]*top.xpbar[iw])*top.vzbar[iw]/
           top.xrms[iw]**2)
  vxms = vx - xx*slope
  hmin = min(vxms)
  hmax = max(vxms)

  # --- Bin up the data onto a 2-D grid
  xvxmesh = zeros((ngridx+1,ngridy+1),'d')
  getpsgrd(top.nps[js],xx,vx,ngridx,ngridy,xvxmesh,wmin,wmax,hmin,hmax,
           top.zwindows[0,iw],top.zwindows[1,iw],zz,vz,slope)
  xvxmesh = xvxmesh/maxnd(xvxmesh)
  # --- Scatter the resulting normalized density back to the particle
  # --- locations.
  dd = scatter_2d(xvxmesh,wmin,hmin,wmax,hmax,xx,vxms)

  # --- Select particles that are in the regions where the threshold is met.
  xx = compress(greater(dd,threshold),xx)
  vx = compress(greater(dd,threshold),vx)

  # --- Calculate emittance from the remaining particles.
  txe = ((ave(xx**2) - ave(xx)**2)*(ave(vx**2) - ave(vx)**2) -
         (ave(xx*vx) - ave(xx)*ave(vx))**2)
  epsnx = 4.*sqrt(txe)/top.clight

  # --- Repeat for y
  wmin = min(yy)
  wmax = max(yy)
  slope = ((top.yypbar[iw]-top.ybar[iw]*top.ypbar[iw])*top.vzbar[iw]/
           top.yrms[iw]**2)
  vyms = vy - yy*slope
  hmin = min(vyms)
  hmax = max(vyms)
  yvxmesh = zeros((ngridx+1,ngridy+1),'d')
  getpsgrd(top.nps[js],yy,vy,ngridx,ngridy,yvxmesh,wmin,wmax,hmin,hmax,
           top.zwindows[0,iw],top.zwindows[1,iw],zz,vz,slope)
  yvxmesh = yvxmesh/maxnd(yvxmesh)
  dd = scatter_2d(yvxmesh,wmin,hmin,wmax,hmax,yy,vyms)
  yy = compress(greater(dd,threshold),yy)
  vy = compress(greater(dd,threshold),vy)
  tye = ((ave(yy**2) - ave(yy)**2)*(ave(vy**2) - ave(vy)**2) -
         (ave(yy*vy) - ave(yy)*ave(vy))**2)
  epsny = 4.*sqrt(tye)/top.clight

  # --- return the resulting emittances.
  return (epsnx,epsny)








