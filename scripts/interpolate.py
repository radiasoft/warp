from warp import *
interepolate_version = "$Id: interpolate.py,v 1.2 2001/07/09 20:39:40 dave Exp $"

def interpolate(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,
                serial=1):
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  xx = take(top.xp,ii)
  yy = take(top.yp,ii)
  zz = take(top.zp,ii)
  gg = take(top.gaminv,ii)
  vx = take(top.uxp,ii)*gg
  vy = take(top.uyp,ii)*gg
  vz = take(top.uzp,ii)*gg
  if zl and zu:
    zcent = 0.5*(zl + zu)
  elif iz:
    if npes > 0: zcent = top.zmslmin[0] + iz*w3d.dz
    else: zcent = w3d.zmmin + iz*w3d.dz
  else:
    zcent = top.zbeam + 0.5*(top.zwindows[0,iw] + top.zwindows[1,iw])
  delt = (zcent - zz)/vz
  xx = xx + vx*delt
  yy = yy + vy*delt
  zz = zz + vz*delt
  return (xx,yy,zz,vx,vy,vz,gg)


# --- Get a value at a specified location using linear interpolation
def interpolate3d(x,y,z):
  ix = int((x - w3d.xmmin)/w3d.dx)
  iy = int((y - w3d.ymmin)/w3d.dy)
  iz = int((z - w3d.zmmin)/w3d.dz)
  wx = (x - w3d.xmmin)/w3d.dx - ix
  wy = (y - w3d.ymmin)/w3d.dy - iy
  wz = (z - w3d.zmmin)/w3d.dz - iz
  return w3d.phi[ix  ,iy  ,iz+1  ]*(1.-wx)*(1.-wy)*(1.-wz) + \
         w3d.phi[ix+1,iy  ,iz+1  ]*    wx *(1.-wy)*(1.-wz) + \
         w3d.phi[ix  ,iy+1,iz+1  ]*(1.-wx)*    wy *(1.-wz) + \
         w3d.phi[ix+1,iy+1,iz+1  ]*    wx *    wy *(1.-wz) + \
         w3d.phi[ix  ,iy  ,iz+1+1]*(1.-wx)*(1.-wy)*    wz  + \
         w3d.phi[ix+1,iy  ,iz+1+1]*    wx *(1.-wy)*    wz  + \
         w3d.phi[ix  ,iy+1,iz+1+1]*(1.-wx)*    wy *    wz  + \
         w3d.phi[ix+1,iy+1,iz+1+1]*    wx *    wy *    wz

