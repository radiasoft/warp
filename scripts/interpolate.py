from warp import *
interepolate_version = "$Id: interpolate.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

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

