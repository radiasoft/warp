from warp import *
sortwxy_version = "$Id: sortwxy.py,v 1.5 2009/02/05 18:18:40 dave Exp $"
# Sorts the particles for better cache use for the gather/scatter.

def sortwxy(pgroup=None):
  if pgroup is None: pgroup = top.pgroup
  for j in range(pgroup.ns):
    xx = getx(js=j,gather=0)
    yy = gety(js=j,gather=0)
    zz = getz(js=j,gather=0)
    ux = getux(js=j,gather=0)
    uy = getuy(js=j,gather=0)
    uz = getuz(js=j,gather=0)
    gi = getgaminv(js=j,gather=0)
    ix = (abs(xx - w3d.xmmin)/w3d.dx).astype(long)
    iy = (abs(yy - w3d.ymmin)/w3d.dy).astype(long)
    ixy = ix + iy*(w3d.nx+1)
    ii = argsort(ixy)
    pgroup.nps[j] = len(ix)
    i1 = pgroup.ins[j] - 1
    i2 = pgroup.ins[j] + pgroup.nps[j] - 1
    pgroup.xp[i1:i2] = take(xx,ii)
    pgroup.yp[i1:i2] = take(yy,ii)
    pgroup.zp[i1:i2] = take(zz,ii)
    pgroup.uxp[i1:i2] = take(ux,ii)
    pgroup.uyp[i1:i2] = take(uy,ii)
    pgroup.uzp[i1:i2] = take(uz,ii)
    pgroup.gaminv[i1:i2] = take(gi,ii)

def sort3d(pgroup=None):
  if pgroup is None: pgroup = top.pgroup
  for j in range(pgroup.ns):
    xx = getx(js=j,gather=0)
    yy = gety(js=j,gather=0)
    zz = getz(js=j,gather=0)
    ux = getux(js=j,gather=0)
    uy = getuy(js=j,gather=0)
    uz = getuz(js=j,gather=0)
    gi = getgaminv(js=j,gather=0)
    ix = (abs(xx - w3d.xmmin)/w3d.dx).astype(long)
    iy = (abs(yy - w3d.ymmin)/w3d.dy).astype(long)
    iz = (abs(zz - w3d.zmmin)/w3d.dz).astype(long)
    ixy = ix + iy*(w3d.nx+1) + iz*(w3d.nx+1)*(w3d.ny+1)
    ii = argsort(ixy)
    pgroup.nps[j] = len(ix)
    i1 = pgroup.ins[j] - 1
    i2 = pgroup.ins[j] + pgroup.nps[j] - 1
    pgroup.xp[i1:i2] = take(xx,ii)
    pgroup.yp[i1:i2] = take(yy,ii)
    pgroup.zp[i1:i2] = take(zz,ii)
    pgroup.uxp[i1:i2] = take(ux,ii)
    pgroup.uyp[i1:i2] = take(uy,ii)
    pgroup.uzp[i1:i2] = take(uz,ii)
    pgroup.gaminv[i1:i2] = take(gi,ii)

