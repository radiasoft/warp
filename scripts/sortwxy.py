from warp import *
sortwxy_version = "$Id: sortwxy.py,v 1.3 2002/07/10 17:57:25 dave Exp $"
# Sorts the particles for better cache use for the gather/scatter.

def sortwxy():
  for j in range(top.ns):
    xx = getx(js=j,gather=0)
    yy = gety(js=j,gather=0)
    zz = getz(js=j,gather=0)
    ux = getux(js=j,gather=0)
    uy = getuy(js=j,gather=0)
    uz = getuz(js=j,gather=0)
    gi = getgaminv(js=j,gather=0)
    ix = (abs(xx - w3d.xmmin)/w3d.dx).astype(Int)
    iy = (abs(yy - w3d.ymmin)/w3d.dy).astype(Int)
    ixy = ix + iy*(w3d.nx+1)
    ii = argsort(ixy)
    top.nps[j] = len(ix)
    i1 = top.ins[j] - 1
    i2 = top.ins[j] + top.nps[j] - 1
    top.xp[i1:i2] = take(xx,ii)
    top.yp[i1:i2] = take(yy,ii)
    top.zp[i1:i2] = take(zz,ii)
    top.uxp[i1:i2] = take(ux,ii)
    top.uyp[i1:i2] = take(uy,ii)
    top.uzp[i1:i2] = take(uz,ii)
    top.gaminv[i1:i2] = take(gi,ii)

def sort3d():
  for j in range(top.ns):
    xx = getx(js=j,gather=0)
    yy = gety(js=j,gather=0)
    zz = getz(js=j,gather=0)
    ux = getux(js=j,gather=0)
    uy = getuy(js=j,gather=0)
    uz = getuz(js=j,gather=0)
    gi = getgaminv(js=j,gather=0)
    ix = (abs(xx - w3d.xmmin)/w3d.dx).astype(Int)
    iy = (abs(yy - w3d.ymmin)/w3d.dy).astype(Int)
    iz = (abs(zz - w3d.zmmin)/w3d.dz).astype(Int)
    ixy = ix + iy*(w3d.nx+1) + iz*(w3d.nx+1)*(w3d.ny+1)
    ii = argsort(ixy)
    top.nps[j] = len(ix)
    i1 = top.ins[j] - 1
    i2 = top.ins[j] + top.nps[j] - 1
    top.xp[i1:i2] = take(xx,ii)
    top.yp[i1:i2] = take(yy,ii)
    top.zp[i1:i2] = take(zz,ii)
    top.uxp[i1:i2] = take(ux,ii)
    top.uyp[i1:i2] = take(uy,ii)
    top.uzp[i1:i2] = take(uz,ii)
    top.gaminv[i1:i2] = take(gi,ii)

