from warp import *
sortwxy_version = "$Id: sortwxy.py,v 1.1 2001/06/18 20:43:55 dave Exp $"
# Sorts the particles for better cache use for the gather/scatter.

def sortwxy():
  for j in range(top.ns):
    js = j + 1
    xx = getx(js=js)
    yy = gety(js=js)
    zz = getz(js=js)
    ux = getux(js=js)
    uy = getuy(js=js)
    uz = getuz(js=js)
    gi = getgaminv(js=js)
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

