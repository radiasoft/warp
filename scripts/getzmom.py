from warp import *
getzmom_version = "$Id: getzmom.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

def zmmnt(itask=0):
  uxpo = top.uxp
  uypo = top.uyp
  uzpo = top.uzp

  # If changed, fix Win_Moments arrays
  if (len(top.pnumz) != top.nzmmnt+1): gchange ("Z_Moments")

  # If itask is greater than zero, assume that the intent was to do the
  # calculation regardless of the value of laccumulate_zmoments.  Save
  # the value of laccumulate_zmoments and then set it to false.
  if (itask > 0):
    l1 = top.laccumulate_zmoments
    top.laccumulate_zmoments = false

  # Zero out the moments arrays
  if (itask == 0 or itask == 1):
    getzmmnt(1,top.xp,top.yp,top.zp,top.uxp,top.uyp,top.uzp,top.gaminv,
             top.sq,top.sm,top.sw,top.dt/2.,1,
             top.nplive,uxpo,uypo,uzpo,1,top.ns)

  # Calculate the moments
  if (itask == 0 or itask == 2):
    for js in xrange(top.ns):
      for ipmin in xrange(top.ins[js]-1,top.ins[js]+top.nps[js]-1,256):
         ip = min(256, top.ins[js]+top.nps[js]-ipmin-1)
         getzmmnt(ip,top.xp[ipmin:],top.yp[ipmin:],top.zp[ipmin:],
                  top.uxp[ipmin:],top.uyp[ipmin:],top.uzp[ipmin:],
                  top.gaminv[ipmin:],top.sq[js],top.sm[js],top.sw[js],
                  top.dt,2,top.nplive,
                  uxpo[ipmin:],uypo[ipmin:],uzpo[ipmin:],js,top.ns)

# Do final summing and averaging of the moments
  if (itask == 0 or itask == 3):
    getzmmnt(1,top.xp,top.yp,top.zp,top.uxp,top.uyp,top.uzp,top.gaminv,
             top.sq,top.sm,top.sw,top.dt/2.,3,
             top.nplive,uxpo,uypo,uzpo,1,top.ns)

# Restore the value of laccumulate_zmoments
  if (itask > 0):
    top.laccumulate_zmoments = l1

