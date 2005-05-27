"""Fixes beam so that it exactly agress with the specified beam paramters
"""
from warp import *
fixwxy_version = "$Id: fixwxy.py,v 1.8 2005/05/27 22:40:03 dave Exp $"

# --- Fixes 1st and 2nd moments
def fixwxy2(a=None,b=None,ap=None,bp=None,x=None,y=None,xp=None,yp=None, 
            emitx=None,emity=None):
  if a  is None: a  = top.a0
  if ap is None: ap = top.ap0
  if b  is None: b  = top.b0
  if bp is None: bp = top.bp0
  if x  is None: x  = top.x0
  if xp is None: xp = top.xp0
  if y  is None: y  = top.y0
  if yp is None: yp = top.yp0
  # --- x-emittance 
  if emitx is None: 
    if top.emitn == 0. and top.emitnx == 0.:
      # --- set from unnormalized emittance  
      if top.emitx == 0.: 
        emitx = top.emit
      else:
        emitx = top.emitx 
    else: 
      # --- set from normalized emittance 
      if top.emitnx == 0.: 
        emitx = top.emitn*top.clight/top.vbeam
      else: 
        emitx = top.emitnx*top.clight/top.vbeam
  # --- y-emittance 
  if emity is None: 
    if top.emitn == 0. and top.emitny == 0.:
      # --- set from unnormalized emittance  
      if top.emity == 0.: 
        emity = top.emit
      else: 
        emity = top.emity 
    else: 
      # --- set from normalized emittance 
      if top.emitny == 0.: 
        emity = top.emitn*top.clight/top.vbeam
      else: 
        emity = top.emitny*top.clight/top.vbeam

  vx = xp*top.vbeam
  vy = yp*top.vbeam

  # --- Fix the beam center and size to match exactly to the input
  top.xp[:] = where(not_equal(top.uzp, 0.),top.xp - top.xbar[0,-1] + x,0.)
  top.yp[:] = where(not_equal(top.uzp, 0.),top.yp - top.ybar[0,-1] + y,0.)
  top.xp[:] = top.xp*a/(2.*top.xrms[0,-1])
  top.yp[:] = top.yp*b/(2.*top.yrms[0,-1])

  # --- Fix velocity to match emittance and envelope angle
  # --- First, remove any average velocity
  if top.lrelativ: gi = top.gaminv
  else:            gi = 1.
  top.uxp[:] = where(not_equal(top.uzp, 0.),top.uxp - top.vxbar[0,-1]/gi,0.)
  top.uyp[:] = where(not_equal(top.uzp, 0.),top.uyp - top.vybar[0,-1]/gi,0.)
  # --- Then remove any coherent velocity
  top.xxpbar[0] = ave(getx()*getxp())
  top.yypbar[0] = ave(gety()*getyp())
  slopex = top.xxpbar[0,-1]/(a/2.)**2*top.vbeam
  slopey = top.yypbar[0,-1]/(b/2.)**2*top.vbeam
  top.uxp[:] = where(not_equal(top.uzp, 0.),top.uxp - slopex*top.xp,0.)
  top.uyp[:] = where(not_equal(top.uzp, 0.),top.uyp - slopey*top.yp,0.)
  # --- Scale to get correct thermal spread
  top.uxp[:] = top.uxp[:]*(emitx/a)/(top.epsx[0,-1]/(2.*top.xrms[0,-1]))
  top.uyp[:] = top.uyp[:]*(emity/b)/(top.epsy[0,-1]/(2.*top.yrms[0,-1]))
  # --- Added back in specified coherent velocity and average
  slopex = ap/a*top.vbeam
  slopey = bp/b*top.vbeam
  top.uxp[:]=where(not_equal(top.uzp,0.),top.uxp+slopex*top.xp+vx,0.)
  top.uyp[:]=where(not_equal(top.uzp,0.),top.uyp+slopey*top.yp+vy,0.)

  # --- Now fix up rho and diagnostics
  w3d.rho = 0.
  loadrho()
  import getzmom
  getzmom.zmmnt()
  srhoax3d()
  rhodia3d()
  top.jhist = -1
  savehist(0.)


# --- Fixes only 1st moments (beam centroid)
def fixwxy1(x=None,y=None,xp=None,yp=None,xerr=0.,yerr=0.,xperr=0.,yperr=0.,
            replacehist=1):
  if x  is None: x  = top.x0
  if xp is None: xp = top.xp0
  if y  is None: y  = top.y0
  if yp is None: yp = top.yp0

  rr = sqrt(ranf())
  tt = ranf()*2.*pi
  xx = x + xerr*rr*cos(tt)
  yy = y + xerr*rr*sin(tt)
  rr = sqrt(ranf())
  tt = ranf()*2.*pi
  vx = (xp + xperr*rr*cos(tt))*top.vbeam
  vy = (yp + yperr*rr*sin(tt))*top.vbeam

  # --- Fix the beam center to match exactly to the input
  top.xp[:] = where(not_equal(top.uzp, 0.),top.xp - top.xbar[0,-1] + xx,0.)
  top.yp[:] = where(not_equal(top.uzp, 0.),top.yp - top.ybar[0,-1] + yy,0.)

  # --- Fix velocity to match envelope angle
  top.uxp[:] = where(not_equal(top.uzp, 0.),top.uxp - top.vxbar[0,-1] + vx,0.)
  top.uyp[:] = where(not_equal(top.uzp, 0.),top.uyp - top.vybar[0,-1] + vy,0.)

  # --- Now fix up rho and diagnostics
  w3d.rho = 0.
  loadrho()
  import getzmom
  getzmom.zmmnt()
  srhoax3d()
  rhodia3d()
  if replacehist: top.jhist = top.jhist - 1
  savehist(top.time)

