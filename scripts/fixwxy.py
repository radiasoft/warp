from warp import *
fixwxy_version = "$Id: fixwxy.py,v 1.2 2001/03/01 22:23:37 dave Exp $"
# Fixes beam so that it exactly agress with the specified beam paramters

# --- Fixes 1st and 2nd moments
def fixwxy2(a=None,b=None,ap=None,bp=None,x=None,y=None,xp=None,yp=None):
  if not a: a = top.a0
  if not ap: ap = top.ap0
  if not b: b = top.b0
  if not bp: bp = top.bp0
  if not x: x = top.x0
  if not xp: xp = top.xp0
  if not y: y = top.y0
  if not yp: yp = top.yp0
  vx = xp*top.vbeam
  vy = yp*top.vbeam

  # --- Fix the beam center and size to match exactly to the input
  top.xp[:] = where(not_equal(top.uzp, 0.),top.xp - top.xbar[0] + x,0.)
  top.yp[:] = where(not_equal(top.uzp, 0.),top.yp - top.ybar[0] + y,0.)
  top.xp[:] = top.xp*a/(2.*top.xrms[0])
  top.yp[:] = top.yp*b/(2.*top.yrms[0])

  # --- Fix velocity to match emittance and envelope angle
  # --- First, remove any average velocity
  top.uxp[:] = where(not_equal(top.uzp, 0.),top.uxp - top.vxbar[0],0.)
  top.uyp[:] = where(not_equal(top.uzp, 0.),top.uyp - top.vybar[0],0.)
  # --- Then remove any coherent velocity
  slopex = top.xxpbar[0]/top.xsqbar[0]*top.vbeam
  slopey = top.yypbar[0]/top.ysqbar[0]*top.vbeam
  top.uxp[:] = where(not_equal(top.uzp, 0.),top.uxp - slopex*top.xp,0.)
  top.uyp[:] = where(not_equal(top.uzp, 0.),top.uyp - slopey*top.yp,0.)
  # --- Scale to get correct thermal spread
  if top.emitn != 0:
    emit = top.emitn*top.clight/top.vbeam
  else:
    emit = top.emit
  top.uxp[:] = top.uxp[:]*(emit/top.a0)/(top.epsx[0]/(2.*top.xrms[0]))
  top.uyp[:] = top.uyp[:]*(emit/top.b0)/(top.epsy[0]/(2.*top.yrms[0]))
  # --- Added back in specified coherent velocity and average
  slopex = top.ap0/top.a0*top.vbeam
  slopey = top.bp0/top.b0*top.vbeam
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
  if not x: x = top.x0
  if not xp: xp = top.xp0
  if not y: y = top.y0
  if not yp: yp = top.yp0
  rr = sqrt(ranf())
  tt = ranf()*2.*pi
  xx = x + xerr*rr*cos(tt)
  yy = y + xerr*rr*sin(tt)
  rr = sqrt(ranf())
  tt = ranf()*2.*pi
  vx = (xp + xperr*rr*cos(tt))*top.vbeam
  vy = (yp + yperr*rr*sin(tt))*top.vbeam

  # --- Fix the beam center to match exactly to the input
  top.xp[:] = where(not_equal(top.uzp, 0.),top.xp - top.xbar[0] + xx,0.)
  top.yp[:] = where(not_equal(top.uzp, 0.),top.yp - top.ybar[0] + yy,0.)

  # --- Fix velocity to match envelope angle
  top.uxp[:] = where(not_equal(top.uzp, 0.),top.uxp - top.vxbar[0] + vx,0.)
  top.uyp[:] = where(not_equal(top.uzp, 0.),top.uyp - top.vybar[0] + vy,0.)

  # --- Now fix up rho and diagnostics
  w3d.rho = 0.
  loadrho()
  import getzmom
  getzmom.zmmnt()
  srhoax3d()
  rhodia3d()
  if replacehist: top.jhist = top.jhist - 1
  savehist(top.time)

