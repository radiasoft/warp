from warp import *
errorcheck_version = "$Id:"

def errorcheckdoc():
  print "errorcheck checks for consistency errors and other possible mistakes"

def errorcheck():
  """
Checks for consistency errors and other possible mistakes.
It is not quaranteed to find all mistakes.
  """
  # --- Check for use of symmetry
  if w3d.l4symtry:
    ok = 1
    if top.x0 != 0: ok = 0
    if top.y0 != 0: ok = 0
    if top.xp0 != 0: ok = 0
    if top.yp0 != 0: ok = 0
    if max(abs(top.x0_s)) != 0: ok = 0
    if max(abs(top.y0_s)) != 0: ok = 0
    if max(abs(top.xp0_s)) != 0: ok = 0
    if max(abs(top.yp0_s)) != 0: ok = 0
    if max(abs(top.xcent_s)) != 0: ok = 0
    if max(abs(top.ycent_s)) != 0: ok = 0
    if max(abs(top.xpcent_s)) != 0: ok = 0
    if max(abs(top.ypcent_s)) != 0: ok = 0
    if not ok:
      raise "ERROR: Beam is being offset with four-fold symmetry turned on"

  if w3d.l2symtry:
    ok = 1
    if top.y0 != 0: ok = 0
    if top.yp0 != 0: ok = 0
    if max(abs(top.y0_s)) != 0: ok = 0
    if max(abs(top.yp0_s)) != 0: ok = 0
    if max(abs(top.ycent_s)) != 0: ok = 0
    if max(abs(top.ypcent_s)) != 0: ok = 0
    if not ok:
      raise "ERROR: Beam is being offset in y with two-fold symmetry turned on"

  # --- Make sure the envelope doesn't extend outside the grid when it is
  # --- being used to load particles.
  if env.nenv > 0:
    izl = (top.zimin - env.zl)/env.dzenv
    izr = (top.zimax - env.zl)/env.dzenv
    if izl < 0 or izr > 0:
      raise "ERROR: beam axial extent extends beyond range of envelope calculation"
    if (max(env.aenv) > w3d.xmmax or min(env.aenv) < w3d.xmmin or
    if  max(env.benv) > w3d.ymmax or min(env.benv) < w3d.ymmin):




