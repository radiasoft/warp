from warp import *
errorcheck_version = "$Id:"

def errorcheckdoc():
  print "errorcheck: runs all checks described below"
  print "checksymmetry: checks use of symmetry"
  print "checkparticleload: checks particle loading"

############################################################################
def errorcheck():
  """
Checks for consistency errors and other possible mistakes.
It is not quaranteed to find all mistakes.
  """
  checksymmetry()
  checkparticleload()

############################################################################
############################################################################
############################################################################
def checksymmetry():
  """Check for use of symmetry"""
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

############################################################################
def checkparticleload():
  """Checks for errors in particle loading"""
  # --- Get current package
  currpkg = package()[0]

  # --- Calculate some handy temps
  gridxmax = w3d.xmmax
  gridymax = w3d.ymmax
  if w3d.l4symtry:
    gridxmin = -w3d.xmmax
    gridymin = -w3d.ymmax
  elif w3d.l2symtry:
    gridxmin = w3d.xmmin
    gridymin = -w3d.ymmax
  else:
    gridxmin = w3d.xmmin
    gridymin = w3d.ymmin

  # --- Make sure the envelope doesn't extend outside the grid when it is
  # --- being used to load particles.
  if currpkg=='w3d' and env.nenv > 0 and not w3d.cylinder and w3d.nenvofz == 0:
    izl = (top.zimin - env.zl)/env.dzenv
    izr = (top.zimax - env.zl)/env.dzenv
    if izl < 0 or izr > 0:
      raise "ERROR: beam axial extent extends beyond range of envelope calculation"
    if (max(+env.aenv[izl:izr+1]+env.xenv[izl:izr+1]) > gridxmax or
        min(-env.aenv[izl:izr+1]+env.xenv[izl:izr+1]) < gridxmin or
        max(+env.benv[izl:izr+1]+env.yenv[izl:izr+1]) > gridymax or
        min(-env.benv[izl:izr+1]+env.yenv[izl:izr+1]) < gridymin):
      raise "ERROR: transverse extent of envelope is larger than the field grid"

  # --- Make sure user supplied envelope is not bigger than the mesh.
  if currpkg=='w3d' and w3d.nenvofz > 0:
    if (max(+w3d.aofz+w3d.xofz) > gridxmax or
        min(-w3d.aofz+w3d.xofz) < gridxmin or
        max(+w3d.bofz+w3d.yofz) > gridymax or
        min(-w3d.bofz+w3d.yofz) < gridymin):
      raise "ERROR: User specified axially varying envelope extends beyond the field grid"

  # --- Make sure slice particles are loaded within the mesh
  if currpkg=='wxy':
    if (+top.a0 + top.x0 > gridxmax or
        -top.a0 + top.x0 < gridxmin or
        +top.b0 + top.y0 > gridymax or
        -top.b0 + top.y0 < gridymin):
      raise "ERROR: grid size extends beyond the field grid"



