from warp import *
adjustmesh3d_version = "$Id: adjustmesh3d.py,v 1.6 2001/07/11 22:40:16 dave Exp $"

def adjustmesh3ddoc():
  print "adjustmeshz: Adjust the longitudinal length of the mesh."

# -------------------------------------------------------------------------
def adjustmeshz(newlen,dorho=1,dofs=0,keepcentered=0):
  """Adjust the longitudinal length of the mesh.
  - newlen: the new length of the mesh
  - dorho=1: when true, the charge density is recalculated
  - dofs=0: when true, the fieldsolver is called
  - keepcentered=0: when true, the center of the mesh is unchanged, otherwise,
                    zmmin and zmmax are simply scaled by the ratio of the
                    new length to the old length
  """
  # --- Save old grid cell and mesh length
  olddz = w3d.dz
  if not lparallel:
    oldcenter = 0.5*(w3d.zmmin + w3d.zmmax)
  else:
    oldcenter = 0.5*(top.zmslmin[0] + top.zmslmax[-1])
  # --- Set new mesh length by first scaling the min and max
  w3d.dz = newlen/w3d.nzfull
  w3d.zmmin = w3d.zmmin*w3d.dz/olddz
  w3d.zmmax = w3d.zmmax*w3d.dz/olddz
  if lparallel:
    top.zmslmin[:] = top.zmslmin*w3d.dz/olddz
    top.zmslmax[:] = top.zmslmax*w3d.dz/olddz
    top.zpslmin[:] = top.zpslmin*w3d.dz/olddz
    top.zpslmax[:] = top.zpslmax*w3d.dz/olddz
  # --- If requested, recenter the mesh about its old center.
  if keepcentered:
    if not lparallel:
      newcenter = 0.5*(w3d.zmmin + w3d.zmmax)
    else:
      newcenter = 0.5*(top.zmslmin[0] + top.zmslmax[-1])
    w3d.zmmin = w3d.zmmin + (oldcenter - newcenter)
    w3d.zmmax = w3d.zmmax + (oldcenter - newcenter)
    if lparallel:
      top.zmslmin[:] = top.zmslmin + (oldcenter - newcenter)
      top.zmslmax[:] = top.zmslmax + (oldcenter - newcenter)
      top.zpslmin[:] = top.zpslmin + (oldcenter - newcenter)
      top.zpslmax[:] = top.zpslmax + (oldcenter - newcenter)
  # --- Recalculate zmesh
  w3d.zmesh[:] = w3d.zmmin + iota(0,w3d.nz)*w3d.dz
  # --- Adjust all of the axial meshes
  if top.nzl == w3d.nz:
    top.dzl = w3d.dz
    top.dzli = 1./w3d.dz
    top.zlmin = w3d.zmmin
    top.zlmax = w3d.zmmax
    top.zlmesh[:] = top.zlmin + iota(0,top.nzl)*top.dzl
    setlatt()
  if top.nzzarr == w3d.nz:
    top.dzz = w3d.dz
    top.dzzi = 1./w3d.dz
    top.zzmin = w3d.zmmin
    top.zzmax = w3d.zmmax
    top.zplmesh[:] = top.zzmin + iota(0,top.nzzarr)*top.dzz
  if top.nzmmnt == w3d.nz:
    top.dzm = w3d.dz
    top.dzmi = 1./w3d.dz
    top.zmmntmin = w3d.zmmin
    top.zmmntmax = w3d.zmmax
    top.zmntmesh[:] = top.zmmntmin + iota(0,top.nzmmnt)*top.dzm
  # --- Rearrange the particles
  if npes > 0:
    reorgparticles()
  else:
    zpartbnd(w3d.zmmax,w3d.zmmin,w3d.dz,top.zgrid)
  # --- Reset field solve parameters (kzsq)
  if top.fstype >= 0:
    fstypesave = top.fstype
    top.fstype = 0
    fieldsol(1)
    top.fstype = fstypesave
  # --- Redeposit the charge density
  if dorho:
    w3d.rho = 0.
    loadrho()
  # --- Now ready for next field solve
  if dofs:
    fieldsol(-1)

# -------------------------------------------------------------------------
def adjustmeshxy(newsize,dorho=1,dofs=0,keepcentered=0):
  """Adjust the longitudinal length of the mesh.
  - newsize: the new transverse size of the mesh
  - dorho=1: when true, the charge density is recalculated
  - dofs=0: when true, the fieldsolver is called
  - keepcentered=0: when true, the center of the mesh is unchanged, otherwise,
                    zmmin and zmmax are simply scaled by the ratio of the
                    new length to the old length
  """
  # --- Save old grid cell and mesh length
  olddx = w3d.dx
  olddy = w3d.dy
  oldcenterx = 0.5*(w3d.xmmin + w3d.xmmax)
  oldcentery = 0.5*(w3d.ymmin + w3d.ymmax)
  # --- Set new mesh length by first scaling the min and max
  w3d.dx = newsize/w3d.nx
  w3d.dy = newsize/w3d.ny
  w3d.xmmin = w3d.xmmin*w3d.dx/olddx
  w3d.xmmax = w3d.xmmax*w3d.dx/olddx
  w3d.ymmin = w3d.ymmin*w3d.dy/olddy
  w3d.ymmax = w3d.ymmax*w3d.dy/olddy
  # --- If requested, recenter the mesh about its old center.
  if keepcentered:
    newcenterx = 0.5*(w3d.xmmin + w3d.xmmax)
    newcentery = 0.5*(w3d.ymmin + w3d.ymmax)
    w3d.xmmin = w3d.xmmin + (oldcenterx - newcenterx)
    w3d.xmmax = w3d.xmmax + (oldcenterx - newcenterx)
    w3d.ymmin = w3d.ymmin + (oldcentery - newcentery)
    w3d.ymmax = w3d.ymmax + (oldcentery - newcentery)
  # --- Recalculate mesh
  w3d.xmesh[:] = w3d.xmmin + iota(0,w3d.nx)*w3d.dx
  w3d.ymesh[:] = w3d.ymmin + iota(0,w3d.ny)*w3d.dy
  # --- Recheck particle boundary conditions
  ###
  # --- Reset field solve parameters (kxsq and kysq)
  if top.fstype >= 0:
    fstypesave = top.fstype
    top.fstype = 0
    fieldsol(1)
    top.fstype = fstypesave
  # --- Redeposit the charge density
  if dorho:
    w3d.rho = 0.
    loadrho()
  # --- Now ready for next field solve
  if dofs:
    fieldsol(-1)

