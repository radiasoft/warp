"""
Routines for adjusting the mesh for the 3-D self field
resizemesh: Change size of mesh
adjustmeshz: Adjust the longitudinal length of the mesh.
adjustmeshxy: Adjust the longitudinal length of the mesh.
"""
from warp import *
adjustmesh3d_version = "$Id: adjustmesh3d.py,v 1.10 2003/03/06 00:28:35 jlvay Exp $"

def adjustmesh3ddoc():
  import adjustmesh3d
  print adjustmesh3d.__doc__


# -------------------------------------------------------------------------
def resizemesh(nx=None,ny=None,nz=None,lloadrho=1,lfieldsol=1,
               linj=0,lzmom=0,lzarray=0,setobjects=None):
  """
Changes the number of grid points in the mesh.
Warning - this does not yet work in parallel
  """
  # --- Todo for parallel ...
  # ---  reset izextra if needed
  # ---  recalculate domain decomposition

  # --- Set defaults to original values
  if nx is None: nx = w3d.nx
  if ny is None: ny = w3d.ny
  if nz is None: nz = w3d.nz

  # --- If nothing changes, then just return
  if nx == w3d.nx and ny == w3d.ny and nz == w3d.nz: return

  if(w3d.nx>0):
    rx = float(nx)/float(w3d.nx)
  else:
    rx = 1.
  if(w3d.ny>0):
    ry = float(ny)/float(w3d.ny)
  else:
    ry = 1.
  rz = float(nz)/float(w3d.nz)

  # --- Set scalars
  w3d.nx = nx
  w3d.ny = ny
  w3d.nz = nz
  w3d.nzfull = w3d.nz
  w3d.izfsmax = w3d.nz
  w3d.nmxy  = max(w3d.nx,w3d.ny)
  w3d.nmxyz = max(w3d.nx,w3d.ny,w3d.nzfull)
  w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
  w3d.dz = (w3d.zmmax - w3d.zmmin)/w3d.nz

  # --- Reallocate the fields
  try:
    gallot("SelfFieldGrid3d")
  except:
    gallot("Fields3d")
  if w3d.solvergeom is w3d.RZgeom:
    frz.del_base()
    frz.init_base(w3d.nx,w3d.nz,w3d.dx,w3d.dz,w3d.xmmin,w3d.zmmin,false)

  # --- Calculate the mesh points
  w3d.xmesh[:] = w3d.xmmin + arange(w3d.nx+1)*w3d.dx
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.ymesh[:] = w3d.ymmin + arange(w3d.ny+1)*w3d.dy
  w3d.zmesh[:] = w3d.zmmin + arange(w3d.nz+1)*w3d.dz

  # --- Find the grid axis
  w3d.ix_axis = nint(-w3d.xmmin/w3d.dx)
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.iy_axis = nint(-w3d.ymmin/w3d.dy)
  w3d.iz_axis = nint(-w3d.zmmin/w3d.dz)

  # --- if requested, resize injection arrays accordingly
  if(linj):
    w3d.inj_nx = w3d.inj_nx*rx
    w3d.inj_ny = w3d.inj_ny*ry
    w3d.inj_nz = w3d.inj_nz*rz
    w3d.inj_dx = w3d.inj_dx/rx 
    w3d.inj_dy = w3d.inj_dy/ry 
    w3d.inj_dz = w3d.inj_dz/rz 
    injctint()

  # --- if requested, resize Z_Moments arrays accordingly 
  if(lzmom):
    top.nzmmnt = top.nzmmnt*rz
    gchange('Z_Moments')
    top.dzm = top.dzm/rz
    top.dzmi = 1./top.dzm
    for k in range(0,top.nzmmnt+1):
        top.zmntmesh[k] = top.zmmntmin + k*top.dzm

  # --- if requested, resize Z_arrays accordingly
  if(lzarray):
    top.nzzarr = top.nzzarr*rz
    gchange('Z_arrays')
    top.dzz = top.dzz/rz
    top.dzzi = 1./top.dzz
    for k in range(0,top.nzzarr+1):
        top.zplmesh[k] = top.zzmin + k*top.dzz
        top.prwallz[k] = top.prwall
        top.prwallxz[k] = top.prwallx
        top.prwallyz[k] = top.prwally
        top.prwelips[k] = 1.

  # --- Re-initialize any field solve parameters
  fieldsol(1)

  # --- Call subroutine for setting objects if provided
  if setobjects is not None:
    setobjects()

  # --- If requested, reload rho
  if lloadrho:
    w3d.rho = 0
    loadrho()

  # --- If requested, calculate the new fields
  if lfieldsol:
    fieldsol(-1)


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
  # --- Reset the global values of the mesh extent
  w3d.zmminglobal = w3d.zmmin
  w3d.zmmaxglobal = w3d.zmmax
  if lparallel:
    w3d.zmminglobal = top.zmslmin[0]
    w3d.zmmaxglobal = top.zmslmax[-1]
  # --- Recalculate zmesh
  w3d.zmesh[:] = w3d.zmmin + iota(0,w3d.nz)*w3d.dz
  # --- Adjust all of the axial meshes
  if top.nzl == w3d.nzfull:
    top.dzl = w3d.dz
    top.dzli = 1./w3d.dz
    top.zlmin = w3d.zmminglobal
    top.zlmax = w3d.zmmaxglobal
    top.zlmesh[:] = top.zlmin + iota(0,top.nzl)*top.dzl
    setlatt()
  if top.nzzarr == w3d.nzfull:
    top.dzz = w3d.dz
    top.dzzi = 1./w3d.dz
    top.zzmin = w3d.zmminglobal
    top.zzmax = w3d.zmmaxglobal
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

