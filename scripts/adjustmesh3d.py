"""
Routines for adjusting the mesh for the 3-D self field
resizemesh: Change size of mesh
adjustmeshz: Adjust the longitudinal length of the mesh.
adjustmeshxy: Adjust the longitudinal length of the mesh.
"""
from warp import *
adjustmesh3d_version = "$Id: adjustmesh3d.py,v 1.32 2011/11/07 23:01:05 grote Exp $"

def adjustmesh3ddoc():
  import adjustmesh3d
  print adjustmesh3d.__doc__

def resizeZ_arrays(zzmin=None,zzmax=None,nzzarr=None):
  if zzmin is not None: top.zzmin = zzmin
  if zzmax is not None: top.zzmax = zzmax
  if nzzarr is not None: top.nzzarr = nzzarr
  top.dzz = (top.zzmax - top.zzmin)/top.nzzarr
  top.dzzi = 1./top.dzz
  gchange('Z_arrays')
  top.zplmesh = span(top.zzmin,top.zzmax,top.nzzarr+1)

def resizeLatticeInternal(zlmin=None,zlmax=None,nzlmax=None):
  if zlmin is not None: top.zlmin = zlmin
  if zlmax is not None: top.zlmax = zlmax
  if nzlmax is not None: top.nzlmax = nzlmax
  top.nzl = top.nzlmax
  top.dzl = (top.zlmax - top.zlmin)/top.nzlmax
  top.dzli = 1./top.dzl
  gchange("LatticeInternal")
  top.zlmesh = span(top.zlmin,top.zlmax,top.nzlmax+1)
  setlatt()

def resizeZ_Momments(zmmntmin=None,zmmntmax=None,nzmmnt=None):
  if zmmntmin is not None: top.zmmntmin = zmmntmin
  if zmmntmax is not None: top.zmmntmax = zmmntmax
  if nzmmnt is not None: top.nzmmnt = nzmmnt
  top.dzm = (top.zmmntmax - top.zmmntmin)/top.nzmmnt
  top.dzmi = 1./top.dzm
  gchange("Z_Moments")
  top.zmntmesh = span(top.zmmntmin,top.zmmntmax,top.nzmmnt+1)

def resizetopmeshes(zmmin=None,zmmax=None,nz=None):
  if zmmin is None: zmmin = w3d.zmmin
  if zmmax is None: zmmax = w3d.zmmax
  if nz is None: nz = w3d.nz
  resizeZ_arrays(zmmin,zmmax,nz)
  resizeLatticeInternal(zmmin,zmmax,nz)
  resizeZ_Momments(zmmin,zmmax,nz)

# -------------------------------------------------------------------------
def resizemesh(nx=None,ny=None,nz=None,lloadrho=True,lfieldsol=True,
               linj=False,lzmom=False,lzarray=False,llattice=False,conductors=None):
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
  w3d.nzlocal = w3d.nz
  w3d.nmxy  = max(w3d.nx,w3d.ny)
  w3d.nmxyz = max(w3d.nx,w3d.ny,w3d.nz)
  w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
  w3d.dz = (w3d.zmmax - w3d.zmmin)/w3d.nz

  # --- Reallocate the fields
  try:
    gallot("SelfFieldGrid3d")
  except:
    gallot("Fields3d")
  setupFields3dParticles()
  if w3d.solvergeom is w3d.RZgeom:
    frz.del_base()
    frz.init_base(w3d.nx,w3d.nz,w3d.dx,w3d.dz,w3d.xmmin,w3d.zmminlocal)

  # --- Calculate the mesh points
  w3d.xmesh[:] = w3d.xmmin + arange(w3d.nx+1)*w3d.dx
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.ymesh[:] = w3d.ymmin + arange(w3d.ny+1)*w3d.dy
  w3d.zmesh[:] = w3d.zmmin + arange(w3d.nz+1)*w3d.dz
  w3d.zmeshlocal[:] = w3d.zmminlocal + arange(w3d.nzlocal+1)*w3d.dz

  # --- Find the grid axis
  w3d.ix_axis = nint(-w3d.xmmin/w3d.dx)
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.iy_axis = nint(-w3d.ymmin/w3d.dy)
  w3d.iz_axis = nint(-w3d.zmmin/w3d.dz)

  # --- if requested, resize injection arrays accordingly
  if(linj):
    w3d.inj_nx = nint(w3d.inj_nx*rx)
    w3d.inj_ny = nint(w3d.inj_ny*ry)
    w3d.inj_nz = nint(w3d.inj_nz*rz)
    w3d.inj_dx = w3d.inj_dx/rx 
    w3d.inj_dy = w3d.inj_dy/ry 
    w3d.inj_dz = w3d.inj_dz/rz 
    injctint()

  # --- resize various things if requested
  if lzmom: resizeZ_Momments(zmmntmin=w3d.zmmin,zmmntmax=w3d.zmmax,nzmmnt=w3d.nz)
  if lzarray: resizeZ_arrays(zzmin=w3d.zmmin,zzmax=w3d.zmmax,nzzarr=w3d.nz)
  if llattice: resizeLatticeInternal(zlmin=w3d.zmmin,zlmax=w3d.zmmax,nzlmax=w3d.nz)

  # --- Re-initialize any field solve parameters
  fieldsol(1)

  # --- Reinstall conductors
  if conductors is not None:
    for conductor in conductors:
        conductor.install()

  # --- If requested, reload rho
  if lloadrho:
    w3d.rho = 0
    loadrho()

  # --- If requested, calculate the new fields
  if lfieldsol:
    fieldsol(-1,lbeforefs=1,lafterfs=1)

# -------------------------------------------------------------------------
def resizemeshxy(nx=None,ny=None,xmmin=None,xmmax=None,ymmin=None,ymmax=None,
                 lloadrho=1,lfieldsol=1,setobjects=None):
  """
Resizes the transverse size of the mesh
 - nx, ny: The new values of the number of grid points. If not given, default
           to current values.
 - xmmin,xmmax,ymmin,ymmax: The new values of the extent of the mesh. If not
                            given, default to current values.
 - lloadrho=1: Flag specifying whether charge density is to be reloaded.
 - llfieldsol=1: Flag specifying whether the field should be recalculated.
 - setobjects: Optional function to regenerate data for conductors.
  """
  # --- Set defaults to original values
  if nx is None: nx = w3d.nx
  if ny is None: ny = w3d.ny
  if xmmin is None: xmmin = w3d.xmmin
  if xmmax is None: xmmax = w3d.xmmax
  if ymmin is None: ymmin = w3d.ymmin
  if ymmax is None: ymmax = w3d.ymmax

  # --- If nothing changes, then just return
  if (nx == w3d.nx and ny == w3d.ny and
      xmmin==w3d.xmmin and xmmax==w3d.xmmax and
      ymmin==w3d.ymmin and ymmax==w3d.ymmax): return

  # --- Make sure nx and ny are reasonable
  assert nx>0 and ny>0,"nx and ny must be greater than zero"
  assert xmmax>xmmin and ymmax>=ymmin,"max's must be greater than min's"

  # --- Set scalars
  w3d.nx = nx
  w3d.ny = ny
  w3d.xmmin = xmmin
  w3d.xmmax = xmmax
  w3d.ymmin = ymmin
  w3d.ymmax = ymmax
  if w3d.solvergeom==w3d.XYZgeom:
    if w3d.l2symtry:
      w3d.ymmin = 0.
    elif w3d.l4symtry:
      w3d.xmmin = 0.
      w3d.ymmin = 0.
  elif w3d.solvergeom==w3d.XZgeom:
    w3d.ymmin = 0.
    if w3d.l2symtry or w3d.l4symtry: w3d.xmmin = 0.
  elif w3d.solvergeom==w3d.RZgeom or w3d.solvergeom==w3d.Zgeom:
    w3d.xmmin = 0.
    w3d.ymmin = 0.

  w3d.nmxy  = max(w3d.nx,w3d.ny)
  w3d.nmxyz = max(w3d.nx,w3d.ny,w3d.nz)
  w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny

  # --- Reallocate the fields
  gallot("Fields3d")
  setupFields3dParticles()
  if w3d.solvergeom is w3d.RZgeom:
    frz.del_base()
    frz.init_base(w3d.nx,w3d.nz,w3d.dx,w3d.dz,w3d.xmmin,w3d.zmminlocal)

  # --- Calculate the mesh points
  w3d.xmesh[:] = w3d.xmmin + arange(w3d.nx+1)*w3d.dx
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.ymesh[:] = w3d.ymmin + arange(w3d.ny+1)*w3d.dy

  # --- Find the grid axis
  w3d.ix_axis = nint(-w3d.xmmin/w3d.dx)
  if w3d.solvergeom in [w3d.XYZgeom, w3d.AMRgeom]:
    w3d.iy_axis = nint(-w3d.ymmin/w3d.dy)

  # --- Fix selfe if needed
  allocateselfepforparticles(false)

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
    fieldsol(-1,lbeforefs=1,lafterfs=1)


# -------------------------------------------------------------------------
def adjustmeshz(newlen,dorho=1,dofs=0,keepcentered=0,adjustppdecomp=true):
  """Adjust the longitudinal length of the mesh.
  - newlen: the new length of the mesh
  - dorho=1: when true, the charge density is recalculated
  - dofs=0: when true, the fieldsolver is called
  - keepcentered=0: when true, the center of the mesh is unchanged, otherwise,
                    zmmin and zmmax are simply scaled by the ratio of the
                    new length to the old length
  - adjustppdecomp=true: when true, adjust the particle decomposition
                         by the same scaling
  """
  # --- Save old grid cell and mesh length
  olddz = w3d.dz

  # --- Scale everything by the ratio of the new and old grid cell size
  w3d.dz = newlen/w3d.nz

  # --- Shift the center if requested.
  if keepcentered:
    oldcenter = 0.5*(w3d.zmmin + w3d.zmmax)
    newcenter = 0.5*(w3d.zmmin + w3d.zmmax)*w3d.dz/olddz
    delcenter = oldcenter - newcenter
  else:
    delcenter = 0.

  # --- Scale the mins and maxes.
  w3d.zmmin = w3d.zmmin*w3d.dz/olddz + delcenter
  w3d.zmmax = w3d.zmmax*w3d.dz/olddz + delcenter
  w3d.zmminlocal = w3d.zmminlocal*w3d.dz/olddz + delcenter
  w3d.zmmaxlocal = w3d.zmmaxlocal*w3d.dz/olddz + delcenter
  top.fsdecomp.zmin = top.fsdecomp.zmin*w3d.dz/olddz + delcenter
  top.fsdecomp.zmax = top.fsdecomp.zmax*w3d.dz/olddz + delcenter
  if adjustppdecomp:
    # --- Scale the particle extrema the same way.
    top.ppdecomp.zmin = top.ppdecomp.zmin*w3d.dz/olddz + delcenter
    top.ppdecomp.zmax = top.ppdecomp.zmax*w3d.dz/olddz + delcenter
    # --- These are no longer used, but update them anyway
    top.zpslmin = top.zpslmin*w3d.dz/olddz + delcenter
    top.zpslmax = top.zpslmax*w3d.dz/olddz + delcenter
  else:
    # --- The min and max must be adjusted to be within the new grid
    top.ppdecomp.zmin = top.ppdecomp.zmin.clip(w3d.zmmin,w3d.zmmax)
    top.ppdecomp.zmax = top.ppdecomp.zmax.clip(w3d.zmmin,w3d.zmmax)
    # --- iz and nz must be adjusted to be relative to the new zmin
    top.ppdecomp.iz = int((top.ppdecomp.zmin - w3d.zmmin)/w3d.dz)
    top.ppdecomp.nz = int((top.ppdecomp.zmax - w3d.zmmin)/w3d.dz -
                          top.ppdecomp.iz)
    # --- Check for the cases where int rounded down.
    zmax1 = w3d.zmmin + (top.ppdecomp.iz + top.ppdecomp.nz)*w3d.dz
    top.ppdecomp.nz[top.ppdecomp.zmax > zmax1] += 1

  w3d.zmminp = w3d.zmmin + top.ppdecomp.iz[top.ppdecomp.izproc]*w3d.dz
  w3d.zmmaxp = w3d.zmmin + ((top.ppdecomp.iz[top.ppdecomp.izproc] +
                             top.ppdecomp.nz[top.ppdecomp.izproc])*w3d.dz)
  top.zpminlocal = top.ppdecomp.zmin[top.ppdecomp.izproc]
  top.zpmaxlocal = top.ppdecomp.zmax[top.ppdecomp.izproc]

  # --- Recalculate zmesh
  w3d.zmesh[:] = w3d.zmmin + arange(1+w3d.nz)*w3d.dz
  w3d.zmeshlocal[:] = w3d.zmminlocal + arange(1+w3d.nzlocal)*w3d.dz

  # --- Adjust all of the axial meshes
  if top.nzl == w3d.nz:
    top.dzl = w3d.dz
    top.dzli = 1./w3d.dz
    top.zlmin = w3d.zmmin
    top.zlmax = w3d.zmmax
    top.zlmesh[:] = top.zlmin + arange(1+top.nzl)*top.dzl
    setlatt()
  if top.nzzarr == w3d.nz:
    top.dzz = w3d.dz
    top.dzzi = 1./w3d.dz
    top.zzmin = w3d.zmmin
    top.zzmax = w3d.zmmax
    top.zplmesh[:] = top.zzmin + arange(1+top.nzzarr)*top.dzz
  if top.nzmmnt == w3d.nz:
    top.dzm = w3d.dz
    top.dzmi = 1./w3d.dz
    top.zmmntmin = w3d.zmmin
    top.zmmntmax = w3d.zmmax
    top.zmntmesh[:] = top.zmmntmin + arange(1+top.nzmmnt)*top.dzm

  # --- Rearrange the particles
  particlegridboundaries3d(top.pgroup,-1)
  for js in xrange(top.ns):
    processlostpart(top.pgroup,js+1,top.clearlostpart,top.time,top.zbeam)

  # --- Reset field solve parameters (kzsq)
  if 0 <= top.fstype <= 4:
    vp3d(1)

  # --- Redeposit the charge density
  if dorho:
    loadrho()

  # --- Now ready for next field solve
  if dofs:
    fieldsolve(-1,lbeforefs=1,lafterfs=1)

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
    fieldsol(-1,lbeforefs=1,lafterfs=1)

