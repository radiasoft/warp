from warp import *

loadbalance_version = "$Id: loadbalance.py,v 1.1 2001/07/12 00:28:38 dave Exp $"

def loadbalancedoc():
  print """
Various routines for doing loading balancing for the parallel version
setparticledomains: Applies decomposition given a list of domain sizes
loadbalanceparticles: Load balances the particles based on pnumz
  """

#########################################################################
def setparticledomains(zslave,lloadrho=1,dofs=1):
  """
Sets the particles domains from the input, zslave, in the same way as done
with top.zslave during the generate. This is only meant to be used after
that has already been done.
 - zslave: list of sizes of domains - does not need any normalization
 - lloadrho=1: when true, the charge density is redoposited
 - dofs=1: when true, the fields are recalculated
  """
  # --- It is assumed that the user supplied decomposition is specified
  # --- in the array zslave, which is an unscaled weighting of the z-ranges
  # --- of the particles for each processor.

  # --- All values of zslave must be > 0.
  assert min(zslave) > 0.,"The length of all particle domains must be > 0."

  # --- Get sum of zslave to allow proper scaling.
  sumzslave = sum(zslave)

  # --- Set domain of each processor.
  zlast = top.zmslmin[0]
  for i in range(npes):
    top.zpslmin[i] = zlast
    top.zpslmax[i] = zlast+zslave[i]/sumzslave*(top.zmslmax[-1]-top.zmslmin[0])
    zlast = zpslmax[i]

  # --- This is only needed to avoid problems from round off in the
  # --- accumulation. From the loop above, zpslmax[-1] will
  # --- not be exactly the same as zmmax due to roundoff.
  top.zpslmax[-1] = top.zmslmax[-1]

  # --- Set iz and nz. This is done so that zmesh[izpslave] < zpslmin, and
  # --- zmesh[izpslave+nzpslave] > zpslmax.
  for i in range(npes):
    top.izpslave[i] = int((top.zpslmin[i] - top.zmslmin[0])/w3d.dz)
    top.nzpslave[i] = int((top.zpslmax[i] - top.zmslmin[0])/w3d.dz) - \
                      top.izpslave[i] + 1

  # --- Make sure that the last processors doesn't have grid cells
  # --- sticking out the end.
  top.nzpslave[-1] = w3d.nzfull - top.izpslave[-1]

  #---------------------------------------------------------------------------
  # --- Now set the axial extent of each slaves domain to include
  # --- both the particle and field solve domain.
  for i in range(npes):
    top.izslave[i] = min(top.izpslave[i],top.izfsslave[i])
    top.nzslave[i] = max(top.izpslave[i] + top.nzpslave[i], \
                         top.izfsslave[i] + top.nzfsslave[i]) - top.izslave[i]
    top.zmslmin[i] = top.izslave[i]*w3d.dz + top.zmslmin[0]
    top.zmslmax[i] = (top.izslave[i] + top.nzslave[i])*w3d.dz + top.zmslmin[0]

  #---------------------------------------------------------------------------
  # --- Reset local values
  w3d.nz     = top.nzslave[me]
  top.nzzarr = top.nzpslave[me]
  top.nzl    = top.nzpslave[me]
  top.nzlmax = top.nzpslave[me]
  top.nzmmnt = top.nzpslave[me]
  zpmin = top.zmslmin[0] + top.izpslave[me]*w3d.dz
  zpmax = (top.izpslave[me]+top.nzpslave[me])*w3d.dz + top.zmslmin[0]
  top.zzmin = zpmin
  top.zzmax = zpmax
  top.zlmin = zpmin
  top.zlmax = zpmax
  top.zmmntmin = zpmin
  top.zmmntmax = zpmax
  w3d.izfsmin = top.izfsslave[me] - top.izslave[me]
  w3d.izfsmax = top.izfsmin + top.nzfsslave[me]
  w3d.zmmin = top.zmslmin[me]
  w3d.zmmax = top.zmslmax[me]

  # --- Change the alocation of everything effected are reset the meshes.
  gchange("Fields3d")
  gchange("Z_arrays")
  gchange("LatticeInternal")
  gchange("Z_Moments")
  w3d.zmesh[:] = w3d.zmmin + iota(0,w3d.nz)*w3d.dz
  top.zplmesh[:] = top.zzmin + iota(0,top.nzzarr)*top.dzz
  top.zlmesh[:] = top.zlmin + iota(0,top.nzl)*top.dzl
  top.zmntmesh[:] = top.zmmntmin + iota(0,top.nzmmnt)*top.dzm
  
  # --- Reorganize the particles
  reorgparticles()

  # --- Do some additional work if requested
  if lloadrho: loadrho()
  if dofs: fieldsol(0)


#########################################################################
def loadbalanceparticles(lloadrho=1,dofs=1):
  """
Load balances the particles as evenly as possible. The load balancing is
based off of the data in top.pnumz which of course must already have
been calculated. The number density is assumed to vary linearly between
grid points.
 - lloadrho=1: when true, the charge density is redoposited
 - dofs=1: when true, the fields are recalculated
  """
  # --- Gather pnumz. The commented out line does not work since there
  # --- maybe a complicated series of overlap among the processors.
  #pnumz = gatherarray(top.pnumz
  pnumz = zeros(w3d.nzfull+1,'d')
  for iz in range(0,w3d.nzfull+1):
    pe = convertiztope(iz)
    pnumz[iz] = mpi.bcast(top.pnumz[iz-top.izpslave[me]],pe)

  # --- Integrate pnumz, assuming linear variation between grid points
  np = 0.5*pnumz[0] + sum(pnumz[1:-1]) + 0.5*pnumz[-1]
  npperpe = 1.*np/npes

  zslave = zeros(npes,'d')
  iz = 0
  delta = 0.
  for ip in range(npes):
    npint = 0.
    npnext = pnumz[iz  ]*((1.-delta)+0.5*(delta**2-1.)) + \
             pnumz[iz+1]*0.5*(1. - delta**2)
    # --- Get the remaining bit from the last cell if it is not too much.
    if npnext < npperpe:
      zslave[ip] = 1. - delta
      iz = iz + 1
      delta = 0.
      npint = npnext
    # --- Keep adding cells until the number per processor is reached.
    while npint + 0.5*(pnumz[iz]+pnumz[iz+1]) < npperpe:
      zslave[ip] = zslave[ip] + 1.
      delta = 0.
      npint = npint + 0.5*(pnumz[iz]+pnumz[iz+1])
      iz = iz + 1
      if iz == w3d.nzfull: break
    if iz == w3d.nzfull: break
    # --- Add the last little bit to get to exactly npperpe.
    if pnumz[iz] != pnumz[iz+1]:
      a = 0.5*pnumz[iz] - 0.5*pnumz[iz+1]
      b = -pnumz[iz]
      c = pnumz[iz]*(delta - 0.5*delta**2) + 0.5*pnumz[iz+1]*delta**2 + \
          npperpe - npint
      delta = (-b + sqrt(b**2 - 4.*a*c))/(2.*a)
    else:
      b = -pnumz[iz]
      c = pnumz[iz]*(delta - 0.5*delta**2) + 0.5*pnumz[iz+1]*delta**2 + \
          npperpe - npint
      delta = -c/b
    zslave[ip] = zslave[ip] + delta

  # --- Apply the new domain decomposition.
  setparticledomains(zslave,lloadrho=lloadrho,dofs=dofs):

