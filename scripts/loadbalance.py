from warp import *

#!#!#!#!#!#!#!#!#!#!#!#!#!#
##!#!#  TODO   #!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#
# realign the z-moments histories data

loadbalance_version = "$Id: loadbalance.py,v 1.5 2001/08/09 22:07:40 dave Exp $"

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
 - lloadrho=1: when true, the charge density is redeposited
 - dofs=1: when true, the fields are recalculated
  """
  # --- It is assumed that the user supplied decomposition is specified
  # --- in the array zslave, which is an unscaled weighting of the z-ranges
  # --- of the particles for each processor.

  # --- All values of zslave must be > 0.
  assert min(zslave) > 0.,"The length of all particle domains must be > 0."

  # --- Save some data which will need to be redistributed
  eearsofz = _gatherallzarray(top.eearsofz)
  prwallz  = _gatherallzarray(top.prwallz)
  prwallxz = _gatherallzarray(top.prwallxz)
  prwallyz = _gatherallzarray(top.prwallyz)
  prwelips = _gatherallzarray(top.prwelips)
  lostpars = _gatherallzarray(top.lostpars,'i')
  phi = w3d.phi + 0.
  rho = w3d.rho + 0.

  # --- Broadcast the window moments to all processors and gather lab window
  # --- data onto PE0. This is needed since the processors which own the
  # --- windows may change.
  getwin_moments()
  gethist()
  getlabmoments()

  # --- Save the current extent of the grid. This is used to correct the
  # --- z location of any conductor points for the field-solver.
  oldiz = top.izslave[me]
  oldnz = top.nzslave[me]

  # --- Get sum of zslave to allow proper scaling.
  sumzslave = sum(zslave)

  # --- Set domain of each processor.
  zlast = top.zmslmin[0]
  for i in range(npes):
    top.zpslmin[i] = zlast
    top.zpslmax[i] = zlast+zslave[i]/sumzslave*(top.zmslmax[-1]-top.zmslmin[0])
    zlast = top.zpslmax[i]

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
  w3d.izfsmax = w3d.izfsmin + top.nzfsslave[me]
  w3d.zmmin = top.zmslmin[me]
  w3d.zmmax = top.zmslmax[me]
  if top.fstype in [3,7] or f3d.nsorerr > 0:
    # --- The additional check of nsorerr>0 is done since in some cases,
    # --- the field solver may be turned off when the load balancing is done
    # --- but used otherwise.
    f3d.nzpsor = w3d.nz
    f3d.nsorerr = (w3d.nx+1)*(w3d.ny+1)*(w3d.nz+1)/(.9*w3d.nx-1)

  # --- Change the alocation of everything effected are reset the meshes.
  gchange("Fields3d")
  gchange("Z_arrays")
  gchange("LatticeInternal")
  gchange("Z_Moments")
  if top.fstype in [3,7] or f3d.nsorerr > 0: gchange("PSOR3d")
  w3d.zmesh[:] = w3d.zmmin + iota(0,w3d.nz)*w3d.dz
  top.zplmesh[:] = top.zzmin + iota(0,top.nzzarr)*top.dzz
  top.zlmesh[:] = top.zlmin + iota(0,top.nzl)*top.dzl
  top.zmntmesh[:] = top.zmmntmin + iota(0,top.nzmmnt)*top.dzm
  
  # --- Reset the lattice
  resetlat()
  setlatt()

  # --- Reorganize the particles
  reorgparticles()

  # --- Shift the existing charge density and phi
  izstart = max(oldiz,top.izslave[me])
  izend = min(oldiz+oldnz,top.izslave[me]+top.nzslave[me])
  newiz1 = izstart - top.izslave[me]
  newiz2 = izend - top.izslave[me] + 1
  oldiz1 = izstart - oldiz
  oldiz2 = izend - oldiz + 1
  w3d.phi[:,:,newiz1+1:newiz2+1] = phi[:,:,oldiz1+1:oldiz2+1]
  w3d.rho[:,:,newiz1:newiz2] = rho[:,:,oldiz1:oldiz2]

  # --- Restore some data which needed to be redistributed
  top.eearsofz[:] = _scatterallzarray(eearsofz)
  top.prwallz[:]  = _scatterallzarray(prwallz)
  top.prwallxz[:] = _scatterallzarray(prwallxz)
  top.prwallyz[:] = _scatterallzarray(prwallyz)
  top.prwelips[:] = _scatterallzarray(prwelips)
  top.lostpars[:] = _scatterallzarray(lostpars)

  # --- Correct the locations of conductor points for the field-solver.
  newiz = top.izslave[me]
  if f3d.ncond > 0:
    f3d.izcond[:f3d.ncond] = f3d.izcond[:f3d.ncond] + oldiz - newiz
  if f3d.necndbdy > 0:
    f3d.iecndz[:f3d.necndbdy] = f3d.iecndz[:f3d.necndbdy] + oldiz - newiz
  if f3d.nocndbdy > 0:
    f3d.iocndz[:f3d.nocndbdy] = f3d.iocndz[:f3d.nocndbdy] + oldiz - newiz
  cleanconductors()

  # --- Correct location of injection source.
  if top.inject > 0:
    w3d.inj_grid[:,:] = w3d.inj_grid + oldiz - newiz

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
  pnumz = gatherallzarray(top.pnumz)
# pnumz = zeros(w3d.nzfull+1,'d')
# for iz in range(0,w3d.nzfull+1):
#   pe = convertiztope(iz)
#   if me == pe: nn = top.pnumz[iz-top.izpslave[me]]
#   else:        nn = 0
#   pnumz[iz] = mpi.bcast(nn,pe)

  # --- Convert the number of particles to a decomposition
  zslave = decompose(pnumz,npes)

  # --- Apply the new domain decomposition.
  setparticledomains(zslave,lloadrho=lloadrho,dofs=dofs)

#########################################################################
def loadbalancesor():
  """
Load balance the SOR field solver based off of the current timings. This is
needed since some processors may have more conductor points than others.
  """
  # --- Gather the field solve timings
  weight = gatherarray(top.fstime)
  weight = broadcast(weight)

  # --- Divide by nzfsslave to get the relative timings per nz
  weight = weight/(top.nzfsslave-1)

  # --- Convert to a decomposition and scale to w3d.nzfull
  zslave = decompose(weight,npes)
  zslave = zslave/max(zslave)*w3d.nzfull

  # --- Save the old values
  oldiz = top.izslave + 0
  oldnz = top.nzslave + 0
  oldizfs = top.izfsslave + 0
  oldnzfs = top.nzfsslave + 0

  # --- Set domain of each processor.
  zlast = 0
  for i in range(npes):
    top.izfsslave[i] = zlast
    top.nzfsslave[i] = nint(zslave[i]) + 1
    zlast = top.izfsslave[i] + top.nzfsslave[i] - 1
  top.nzfsslave[-1] = w3d.nzfull - w3d.izfsslave[i]

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
  w3d.izfsmin = top.izfsslave[me] - top.izslave[me]
  w3d.izfsmax = w3d.izfsmin + top.nzfsslave[me]
  w3d.zmmin = top.zmslmin[me]
  w3d.zmmax = top.zmslmax[me]
  if top.fstype in [3,7] or f3d.nsorerr > 0:
    # --- The additional check of nsorerr>0 is done since in some cases,
    # --- the field solver may be turned off when the load balancing is done
    # --- but used otherwise.
    f3d.nzpsor = w3d.nz
    f3d.nsorerr = (w3d.nx+1)*(w3d.ny+1)*(w3d.nz+1)/(.9*w3d.nx-1)

  # --- Change the alocation of everything effected are reset the meshes.
  gchange("Fields3d")
  if top.fstype in [3,7] or f3d.nsorerr > 0: gchange("PSOR3d")
  w3d.zmesh[:] = w3d.zmmin + iota(0,w3d.nz)*w3d.dz
  
  # --- Shift the existing charge density and phi
  izstart = max(oldiz,top.izslave[me])
  izend = min(oldiz+oldnz,top.izslave[me]+top.nzslave[me])
  newiz1 = izstart - top.izslave[me]
  newiz2 = izend - top.izslave[me] + 1
  oldiz1 = izstart - oldiz
  oldiz2 = izend - oldiz + 1
  w3d.phi[:,:,newiz1+1:newiz2+1] = phi[:,:,oldiz1+1:oldiz2+1]
  w3d.rho[:,:,newiz1:newiz2] = rho[:,:,oldiz1:oldiz2]

  # --- Correct the locations of conductor points for the field-solver.
  newiz = top.izslave
  newnz = top.nzslave
  newizfs = top.izfsslave
  newnzfs = top.nzfsslave
  reorgconductors(oldiz,oldnz,oldizfs,oldnzfs,
                  newiz,newnz,newizfs,newnzfs)

  # --- Correct location of injection source.
  if top.inject > 0:
    w3d.inj_grid[:,:] = w3d.inj_grid + oldiz - newiz

#########################################################################
def decompose(weight,npes):
  """
Converts a weight into the size of the domains.
 - weight: array of relative weights of the work done by each processor
 - npes: number of processors
Returns an array of the same length which is which the relative length of each
of the domains.
  """
  # --- Integrate weight, assuming linear variation between grid points
  nz = len(weight)
  np = 0.5*weight[0] + sum(weight[1:-1]) + 0.5*weight[-1]
  npperpe = 1.*np/npes

  zslave = zeros(npes,'d')
  iz = 0
  delta = 0.
  for ip in range(npes-1):
    npint = 0.
    npnext = weight[iz  ]*((1.-delta)+0.5*(delta**2-1.)) + \
             weight[iz+1]*0.5*(1. - delta**2)
    # --- Get the remaining bit from the last cell if it is not too much.
    if npnext < npperpe:
      zslave[ip] = 1. - delta
      iz = iz + 1
      delta = 0.
      npint = npnext
    # --- Keep adding cells until the number per processor is reached.
    while npint + 0.5*(weight[iz]+weight[iz+1]) < npperpe:
      zslave[ip] = zslave[ip] + 1.
      delta = 0.
      npint = npint + 0.5*(weight[iz]+weight[iz+1])
      iz = iz + 1
      if iz == nz: break
    if iz == nz: break
    # --- Add the last little bit to get to exactly npperpe.
    a = 0.5*weight[iz] - 0.5*weight[iz+1]
    b = weight[iz]
    c = weight[iz]*(delta - 0.5*delta**2) + 0.5*weight[iz+1]*delta**2 + \
        npperpe - npint
    if b != 0.:
      delta = 2.*c/(sqrt(b**2 - 4.*a*c) + b)
    else:
      delta = sqrt(-c/a)
    zslave[ip] = zslave[ip] + delta

  # --- The last processor gets everything left over
  zslave[-1] = nz - sum(zslave)

  return zslave

#########################################################################
#########################################################################
# --- These are the messy routines for reorganizing the conductor data
def reorgconductors(oldiz,oldnz,oldizfs,oldnzfs,
                    newiz,newnz,newizfs,newnzfs):
  if globalsum(f3d.ncond) > 0:
    # --- Make things easier to deal with by ensuring that all arrays
    # --- are allocated.
    f3d.ncondmax = f3d.ncond + 1
    gchange("PSOR3d")

    # --- Shift the data to be relative to the global system
    f3d.izcond[:] = f3d.izcond[:] + oldiz[me]

    # --- Do the work
    results = _reorgconductorarrays([f3d.ixcond[:f3d.ncond], \
                                     f3d.iycond[:f3d.ncond], \
                                     f3d.izcond[:f3d.ncond], \
                                     f3d.condvolt[:f3d.ncond]], \
                                    f3d.izcond[:f3d.ncond]+0, \
                                    oldiz,oldnz,oldizfs,oldnzfs, \
                                    newiz,newnz,newizfs,newnzfs)

    # --- Change array sizes and copy the data, localizing it.
    f3d.ncond = len(results[0])
    f3d.ncondmax = f3d.ncond
    gchange("PSOR3d")
    if f3d.ncond > 0:
      f3d.ixcond[:] = results[0]
      f3d.iycond[:] = results[1]
      f3d.izcond[:] = results[2] - newiz[me]
      f3d.condvolt[:] = results[3]

  if globalsum(f3d.necndbdy) > 0:
    # --- Make things easier to deal with by ensuring that all arrays
    # --- are allocated.
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy) + 1
    gchange("PSOR3d")

    # --- Shift the data to be relative to the global system
    f3d.iecndz[:] = f3d.iecndz[:] + oldiz[me]

    # --- Do the work
    results = _reorgconductorarrays([f3d.iecndx[:f3d.necndbdy], \
                                     f3d.iecndy[:f3d.necndbdy], \
                                     f3d.iecndz[:f3d.necndbdy], \
                                     f3d.ecdelmx[:f3d.necndbdy], \
                                     f3d.ecdelmy[:f3d.necndbdy], \
                                     f3d.ecdelmz[:f3d.necndbdy], \
                                     f3d.ecdelpx[:f3d.necndbdy], \
                                     f3d.ecdelpy[:f3d.necndbdy], \
                                     f3d.ecdelpz[:f3d.necndbdy], \
                                     f3d.ecvolt[:f3d.necndbdy]], \
                                    f3d.iecndz[:f3d.necndbdy]+0, \
                                    oldiz,oldnz,oldizfs,oldnzfs, \
                                    newiz,newnz,newizfs,newnzfs)

    # --- Change array sizes and copy the data, localizing it.
    f3d.necndbdy = len(results[0])
    if f3d.necndbdy > f3d.ncndmax:
      f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
      gchange("PSOR3d")
    if f3d.necndbdy > 0:
      f3d.iecndx[:f3d.necndbdy] = results[0]
      f3d.iecndy[:f3d.necndbdy] = results[1]
      f3d.iecndz[:f3d.necndbdy] = results[2] - newiz[me]
      f3d.ecdelmx[:f3d.necndbdy] = results[3]
      f3d.ecdelmy[:f3d.necndbdy] = results[4]
      f3d.ecdelmz[:f3d.necndbdy] = results[5]
      f3d.ecdelpx[:f3d.necndbdy] = results[6]
      f3d.ecdelpy[:f3d.necndbdy] = results[7]
      f3d.ecdelpz[:f3d.necndbdy] = results[8]
      f3d.ecvolt[:f3d.necndbdy] = results[9]

  if globalsum(f3d.nocndbdy) > 0:
    # --- Make things easier to deal with by ensuring that all arrays
    # --- are allocated.
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy) + 1
    gchange("PSOR3d")

    # --- Shift the data to be relative to the global system
    f3d.iocndz[:] = f3d.iocndz[:] + oldiz[me]

    # --- Do the work
    results = _reorgconductorarrays([f3d.iocndx[:f3d.nocndbdy], \
                                     f3d.iocndy[:f3d.nocndbdy], \
                                     f3d.iocndz[:f3d.nocndbdy], \
                                     f3d.ocdelmx[:f3d.nocndbdy], \
                                     f3d.ocdelmy[:f3d.nocndbdy], \
                                     f3d.ocdelmz[:f3d.nocndbdy], \
                                     f3d.ocdelpx[:f3d.nocndbdy], \
                                     f3d.ocdelpy[:f3d.nocndbdy], \
                                     f3d.ocdelpz[:f3d.nocndbdy], \
                                     f3d.ocvolt[:f3d.nocndbdy]], \
                                    f3d.iocndz[:f3d.nocndbdy]+0, \
                                    oldiz,oldnz,oldizfs,oldnzfs, \
                                    newiz,newnz,newizfs,newnzfs)

    # --- Change array sizes and copy the data, localizing it.
    f3d.nocndbdy = len(results[0])
    if f3d.nocndbdy > f3d.ncndmax:
      f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
      gchange("PSOR3d")
    if f3d.nocndbdy > 0:
      f3d.iocndx[:f3d.nocndbdy] = results[0]
      f3d.iocndy[:f3d.nocndbdy] = results[1]
      f3d.iocndz[:f3d.nocndbdy] = results[2] - newiz[me]
      f3d.ocdelmx[:f3d.nocndbdy] = results[3]
      f3d.ocdelmy[:f3d.nocndbdy] = results[4]
      f3d.ocdelmz[:f3d.nocndbdy] = results[5]
      f3d.ocdelpx[:f3d.nocndbdy] = results[6]
      f3d.ocdelpy[:f3d.nocndbdy] = results[7]
      f3d.ocdelpz[:f3d.nocndbdy] = results[8]
      f3d.ocvolt[:f3d.nocndbdy] = results[9]

#-------------------------------------------------------------------------
def _reorgconductorarrays(arrays,z,oldiz,oldnz,oldizfs,oldnzfs,
                                   newiz,newnz,newizfs,newnzfs):
  # --- Create list to save the incoming data in.
  results = len(arrays)*[[]]

 # --- Loop over global extent of grid, gathering data from other processors
 for iz in range(0,w3d.nzfull+1):
   # --- Get the processor which "owns" the data, relative to the old
   # --- grid extents.
   pe = compress(logical_and(less_equal(oldizfs,iz),
                 less_equal(iz,oldizfs+oldnzfs)),arange(npes))[-1]
   # --- Check if data at this iz is needed locally
   if not (oldiz[me] <= iz <= oldiz[me]+oldnz[me]) and \
          (newizfs[me] <= iz <= newizfs[me]+newnzfs[me]):
     for i in range(len(arrays)):
       results[i] = results[i] + list(getarray(pe,0,me))
   elif me == pe:
     # --- Loop over processors to check which ones need data
     for ip in range(npes):
       if not (oldiz[ip] <= iz <= oldiz[ip]+oldnz[ip]) and \
              (newizfs[ip] <= iz <= newizfs[ip]+newnzfs[ip]):
         ii = compress(equal(iz,z),arange(len(arrays[0])))
         for i in range(len(arrays)):
           temp = getarray(me,take(arrays[i],ii),ip)

  # --- Make sure all processors are done before continuing
  barrier()

  # --- Gather any data that is still local
  for iz in range(newizfs[me],newizfs[me]+newnzfs[me]+1):
    ii = compress(equal(iz,z),arange(len(arrays[0])))
    for i in range(len(arrays)):
      results[i] = results[i] + list(take(arrays[i],ii))

  for i in range(len(results)): results[i] = array(results[i])
  return results



