from warp import *

#!#!#!#!#!#!#!#!#!#!#!#!#!#
##!#!#  TODO   #!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#
# realign the z-moments histories data

loadbalance_version = "$Id: loadbalance.py,v 1.24 2002/07/18 01:05:24 dave Exp $"

def loadbalancedoc():
  print """
Various routines for doing loading balancing for the parallel version
setparticledomains: Applies decomposition given a list of domain sizes
loadbalanceparticles: Load balances the particles based on pnumz
loadbalancesor: Load balances the SOR solver, balancing the total work in
                the solver including the work specifying the conductors.
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
  if not lparallel: return
  # --- It is assumed that the user supplied decomposition is specified
  # --- in the array zslave, which is an unscaled weighting of the z-ranges
  # --- of the particles for each processor.

  # --- All values of zslave must be > 0.
  assert min(zslave) > 0.,"The length of all particle domains must be > 0."

  # --- Save some data which will need to be redistributed
  eearsofz = gatherallzarray(top.eearsofz)
  prwallz  = gatherallzarray(top.prwallz)
  prwallxz = gatherallzarray(top.prwallxz)
  prwallyz = gatherallzarray(top.prwallyz)
  prwelips = gatherallzarray(top.prwelips)
  lostpars = gatherallzarray(top.lostpars)
  if(w3d.solvergeom == w3d.XYZgeom):
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

  # --- Make sure that the last processor doesn't have grid cells
  # --- sticking out the end.
  top.nzpslave[-1] = w3d.nzfull - top.izpslave[-1]

  # --- Adjust the Z data
  _adjustz()

  # --- Reorganize the particles
  reorgparticles()

  # --- Shift the existing charge density and phi
  if(w3d.solvergeom == w3d.XYZgeom):
    izstart = max(oldiz,top.izslave[me])
    izend = min(oldiz+oldnz,top.izslave[me]+top.nzslave[me])
    newiz1 = izstart - top.izslave[me]
    newiz2 = izend - top.izslave[me] + 1
    oldiz1 = izstart - oldiz
    oldiz2 = izend - oldiz + 1
    w3d.phi[:,:,newiz1+1:newiz2+1] = phi[:,:,oldiz1+1:oldiz2+1]
    w3d.rho[:,:,newiz1:newiz2] = rho[:,:,oldiz1:oldiz2]

  # --- Restore some data which needed to be redistributed
  top.eearsofz[:] = scatterallzarray(eearsofz)
  top.prwallz[:]  = scatterallzarray(prwallz)
  top.prwallxz[:] = scatterallzarray(prwallxz)
  top.prwallyz[:] = scatterallzarray(prwallyz)
  top.prwelips[:] = scatterallzarray(prwelips)
  top.lostpars[:] = scatterallzarray(lostpars)

  # --- Correct the locations of conductor points for the field-solver.
  if(w3d.solvergeom == w3d.XYZgeom):
    if top.fstype == 3:
      newiz = top.izslave[me]
      if f3d.ncond > 0:
        f3d.izcond[:f3d.ncond] = f3d.izcond[:f3d.ncond] + oldiz - newiz
      if f3d.necndbdy > 0:
        f3d.iecndz[:f3d.necndbdy] = f3d.iecndz[:f3d.necndbdy] + oldiz - newiz
      if f3d.nocndbdy > 0:
        f3d.iocndz[:f3d.nocndbdy] = f3d.iocndz[:f3d.nocndbdy] + oldiz - newiz
      cleanconductors()

  # --- Correct location of injection source.
  if(w3d.solvergeom == w3d.XYZgeom):
    if top.inject > 0:
      newiz = top.izslave[me]
      w3d.inj_grid[:,:,:] = w3d.inj_grid + oldiz - newiz
  else:
    gchange_rhop_phip_rz()
    
  # --- Do some additional work if requested
  if lloadrho: loadrho()
  if dofs: fieldsol(0)


#########################################################################
def loadbalanceparticles(lloadrho=1,dofs=1,spread=1.):
  """
Load balances the particles as evenly as possible. The load balancing is
based off of the data in top.pnumz which of course must already have
been calculated. The number density is assumed to vary linearly between
grid points.
 - lloadrho=1: when true, the charge density is redoposited
 - dofs=1: when true, the fields are recalculated
 - spread=1.: fraction of processors to spread the work among
  """
  if not lparallel: return
  # --- Gather pnumz. Return if there is no data.
  pnumz = gatherallzarray(top.pnumz)
  if max(pnumz) == 0.: return

  # --- Add fictitious data so that actual work is spread only to the
  # --- requested fraction of the processors.
  assert (0. < spread <= 1.),"spread must be between 0 and 1 or 1."
  avepnumz = ave(pnumz)
  pnumz[:] = pnumz + avepnumz*(1./spread - 1.)

  # --- Convert the number of particles to a decomposition
  zslave = decompose(pnumz,npes)

  # --- Apply the new domain decomposition.
  setparticledomains(zslave,lloadrho=lloadrho,dofs=dofs)

#########################################################################
def loadbalancesor(sgweight=7.0,condweight=2.0):
  """
Load balance the SOR field solver based off of the current timings. This is
needed since some processors may have more conductor points than others.
 - sgweight=7.: weight (in timing) of subgrid points relative to weight of a
                grid cell
 - condweight=2.: weight (in timing) of a conductor points relative to weight
                  of a grid cell
  """
  if not lparallel: return
  # --- Save the old values
  oldiz = top.izslave + 0
  oldnz = top.nzslave + 0
  oldizfs = top.izfsslave + 0
  oldnzfs = top.nzfsslave + 0
  oldphi = w3d.phi + 0.
  oldrho = w3d.rho + 0.

  # --- Make sure that the conductor arrays are allocated.
  if f3d.ncondmax == 0: f3d.ncondmax = 1
  if f3d.ncndmax == 0: f3d.ncndmax = 1
  gchange("Conductor3d")
    
  # --- Gather the field solve weights. For each z plane, sum the number of
  # --- grid cells, subgrid points, and conductor points, appropriately
  # --- weighted.
  weight = zeros(top.nzfsslave[me]+1,'d')
  for iz in iota(w3d.izfsmin,w3d.izfsmax):
    nec = len(nonzero(logical_not(f3d.iecndz[:f3d.necndbdy]-iz)))
    noc = len(nonzero(logical_not(f3d.iocndz[:f3d.nocndbdy]-iz)))
    nc  = len(nonzero(logical_not(f3d.izcond[:f3d.ncond]-iz)))
    weight[iz-w3d.izfsmin] = (w3d.nx+1)*(w3d.ny+1) + \
                             sgweight*(nec + noc) + \
                             condweight*nc
  weight = gatherallzfsarray(weight)

  # --- Convert to a decomposition
  zslave = decompose(weight,npes)

  # --- Set domain of each processor.
  # --- This coding ensures that all of the processors have nzfsslave at least
  # --- 2 or greater and that the last processor isn't left with too few cells.
  # --- In cases where all of the processors have values only slightly
  # --- greater than 2, the last processor will likely end up with too few
  # --- cells because of the accumulation of rounding up. Find the processors
  # --- which have the largest amount of roundup and take away one of their
  # --- grid cells until the last processor has enough. Do the same for the
  # --- case where the last processor has too many.
  zlast = 0
  for i in range(npes):
    top.izfsslave[i] = zlast
    top.nzfsslave[i] = max(nint(zslave[i]) + 1,2)
    zlast = top.izfsslave[i] + top.nzfsslave[i] - 1
  top.nzfsslave[-1] = w3d.nzfull - top.izfsslave[-1]
  while (zslave[-1]-top.nzfsslave[-1]) > max(zslave[:-1]-top.nzfsslave[:-1]) \
         or top.nzfsslave[-1] < 2:
    i = argmax(where(greater(top.nzfsslave[:-1],2),
                     top.nzfsslave[:-1]-zslave[:-1],-10000.))
    top.nzfsslave[i] = top.nzfsslave[i] - 1
    top.izfsslave[i+1:] = top.izfsslave[i+1:] - 1
    top.nzfsslave[-1] = w3d.nzfull - top.izfsslave[-1]
  while (top.nzfsslave[-1]-zslave[-1]) > max(top.nzfsslave[:-1]-zslave[:-1]):
    i = argmax(zslave[:-1]-top.nzfsslave[:-1])
    top.nzfsslave[i] = top.nzfsslave[i] + 1
    top.izfsslave[i+1:] = top.izfsslave[i+1:] + 1
    top.nzfsslave[-1] = w3d.nzfull - top.izfsslave[-1]

  # --- Adjust the Z data
  _adjustz()

  # --- Shift the existing charge density and phi
  izstart = max(oldiz[me],top.izslave[me])
  izend = min(oldiz[me]+oldnz[me],top.izslave[me]+top.nzslave[me])
  newiz1 = izstart - top.izslave[me]
  newiz2 = izend - top.izslave[me] + 1
  oldiz1 = izstart - oldiz[me]
  oldiz2 = izend - oldiz[me] + 1
  w3d.phi[:,:,newiz1+1:newiz2+1] = oldphi[:,:,oldiz1+1:oldiz2+1]
  w3d.rho[:,:,newiz1:newiz2] = oldrho[:,:,oldiz1:oldiz2]

  # --- Correct the locations of conductor points for the field-solver.
  newiz = top.izslave
  newnz = top.nzslave
  newizfs = top.izfsslave
  newnzfs = top.nzfsslave
  reorgconductors(oldiz,oldnz,oldizfs,oldnzfs,
                  newiz,newnz,newizfs,newnzfs)

  # --- Correct location of injection source.
  if top.inject > 0:
    w3d.inj_grid[:,:,:] = w3d.inj_grid + oldiz[me] - newiz[me]

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
  nz = len(weight) - 1
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
      if iz == nz+1: break
    if iz == nz+1: break
    # --- Add the last little bit to get to exactly npperpe.
    delta1 = delta
    a = 0.5*weight[iz] - 0.5*weight[iz+1]
    b = weight[iz]
    c = weight[iz]*(delta1 - 0.5*delta1**2) + 0.5*weight[iz+1]*delta1**2 + \
        npperpe - npint
    if b != 0.:
      delta = 2.*c/(sqrt(b**2 - 4.*a*c) + b)
    else:
      delta = sqrt(-c/a)
    #npint = npint + weight[iz]*((delta-delta1) + 0.5*(delta1**2-delta**2)) + \
    #                weight[iz+1]*0.5*(delta**2 - delta1**2)
    zslave[ip] = zslave[ip] + delta - delta1

  # --- The last processor gets everything left over
  zslave[-1] = nz - sum(zslave)

  return zslave


#########################################################################
def _adjustz():

  #---------------------------------------------------------------------------
  # --- If space-charge limited injection is turned on, then add additional
  # --- grid cells to the first processor so it will have the injection
  # --- surface completely covered.
 #if(w3d.solvergeom == w3d.XYZgeom):
 # if top.inject > 1:
 #  zinjmax = top.ainject**2/(top.rinject+sqrt(top.rinject**2-top.ainject**2))
 #  nzinj = int(max((top.zinject+zinjmax-top.zmslmin[0])/w3d.dz+top.inj_d)) + 1
 #  top.nzpslave[0] = max(top.nzpslave[0],nzinj)

  #---------------------------------------------------------------------------
  # --- Set the axial extent of each slaves domain to include
  # --- both the particle and field solve domain.
  if(w3d.solvergeom == w3d.XYZgeom):
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

  # --- Change the alocation of everything effected are reset the meshes.
  gchange("Fields3d")
  gchange("Z_arrays")
  gchange("LatticeInternal")
  gchange("Z_Moments")
  gchange("Hist")
  w3d.zmesh[:] = w3d.zmmin + iota(0,w3d.nz)*w3d.dz
  top.zplmesh[:] = top.zzmin + iota(0,top.nzzarr)*top.dzz
  top.zlmesh[:] = top.zlmin + iota(0,top.nzl)*top.dzl
  top.zmntmesh[:] = top.zmmntmin + iota(0,top.nzmmnt)*top.dzm
  
  # --- Reset the lattice
  setlatt()

#########################################################################
#########################################################################
# --- These are the messy routines for reorganizing the conductor data
def reorgconductors(oldiz,oldnz,oldizfs,oldnzfs,
                    newiz,newnz,newizfs,newnzfs):
  if globalsum(f3d.ncond) > 0:
    # --- Make things easier to deal with by ensuring that all arrays
    # --- are allocated.
    f3d.ncondmax = f3d.ncond + 1
    gchange("Conductor3d")

    # --- Shift the data to be relative to the global system
    f3d.izcond[:] = f3d.izcond[:] + oldiz[me]

    # --- Do the work
    results = _reorgconductorarrays([f3d.ixcond[:f3d.ncond], \
                                     f3d.iycond[:f3d.ncond], \
                                     f3d.izcond[:f3d.ncond], \
                                     f3d.condvolt[:f3d.ncond], \
                                     f3d.condnumb[:f3d.ncond]], \
                                    f3d.izcond[:f3d.ncond]+0, \
                                    oldiz,oldnz,oldizfs,oldnzfs, \
                                    newiz,newnz,newizfs,newnzfs)

    # --- Change array sizes and copy the data, localizing it.
    f3d.ncond = len(results[0])
    f3d.ncondmax = f3d.ncond
    gchange("Conductor3d")
    if f3d.ncond > 0:
      f3d.ixcond[:] = results[0]
      f3d.iycond[:] = results[1]
      f3d.izcond[:] = results[2] - newiz[me]
      f3d.condvolt[:] = results[3]
      f3d.condnumb[:] = results[4]

  if globalsum(f3d.necndbdy) > 0:
    # --- Make things easier to deal with by ensuring that all arrays
    # --- are allocated.
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy) + 1
    gchange("Conductor3d")

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
                                     f3d.ecvolt[:f3d.necndbdy], \
                                     f3d.ecnumb[:f3d.necndbdy]], \
                                    f3d.iecndz[:f3d.necndbdy]+0, \
                                    oldiz,oldnz,oldizfs,oldnzfs, \
                                    newiz,newnz,newizfs,newnzfs)

    # --- Change array sizes and copy the data, localizing it.
    f3d.necndbdy = len(results[0])
    if f3d.necndbdy > f3d.ncndmax:
      f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
      gchange("Conductor3d")
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
      f3d.ecnumb[:f3d.necndbdy] = results[10]

  if globalsum(f3d.nocndbdy) > 0:
    # --- Make things easier to deal with by ensuring that all arrays
    # --- are allocated.
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy) + 1
    gchange("Conductor3d")

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
                                     f3d.ocvolt[:f3d.nocndbdy], \
                                     f3d.ocnumb[:f3d.nocndbdy]], \
                                    f3d.iocndz[:f3d.nocndbdy]+0, \
                                    oldiz,oldnz,oldizfs,oldnzfs, \
                                    newiz,newnz,newizfs,newnzfs)

    # --- Change array sizes and copy the data, localizing it.
    f3d.nocndbdy = len(results[0])
    if f3d.nocndbdy > f3d.ncndmax:
      f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
      gchange("Conductor3d")
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
      f3d.ocnumb[:f3d.nocndbdy] = results[10]

#-------------------------------------------------------------------------
def _reorgconductorarrays(arrays,z,oldiz,oldnz,oldizfs,oldnzfs,
                                   newiz,newnz,newizfs,newnzfs):
  # --- Create list to save the incoming data in.
  results = len(arrays)*[[]]

  # --- Loop over global extent of grid, gathering data that is needed
  for iz in range(0,w3d.nzfull+1):

    # --- If me has the data then get the indices of it.
    if (oldizfs[me] <= iz <= oldizfs[me]+oldnzfs[me]):
      ii = compress(equal(iz,z),arange(len(z)))

      # --- If the data is needed by me, just copy it.
      if (newizfs[me] <= iz <= newizfs[me]+newnzfs[me]):
        for i in range(len(arrays)):
          results[i] = results[i] + list(take(arrays[i],ii))

    # --- Get the processor which "owns" the data, relative to the old
    # --- grid extents.
    pe = compress(logical_and(less_equal(oldizfs,iz),
                  less_equal(iz,oldizfs+oldnzfs)),arange(npes))[-1]

    if me == pe:
      # --- Loop over processors to check which ones need data
      # --- If the data is needed by others, then send it.
      for ip in range(npes):
        if not (oldizfs[ip] <= iz <= oldizfs[ip]+oldnzfs[ip]) and \
               (newizfs[ip] <= iz <= newizfs[ip]+newnzfs[ip]):
          for i in range(len(arrays)):
            temp = getarray(me,take(arrays[i],ii),ip)

    else:
      # --- Get the data that is needed from other procesors.
      if not (oldizfs[me] <= iz <= oldizfs[me]+oldnzfs[me]) and \
             (newizfs[me] <= iz <= newizfs[me]+newnzfs[me]):
        for i in range(len(arrays)):
          results[i] = results[i] + list(getarray(pe,0,me))

  # --- Make sure all processors are done before continuing
  barrier()

  for i in range(len(results)): results[i] = array(results[i])
  return results



