"""
Various routines for doing loading balancing for the parallel version
LoadBalancer: class wrapping particle load balancing. Sets up automatic
              periodic load balancing of particles.
setparticledomains: Applies decomposition given a list of domain sizes
loadbalanceparticles: Load balances the particles based on pnumz
loadbalancesor: Load balances the SOR solver, balancing the total work in
                the solver including the work specifying the conductors.
"""
from warp import *

loadbalance_version = "$Id: loadbalance.py,v 1.37 2005/03/17 21:40:29 dave Exp $"

def loadbalancedoc():
  import loadbalance
  print loadbalance.__doc__

#########################################################################
#########################################################################
class LoadBalancer:
  """
Installs load balancer.
Creation arguments:
 - when: dictionary of when to do the load balancing. Keys are time step
         numbers, values are frequency of loadbalancing when top.it is less
         than key. Default is {10:1,100:10,1000000:20}
 - padright: Amount of space added to right end of grid. When not specified,
             it is product of max(vz)*top.dt*2 and the number of steps between
             load balances.
 - padleft=0: Amount of space added to left end of grid.
 - doloadrho=0: Specifies whether the charge density is recalculated
 - dofs=0: Specifies whether the fields are recalculated
 - verbose=0: Prints output
Note, if particles on cover a few grid cells, then distribution is
recalculated on a finer mesh to give better balancing.
  """
  def __init__(self,padright=None,padleft=0.,when=None,doitnow=0,
               doloadrho=0,dofs=0,verbose=0):
    if not lparallel: return
    if when is None:
      self.when = {10:1,100:10,1000000:20}
    else:
      self.when = when
    self.padright = padright
    self.padleft = padleft
    self.doloadrho = doloadrho
    self.dofs = dofs
    self.verbose = verbose
    if doitnow: self.doloadbalance()
    installafterstep(self.doloadbalance)

  def doloadbalance(self,lforce=0,doloadrho=None,dofs=None,reorg=None):
    if not lparallel: return

    # --- Check if there are any particles anywhere
    if getn(jslist=-1) == 0: return

    if not top.lmoments or top.laccumulate_zmoments:
      # --- If the moments were not calculated, then top.zminp and top.zmaxp
      # --- are not reliable and so need to be calculated.
      zz = getz(jslist=-1,gather=0)
      if len(zz) > 0:
        zminp = min(zz)
        zmaxp = max(zz)
      else:
        zminp = +largepos
        zmaxp = -largepos
      zminp = globalmin(zminp)
      zmaxp = globalmax(zmaxp)
    else:
      # --- Otherwise, use the values already calculated.
      zz = None
      zminp = top.zminp[-1]
      zmaxp = top.zmaxp[-1]

    # --- Special check when injection is turned on
    if top.inject:
      # --- Make sure that all of the injection sources are included.
      # --- These are crude estimates of the min and max when xpinject
      # --- and ypinject are nonzero.
      rinj = sqrt(top.ainject**2 + top.binject**2)
      rpinj = sqrt(top.xpinject**2 + top.ypinject**2)
      zinjectmin = max(w3d.zmminglobal,min(top.zinject - rinj*rpinj))
      zinjectmax = min(w3d.zmmaxglobal,max(top.zinject + rinj*rpinj))
      zminp = minimum(zinjectmin,zminp)
      zmaxp = maximum(zinjectmax,zmaxp)

    # --- Check if rightmost particle is close to edge of last processor
    # --- If so, then force a reloadbalance.
    if top.zpslmax[-1] < w3d.zmmaxglobal-w3d.dz:
      if zmaxp > top.zpslmax[-1]-2*w3d.dz + top.zbeam:
        lforce = true
        if self.verbose:
          print "Load balancing since particles near right end of mesh ",top.zpslmax[-1],w3d.zmmaxglobal,zmaxp,top.zpslmax[-1]-2*w3d.dz

    # --- Find frequency of load balancing
    ii = max(self.when.values())
    for key,value in self.when.items():
      if top.it < key: ii = min(ii,value)

    # --- Just return is load balancing not done now.
    if not lforce and (top.it%ii) != 0: return

    if (top.it%ii) == 0 and self.verbose:
      print "Load balancing based on frequency"

    # --- On step zero, a complete reorganization is done so the reorg flag
    # --- is set to true to use the particle sorter which is more efficient
    # --- in that case.
    if reorg is None: reorg = (top.it==1)

    if doloadrho is None: doloadrho = self.doloadrho
    if dofs is None: dofs = self.dofs

    if ((zmaxp - zminp)/w3d.dz < 10 or not top.lmoments or
        top.laccumulate_zmoments):
      # --- If the particles only extend over a few grid cells, recalculate
      # --- the distribution on a finer grid to get better loadbalancing.
      # --- Also, calculate the distribution if the moments were not
      # --- calculated on the most recent step.
      pnumz = zeros(1001,'d')
      zmin = max(w3d.zmminglobal+top.zbeam,zminp-w3d.dz)
      zmax = min(w3d.zmmaxglobal+top.zbeam,zmaxp+w3d.dz)
      if zz is None:
        # --- If zz was not fetched above, then get it here.
        zz = getz(jslist=-1,gather=0)
      setgrid1d(len(zz),zz,1000,pnumz,zmin,zmax)
      pnumz = parallelsum(pnumz)
      dz = (zmax - zmin)/1000.
    else:
      # --- Otherwise use the already calculated z-moment
      pnumz = top.pnumz[:,-1]
      dz = w3d.dz

    # --- Calculate the right hand side padding.
    if self.padright is None:
      if not top.lmoments:
        vz = getvz(jslist=-1,gather=0)
        if len(vz) > 0: vzmaxp = max(vz)
        else:           vzmaxp = -largepos
        del vz
        vzmaxp = globalmax(vzmaxp)
      else:
        vzmaxp = top.vzmaxp
      if vzmaxp > 0.: padright = vzmaxp*top.dt*ii*2
      else:           padright = ii*w3d.dz
    else:             padright = self.padright
    padright = max(padright,self.padright)
    if self.verbose:
      print "Load balancing padright = ",padright

    loadbalanceparticles(doloadrho=doloadrho,dofs=dofs,
                         padright=padright,padleft=self.padleft,
                         reorg=reorg,pnumz=pnumz,zmin=w3d.zmminglobal,dz=dz,
                         zminp=zminp,zmaxp=zmaxp,verbose=self.verbose)
    getphiforparticles()

#########################################################################
#########################################################################
def setparticledomains(zslave,doloadrho=1,dofs=1,padleft=0.,padright=0.,reorg=0):
  """
Sets the particles domains from the input, zslave, in the same way as done
with top.zslave during the generate. This is only meant to be used after
that has already been done.
 - zslave: list of starting locations of the domains in grid cell units
 - doloadrho=1: when true, the charge density is redeposited
 - dofs=1: when true, the fields are recalculated
 - padleft=0, padright=0: extra space added on to leftmost and rightmost
                          domains (up to edge of system) in units of meters
 - reorg=0: when true, call reorg_particles which is fastest when particles
            are to be shifted across multiple processors, otherwise use
            zpartbnd which is fastest when particles are to be shifted only to
            nearest neighbors.
  """
  if not lparallel: return
  # --- It is assumed that the user supplied decomposition is specified
  # --- in the array zslave.

  # --- All values of zslave must be > 0.
  assert min(zslave[1:]-zslave[:-1]) > 0.,"The length of all particle domains must be > 0."

  # --- Set domain of each processor.
  for i in range(npes):
    top.zpslmin[i] = w3d.zmminglobal + zslave[i]*w3d.dz
    top.zpslmax[i] = w3d.zmminglobal + zslave[i+1]*w3d.dz

  top.zpslmin[0] = max(w3d.zmminglobal,top.zpslmin[0] - padleft)
  top.zpslmax[-1] = min(w3d.zmmaxglobal,top.zpslmax[-1] + padright)

  # --- Set iz and nz. This is done so that zmesh[izpslave] < zpslmin, and
  # --- zmesh[izpslave+nzpslave] > zpslmax.
  top.izpslave[:] = int((top.zpslmin - w3d.zmminglobal)/w3d.dz)
  top.nzpslave[:] = int((top.zpslmax - w3d.zmminglobal)/w3d.dz)-top.izpslave+1

  # --- Make sure that the last processor doesn't have grid cells
  # --- sticking out the end.
  if top.izpslave[-1]+top.nzpslave[-1] > w3d.nzfull:
    top.nzpslave[-1] = w3d.nzfull - top.izpslave[-1]

  # --- Adjust the Z data
  #_adjustz()

  # --- Reorganize the particles
  if reorg:
    reorgparticles()
  else:
    zpartbnd(w3d.zmmax,w3d.zmmin,w3d.dz,top.zgrid)

  # --- Update sizes of arrays for particles
  if(w3d.solvergeom == w3d.XYZgeom):
    w3d.nzp = top.nzpslave[me]
    w3d.zmminp = w3d.zmminglobal + top.izpslave[me]*w3d.dz
    w3d.zmmaxp = w3d.zmminp + w3d.nzp*w3d.dz
    gchange("Fields3dParticles")
  else:
    gchange_rhop_phip_rz()

  # --- Do some additional work if requested
  if doloadrho: loadrho()
  if dofs: fieldsol(0)


#########################################################################
def loadbalanceparticles(doloadrho=1,dofs=1,spread=1.,padleft=0.,padright=0.,
                         reorg=0,pnumz=None,zmin=None,dz=None,
                         zminp=None,zmaxp=None,verbose=0):
  """
Load balances the particles as evenly as possible. The load balancing is
based off of the data in top.pnumz which of course must already have
been calculated. The number density is assumed to vary linearly between
grid points.
 - doloadrho=1: when true, the charge density is redoposited
 - dofs=1: when true, the fields are recalculated
 - spread=1.: fraction of processors to spread the work among
 - padleft=0, padright=0: extra space added on to leftmost and rightmost
                          domains (up to edge of system) in units of meters
 - reorg=0: when true, call reorg_particles which is fastest when particles
            are to be shifted across multiple processors, otherwise use
            zpartbnd which is fastest when particles are to be shifted only to
            nearest neighbors.
 - pnumz=top.pznum: the particle distribution to base the load balancing on
 - zmin=None: optional z-minimum of the grid
 - dz=None: optional grid cell size
 - zminp,zmaxp=None: optional min and max of the region that must be included
                     in the decomposition
 - verbose=0: when true, prints out timing information
  """
  if not lparallel: return
  starttime = wtime()

  if pnumz is None:
    # --- Gather pnumz. Return if there is no data.
    pnumz = top.pnumz[:,-1]
    if max(pnumz) == 0.: return

  # --- Add fictitious data so that actual work is spread only to the
  # --- requested fraction of the processors.
  assert (0. < spread <= 1.),"spread must be between 0 and 1 or 1."
  avepnumz = ave(pnumz)
  pnumz[:] = pnumz + avepnumz*(1./spread - 1.)

  # --- Convert the number of particles to a decomposition
  zslave = decompose(pnumz,npes,lfullcoverage=0)

  # --- Scale to specified grid if zmin and/or dz input
  if dz is not None: zslave = zslave*dz/w3d.dz
  if zmin is not None: zslave = zslave + zmin
  if zminp is not None: zslave[0] = min(zslave[0],zminp)
  if zmaxp is not None: zslave[-1] = max(zslave[-1],zmaxp)

  # --- Apply the new domain decomposition.
  setparticledomains(zslave,doloadrho=doloadrho,dofs=dofs,
                     padleft=padleft,padright=padright,reorg=reorg)
  endtime = wtime()
  if verbose: print "Load balance time = ",endtime - starttime

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
  zslave = decompose(weight,npes,lfullcoverage=1)

  # --- Set domain of each processor.
  # --- This coding ensures that all of the processors have nzfsslave at least
  # --- 2 or greater and that the last processor isn't left with too few cells.
  # --- In cases where all of the processors have values only slightly
  # --- greater than 2, the last processor will likely end up with too few
  # --- cells because of the accumulation of rounding up. Find the processors
  # --- which have the largest amount of roundup and take away one of their
  # --- grid cells until the last processor has enough. Do the same for the
  # --- case where the last processor has too many.
  top.izfsslave[0] = 0
  top.nzfsslave[0] = max(nint(zslave[1]) + 1,2)
  for i in range(1,npes):
    top.izfsslave[i] = top.izfsslave[i-1] + top.nzfsslave[i-1] - 1
    top.nzfsslave[i] = max(nint(zslave[i+1]-zslave[i]) + 1,2)
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
def decompose(weight,npes,lfullcoverage=0):
  """
Converts a weight into the size of the domains.
 - weight: array of relative weights of the work done by each processor
 - npes: number of processors
 - lfullcoverage=0: when true, the domains cover the full extent of
                    the system
Returns an array of the same length which is the relative length of each
of the domains.
  """
  # --- Integrate weight, assuming linear variation between grid points
  nz = len(weight) - 1
  np = 0.5*weight[0] + sum(weight[1:-1]) + 0.5*weight[-1]
  npperpe = 1.*np/npes

  zslave = zeros(npes+1,'d')
  iz = 0
  if not lfullcoverage:
    # --- First first non-zero weight, making sure to check first cell too.
    while weight[iz] == 0. and weight[iz+1] == 0.: iz = iz + 1
  delta = 0.
  zslave[0] = iz
  for ip in xrange(1,npes):
    fract = 0.
    npint = 0.
    npnext = weight[iz  ]*((1.-delta)+0.5*(delta**2-1.)) + \
             weight[iz+1]*0.5*(1. - delta**2)
    # --- Get the remaining bit from the previous cell if it is not too much.
    if npnext < npperpe:
      fract = 1. - delta
      iz = iz + 1
      delta = 0.
      npint = npnext
    # --- Keep adding cells until the number per processor is reached.
    while npint + 0.5*(weight[iz]+weight[iz+1]) < npperpe:
      fract = fract + 1.
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
    fract = fract + delta - delta1
    zslave[ip] = zslave[ip-1] + fract

  # --- Set the end of the last domain
  if not lfullcoverage:
    # --- Find the last place with non-zero weight, and give the last processor
    # --- everything up to that point.
    for ii in xrange(iz,nz):
      if weight[ii] > 0.: zslave[-1] = ii+1
  else:
    zslave[-1] = nz
    
  return zslave

#########################################################################
def _adjustz():

  #---------------------------------------------------------------------------
  # --- Set the axial extent of each slaves domain to include
  # --- both the particle and field solve domain.
  if(w3d.solvergeom == w3d.XYZgeom):
   for i in range(npes):
    top.izslave[i] = top.izfsslave[i]
    top.nzslave[i] = top.nzfsslave[i]
    top.zmslmin[i] = top.izfsslave[i]*w3d.dz + w3d.zmminglobal
    top.zmslmax[i] = (top.izfsslave[i]+top.nzfsslave[i])*w3d.dz+w3d.zmminglobal

  #---------------------------------------------------------------------------
  # --- Reset local values
  w3d.nz     = top.nzfsslave[me]
  zpmin = w3d.zmminglobal + top.izpslave[me]*w3d.dz
  zpmax = (top.izpslave[me]+top.nzpslave[me])*w3d.dz + w3d.zmminglobal
  w3d.izfsmin = 0.
  w3d.izfsmax = top.nzfsslave[me]
  w3d.zmmin = top.zmslmin[me]
  w3d.zmmax = top.zmslmax[me]

  # --- Change the alocation of everything effected are reset the meshes.
  gchange("Fields3d")
  gchange("Z_Moments")
  gchange("Hist")
  w3d.zmesh[:] = w3d.zmmin + iota(0,w3d.nz)*w3d.dz
  
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



