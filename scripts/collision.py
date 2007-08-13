"""Defines classes to handle Coulomb collisions
"""
from warp import *
collision_version = "$Id: collision.py,v 1.1 2007/08/13 22:11:39 dave Exp $"

def collisiondoc():
  import collision
  print collision.__doc__

class LangevinCollisions(object):
  """
Implements a Langevin collision operator as described in
Manheimer, Lampe, Joyce, JCP 138, 563 (1997).
Also, see Rognlien and Cutler, Nuc Fusion 20, 1003 1980.
  """
  # --------------------------------------------------------------------
  def __init__(self,ncint=1,loglambda=20.,epvth=0.95,
                    vthsqinit1=None,vthsqinit2=None,
                    nx=None,ny=None,nz=None,
                    xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None,
                    pboundxy=None,pbound0=None,pboundnz=None,
                    l2symtry=None,l4symtry=None,pbounds=None):

    # --- ncint specifies how often the collision operator is applied.
    self.ncint = ncint

    # --- The standard loglambda collision coefficient.
    self.loglambda = loglambda

    # --- This is a weighting factor which is used to avoid a cooling instability
    # --- See Cohen et al PoP 13, 22705 (2006) after eq.(4).
    self.epvth = epvth

    # --- These are the initial values for the velocity distribution which are
    # --- used to avoid division by small numbers (near the plasma edge for example
    # --- when there is poor statistics) and with epvth.
    # --- For now, these must be supplied. Could the average be used instead of
    # --- the initial value?
    assert (vthsqinit1 is not None),'vthsqinit1 must be supplied'
    assert (vthsqinit2 is not None),'vthsqinit2 must be supplied'
    self.vthsqinit1 = vthsqinit1
    self.vthsqinit2 = vthsqinit2

    # --- Setup grid parameters
    _default = lambda x,d: (x,d)[x is None]
    self.nx = _default(nx,w3d.nx)
    self.ny = _default(ny,w3d.ny)
    self.nz = _default(nz,w3d.nz)
    self.xmin = _default(xmin,w3d.xmmin)
    self.xmax = _default(xmax,w3d.xmmax)
    self.ymin = _default(ymin,w3d.ymmin)
    self.ymax = _default(ymax,w3d.ymmax)
    self.zmin = _default(zmin,w3d.zmminlocal)
    self.zmax = _default(zmax,w3d.zmmaxlocal)
    self.pboundxy = _default(pboundxy,top.pboundxy)
    self.pbound0 = _default(pbound0,top.pbound0)
    self.pboundnz = _default(pboundnz,top.pboundnz)
    self.l2symtry = _default(l2symtry,w3d.l2symtry)
    self.l4symtry = _default(l4symtry,w3d.l4symtry)
    self.pbounds = pbounds
    self.setuppbounds()

    # --- Create a serial list of the interaction pairs.
    self.collisionpairs = []

    # --- Create the dictionary that will hold the test species that
    # --- will collide with each of the field species.
    # --- The keys of fielddict are the field species. Each associated value
    # --- is the list of test species to collide on the field species.
    self.fielddict = {}

    # --- Turn the operator on.
    self.enabled = 0
    self.enable()

  # --------------------------------------------------------------------
  def enable(self):
    "Enable the collision operator"
    if not self.enabled:
      self.enabled = 1
      installafterstep(self.docollisions)

  # --------------------------------------------------------------------
  def disable(self):
    "Disable the collision operator"
    if self.enabled:
      self.enabled = 0
      uninstallafterstep(self.docollisions)

  # --------------------------------------------------------------------
  def setuppbounds(self):
    # --- Note that this is more or less directly copied from fieldsolver
    if self.pbounds is not None: return
    self.pbounds = zeros(6,'l')
    self.pbounds[0] = self.pboundxy
    self.pbounds[1] = self.pboundxy
    self.pbounds[2] = self.pboundxy
    self.pbounds[3] = self.pboundxy
    self.pbounds[4] = self.pbound0
    self.pbounds[5] = self.pboundnz
    if self.l2symtry:
      self.pbounds[2] = reflect
      if self.pboundxy == periodic: self.pbounds[3] = reflect
      self.ymmin = 0.
    elif self.l4symtry:
      self.pbounds[0] = reflect
      self.pbounds[2] = reflect
      if self.pboundxy == periodic: self.pbounds[1] = reflect
      if self.pboundxy == periodic: self.pbounds[3] = reflect
      self.xmmin = 0.
      self.ymmin = 0.
    if self.solvergeom == w3d.RZgeom:
      self.pbounds[0] = reflect
      self.pbounds[2] = reflect
      self.pbounds[3] = reflect
      if self.xmmin < 0.: self.xmmin = 0.
    elif self.solvergeom == w3d.XZgeom:
      self.pbounds[2] = reflect
      self.pbounds[3] = reflect

  # --------------------------------------------------------------------
  def addpair(self,testspecies,fieldspecies=None,mutual=1):
    """
Add pairs of species that will collide against each other. The testspecies
is the one that is affected by colliding against the fieldspecies. If no
fieldspecies is given, then the collisions of the testspecies are against
itself. If the mutual flag is set, then the reverse collision also happens,
the fieldspecies is affected by collision against the testspecies.
 - testspecies:
 - fieldspecies=None:
 - mutual=1:
    """
    # --- If fieldspecies is not specified, then the species is
    # --- colliding with itself.
    if fieldspecies is None: fieldspecies = testspecies

    # --- If the fieldspecies and testspecies are the same, then
    # --- turn off mutual (since it would be redundant).
    if fieldspecies == testspecies: mutual = 0

    # --- Add the collision pair
    self.processpair([testspecies,fieldspecies])

    # --- If mutual, then add the pair with the roles reversed.
    if mutual:
      self.processpair([fieldspecies,testspecies])

  # --------------------------------------------------------------------
  def processpair(self,pair):
    # --- Add the pair to the serial list of pairs.
    self.collisionpairs.append(pair)

    # --- Get the list of test species associated with the field speices.
    # --- Create a new list in the dictionary is there is not one already.
    testlist = self.fielddict.setdefault(pair[1],[])

    # --- Add the test species if it is not already in the list.
    if pair[0] not in testlist: testlist.append(pair[0])

  # --------------------------------------------------------------------
  def handlegridboundaries(self,grid):
    if self.pbounds[0] == 1: grid[0,:,:,...] *= 2.
    if self.pbounds[0] == 2: grid[0,:,:,...] += grid[-1,:,:,...]
    if self.pbounds[1] == 1: grid[-1,:,:,...] *= 2.
    if self.pbounds[1] == 2: grid[-1,:,:,...] = grid[0,:,:,...]

    if self.pbounds[2] == 1: grid[:,0,:,...] *= 2.
    if self.pbounds[2] == 2: grid[:,0,:,...] += grid[:,-1,:,...]
    if self.pbounds[3] == 1: grid[:,-1,:,...] *= 2.
    if self.pbounds[3] == 2: grid[:,-1,:,...] = grid[:,0,:,...]

    if self.pbounds[4] == 1: grid[:,:,0,...] *= 2.
    if self.pbounds[4] == 2: grid[:,:,0,...] += grid[:,:,-1,...]
    if self.pbounds[5] == 1: grid[:,:,-1,...] *= 2.
    if self.pbounds[5] == 2: grid[:,:,-1,...] = grid[:,:,0,...]

  # --------------------------------------------------------------------
  def docollisions(self):

    # --- Only do the collisions every ncint steps.
    if top.it%self.ncint > 0: return

    # --- Loop over the field species, colliding all of the test species
    # --- against it. This way, the averages are only calculated once for
    # --- each field species.
    for field,testspecies in self.fielddict.items():

      # --- Get the particle data of the field species. Note that the
      # --- collisions can always be done locally.
      fnp = getn(js=field,gather=0)
      x = getx(js=field,gather=0)
      y = gety(js=field,gather=0)
      z = getz(js=field,gather=0)
      ux = getux(js=field,gather=0)
      uy = getuy(js=field,gather=0)
      uz = getuz(js=field,gather=0)
      if top.wpid > 0: w = getpid(js=field,id=top.wpid-1)
      else:            w = ones(fnp,'d')

      # --- Create the various arrays
      self.densitygrid = fzeros((1+nx,1+ny,1+nz),'d')
      self.velocitygrid = fzeros((1+nx,1+ny,1+nz,3),'d')
      self.temperaturegrid = fzeros((1+nx,1+ny,1+nz),'d')

      # --- Calculate the velocity averages, as well as the density.
      deposgrid3dvect(1,fnp,x,y,z,ux,uy,uz,w,
                      self.nx,self.ny,self.nz,
                      self.velocitygrid,self.densitygrid,
                      self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)
      self.handlegridboundaries(self.densitygrid)
      self.handlegridboundaries(self.velocitygrid)
      gridcount = where(self.densitygrid > 0.,self.densitygrid,1.)
      self.velocitygrid /= gridcount[...,NewAxis]

      # --- Fetch the average velocity for the field particles, which is used
      # --- to calculate the temperature.
      uxbar = zeros(fnp,'d')
      uybar = zeros(fnp,'d')
      uzbar = zeros(fnp,'d')
      getgrid3d(fnp,x,y,z,uxbar,self.nx,self.ny,self.nz,self.velocitygrid[...,0],
                self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                l2symtry,l4symtry)
      getgrid3d(fnp,x,y,z,uybar,self.nx,self.ny,self.nz,self.velocitygrid[...,1],
                self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                l2symtry,l4symtry)
      getgrid3d(fnp,x,y,z,uzbar,self.nx,self.ny,self.nz,self.velocitygrid[...,2],
                self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                l2symtry,l4symtry)

      # --- Calculate and fetch the average temperature.
      temperature = ((ux - uxbar)**2 + (uy - uybar)**2 + (uz - uzbar)**2)/3.
      junk = fzeros(self.densitygrid,'d') # --- count won't include w so is not used
      # --- Note that the temperature is weighted to be consistent with gridcount
      # --- as calculated above.
      deposgrid3d(1,fnp,x,y,z,w*temperature,
                  self.nx,self.ny,self.nz,self.temperaturegrid,junk,
                  self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)
      self.handlegridboundaries(self.temperaturegrid)
      self.temperaturegrid /= gridcount

      # --- Now, loop over the test species, colliding each one against the field.
      for test in testspecies:

        # --- Get the particle data of the test species. Note that the
        # --- collisions can always be done locally.
        tnp = getn(js=test,gather=0)
        x = getx(js=test,gather=0)
        y = gety(js=test,gather=0)
        z = getz(js=test,gather=0)
        ux = getux(js=test,gather=0)
        uy = getuy(js=test,gather=0)
        uz = getuz(js=test,gather=0)

        # --- Fetch the average velocity for the test particles
        density = zeros(tnp,'d')
        uxbar = zeros(tnp,'d')
        uybar = zeros(tnp,'d')
        uzbar = zeros(tnp,'d')
        getgrid3d(tnp,x,y,z,density,self.nx,self.ny,self.nz,self.densitygrid,
                  self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                  l2symtry,l4symtry)
        getgrid3d(tnp,x,y,z,uxbar,self.nx,self.ny,self.nz,self.velocitygrid[...,0],
                  self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                  l2symtry,l4symtry)
        getgrid3d(tnp,x,y,z,uybar,self.nx,self.ny,self.nz,self.velocitygrid[...,1],
                  self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                  l2symtry,l4symtry)
        getgrid3d(tnp,x,y,z,uzbar,self.nx,self.ny,self.nz,self.velocitygrid[...,2],
                  self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                  l2symtry,l4symtry)

        # --- Get the temperature, and convert to thermal velocity.
        temperature = ((ux - uxbar)**2 + (uy - uybar)**2 + (uz - uzbar)**2)/3.
        getgrid3d(tnp,x,y,z,temperature,self.nx,self.ny,self.nz,self.temperaturegrid,
                  self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
                  l2symtry,l4symtry)
        vthsq = sqrt(temperature)

        # --- Now, apply the operator.
        langevincollisions2d(test==field,
                             tnp,ux,uy,uz,uxbar,uybar,uzbar,
                             density,vthsq,
                             top.pgroup.sq[test],top.pgroup.sq[field],
                             top.pgroup.sm[test],top.pgroup.sm[field],
                             vthsqinit1,vthsqinit2,
                             top.dt*self.ncint,self.loglambda,self.epvth)

