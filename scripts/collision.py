"""Defines classes to handle Coulomb collisions
"""
from warp import *
import time
collision_version = "$Id: collision.py,v 1.5 2008/04/24 21:33:01 dave Exp $"

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
  def __init__(self,ncint=1,loglambda=None,epvth=0.95,
                    nx=None,ny=None,nz=None,
                    xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None,
                    pboundxy=None,pbound0=None,pboundnz=None,
                    l2symtry=None,l4symtry=None,pbounds=None,
                    geometry='XYZgeom'):

    # --- ncint specifies how often the collision operator is applied.
    self.ncint = ncint

    # --- The standard loglambda collision coefficient. If it is not supplied,
    # --- then it is calculated automatically from averaged values.
    self.loglambda = loglambda

    # --- This is a weighting factor which is used to avoid a cooling instability
    # --- See Cohen et al PoP 13, 22705 (2006) after eq.(4).
    self.epvth = epvth

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
    if self.nx > 0: self.dx = (self.xmax - self.xmin)/self.nx
    else:           self.dx = 1.
    if self.ny > 0: self.dy = (self.ymax - self.ymin)/self.ny
    else:           self.dy = 1.
    if self.nz > 0: self.dz = (self.zmax - self.zmin)/self.nz
    else:           self.dz = 1.
    self.pboundxy = _default(pboundxy,top.pboundxy)
    self.pbound0 = _default(pbound0,top.pbound0)
    self.pboundnz = _default(pboundnz,top.pboundnz)
    self.l2symtry = _default(l2symtry,w3d.l2symtry)
    self.l4symtry = _default(l4symtry,w3d.l4symtry)
    self.pbounds = pbounds
    self.geometry = geometry
    self.setuppbounds()

    # --- Create a serial list of the interaction pairs.
    self.collisionpairs = []

    # --- Create the dictionary that will hold the test species that
    # --- will collide with each of the field species.
    # --- The keys of fielddict are the field species. Each associated value
    # --- is the list of test species to collide on the field species.
    self.fielddict = {}

    # --- This holds the initial values for the velocity distributions which are
    # --- used to avoid division by small numbers (near the plasma edge for example
    # --- when there is poor statistics) and with epvth.
    # --- They are calculated the first time they are needed.
    self.vthsqinit = {}

    # --- Turn the operator on.
    self.enabled = 0
    self.enable()
    self.timer = 0.

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
    if self.geometry == w3d.RZgeom:
      self.pbounds[0] = reflect
      self.pbounds[2] = reflect
      self.pbounds[3] = reflect
      if self.xmmin < 0.: self.xmmin = 0.
    elif self.geometry == w3d.XZgeom:
      self.pbounds[2] = reflect
      self.pbounds[3] = reflect

  # --------------------------------------------------------------------
  def addpair(self,testspecies,fieldspecies=None,mutual=1,loglambda=None):
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
    if isinstance(testspecies,Species): testspecies = testspecies.jslist[0]
    if isinstance(fieldspecies,Species): fieldspecies = fieldspecies.jslist[0]
    if loglambda is None: loglambda = self.loglambda

    # --- If fieldspecies is not specified, then the species is
    # --- colliding with itself.
    if fieldspecies is None: fieldspecies = testspecies

    # --- If the fieldspecies and testspecies are the same, then
    # --- turn off mutual (since it would be redundant).
    if fieldspecies == testspecies: mutual = 0

    # --- Add the collision pair
    self.processpair([testspecies,fieldspecies],loglambda=loglambda)

    # --- If mutual, then add the pair with the roles reversed.
    if mutual:
      self.processpair([fieldspecies,testspecies],loglambda=loglambda)

  # --------------------------------------------------------------------
  def processpair(self,pair,loglambda=None):
    # --- Add the pair to the serial list of pairs.
    self.collisionpairs.append(pair)

    # --- Get the list of test species associated with the field speices.
    # --- Create a new list in the dictionary is there is not one already.
    testlist = self.fielddict.setdefault(pair[1],[])

    # --- Add the test species if it is not already in the list.
    if pair[0] not in testlist: testlist.append([pair[0],loglambda])

  # --------------------------------------------------------------------
  def handlegridboundaries(self,grid):
    if self.nx > 0:
      if self.pbounds[0] == 1: grid[0,:,:,...] *= 2.
      if self.pbounds[0] == 2: grid[0,:,:,...] += grid[-1,:,:,...]
      if self.pbounds[1] == 1: grid[-1,:,:,...] *= 2.
      if self.pbounds[1] == 2: grid[-1,:,:,...] = grid[0,:,:,...]

    if self.ny > 0:
      if self.pbounds[2] == 1: grid[:,0,:,...] *= 2.
      if self.pbounds[2] == 2: grid[:,0,:,...] += grid[:,-1,:,...]
      if self.pbounds[3] == 1: grid[:,-1,:,...] *= 2.
      if self.pbounds[3] == 2: grid[:,-1,:,...] = grid[:,0,:,...]

    if self.nz > 0:
      if self.pbounds[4] == 1: grid[:,:,0,...] *= 2.
      if self.pbounds[4] == 2: grid[:,:,0,...] += grid[:,:,-1,...]
      if self.pbounds[5] == 1: grid[:,:,-1,...] *= 2.
      if self.pbounds[5] == 2: grid[:,:,-1,...] = grid[:,:,0,...]

  # --------------------------------------------------------------------
  def deposgrid3dvect(self,np,x,y,z,ux,uy,uz,w,velocitygrid,densitygrid):
    deposgrid3dvect(1,np,x,y,z,ux,uy,uz,w,
                    self.nx,self.ny,self.nz,
                    velocitygrid,densitygrid,
                    self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)

  def getgrid3d(self,np,x,y,z,data,grid):
    getgrid3d(np,x,y,z,data,self.nx,self.ny,self.nz,grid,
              self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax,
              self.l2symtry,self.l4symtry)

  def deposgrid3d(self,np,x,y,z,data,grid,count):
      deposgrid3d(1,np,x,y,z,data,self.nx,self.ny,self.nz,grid,count,
                  self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)

  # --------------------------------------------------------------------
  def getvthsq(self,species):

    # --- Get the particle data of the species.
    np = getn(js=species,gather=0)
    x = getx(js=species,gather=0)
    y = gety(js=species,gather=0)
    z = getz(js=species,gather=0)
    ux = getux(js=species,gather=0)
    uy = getuy(js=species,gather=0)
    uz = getuz(js=species,gather=0)
    if top.wpid > 0: w = getpid(js=species,id=top.wpid-1)
    else:            w = ones(np,'d')

    # --- Create the various arrays
    densitygrid = fzeros((1+self.nx,1+self.ny,1+self.nz),'d')
    velocitygrid = fzeros((1+self.nx,1+self.ny,1+self.nz,3),'d')

    # --- Calculate the velocity averages, as well as the density.
    self.deposgrid3dvect(np,x,y,z,ux,uy,uz,w,velocitygrid,densitygrid)
    self.handlegridboundaries(densitygrid)
    self.handlegridboundaries(velocitygrid)
    gridcount = where(densitygrid > 0.,densitygrid,1.)
    velocitygrid /= gridcount[...,NewAxis]
    densitygrid *= (top.pgroup.sw[species]/(self.dx*self.dy*self.dz))

    # --- Fetch the average velocity for the field particles, which is used
    # --- to calculate vthermal**2.
    uxbar = zeros(np,'d')
    uybar = zeros(np,'d')
    uzbar = zeros(np,'d')
    self.getgrid3d(np,x,y,z,uxbar,velocitygrid[...,0])
    self.getgrid3d(np,x,y,z,uybar,velocitygrid[...,1])
    self.getgrid3d(np,x,y,z,uzbar,velocitygrid[...,2])

    # --- Calculate the initial vthermal**2
    vthsq = ave(((ux - uxbar)**2 + (uy - uybar)**2 + (uz - uzbar)**2)/3.)
    return vthsq

  # --------------------------------------------------------------------
  def docollisions(self):

    # --- Only do the collisions every ncint steps.
    if top.it%self.ncint > 0: return

    starttime = time.clock()

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
      self.densitygrid = fzeros((1+self.nx,1+self.ny,1+self.nz),'d')
      self.velocitygrid = fzeros((1+self.nx,1+self.ny,1+self.nz,3),'d')
      self.vthsqgrid = fzeros((1+self.nx,1+self.ny,1+self.nz),'d')

      # --- Calculate the velocity averages, as well as the density.
      self.deposgrid3dvect(fnp,x,y,z,ux,uy,uz,w,self.velocitygrid,self.densitygrid)
      self.handlegridboundaries(self.densitygrid)
      self.handlegridboundaries(self.velocitygrid)
      gridcount = where(self.densitygrid > 0.,self.densitygrid,1.)
      self.velocitygrid /= gridcount[...,NewAxis]
      self.densitygrid *= (top.pgroup.sw[field]/(self.dx*self.dy*self.dz))

      # --- Fetch the average velocity for the field particles, which is used
      # --- to calculate vthermal**2.
      uxbar = zeros(fnp,'d')
      uybar = zeros(fnp,'d')
      uzbar = zeros(fnp,'d')
      self.getgrid3d(fnp,x,y,z,uxbar,self.velocitygrid[...,0])
      self.getgrid3d(fnp,x,y,z,uybar,self.velocitygrid[...,1])
      self.getgrid3d(fnp,x,y,z,uzbar,self.velocitygrid[...,2])

      # --- Calculate and fetch the average vthermal**2
      vthsq = ((ux - uxbar)**2 + (uy - uybar)**2 + (uz - uzbar)**2)/3.
      junk = fzeros(self.densitygrid.shape,'d') # count won't include w so is not used
      # --- Note that the vthermal**2 is weighted to be consistent with gridcount
      # --- as calculated above.
      self.deposgrid3d(fnp,x,y,z,w*vthsq,self.vthsqgrid,junk)
      self.handlegridboundaries(self.vthsqgrid)
      self.vthsqgrid /= gridcount

      # --- Calculate the initial vthermal**2 if needed.
      if field not in self.vthsqinit:
        self.vthsqinit[field] = ave(vthsq)

      # --- Now, loop over the test species, colliding each one against the field.
      for test,loglambda in testspecies:

        # --- Get the particle data of the test species. Note that the
        # --- collisions can always be done locally.
        tnp = getn(js=test,gather=0)
        x = getx(js=test,gather=0)
        y = gety(js=test,gather=0)
        z = getz(js=test,gather=0)
        ux = getux(js=test,gather=0)
        uy = getuy(js=test,gather=0)
        uz = getuz(js=test,gather=0)

        # --- Fetch the average velocity of the field particles at the test
        # --- particle locations.
        density = zeros(tnp,'d')
        uxbar = zeros(tnp,'d')
        uybar = zeros(tnp,'d')
        uzbar = zeros(tnp,'d')
        self.getgrid3d(tnp,x,y,z,density,self.densitygrid)
        self.getgrid3d(tnp,x,y,z,uxbar,self.velocitygrid[...,0])
        self.getgrid3d(tnp,x,y,z,uybar,self.velocitygrid[...,1])
        self.getgrid3d(tnp,x,y,z,uzbar,self.velocitygrid[...,2])

        # --- Calculate the initial vthermal**2 if needed. This is done in
        # --- a separate routine which calculates and subtracts the average
        # --- velocity of the test speices. Note that the uxbar etc
        # --- calculated above are from the field species and is not the
        # --- correct thing to subtract.
        if test not in self.vthsqinit:
          self.vthsqinit[test] = self.getvthsq(test)

        # --- Get the vthermal**2 of the field speices at the test particle
        # --- locations.
        vthsqfield = zeros(tnp,'d')
        self.getgrid3d(tnp,x,y,z,vthsqfield,self.vthsqgrid)

        # --- Calculate log(lambda) if needed.
        if loglambda is None:
          mf = top.pgroup.sm[field]
          fieldiselectron = (abs((mf-emass)/emass) < .5)
          mt = top.pgroup.sm[test]
          testiselectron = (abs((mt-emass)/emass) < .5)
          if testiselectron and fieldiselectron:
            # --- Use the equations from the NRL plasma formulary for
            # --- electron-electron collisions
            Te = vthsqfield*mf/jperev
            loglambda = (23.5 - log(sqrt(density*1.e-6)*Te**(-5./4.)) - 
                         sqrt(1.e-5 + (log(Te)-2.)**2/16.))
          else:
            q2 = (abs(top.pgroup.sq[test]*top.pgroup.sq[field])/echarge**2)
            qto3 = q2**(3./2.)
            if test == field:
              Te = vthsqfield*mf/jperev
            else:
              Te = self.getvthsq(test)*mt/jperev
            T3 = Te**(-3./2.)
            # --- This expression was grabbed from the scatterParticleGroup
            # --- routine of LSP.
            loglambda = 23.0 - log(1.0 + sqrt(2.*density*1.e-6)*qto3*T3)

        #   zpc = q2
        #   rpc = top.pgroup.sm[test]/top.pgroup.sm[field]
        #   mpc = (top.pgroup.sm[field]/
        #          (top.pgroup.sm[test]+top.pgroup.sm[field]))

        #   loglambda = 23.0 -
        #            log(1. + zpc*sqrt(M*zpc/Tm + N*zpc/Tn)/
        #                     ((Tm + Tn*rpc)*mpc + sem[m]*mpc*0.3333*vmmn))


            """
          mu = (top.pgroup.sm[test]*top.pgroup.sm[field]/
                (top.pgroup.sm[test] + top.pgroup.sm[field]))
          vthsqave = ave(vthsqfield)
          if test != field:
            # --- This averages the vthermal**2 of the test and field species.
            # --- Is this a good thing to do? Does it matter?
            vthsqtest = ave(((ux - uxbar)**2 + (uy - uybar)**2 + (uz - uzbar)**2)/3.)
            vthsqave = 0.5*(vthsqave + vthsqtest)
          b0 = (abs(top.pgroup.sq[test]*top.pgroup.sq[field])/
                (4.*pi*eps0*mu*vthsqave))
          # --- Note that that density used here is of the field species.
          # --- A better way may be an average of the field and test species,
          # --- but the density of the test species is not otherwise calculated
          # --- here.
          omegape = sqrt(ave(density)*echarge**2/(emass*eps0))
          lambdadb = sqrt(vthsqave)/omegape
          loglambda = log(lambdadb/b0)
            """
        else:
          loglambda = loglambda*ones(tnp,'d')

        # --- Now, apply the operator.
        langevincollisions3d(test==field,
                             tnp,ux,uy,uz,uxbar,uybar,uzbar,
                             density,vthsqfield,
                             top.pgroup.sq[test],top.pgroup.sq[field],
                             top.pgroup.sm[test],top.pgroup.sm[field],
                             self.vthsqinit[test],self.vthsqinit[field],
                             top.dt*self.ncint,loglambda,self.epvth)

    endtime = time.clock()
    self.timer += (endtime - starttime)

