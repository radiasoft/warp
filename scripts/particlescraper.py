"""
ParticleScraper: class for creating particle scraping
"""
from warp import *
from generateconductors import *

particlescraper_version = "$Id: particlescraper.py,v 1.11 2004/02/11 21:49:01 dave Exp $"
def particlescraperdoc():
  import particlescraper
  print particlescraper.__doc__


class ParticleScraper:
  """
Class for creating particle scraper for conductors
 - conductors: a conductor or list of conductors which act as particle scrapers
               Note that each conductor MUST have a unique id.
 - lsavecondid: when true, the id of the conductor where the particle is
                lost is save. The id is saved in the array top.pidlost[:,-1].
 - install=1: flag whether or not to install the scraper so that the scraping
              automatically happens every time step.
After an instance is created, additional conductors can be added by calling
the method registerconductors which takes either a conductor or a list of
conductors are an argument.
  """
  _scrapeid = None
  def __init__(self,conductors,lsavecondid=0,install=1): 
    # --- Create grid
    self.grid = Grid()
    # --- register any initial conductors
    self.conductors = []
    self.registerconductors(conductors)
    # --- Add extra column(s) to top.pid for work space. Only do this if this
    # --- is the first time an instance of this class is created or if the
    # --- the extra column(s) had been removed. Otherwise use the already
    # --- allocated column(s) in top.pid.
    if (ParticleScraper._scrapeid is None or
        top.npid-1 < ParticleScraper._scrapeid or
        (w3d.solvergeom != w3d.XYZgeom and
         top.npid-1 < ParticleScraper._scrapeid+1)):
      ParticleScraper._scrapeid = top.npid
      self.npid = top.npid
      top.npid = top.npid + 1
      if w3d.solvergeom != w3d.XYZgeom:
        self.pwork = top.npid
        top.npid = top.npid + 1
      top.npidmax = max(top.npidmax,top.npid)
      top.npidlostmax = max(top.npidlostmax,top.npidmax)
    else:
      self.npid = _scrapeid
      if w3d.solvergeom != w3d.XYZgeom:
        self.pwork = _scrapeid + 1
    # --- Make sure that npmaxi is set
    top.npmaxi = max(top.npmax,2)
    gchange("Particles")
    gchange("LostParticles")
    # --- If the conductor id where particles are lost is being saved,
    # --- need to turn on saving of lost particles.
    self.lsavecondid = lsavecondid
    if lsavecondid:
      top.lsavelostpart = true
    # --- Install the call to scrape particles if requested
    if install: self.installscraper()

  def installscraper(self):
    # --- Install the call to scrape particles
    installparticlescraper(self.scrapeall)

  def registerconductors(self,newconductors):
    if type(newconductors) is not ListType: newconductors = [newconductors]
    for c in newconductors:
      self.conductors.append(c)
      self.grid.getisinside(c)

  def scrapeall(self,clear=0):
    for js in xrange(top.ns):
      self.scrape(js)
      if clear or self.lsavecondid:
        processlostpart(js+1,top.clearlostpart,top.time,top.zbeam)
      if self.lsavecondid:
        self.savecondid(js)

  def scrape(self,js):
    dx = self.grid.dx
    dy = self.grid.dy
    dz = self.grid.dz
    nx = self.grid.nx
    ny = self.grid.ny
    nz = self.grid.nz
    xmin = self.grid.xmin
    xmax = self.grid.xmax
    ymin = self.grid.ymin
    ymax = self.grid.ymax
    zmin = self.grid.zmin
    zmax = self.grid.zmax
    isinside = self.grid.isinside

    # --- First, find any particles near a conductor
    i1 = top.ins[js] - 1
    i2 = top.ins[js] + top.nps[js] - 1
    top.pid[i1:i2,self.npid] = 0.
    xx = top.xp[i1:i2]
    yy = top.yp[i1:i2]
    zz = top.zp[i1:i2]
    if w3d.solvergeom == w3d.XYZgeom:
      getgrid3d(top.nps[js],xx,yy,zz,top.pid[i1:i2,self.npid],
                nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                w3d.l2symtry,w3d.l4symtry)
    else:
      top.pid[i1:i2,self.pwork] = sqrt(xx**2 + yy**2)
      xx = top.pid[i1:i2,self.pwork]
      yy = zeros(len(xx),'d')
      getgrid2d(top.nps[js],xx,zz,top.pid[i1:i2,self.npid],
                nx,nz,isinside,xmin,xmax,zmin,zmax)
    iscrape = compress(top.pid[i1:i2,self.npid]>0.,arange(i1,i2))
    if len(iscrape) == 0: return

    # --- Duplicate the particle list eight times, once for each corner.
    # --- The direction in which the corners lie depends on the symmetry
    # --- and the sign of the coordinate.
    iscrape = repeat(iscrape,8)
    nn = len(iscrape)
    xx = take(xx,iscrape-i1)
    yy = take(yy,iscrape-i1)
    zz = take(zz,iscrape-i1)
    if w3d.l4symtry: sx = sign(ones(nn),xx)
    else: sx = 1.
    if w3d.l4symtry or w3d.l2symtry: sy = sign(ones(nn),yy)
    else: sy = 1.
    xg = xmin+int(abs(xx-xmin)/dx)*dx + array(nn/8*[0.,dx,0.,dx,0.,dx,0.,dx])
    yg = ymin+int(abs(yy-ymin)/dy)*dy + array(nn/8*[0.,0.,dy,dy,0.,0.,dy,dy])
    zg = zmin+int(abs(zz-zmin)/dz)*dz + array(nn/8*[0.,0.,0.,0.,dz,dz,dz,dz])
    pp = zeros(nn,'d')

    # --- Get conductor id that particles are near
    if w3d.solvergeom == w3d.XYZgeom:
      getgridngp3d(nn,xg,yg,zg,pp,
                   nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                   w3d.l2symtry,w3d.l4symtry)
    else:
      getgridngp2d(nn,xg,zg,pp,nx,nz,isinside,xmin,xmax,zmin,zmax)

    # --- Loop over the conductors, removing particles inside of each.
    for c in self.conductors:
      ii = compress(pp == c.condid,arange(nn))
      if len(ii) == 0: continue
      xc = take(xx,ii)
      yc = take(yy,ii)
      zc = take(zz,ii)
      ic = take(iscrape,ii)
      ic = compress(c.isinside(xc,yc,zc).isinside,ic)
      # --- For particles which are inside, set gaminv to 0, the lost particle
      # --- flag
      put(top.gaminv,ic,0.)


  def savecondid(self,js):
    # --- First make sure there is extra space in the pidlost array.
    if top.npidlostmax < top.npidmax+1:
      top.npidlostmax = top.npidmax + 1
      gchange("LostParticles")

    # --- Much of this code is duplicated from scrape above so if it changes,
    # --- this should change as well.
    dx = self.grid.dx
    dy = self.grid.dy
    dz = self.grid.dz
    nx = self.grid.nx
    ny = self.grid.ny
    nz = self.grid.nz
    xmin = self.grid.xmin
    xmax = self.grid.xmax
    ymin = self.grid.ymin
    ymax = self.grid.ymax
    zmin = self.grid.zmin
    zmax = self.grid.zmax
    isinside = self.grid.isinside

    i1 = top.inslost[js] - 1
    i2 = top.inslost[js] + top.npslost[js] - 1
    xx = top.xplost[i1:i2]
    yy = top.yplost[i1:i2]
    zz = top.zplost[i1:i2]

    # --- Get the indices of all lost particles that havn't been localized
    # --- to a conductor.
    iscrape = compress(top.pidlost[i1:i2,-1]==0,arange(i1,i2))

    # --- Duplicate the particle list eight times, once for each corner.
    # --- The direction in which the corners lie depends on the symmetry
    # --- and the sign of the coordinate.
    iscrape = repeat(iscrape,8)
    nn = len(iscrape)
    xx = take(xx,iscrape-i1)
    yy = take(yy,iscrape-i1)
    zz = take(zz,iscrape-i1)
    if w3d.l4symtry: sx = sign(ones(nn),xx)
    else: sx = 1.
    if w3d.l4symtry or w3d.l2symtry: sy = sign(ones(nn),yy)
    else: sy = 1.
    xg = xmin+int(abs(xx-xmin)/dx)*dx + array(nn/8*[0.,dx,0.,dx,0.,dx,0.,dx])
    yg = ymin+int(abs(yy-ymin)/dy)*dy + array(nn/8*[0.,0.,dy,dy,0.,0.,dy,dy])
    zg = zmin+int(abs(zz-zmin)/dz)*dz + array(nn/8*[0.,0.,0.,0.,dz,dz,dz,dz])
    pp = zeros(nn,'d')

    # --- Get conductor id that particles are near
    if w3d.solvergeom == w3d.XYZgeom:
      getgridngp3d(nn,xg,yg,zg,pp,
                   nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                   w3d.l2symtry,w3d.l4symtry)
    else:
      getgridngp2d(nn,xg,zg,pp,nx,nz,isinside,xmin,xmax,zmin,zmax)

    # --- Loop over the conductors, removing particles inside of each.
    for c in self.conductors:
      ii = compress(pp == c.condid,arange(nn))
      if len(ii) == 0: continue
      xc = take(xx,ii)
      yc = take(yy,ii)
      zc = take(zz,ii)
      ic = take(iscrape,ii)
      ic = compress(c.isinside(xc,yc,zc).isinside,ic)
      # --- For particles which are inside, set pid to the id of the conductor
      # --- where the particle is lost.
      put(top.pidlost[:,-1],ic,c.condid)

