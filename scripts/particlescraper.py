"""
ParticleScraper: class for creating particle scraping
"""
from warp import *
from generateconductors import *

particlescraper_version = "$Id: particlescraper.py,v 1.7 2003/10/10 00:20:58 dave Exp $"
def particlescraperdoc():
  import particlescraper
  print particlescraper.__doc__

class ParticleScraper:
  """
Class for creating particle scraper for conductors
 - conductors: a conductor or list of conductors which act as particle scrapers
 - install=1: flag whether or not to install the scraper so that the scraping
              automatically happens every time step.
  """
  def __init__(self,conductors,install=1): 
    # --- Create grid
    self.grid = Grid()
    # --- register any initial conductors
    self.conductors = {}
    self.registerconductors(conductors)
    # --- Add extra column(s) to top.pid for work space.
    self.npid = top.npid
    top.npid = top.npid + 1
    if w3d.solvergeom != w3d.XYZgeom:
      self.pwork = top.npid
      top.npid = top.npid + 1
    top.npidmax = max(top.npidmax,top.npid)
    # --- Make sure that npmaxi is set
    top.npmaxi = max(top.npmax,2)
    gchange("Particles")
    gchange("LostParticles")
    # --- Install the call to scrape particles if requested
    if install: self.installscraper()

  def installscraper(self):
    # --- Install the call to scrape particles
    installparticlescraper(self.scrapeall)

  def registerconductors(self,conductors):
    if type(conductors) is not ListType: conductors = [conductors]
    for c in conductors:
      self.conductors[c.condid] = c
      self.grid.getisinside(c)

  def scrapeall(self,clear=0):
    for js in xrange(top.ns):
      self.scrape(js)
      if clear: processlostpart(js+1,top.clearlostpart,top.time,top.zbeam)

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
    xg = int(abs(xx-xmin)/dx)*dx + array(nn/8*[0.,dx,0.,dx,0.,dx,0.,dx])
    yg = int(abs(yy-ymin)/dy)*dy + array(nn/8*[0.,0.,dy,dy,0.,0.,dy,dy])
    zg = int(abs(zz-zmin)/dz)*dz + array(nn/8*[0.,0.,0.,0.,dz,dz,dz,dz])
    pp = zeros(nn,'d')

    # --- Get conductor id that particles are near
    if w3d.solvergeom == w3d.XYZgeom:
      getgridngp3d(nn,xg,yg,zg,pp,
                   nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                   w3d.l2symtry,w3d.l4symtry)
    else:
      getgridngp2d(nn,xg,zg,pp,nx,nz,isinside,xmin,xmax,zmin,zmax)

    # --- Loop over the conductors, removing particles inside of each.
    for cid,c in self.conductors.iteritems():
      ii = compress(pp == cid,arange(nn))
      if len(ii) == 0: continue
      xc = take(xx,ii)
      yc = take(yy,ii)
      zc = take(zz,ii)
      ic = take(iscrape,ii)
      ic = compress(c.isinside(xc,yc,zc).isinside,ic)
      # --- For particles which are inside, set gaminv to 0, the lost particle
      # --- flag
      put(top.gaminv,ic,0.)

