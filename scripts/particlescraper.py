"""
ParticleScraper: class for creating particle scraping
"""
from warp import *
from generateconductors import *

particlescraper_version = "$Id: particlescraper.py,v 1.17 2004/06/03 20:57:53 dave Exp $"
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
 - lsaveintercept: when true, the location and surface normal where the
                   particle intercepted the conductor surface is calculated.
                   The location is overwritten onto the xplost, yplost, and
                   zplost arrays. The angles describing the surface normal are
                   put into pidlost[:,-3] and pidlost[:,-2]. The spherical
                   coordinate angle theta is in -3, and phi is in -2.
                   The time at which the particles are lost is put into
                   pidlost[:,-4].
                   Note that the condid where the particle is lost is also
                   saved in pidlost[:,-1].
 - mglevel=0: Coarsening level for index grid which is used to determine
              which conductors particles are near. This grid is a full size,
              3d (or 2d) array and can require a not insignificant amount of
              time to compute. If it is expected that few particles will be
              lost, using a coarser grid can substantially reduce the memory
              and time requirements for this grid. However, a coarser grid
              will mean that more particles will be flagged as being near the
              conductors, and the detailed check required for these particles
              will be more expensive.  A trade off is required that can really
              only be optimized empirically. A value of 0, 1, or 2 are
              probably optimal.
 - install=1: flag whether or not to install the scraper so that the scraping
              automatically happens every time step.
After an instance is created, additional conductors can be added by calling
the method registerconductors which takes either a conductor or a list of
conductors are an argument.
  """
  def __init__(self,conductors,lsavecondid=0,lsaveintercept=0,mglevel=0,
                    install=1): 
    self.mglevel = mglevel
    # --- Don't create the grid until it is needed.
    self.grid = None
    # --- register any initial conductors
    self.conductors = []
    self.registerconductors(conductors)
    ## --- Make sure that npmaxi is set
    #top.npmaxi = max(top.npmax,2)
    #gchange("Particles")
    gchange("LostParticles")
    # --- If the conductor id where particles are lost is being saved,
    # --- need to turn on saving of lost particles.
    self.lsaveintercept = lsaveintercept
    self.lsavecondid = lsavecondid or lsaveintercept
    if self.lsavecondid:
      top.lsavelostpart = true
    # --- Install the call to scrape particles if requested
    if install: self.installscraper()

  def installscraper(self):
    # --- Install the call to scrape particles
    installparticlescraper(self.scrapeall)

  def disable(self):
    if isinstalledparticlescraper(self.scrapeall):
      uninstallparticlescraper(self.scrapeall)

  def registerconductors(self,newconductors):
    if self.grid is None: self.grid = Grid()
    if type(newconductors) is not ListType: newconductors = [newconductors]
    for c in newconductors:
      self.conductors.append(c)
      self.grid.getisinside(c,mglevel=self.mglevel)

  def unregisterconductors(self,conductor,nooverlap=0):
    self.conductors.remove(conductor)
    if not nooverlap:
      # --- This is horribly inefficient!!!
      self.grid.resetgrid()
      for c in self.conductors:
        self.grid.getisinside(c,mglevel=self.mglevel)
    else:
      self.grid.removeisinside(conductor)

  def scrapeall(self,clear=0):
    if len(self.conductors) == 0: return
    for js in xrange(top.ns):
      self.scrape(js)
      if clear or self.lsavecondid:
        processlostpart(js+1,top.clearlostpart,top.time,top.zbeam)
      if self.lsavecondid:
        self.savecondid(js)

  def scrape(self,js):
    dx,dy,dz,nx,ny,nz,iz = self.grid.getmeshsize(self.mglevel)
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
    xx = top.xp[i1:i2]
    yy = top.yp[i1:i2]
    zz = top.zp[i1:i2]
    pp = zeros(top.nps[js],'d')
    if w3d.solvergeom == w3d.XYZgeom:
      getgrid3d(top.nps[js],xx,yy,zz,pp,
                nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                w3d.l2symtry,w3d.l4symtry)
    else:
      xx = sqrt(xx**2 + yy**2)
      yy = zeros(len(xx),'d')
      getgrid2d(top.nps[js],xx,zz,pp,nx,nz,isinside,xmin,xmax,zmin,zmax)

    iscrape = compress(pp>0.,arange(i1,i2))
    if len(iscrape) == 0: return

    # --- Duplicate the particle list eight times, once for each corner.
    iscrape = repeat(iscrape,8)
    nn = len(iscrape)
    xx = take(xx,iscrape-i1)
    yy = take(yy,iscrape-i1)
    zz = take(zz,iscrape-i1)
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

    # --- Just return if there are no lost particles.
    if top.npslost[js] == 0: return

    # --- First make sure there is extra space in the pidlost array.
    pidspace = 1
    if self.lsaveintercept: pidspace = 4
    if top.npidlostmax < top.npidmax+pidspace:
      top.npidlostmax = top.npidmax + pidspace
      gchange("LostParticles")

    # --- Much of this code is duplicated from scrape above so if it changes,
    # --- this should change as well.
    dx,dy,dz,nx,ny,nz,iz = self.grid.getmeshsize(self.mglevel)
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

    if w3d.solvergeom == w3d.RZgeom:
      xx = sqrt(xx**2 + yy**2)
      yy = zeros(len(xx),'d')

    # --- Get the indices of all lost particles that havn't been localized
    # --- to a conductor.
    iscrape = compress(top.pidlost[i1:i2,-1]==0,arange(i1,i2))

    # --- Duplicate the particle list eight times, once for each corner.
    iscrape = repeat(iscrape,8)
    nn = len(iscrape)
    x8 = take(xx,iscrape-i1)
    y8 = take(yy,iscrape-i1)
    z8 = take(zz,iscrape-i1)
    xg = xmin+int(abs(x8-xmin)/dx)*dx + array(nn/8*[0.,dx,0.,dx,0.,dx,0.,dx])
    yg = ymin+int(abs(y8-ymin)/dy)*dy + array(nn/8*[0.,0.,dy,dy,0.,0.,dy,dy])
    zg = zmin+int(abs(z8-zmin)/dz)*dz + array(nn/8*[0.,0.,0.,0.,dz,dz,dz,dz])
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
      xc = take(x8,ii)
      yc = take(y8,ii)
      zc = take(z8,ii)
      ic = take(iscrape,ii)
      ic = compress(c.isinside(xc,yc,zc).isinside,ic)
      # --- For particles which are inside, set pid to the id of the conductor
      # --- where the particle is lost.
      put(top.pidlost[:,-1],ic,c.condid)
      # --- Save location and surface normal where particle intercepted the
      # --- conductor.
      if self.lsaveintercept:
        xc = take(xx,ic-i1)
        yc = take(yy,ic-i1)
        zc = take(zz,ic-i1)
        if top.lrelativ:
          vx = take(top.uxplost[i1:i2],ic-i1)*take(top.gaminvlost[i1:i2],ic-i1)
          vy = take(top.uyplost[i1:i2],ic-i1)*take(top.gaminvlost[i1:i2],ic-i1)
          vz = take(top.uzplost[i1:i2],ic-i1)*take(top.gaminvlost[i1:i2],ic-i1)
        else:
          vx = take(top.uxplost[i1:i2],ic-i1)
          vy = take(top.uyplost[i1:i2],ic-i1)
          vz = take(top.uzplost[i1:i2],ic-i1)
        intercept = c.intercept(xc,yc,zc,vx,vy,vz)
        put(top.xplost,ic,intercept.xi)
        put(top.yplost,ic,intercept.yi)
        put(top.zplost,ic,intercept.zi)
        put(top.pidlost[:,-3],ic,intercept.itheta)
        put(top.pidlost[:,-2],ic,intercept.iphi)
        dt = (sqrt((xc-xi)**2 + (yc-yi)**2 + (zc-zi)**2)/
              dvnz(sqrt(vx**2 + vy**2 + vz**2)))
        put(top.pidlost[:,-4],ic,top.time - dt)



