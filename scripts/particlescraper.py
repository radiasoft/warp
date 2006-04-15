"""
ParticleScraper: class for creating particle scraping
"""
from warp import *
from generateconductors import *
import timing as t

particlescraper_version = "$Id: particlescraper.py,v 1.39 2006/04/15 00:13:37 dave Exp $"
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
 - lcollectlpdata: When true, the lost particles statistics will be collected for 
                   each conductor in the list lostparticles_data (Assembly class).
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
 - grid=None: A instance of the Grid class can be supplied, allowing control
              over the region where the scraping is done and the resolution
              of the scraping data.
After an instance is created, additional conductors can be added by calling
the method registerconductors which takes either a conductor or a list of
conductors are an argument.
  """
  def __init__(self,conductors,lsavecondid=0,lsaveintercept=0,lcollectlpdata=0,
                    mglevel=0,aura=0.,install=1,grid=None): 
    self.mglevel = mglevel
    self.aura = aura
    # --- Remember if the user specified the grid.
    self.usergrid = (grid is not None)
    # --- Don't create the grid until it is needed.
    self.grid = grid
    # --- register any initial conductors
    self.conductors = []
    self.registerconductors(conductors)
    # --- Allocate arrays for lost particles
    gchange("LostParticles")
    # --- If the conductor id where particles are lost is being saved,
    # --- need to turn on saving of lost particles.
    self.lsaveintercept = lsaveintercept
    self.lsavecondid = lsavecondid or lsaveintercept or lcollectlpdata
    self.lcollectlpdata = lcollectlpdata
    if self.lsavecondid:
      top.lsavelostpart = true
    if self.lsaveintercept:
      # --- Note that nextpid returns numbers based on 1 based indexing
      self.xoldpid=nextpid() - 1
      self.yoldpid=nextpid() - 1
      self.zoldpid=nextpid() - 1
      gchangeparticles()
      installbeforestep(self.saveoldpos)
    self.l_print_timing=0
    # --- Install the call to scrape particles if requested
    if install: self.installscraper()

  def installscraper(self):
    # --- Install the call to scrape particles
    installparticlescraper(self.scrapeall)

  def disable(self):
    if isinstalledparticlescraper(self.scrapeall):
      uninstallparticlescraper(self.scrapeall)

  def __setstate__(self,dict):
    # --- This is called when the instance is unpickled.
    # --- WARNING!!! When an instance in unpickled, the conductors referrenced
    # --- will be copies of the original conductors passed in - the restarted
    # --- run will have two copies, the original one (if it still exists) and
    # --- a copy refered to by the instances conductor list.
    self.__dict__.update(dict)
    self.installscraper()

  def registerconductors(self,newconductors):
    self.updategrid()
    if type(newconductors) is not ListType: newconductors = [newconductors]
    for c in newconductors:
      assert c.condid != 0,"The conductor id must be nonzero in order for the particle scraping to work."
      self.conductors.append(c)
      self.grid.getisinside(c,mglevel=self.mglevel,aura=self.aura)

  def unregisterconductors(self,conductor,nooverlap=0):
    self.conductors.remove(conductor)
    if not nooverlap:
      # --- This is horribly inefficient!!!
      self.grid.resetgrid()
      self.updateconductors()
    else:
      self.grid.removeisinside(conductor)
      
  def updategrid(self,lforce=0):
    """Update the grid to match any changes to the underlying grid, for example
after load balancing."""
    if self.grid is None: lforce = 1
    if self.usergrid and not lforce: return
    if lparallel: nz = top.nzpslave[me]
    else:         nz = w3d.nz
    if (not lforce and (self.grid.nx == w3d.nx and
                        self.grid.ny == w3d.ny and
                        self.grid.nz ==     nz and
                        self.grid.xmmin == w3d.xmmin and
                        self.grid.xmmax == w3d.xmmax and
                        self.grid.ymmin == w3d.ymmin and
                        self.grid.ymmax == w3d.ymmax and
                        self.grid.zmmin == w3d.zmminglobal and
                        self.grid.zmmax == w3d.zmmaxglobal and
                        self.grid.izslave[me] == top.izpslave[me])): return
    # --- Note that copies of the slave arrays are passed in.
    # --- The arrays in top may be changed the next time loadbalancing is
    # --- done, but the arrays in self.grid should not be changed. Instead,
    # --- a whole new grid is created.
    self.grid = Grid(nz=nz,
                     izslave=top.izpslave.copy(),nzslave=top.nzpslave.copy())
    self.updateconductors()

  def updateconductors(self):
    for c in self.conductors:
      self.grid.getisinside(c,mglevel=self.mglevel,aura=self.aura)

  def saveoldpos(self):
    for js in xrange(top.ns):
      if top.nps[js]>0:
        i1 = top.ins[js] - 1
        i2 = top.ins[js] + top.nps[js] - 1
        top.pid[i1:i2,self.xoldpid]=top.xp[i1:i2].copy()
        top.pid[i1:i2,self.yoldpid]=top.yp[i1:i2].copy()
        top.pid[i1:i2,self.zoldpid]=top.zp[i1:i2].copy()
      
  def scrapeall(self,clear=0):
    if len(self.conductors)==0 or parallelsum(sum(top.nps))==0: return
    self.updategrid()
    for js in xrange(top.ns):
      if top.it%top.ndts[js]==0:
        if self.l_print_timing:tt=0.
        if self.l_print_timing:t.start()
        self.scrape(js)
        if self.l_print_timing:t.finish()
        if self.l_print_timing:print js,'scrape',t.milli()
        if self.l_print_timing:t.start()
        if clear or self.lsavecondid:
          processlostpart(top.pgroup,js+1,top.clearlostpart,top.time,top.zbeam)
        if self.l_print_timing:t.finish()
        if self.l_print_timing:print js,'processlosspart',t.milli()
        if self.l_print_timing:t.start()
        if self.lsavecondid:
          self.savecondid(js)
        if self.l_print_timing:t.finish()
        if self.l_print_timing:print js,'savecondid',t.milli()

  def scrape(self,js):
    self.scrape2(js)

  def scrape2(self,js):
    t.start()
    if top.nps[js] == 0: return
    dx,dy,dz,nx,ny,nz,iz = self.grid.getmeshsize(self.mglevel)
    xmin = self.grid.xmin
    xmax = self.grid.xmax
    ymin = self.grid.ymin
    ymax = self.grid.ymax
    zmin = self.grid.zmmin + iz*dz + top.zbeam
    zmax = self.grid.zmmin + (iz+nz)*dz + top.zbeam
    isinside = self.grid.isinside

    # --- First, find any particles near a conductor
    i1 = top.ins[js] - 1
    i2 = top.ins[js] + top.nps[js] - 1
    xx = top.xp[i1:i2]
    yy = top.yp[i1:i2]
#    if js==1:print js,i1,i2,top.zp[i1:i2],top.zbeam
    zz = top.zp[i1:i2]
    pp = zeros(top.nps[js],'d')
    if w3d.solvergeom in [w3d.XYZgeom,w3d.XYZgeomMR]:
      getgrid3d(top.nps[js],xx,yy,zz,pp,
                nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                w3d.l2symtry,w3d.l4symtry)
    elif w3d.solvergeom == w3d.RZgeom:
      xx = sqrt(xx**2 + yy**2)
      yy = zeros(len(xx),'d')
      getgrid2d(top.nps[js],xx,zz,pp,nx,nz,isinside[:,0,:],xmin,xmax,zmin,zmax)
    elif w3d.solvergeom == w3d.XZgeom:
      getgrid2d(top.nps[js],xx,zz,pp,nx,nz,isinside[:,0,:],xmin,xmax,zmin,zmax)
    elif w3d.solvergeom == w3d.XYgeom:
      getgrid2d(top.nps[js],xx,yy,pp,nx,ny,isinside[:,:,0],xmin,xmax,ymin,ymax)
    else:
      raise "The particle scraping only works for XYZ, XY and RZ geometry"

    iscrape = compress(pp>0.,arange(i1,i2))
    if len(iscrape) == 0: return
 
    if w3d.solvergeom in [w3d.XYZgeom,w3d.XYZgeomMR]:
      nd=3
      gdx = [0.,dx,0.,dx,0.,dx,0.,dx]
      gdy = [0.,0.,dy,dy,0.,0.,dy,dy]
      gdz = [0.,0.,0.,0.,dz,dz,dz,dz]
      xx = take(xx,iscrape-i1)
      yy = take(yy,iscrape-i1)
      zz = take(zz,iscrape-i1)
      xg = xmin+int(abs(xx-xmin)/dx)*dx 
      yg = ymin+int(abs(yy-ymin)/dy)*dy 
      zg = zmin+int(abs(zz-zmin)/dz)*dz 
    elif w3d.solvergeom in [w3d.XZgeom,w3d.RZgeom]:
      nd=2
      gdx = [0.,dx,0.,dx]
      gdz = [0.,0.,dz,dz]
      xx = take(xx,iscrape-i1)
      zz = take(zz,iscrape-i1)
      xg = xmin+int(abs(xx-xmin)/dx)*dx 
      zg = zmin+int(abs(zz-zmin)/dz)*dz 
    elif w3d.solvergeom == w3d.XYgeom:
      nd=2
      gdx = [0.,dx,0.,dx]
      gdy = [0.,0.,dy,dy]
      xx = take(xx,iscrape-i1)
      yy = take(yy,iscrape-i1)
      xg = xmin+int(abs(xx-xmin)/dx)*dx 
      yg = ymin+int(abs(yy-ymin)/dy)*dy 
    
    nn = len(iscrape)
    pp = zeros(nn,'d')

    for i in range(2**nd):

      # --- Get conductor id that particles are near
      if w3d.solvergeom in [w3d.XYZgeom,w3d.XYZgeomMR]:
        getgridngp3d(nn,xg+gdx[i],yg+gdy[i],zg+gdz[i],pp,
                     nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,0.,
                     w3d.l2symtry,w3d.l4symtry)
      elif w3d.solvergeom in [w3d.XZgeom,w3d.RZgeom]:
        getgridngp2d(nn,xg+gdx[i],zg+gdz[i],pp,nx,nz,isinside[:,0,:],xmin,xmax,zmin,zmax)
      elif w3d.solvergeom == w3d.XYgeom:
        getgridngp2d(nn,xg+gdx[i],yg+gdy[i],pp,nx,ny,isinside[:,:,0],xmin,xmax,ymin,ymax)

      # --- Loop over the conductors, removing particles inside of each.
      for c in self.conductors:
        ixyz=arange(nn)
        # get indices of particles that scraped in conductor c
        ii = compress(pp == c.condid,ixyz) 
        # if no particle scraped in conductor c, then check next conductor
        if len(ii) == 0: continue          
        # get positions of scraped particles
        xc = take(xx,ii)
        yc = take(yy,ii)
        zc = take(zz,ii)
        # down-select indices of particles that are inside conductor c
        iic = compress(c.isinside(xc,yc,zc).isinside,ii)
        # if no particle are inside conductor c, then check next conductor
        if len(iic) == 0: continue
        # get indices (in lost particles arrays) of particles inside c
        ic = take(iscrape,iic)
        # --- For particles which are inside, set gaminv to 0, the lost particle
        # --- flag
        put(top.gaminv,ic,0.)
        # remove scraped particles from list of lost particles
        put(iscrape,iic,-1)
        iscrape = compress(iscrape>=0,iscrape)        
        nn = len(iscrape)
        if nn == 0: return
        put(ixyz,iic,-1)
        ixyz = compress(ixyz>=0,ixyz)        
        xx = take(xx,ixyz)
        xg = take(xg,ixyz)
        pp = take(pp,ixyz)
        if w3d.solvergeom in [w3d.XYZgeom,w3d.XYZgeomMR,w3d.XYgeom]:
          yy = take(yy,ixyz)
          yg = take(yg,ixyz)
        if w3d.solvergeom in [w3d.XYZgeom,w3d.XYZgeomMR,w3d.XZgeom,w3d.RZgeom]:
          zz = take(zz,ixyz)
          zg = take(zg,ixyz)

  def scrape1(self,js):
    t.start()
    if top.nps[js] == 0: return
    dx,dy,dz,nx,ny,nz,iz = self.grid.getmeshsize(self.mglevel)
    xmin = self.grid.xmin
    xmax = self.grid.xmax
    ymin = self.grid.ymin
    ymax = self.grid.ymax
    zmin = self.grid.zmmin + iz*dz + top.zbeam
    zmax = self.grid.zmmin + (iz+nz)*dz + top.zbeam
    isinside = self.grid.isinside

    # --- First, find any particles near a conductor
    i1 = top.ins[js] - 1
    i2 = top.ins[js] + top.nps[js] - 1
    xx = top.xp[i1:i2]
    yy = top.yp[i1:i2]
#    if js==1:print js,i1,i2,top.zp[i1:i2],top.zbeam
    zz = top.zp[i1:i2]
    pp = zeros(top.nps[js],'d')
    if w3d.solvergeom in [w3d.XYZgeom,w3d.XYZgeomMR]:
      getgrid3d(top.nps[js],xx,yy,zz,pp,
                nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                w3d.l2symtry,w3d.l4symtry)
    elif w3d.solvergeom == w3d.RZgeom:
      xx = sqrt(xx**2 + yy**2)
      yy = zeros(len(xx),'d')
      getgrid2d(top.nps[js],xx,zz,pp,nx,nz,isinside[:,0,:],xmin,xmax,zmin,zmax)
    elif w3d.solvergeom == w3d.XYgeom:
      getgrid2d(top.nps[js],xx,yy,pp,nx,ny,isinside[:,:,0],xmin,xmax,ymin,ymax)
    else:
      raise "The particle scraping only works for XYZ and RZ geometry"

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
    if w3d.solvergeom in [w3d.XYZgeom,w3d.XYZgeomMR]:
      getgridngp3d(nn,xg,yg,zg,pp,
                   nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,0.,
                   w3d.l2symtry,w3d.l4symtry)
    elif w3d.solvergeom == w3d.RZgeom:
      getgridngp2d(nn,xg,zg,pp,nx,nz,isinside[:,0,:],xmin,xmax,zmin,zmax)
    elif w3d.solvergeom == w3d.XYgeom:
      getgridngp2d(nn,xg,yg,pp,nx,ny,isinside[:,:,0],xmin,xmax,ymin,ymax)

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
    if top.npslost[js] == 0: 
      if self.lcollectlpdata:
        for c in self.conductors:
          w=parallelsum(0.)
          if me==0 and w<>0.:
            c.lostparticles_data += [[top.time, 
                                      w*top.sq[js]*top.sw[js],
                                      top.dt,
                                      js]]
      return

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
    zmin = self.grid.zmmin + iz*dz + top.zbeam
    zmax = self.grid.zmmin + (iz+nz)*dz + top.zbeam
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
    if self.lcollectlpdata:iscrape1=iscrape.copy()

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
    if w3d.solvergeom in [w3d.XYZgeom,w3d.XYZgeomMR]:
      getgridngp3d(nn,xg,yg,zg,pp,
                   nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,0.,
                   w3d.l2symtry,w3d.l4symtry)
    elif w3d.solvergeom == w3d.RZgeom:
      getgridngp2d(nn,xg,zg,pp,nx,nz,isinside[:,0,:],xmin,xmax,zmin,zmax)
    elif w3d.solvergeom == w3d.XYgeom:
      getgridngp2d(nn,xg,yg,pp,nx,ny,isinside[:,:,0],xmin,xmax,ymin,ymax)
    else:
      raise "The particle scraping only works for XYZ and RZ geometry"


    if w3d.solvergeom == w3d.RZgeom:
      xx = top.xplost[i1:i2]
      yy = top.yplost[i1:i2]
      x8 = take(xx,iscrape-i1)
      y8 = take(yy,iscrape-i1)

    # --- Loop over the conductors, removing particles inside of each.
    for c in self.conductors:
      ii = compress(pp == c.condid,arange(nn))
      if len(ii) == 0: 
        if self.lcollectlpdata:
          w=parallelsum(0.)
        continue
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
#        if top.lrelativ:
#          vx = take(top.uxplost[i1:i2],ic-i1)*take(top.gaminvlost[i1:i2],ic-i1)
#          vy = take(top.uyplost[i1:i2],ic-i1)*take(top.gaminvlost[i1:i2],ic-i1)
#          vz = take(top.uzplost[i1:i2],ic-i1)*take(top.gaminvlost[i1:i2],ic-i1)
#        else:
#          vx = take(top.uxplost[i1:i2],ic-i1)
#          vy = take(top.uyplost[i1:i2],ic-i1)
#          vz = take(top.uzplost[i1:i2],ic-i1)
        xo = take(top.pidlost[i1:i2,self.xoldpid],ic-i1)
        yo = take(top.pidlost[i1:i2,self.yoldpid],ic-i1)
        zo = take(top.pidlost[i1:i2,self.zoldpid],ic-i1)
        vx = (xc-xo)/top.dt
        vy = (yc-yo)/top.dt
        vz = (zc-zo)/top.dt
        intercept = c.intercept(xc,yc,zc,vx,vy,vz)
        put(top.xplost,ic,intercept.xi)
        put(top.yplost,ic,intercept.yi)
        put(top.zplost,ic,intercept.zi)
        put(top.pidlost[:,-3],ic,intercept.itheta)
        put(top.pidlost[:,-2],ic,intercept.iphi)
        dt = (sqrt((xc - intercept.xi)**2 +
                   (yc - intercept.yi)**2 +
                   (zc - intercept.zi)**2)/
              dvnz(sqrt(vx**2 + vy**2 + vz**2)))
        put(top.pidlost[:,-4],ic,top.time - dt)
      if self.lcollectlpdata:
        pidlostcondid = take(top.pidlost[:,-1],iscrape1)
        pidtoconsider = compress(pidlostcondid==c.condid,iscrape1)
        if top.wpid==0:
          w = len(pidtoconsider)
        else:
          w = sum(take(top.pidlost[:,wpid],pidtoconsider))
        w=parallelsum(w)
        if me==0:
          c.lostparticles_data += [[top.time, 
                                    w*top.sq[js]*top.sw[js],
                                    top.dt,
                                    js]]

