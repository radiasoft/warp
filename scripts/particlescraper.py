"""
ParticleScraper: class for creating particle scraping
"""
from warp import *
from generateconductors import *
import timing as t

particlescraper_version = "$Id: particlescraper.py,v 1.49 2007/05/21 23:35:22 dave Exp $"
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
 - lrefineintercept: when true, with lsaveintercept, lost particles are advanced
                     from the old positions using a time step small compared to
                     the cyclotron gyroperiod to calculated a better value for
                     the intercept.
 - lrefineallintercept: same as lrefineintercept, but the trajectory of all
                        particles near the conductors are refined rather then 
                        only particles which have already been lost
 - nstepsperorbit=8: number of refined time steps to take when using
                     lrefineintercept.
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
  def __init__(self,conductors,lsavecondid=0,lsaveintercept=0,
                    lrefineintercept=0,lrefineallintercept=0,nstepsperorbit=8,
                    lcollectlpdata=0,mglevel=0,aura=0.,install=1,grid=None): 
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
    self.lsaveintercept = lsaveintercept or lrefineintercept or lrefineallintercept
    self.lrefineintercept = lrefineintercept
    self.lrefineallintercept = lrefineallintercept
    self.nstepsperorbit = nstepsperorbit
    self.lsavecondid = (lsavecondid or lsaveintercept or
                        lrefineintercept or lrefineallintercept or
                        lcollectlpdata)
    self.lcollectlpdata = lcollectlpdata
    if self.lsavecondid:
      top.lsavelostpart = true
    if self.lsaveintercept:
      if not top.lsaveoldpos:top.lsaveoldpos=true
      # --- Note that nextpid returns numbers based on 1 based indexing
      if top.xoldpid==0: top.xoldpid=nextpid()
      if top.yoldpid==0: top.yoldpid=nextpid()
      if top.zoldpid==0: top.zoldpid=nextpid()
      self.xoldpid = top.xoldpid - 1
      self.yoldpid = top.yoldpid - 1
      self.zoldpid = top.zoldpid - 1
      if self.lrefineintercept or self.lrefineallintercept:
        if top.uxoldpid==0: top.uxoldpid=nextpid()
        if top.uyoldpid==0: top.uyoldpid=nextpid()
        if top.uzoldpid==0: top.uzoldpid=nextpid()
        self.uxoldpid = top.uxoldpid - 1
        self.uyoldpid = top.uyoldpid - 1
        self.uzoldpid = top.uzoldpid - 1
      setuppgroup(top.pgroup)
    self.l_print_timing=0
    # --- If the user specified the grid, then add the conductors
    if self.usergrid: self.updateconductors()
    # --- Install the call to scrape particles if requested
    if install: self.installscraper()

  def installscraper(self):
    # --- Install the call to scrape particles
    if not isinstalledparticlescraper(self.scrapeall):
      installparticlescraper(self.scrapeall)

  def disable(self):
    if isinstalledparticlescraper(self.scrapeall):
      uninstallparticlescraper(self.scrapeall)

  def __setstate__(self,dict):
    # --- This is called when the instance is unpickled.
    self.__dict__.update(dict)
    self.installscraper()

  def registerconductors(self,newconductors):
#    self.updategrid()
    if type(newconductors) is not ListType: newconductors = [newconductors]
    for c in newconductors:
      assert c.condid != 0,"The conductor id must be nonzero in order for the particle scraping to work."
      self.conductors.append(c)
#      self.grid.getisinside(c,mglevel=self.mglevel,aura=self.aura)
    self.updategrid()

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
    if top.chdtspid>0:
      if w3d.nxc<>self.grid.nx or w3d.nyc<>self.grid.ny or w3d.nzc<>self.grid.nz:
        w3d.nxc=self.grid.nx
        w3d.nyc=self.grid.ny
        w3d.nzc=self.grid.nz
        gchange('Fields3dParticles')
        sum_neighbors3d(nint(self.grid.isinside),w3d.isnearbycond,self.grid.nx,self.grid.ny,self.grid.nz)

  def updateconductors(self):
    for c in self.conductors:
      self.grid.getisinside(c,mglevel=self.mglevel,aura=self.aura)

  def scrapeall(self,clear=0):
    if len(self.conductors)==0 or parallelsum(sum(top.pgroup.nps))==0: return
    self.updategrid()
    for js in xrange(top.pgroup.ns):
      if top.pgroup.ldts[js]:
        if self.l_print_timing:tt=0.
        if self.l_print_timing:t.start()
        self.scrape(js);
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
    # --- If there are no particles in this species, that nothing needs to be done
    if top.pgroup.nps[js] == 0: return

    # --- Get mesh information into local variables
    dx,dy,dz,nx,ny,nz,iz = self.grid.getmeshsize(self.mglevel)
    xmin = self.grid.xmin
    xmax = self.grid.xmax
    ymin = self.grid.ymin
    ymax = self.grid.ymax
    zmin = self.grid.zmmin + iz*dz + top.zbeam
    zmax = self.grid.zmmin + (iz+nz)*dz + top.zbeam
    isinside = self.grid.isinside

    # --- Get handy references to the particles in the species
    i1 = top.pgroup.ins[js] - 1
    i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
    xx = top.pgroup.xp[i1:i2]
    yy = top.pgroup.yp[i1:i2]
    zz = top.pgroup.zp[i1:i2]
    pp = zeros(top.pgroup.nps[js],'d')
    #if js==1:print js,i1,i2,top.pgroup.zp[i1:i2],top.zbeam

    # --- Find which particles are close to a conductor. This
    # --- interpolates from the isinside grid. The results are
    # --- put into the array pp.
    if w3d.solvergeom in [w3d.XYZgeom]:
      getgrid3d(top.pgroup.nps[js],xx,yy,zz,pp,
                nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,
                w3d.l2symtry,w3d.l4symtry)
    elif w3d.solvergeom == w3d.RZgeom:
      # --- Note that for RZ, the radius is calculated for this, but
      # --- the original particle position is used below.
      rr = sqrt(xx**2 + yy**2)
      getgrid2d(top.pgroup.nps[js],rr,zz,pp,nx,nz,isinside[:,0,:],
                xmin,xmax,zmin,zmax)
    elif w3d.solvergeom == w3d.XZgeom:
      getgrid2d(top.pgroup.nps[js],xx,zz,pp,nx,nz,isinside[:,0,:],
                xmin,xmax,zmin,zmax)
    elif w3d.solvergeom == w3d.XYgeom:
      getgrid2d(top.pgroup.nps[js],xx,yy,pp,nx,ny,isinside[:,:,0],
                xmin,xmax,ymin,ymax)
    else:
      raise "The particle scraping only works for XYZ, XY and RZ geometry"

    # --- Get indices for all of the particles which are close to a
    # --- conductor. If there are none, then immediately return.
    # --- Note, of course, that close may mean inside.
    iclose = compress(pp>0.,arange(i1,i2))
    if len(iclose) == 0: return
 
    # --- Get the positions of particles which are close to a conductor.
    xx = take(xx,iclose-i1)
    yy = take(yy,iclose-i1)
    zz = take(zz,iclose-i1)

    # --- The 'g' lists give the locations of the corners of the grid cell
    # --- relative to the grid location of the particles close to a
    # --- conductor. Also, get those grid locations.
    if w3d.solvergeom in [w3d.XYZgeom]:
      nd = 3
      gdx = [0.,dx,0.,dx,0.,dx,0.,dx]
      gdy = [0.,0.,dy,dy,0.,0.,dy,dy]
      gdz = [0.,0.,0.,0.,dz,dz,dz,dz]
      xg = xmin+int(abs(xx-xmin)/dx)*dx 
      yg = ymin+int(abs(yy-ymin)/dy)*dy 
      zg = zmin+int(abs(zz-zmin)/dz)*dz 
    elif w3d.solvergeom in [w3d.RZgeom]:
      nd = 2
      gdx = [0.,dx,0.,dx]
      gdz = [0.,0.,dz,dz]
      # --- Like above, the radius is calculated in the temporary, but the
      # --- original particle position is used below.
      # --- These two lines calculating rr give the same result, but the second
      # --- is probably faster
      #rr = sqrt(xx**2 + yy**2)
      rr = take(rr(iclose-i1))
      xg = xmin+int(abs(rr-xmin)/dx)*dx 
      zg = zmin+int(abs(zz-zmin)/dz)*dz 
    elif w3d.solvergeom in [w3d.XZgeom]:
      nd = 2
      gdx = [0.,dx,0.,dx]
      gdz = [0.,0.,dz,dz]
      xg = xmin+int(abs(xx-xmin)/dx)*dx 
      zg = zmin+int(abs(zz-zmin)/dz)*dz 
    elif w3d.solvergeom == w3d.XYgeom:
      nd = 2
      gdx = [0.,dx,0.,dx]
      gdy = [0.,0.,dy,dy]
      xg = xmin+int(abs(xx-xmin)/dx)*dx 
      yg = ymin+int(abs(yy-ymin)/dy)*dy 
    
    nn = len(iclose)
    pp = zeros(nn,'d')

    # --- Loop over the corners of the grid cell
    for i in range(2**nd):

      # --- Get id of the conductor that the particles are near
      if w3d.solvergeom in [w3d.XYZgeom]:
        getgridngp3d(nn,xg+gdx[i],yg+gdy[i],zg+gdz[i],pp,
                     nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,0.,
                     w3d.l2symtry,w3d.l4symtry)
      elif w3d.solvergeom in [w3d.XZgeom,w3d.RZgeom]:
        getgridngp2d(nn,xg+gdx[i],zg+gdz[i],pp,nx,nz,isinside[:,0,:],
                     xmin,xmax,zmin,zmax)
      elif w3d.solvergeom == w3d.XYgeom:
        getgridngp2d(nn,xg+gdx[i],yg+gdy[i],pp,nx,ny,isinside[:,:,0],
                     xmin,xmax,ymin,ymax)

      # --- Loop over the conductors, removing particles that are found inside
      # --- of each.
      for c in self.conductors:

        # --- Get indices relative to the temporary arrays.
        # --- Note that iclose is relative to the full particle arrays.
        itempclose=arange(nn)

        # --- Get indices of particles that are close to the conductor
        ii = compress(pp == c.condid,itempclose) 

        # --- If there are no particles close, then skip to the next conductor
        if len(ii) == 0: continue

        # --- Get positions of the particles that are close
        xc = take(xx,ii)
        yc = take(yy,ii)
        zc = take(zz,ii)

        # --- Find the particles that are currently inside and down-select
        # --- the indices. The nint is needed since the quantities is used in
        # --- logical expressions below which require ints.
        currentisinside = nint(c.isinside(xc,yc,zc).isinside)
        iic = compress(currentisinside,ii)
        ic = take(iclose,iic)

        if self.lrefineallintercept:
          # --- Refine whether or not particles are lost by taking small time
          # --- steps, starting from the old position. Note that it is possible
          # --- that particles that were lost could be not lost upon refinement,
          # --- and similarly, particles that were not lost, could be lost upon
          # --- refinement.
          # --- Get the old coordinates of particles that are close.
          iclose1 = take(iclose,ii)
          xo = take(top.pgroup.pid[:,self.xoldpid],iclose1)
          yo = take(top.pgroup.pid[:,self.yoldpid],iclose1)
          zo = take(top.pgroup.pid[:,self.zoldpid],iclose1)
          uxo = take(top.pgroup.pid[:,self.uxoldpid],iclose1)
          uyo = take(top.pgroup.pid[:,self.uyoldpid],iclose1)
          uzo = take(top.pgroup.pid[:,self.uzoldpid],iclose1)

          # --- Get the current fields
          ex = take(top.pgroup.ex,iclose1)
          ey = take(top.pgroup.ey,iclose1)
          ez = take(top.pgroup.ez,iclose1)
          bx = take(top.pgroup.bx,iclose1)
          by = take(top.pgroup.by,iclose1)
          bz = take(top.pgroup.bz,iclose1)

          # --- Create some temporaries
          itime = None
          dt = top.dt*top.pgroup.ndts[js]*top.pgroup.dtscale[js]*ones(len(ii))
          q = top.pgroup.sq[js]
          m = top.pgroup.sm[js]

          # --- Do the refinement calculation. The currentisinside argument controls
          # --- when the current position is replaced by the refined position.
          # --- If the particle is currently lost but in the refined
          # --- calculation is not lost, then the replace the current position
          # --- with that refined position that is not lost. Similarly, if the
          # --- particle is currently not lost, but in the refined calculation
          # --- is lost, then replace the current position with the refined
          # --- position.
          self.refineintercept(c,xc,yc,zc,xo,yo,zo,uxo,uyo,uzo,ex,ey,ez,bx,by,bz,itime,dt,q,m,currentisinside)

          # --- Do the replacements as described above. Note that for lost
          # --- particles, xc,yc,zc hold the positions of the particles one
          # --- small time step into the conductor.
          put(top.pgroup.xp,iclose1,xc)
          put(top.pgroup.yp,iclose1,yc)
          put(top.pgroup.zp,iclose1,zc)
          put(top.pgroup.uxp,iclose1,uxo)
          put(top.pgroup.uyp,iclose1,uyo)
          put(top.pgroup.uzp,iclose1,uzo)

          # --- Determine whether the refined positions are lost.
          refinedisinside = nint(c.isinside(xc,yc,zc).isinside)

          # --- iic lists the particles that are either currently lost or
          # --- lost in the refined calculation. Both types of particles have
          # --- had their positions modified so are included in the list of
          # --- particles that don't need to be further handled.
          iic = compress(currentisinside | refinedisinside,ii)

          # --- ic lists the particles that are lost in the refined
          # --- calculation. These will be scraped.
          ic = take(iclose,compress(refinedisinside,ii))

          # --- Note that the old values of the positions are changed
          # --- only for particles for which the refined calculation
          # --- shows they are lost. This is needed for the interception
          # --- calculation done in savecondid.
          iclose2 = compress(refinedisinside,iclose1)
          if len(iclose2) > 0:
            put(top.pgroup.pid[:,self.xoldpid],iclose2,compress(refinedisinside,xo))
            put(top.pgroup.pid[:,self.yoldpid],iclose2,compress(refinedisinside,yo))
            put(top.pgroup.pid[:,self.zoldpid],iclose2,compress(refinedisinside,zo))

        # --- If no particle are inside the conductor, then skip to the next one
        if len(iic) == 0: continue

        # --- For particles which are inside, set gaminv to 0, the lost
        # --- particle flag
        put(top.pgroup.gaminv,ic,0.)

        # --- Remove the already handled particles, returning if there are no more.
        put(iclose,iic,-1)
        iclose = compress(iclose>=0,iclose)        
        nn = len(iclose)
        if nn == 0: return
        put(itempclose,iic,-1)
        itempclose = compress(itempclose>=0,itempclose)        
        xx = take(xx,itempclose)
        yy = take(yy,itempclose)
        zz = take(zz,itempclose)
        xg = take(xg,itempclose)
        pp = take(pp,itempclose)
        if w3d.solvergeom in [w3d.XYZgeom,w3d.XYgeom]:
          yg = take(yg,itempclose)
        if w3d.solvergeom in [w3d.XYZgeom,w3d.XZgeom,w3d.RZgeom]:
          zg = take(zg,itempclose)

  def savecondid(self,js):
    jsid = top.pgroup.sid[js]

    # --- Just return if there are no lost particles.
    if top.npslost[jsid] == 0: 
      if self.lcollectlpdata:
        # --- If data is being collected, the 0 from this processor must still
        # --- be added to the sum.
        for c in self.conductors:
          # --- This parallelsum coordinates with the ones below.
          w=parallelsum(0.)
          if w<>0.:
            c.lostparticles_data += [[top.time, 
                                      w*top.pgroup.sq[js]*top.pgroup.sw[js],
                                      top.dt,
                                      jsid]]
      return

    # --- First make sure there is extra space in the pidlost array.
    pidspace = 1
    if self.lsaveintercept: pidspace = 4
    if top.npidlost < top.npid+pidspace:
      top.npidlost = top.npid + pidspace
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

    i1 = top.inslost[jsid] - 1
    i2 = top.inslost[jsid] + top.npslost[jsid] - 1
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
    if w3d.solvergeom in [w3d.XYZgeom]:
      getgridngp3d(nn,xg,yg,zg,pp,
                   nx,ny,nz,isinside,xmin,xmax,ymin,ymax,zmin,zmax,0.,
                   w3d.l2symtry,w3d.l4symtry)
    elif w3d.solvergeom == w3d.RZgeom or w3d.solvergeom == w3d.XZgeom:
      getgridngp2d(nn,xg,zg,pp,nx,nz,isinside[:,0,:],xmin,xmax,zmin,zmax)
    elif w3d.solvergeom == w3d.XYgeom:
      getgridngp2d(nn,xg,yg,pp,nx,ny,isinside[:,:,0],xmin,xmax,ymin,ymax)
    else:
      raise "The particle scraping only works for XYZ, XZ, XY and RZ geometry"

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
          # --- This parallelsum coordinates with the other processors
          w=parallelsum(0.)
          if w<>0.:
            c.lostparticles_data += [[top.time, 
                                      w*top.pgroup.sq[js]*top.pgroup.sw[js],
                                      top.dt,
                                      jsid]]
        continue
      xc = take(x8,ii)
      yc = take(y8,ii)
      zc = take(z8,ii)
      ic = take(iscrape,ii)
      ic = compress(c.isinside(xc,yc,zc).isinside,ic)
      if len(ic) == 0: continue
      # --- For particles which are inside, set pid to the id of the conductor
      # --- where the particle is lost.
      put(top.pidlost[:,-1],ic,c.condid)
      # --- Save location and surface normal where particle intercepted the
      # --- conductor.
      if self.lsaveintercept:
        xc = take(xx,ic-i1)
        yc = take(yy,ic-i1)
        zc = take(zz,ic-i1)
        xo = take(top.pidlost[:,self.xoldpid],ic)
        yo = take(top.pidlost[:,self.yoldpid],ic)
        zo = take(top.pidlost[:,self.zoldpid],ic)

        if self.lrefineintercept:
          uxo = take(top.pidlost[:,self.uxoldpid],ic)
          uyo = take(top.pidlost[:,self.uyoldpid],ic)
          uzo = take(top.pidlost[:,self.uzoldpid],ic)
          ex = take(top.exlost,ic)
          ey = take(top.eylost,ic)
          ez = take(top.ezlost,ic)
          bx = take(top.bxlost,ic)
          by = take(top.bylost,ic)
          bz = take(top.bzlost,ic)
          itime = zeros(len(ic),'d')
          dt = top.dt*top.pgroup.ndts[js]*top.pgroup.dtscale[js]*ones(len(ic))
          q = top.pgroup.sq[js]
          m = top.pgroup.sm[js]
          self.refineintercept(c,xc,yc,zc,xo,yo,zo,uxo,uyo,uzo,ex,ey,ez,bx,by,bz,itime,dt,q,m,0)
        else:
          itime = 0.
          dt = top.dt

        if self.lrefineallintercept:
          # --- For this case, the velocities have already been correctly
          # --- calculated and so can be used directly.
          if top.lrelativ:
            vx = take(top.uxplost,ic)*take(top.gaminvlost,ic)
            vy = take(top.uyplost,ic)*take(top.gaminvlost,ic)
            vz = take(top.uzplost,ic)*take(top.gaminvlost,ic)
          else:
            vx = take(top.uxplost,ic)
            vy = take(top.uyplost,ic)
            vz = take(top.uzplost,ic)
        else:
          # --- Otherwise, use an approximate calculation.
          vx = (xc-xo)/dt
          vy = (yc-yo)/dt
          vz = (zc-zo)/dt

        intercept = c.intercept(xc,yc,zc,vx,vy,vz)
        dtintercept = (sqrt((xc - intercept.xi)**2 +
                            (yc - intercept.yi)**2 +
                            (zc - intercept.zi)**2)/
              dvnz(sqrt(vx**2 + vy**2 + vz**2))) + itime

        put(top.xplost,ic,intercept.xi)
        put(top.yplost,ic,intercept.yi)
        put(top.zplost,ic,intercept.zi)

        if self.lrefineintercept:
          # --- Also, reset the velocities
          if top.lrelativ:
            beta = sqrt(vx**2 + vy**2 + vz**2)/clight
            gamma = 1./sqrt((1.-beta)*(1.+beta))
          else:
            gamma = 1.
          put(top.uxplost,ic,vx*gamma)
          put(top.uyplost,ic,vy*gamma)
          put(top.uzplost,ic,vz*gamma)

        # --- Set the angle of incidence and time of interception
        put(top.pidlost[:,-3],ic,intercept.itheta)
        put(top.pidlost[:,-2],ic,intercept.iphi)
        put(top.pidlost[:,-4],ic,top.time - dtintercept)

      if self.lcollectlpdata:
        pidlostcondid = take(top.pidlost[:,-1],iscrape1)
        pidtoconsider = compress(pidlostcondid==c.condid,iscrape1)
        if top.wpid==0:
          w = len(pidtoconsider)
        else:
          w = sum(take(top.pidlost[:,top.wpid],pidtoconsider))
        # --- This parallelsum coordinates with the ones above
        w=parallelsum(w)
        c.lostparticles_data += [[top.time, 
                                  w*top.pgroup.sq[js]*top.pgroup.sw[js],
                                  top.dt,
                                  jsid]]



  def refineintercept(self,c,xc,yc,zc,xo,yo,zo,uxo,uyo,uzo,ex,ey,ez,bx,by,bz,itime,dt,q,m,luserefinedifnotlost):
    """Refine the location of the intercept
c: the conductor
xc,yc,zc: input holding the current particle position
          output holding the point just inside the conductor
xo,yo,zo: input holding the old particle position
          output holding the point just outside the conductor
uxo,uyo,uzo: input holding old velocity, time synchronized with xo,yo,zo
             output holding velocity near time when particle entered the
             conductor
ex,ey,ez,bx,by,bz: fixed E and B fields
itime: output holding time that the particle reached xo,yo,zo
       If input value is None, then the time is not saved
dt: output holding particle time step sizes
q,m: charge and mass of the particles
luserefinedifnotlost: when true, if the refined particle orbit is not lost,
                      then replace replace the current position with the
                      refined position
    """
    # --- Note that this routine is written in a non-pythonic way, where the
    # --- changes are made directly in the input arrays.

    # --- Get the number of particles
    nn = len(xo)

    # --- The cyclotron frequency for each particle
    magB = sqrt(bx**2 + by**2 + bz**2)
    omegac = q/m*magB

    # --- Get the number of steps for each particle. 
    # --- This is set by self.nstepsperorbit which is the number of steps
    # --- per cyclotron orbit.
    isteps = nint(dt*omegac/(2.*pi)*self.nstepsperorbit)

    # --- So that the orbit is always refined, at least a little, set the
    # --- minimum number of steps to be self.nstepsperorbit. This is helpful
    # --- for example if the B field is zero.
    isteps = maximum(isteps,self.nstepsperorbit)

    # --- Get the maximum number of steps needed
    nsteps  = max(isteps)

    # --- Calculate the overcycled step size for each particle
    # --- Note that dtover will be modified in the loop below,
    # --- and dt is used to return the time step size.
    dtover = dt/isteps
    dt[:] = dtover

    # --- Recalculate gaminv. The code could also save the old gaminv, but
    # --- this should be equivalent.
    if top.lrelativ:
      usq = (uxo**2 + uyo**2 + uzo**2)
      gamma = sqrt(1. + usq/clight**2)
      gaminv = 1./gamma
    else:
      gaminv = ones(nn,'d')

    # --- Save the positions. This is needed so that the data can be
    # --- restored if no intercept is found below.
    xcsave = xc.copy()
    ycsave = yc.copy()
    zcsave = zc.copy()
    xosave = xo.copy()
    yosave = yo.copy()
    zosave = zo.copy()
    uxosave = uxo.copy()
    uyosave = uyo.copy()
    uzosave = uzo.copy()

    # --- Get the starting positions of the advance. These should all be
    # --- outside of the conductor. The loop below advances the particles
    # --- using these 'c' arrays.
    xc[:] = xo
    yc[:] = yo
    zc[:] = zo
    if itime is not None:
      itime[:] = 0.

    # --- Do the over cycling loop. Note that all particles are advanced by
    # --- the maximum number of steps, though once a particle goes inside
    # --- the conductor, its time step is set to zero so the coordinates
    # --- don't change anymore.
    # --- One possible optimization is to have the fortran advancing
    # --- routines skip particles that have a zero time step size.
    for it in range(nsteps):

      # --- Do a full split leap-frog advance (with constant E and B fields)
      # --- Note that this does the advance in place, directly changing the
      # --- input arrays.
      bpusht3d(nn,uxo,uyo,uzo,gaminv,bx,by,bz,q,m,dtover,0.5,top.ibpush)
      epusht3d(nn,uxo,uyo,uzo,ex,ey,ez,q,m,dtover,0.5)
      gammaadv(nn,gaminv,uxo,uyo,uzo,top.gamadv,top.lrelativ)
      xpusht3d(nn,xc,yc,zc,uxo,uyo,uzo,gaminv,dtover)
      epusht3d(nn,uxo,uyo,uzo,ex,ey,ez,q,m,dtover,0.5)
      gammaadv(nn,gaminv,uxo,uyo,uzo,top.gamadv,top.lrelativ)
      bpusht3d(nn,uxo,uyo,uzo,gaminv,bx,by,bz,q,m,dtover,0.5,top.ibpush)

      # --- This provides a nice diagnostic for testing
      #plp(yc,xc,marker=circle,color=green)

      # --- Check whether the new positions are inside of the conductor.
      isinside = c.isinside(xc,yc,zc).isinside

      # --- For the particles that are still outside, set the old positions
      # --- to be the updated positions
      xo[:] = where(isinside,xo,xc)
      yo[:] = where(isinside,yo,yc)
      zo[:] = where(isinside,zo,zc)

      # --- If a particle is inside the conductor, then stop advancing it,
      # --- setting its time step size to zero.
      dtover = where(isinside,0.,dtover)

      # --- Now, advance itime. Note that for particles that are now inside,
      # --- itime is not advanced since it is the time just before the particle
      # --- enters the conductor.
      if itime is not None:
        itime[:] = itime + dtover

      # --- Quit the loop if all intercepts have been found.
      if alltrue(dtover==0.): break

    # --- Check for cases where no interception was found. In those cases,
    # --- restore the original data since that will at least not cause
    # --- a code problem. This checks if dtover hasn't been zeroed out
    # --- which means that the particle was never flagged as being inside
    # --- in the loop above.
    if sometrue(dtover > 0.):
      userefined = ((dtover == 0.) | luserefinedifnotlost)
      xc[:] = where(userefined,xc,xcsave)
      yc[:] = where(userefined,yc,ycsave)
      zc[:] = where(userefined,zc,zcsave)
      xo[:] = where(userefined,xo,xosave)
      yo[:] = where(userefined,yo,yosave)
      zo[:] = where(userefined,zo,zosave)
      uxo[:] = where(userefined,uxo,uxosave)
      uyo[:] = where(userefined,uyo,uyosave)
      uzo[:] = where(userefined,uzo,uzosave)
      if itime is not None:
        itime[:] = where(userefined,itime,0.)

