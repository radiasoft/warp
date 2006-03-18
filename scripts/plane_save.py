"""
Contains PlaneSave class that is used to save particle and field data at a
specified z plane. The data is used by PlaneRestore to continue the
simulation. The two simulations are linked together.
"""
from warp import *
plane_save_version = "$Id: plane_save.py,v 1.12 2006/03/18 00:32:03 dave Exp $"

class PlaneSave:
  """
Saves the particle data and phi to a file just after the 
field solve is finished. The positions and velocities are staggered.
It is automatically called every time step after the field solve
Input:
  - zplane: location where simulations are linked together. Units of meters
            relative to the lab frame. Note grid cell nearest zplane is
            the actual location where data is saved.
  - filename=runid.plane: filename where data is stored
  - js: species which are saved. Defaults to all species. Can be single
        integer or a list of integers.
  - allways_save=false: if set to true, particles and potential are saved
                        at all time step. Default is false: saving starts
                        when first particle cross zplane.
  - deltaz: z grid cell size of simulation where the data will be restored.
            Defaults to w3d.dz, must be an integer multiple of w3d.dz. 
  - deltat: time step size of simulation where the data will be restored.
            Defaults to top.dt. top.dt must be an integer multiple of deltat.
  - maxvzdt=deltaz*deltat/dt: maximum distance particles are expected to travel
                              when passing through zplane
  - newfile=0: When true, creates a new file to save data into, otherwise
               append to file if it already exists.

  """

  def __init__(self,zplane,filename=None,js=None,allways_save=false,
                    deltaz=None,deltat=None,maxvzdt=None,newfile=0):

    self.zplane = nint(zplane/w3d.dz)*w3d.dz

    # --- save only if between grid bounds, otherwise raise and exception
    if(self.zplane<w3d.zmminglobal or self.zplane>=w3d.zmmaxglobal):
      raise "The zplane specified is outside the simulation domain"

    if allways_save:
      self.save_this_step = true
    else:
      self.save_this_step = false

    # --- Set distance between saved phi planes
    if deltaz is None: self.deltaz = w3d.dz
    else:              self.deltaz = deltaz
    self.izz = max(1,nint(self.deltaz/w3d.dz))

    # --- Set frequency that data is saved.
    if deltat is None: self.deltat = top.dt
    else:              self.deltat = deltat
    self.itt = max(1,nint(self.deltat/top.dt))

    # --- Maximum distance particles are expected to travel when passing
    # --- through zplane. This is the initial value taken and will be updated
    # --- if faster particles are found.
    if maxvzdt is None: self.maxvzdt = w3d.dz*self.izz*self.itt
    else:               self.maxvzdt = maxvzdt

    # --- initializes list of species
    if type(js) == IntType:
      self.jslist= [js]
    elif js is None:
      self.jslist = range(top.ns)
    else:
      self.jslist = js

    # --- defines useful variables
    self.it = 0

    # --- Set so data is saved to file immediately after a field solve.
    installafterfs(self.saveplane)

    # --- Set the name of the file which will hold the data
    if filename is None:
      self.filename = arraytostr(top.runid)+".plane"
    else:
      self.filename = filename

    # --- The file is only opened and the data written on processor 0
    fileexists = os.access(self.filename,os.F_OK)
    if me == 0 and (newfile or not fileexists):
      self.f = PW.PW(self.filename)

      # --- save plane size and location and time step
      self.f.zplane    = self.zplane
      self.f.nx_plane  = w3d.nx
      self.f.ny_plane  = w3d.ny
      self.f.ixa_plane = w3d.ix_axis
      self.f.iya_plane = w3d.iy_axis
      self.f.xmmin     = w3d.xmmin
      self.f.xmmax     = w3d.xmmax
      self.f.ymmin     = w3d.ymmin
      self.f.ymmax     = w3d.ymmax
      self.f.dt        = top.dt
      self.f.deltaz    = self.deltaz
      self.f.deltat    = self.deltat
      self.f.npid      = top.npid
     
      # --- set sym_plane and write it out
      if (w3d.l4symtry):
        sym_plane = 4
      elif (w3d.l2symtry):
        sym_plane = 2
      else:
        sym_plane = 1
      self.f.sym_plane = sym_plane

      # --- Write out the solver geometry flag
      self.f.solvergeom = w3d.solvergeom

      # --- Write out particle quantities for each species
      for js in self.jslist:
        self.f.write('sq_%d'%js,top.sq[js])
        self.f.write('sm_%d'%js,top.sm[js])
        self.f.write('sw_%d'%js,top.sw[js])

      self.f.set_verbosity(0)
      self.f.close()

  def saveplane(self):
    # --- Only save data at specified frequency
    if (top.it % self.itt) > 0: return

    if me == 0:
      # --- open file in append mode
      self.f = PW.PW(self.filename,'a')

    # --- save for each species
    for js in self.jslist:
        self.saveplanespecies(js)

    # --- phi is saved every time step whether or not there are particles saved
    # --- but only after saving has started.
    if(self.save_this_step):
      # --- get the two planes of phi to be saved
      iz = nint((self.zplane - top.zbeam - w3d.zmminglobal)/w3d.dz)
      #self.phi_save[:,...,0] = getphi(iz=iz-self.izz)
      #self.phi_save[:,...,1] = getphi(iz=iz)
      self.phi_save = transpose(array([transpose(getphi(iz=iz-self.izz)),
                                       transpose(getphi(iz=iz))]))
      if me == 0:
        try:
          self.f.write('phiplane%d'%self.it,self.phi_save)
        except:
          print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print "ERROR: There was an error writing out the particle data"
          print "       This most likely means that an attempt was made"
          print "       to overwrite an existing file."
          print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    if me == 0:
      if(self.save_this_step):
        # --- Write start time and update end time
        if(self.it==1): self.f.tmin = top.time
        self.f.tmax = top.time

      # --- close file
      self.f.set_verbosity(0)
      self.f.close()

  def saveplanespecies(self,js):
    # --- Gather vz of all particles somewhat beyond zplane
    vz = getvz(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
    np = getn(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)

    if np > 0:
      # --- Make sure that the region is large enough to capture all particles
      # --- that may have crossed zplane.
      # --- If there are faster particles, increase maxvzdt appropriately and
      # --- reget vz. Do this until maxvzdt > max(vz)*dt*itt
      maxvz = globalmax(vz)
      while maxvz*top.dt*self.itt > self.maxvzdt:
        self.maxvzdt = maxvz*top.dt*self.itt*1.1
        vz = getvz(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
        maxvz = globalmax(vz)

      # --- Get rest of particle data
      xx = getx(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
      yy = gety(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
      zz = getz(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
      ux = getux(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
      uy = getuy(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
      uz = getuz(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
      gi = getgaminv(js=js,zl=self.zplane,zu=self.zplane+self.maxvzdt)
      if top.npidmax > 0:
        id = []
        for i in range(top.npidmax):
          id.append(getpid(js=js,id=i,zl=self.zplane,zu=self.zplane+self.maxvzdt))
        id = transpose(array(id))

      # --- Recalculate the previous axial location of the particles to check
      # --- if the particles crossed the saving plane.
      old_zz = zz - vz*top.dt*self.itt

      # --- find indices of particles which crossed the plane in the
      # --- last time step
      ii = compress((old_zz<=self.zplane) & (zz>=self.zplane), arange(len(zz)))

      np_save = shape(ii)[0]
      if np_save > 0: self.save_this_step = true
      self.save_this_step = broadcast(self.save_this_step)

      if (self.save_this_step):
        self.it = self.it + 1

      if (np_save > 0 and self.save_this_step and me == 0):
        suffix = '%d_%d'%(self.it,js)
        try:
          self.f.write('xp'+suffix,    take(xx,ii))
          self.f.write('yp'+suffix,    take(yy,ii))
          self.f.write('zp'+suffix,    take(zz,ii))
          self.f.write('uxp'+suffix,   take(ux,ii))
          self.f.write('uyp'+suffix,   take(uy,ii))
          self.f.write('uzp'+suffix,   take(uz,ii))
          self.f.write('gaminv'+suffix,take(gi,ii))
          if top.npidmax > 0:
            self.f.write('pid'+suffix,   take(id,ii))
        except:
          print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print "ERROR: There was an error writing out the particle data"
          print "       This most likely means that an attempt was made"
          print "       to overwrite an existing file."
          print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

