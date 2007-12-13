"""
Contains PlaneSave class that is used to save particle and field data at a
specified z plane. The data is used by PlaneRestore to continue the
simulation. The two simulations are linked together.
"""

__all__ = ['PlaneSave','plane_save_version']

from warp import *
plane_save_version = "$Id: plane_save.py,v 1.18 2007/12/13 02:10:26 dave Exp $"

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
  - newfile=0: When true, creates a new file to save data into, otherwise
               append to file if it already exists.
  - lsavephi=1: When true, saves the potential around the zplane.
  - lsaveparticles=1: When true, save the particles the pass the zplane.

  """

  def __init__(self,zplane,filename=None,js=None,allways_save=false,
                    deltaz=None,deltat=None,newfile=0,
                    lsavephi=1,lsaveparticles=1):

    self.zplane = zplane
    self.js = js
    self.deltat = deltat
    self.deltaz = deltaz
    self.newfile = newfile
    self.lsavephi = lsavephi
    self.lsaveparticles = lsaveparticles

    if allways_save:
      self.save_this_step = true
    else:
      self.save_this_step = false

    # --- defines useful variables
    self.it = 0

    # --- Set so data is saved to file immediately after a field solve.
    installafterfs(self.saveplane)

    # --- Set the name of the file which will hold the data
    self.f = None
    if filename is None:
      self.filename = arraytostr(top.runid)+".plane"
    else:
      self.filename = filename

  def initsaveplane(self):

    # --- Do this here in case that init was called before top.dt was
    # --- initialized.
    if self.deltat is None:
      self.deltat = top.dt

  def initsaveparticles(self):
    if not self.lsaveparticles: return

    # --- initializes list of species
    if self.js is None:
      self.jslist = range(top.ns)
    else:
      try:
        list(self.js)
        self.jslist = self.js
      except TypeError:
        self.jslist= [self.js]

    # --- Create space for the old z position
    self.zoldpid = nextpid() - 1
    setuppgroup(top.pgroup)

    # --- Save the initial old z position
    for js in self.jslist:
      if top.pgroup.nps[js] > 0:
        j1 = top.pgroup.ins[js] - 1
        j2 = j1 + top.pgroup.nps[js]
        top.pgroup.pid[j1:j2,self.zoldpid] = top.pgroup.zp[j1:j2]

  def initsavephi(self):

    # --- Set distance between saved phi planes
    # --- Do this here in case that init was called before w3d.dz was
    # --- initialized.
    if self.deltaz is None:
      self.deltaz = w3d.dz

  def writeinitialdata(self):
    # --- The file is only opened and the data written on processor 0
    if me > 0: return

    fileexists = os.access(self.filename,os.F_OK)
    if self.newfile or not fileexists:
      self.f = PW.PW(self.filename)

      # --- save some initial data, and plane size and location and time step
      self.f.lsavephi   = self.lsavephi
      self.f.lsaveparticles = self.lsaveparticles
      self.f.zplane    = self.zplane
      self.f.dt        = top.dt
      self.f.deltat = self.deltat

      if self.lsaveparticles:
        self.f.npid      = top.npid
        # --- Write out particle quantities for each species
        for js in self.jslist:
          self.f.write('sq_%d'%js,top.pgroup.sq[js])
          self.f.write('sm_%d'%js,top.pgroup.sm[js])
          self.f.write('sw_%d'%js,top.pgroup.sw[js])

      if self.lsavephi:
        # --- Note that the file is already open
        self.f.deltaz    = self.deltaz
        self.f.nx_plane  = w3d.nx
        self.f.ny_plane  = w3d.ny
        self.f.ixa_plane = w3d.ix_axis
        self.f.iya_plane = w3d.iy_axis
        self.f.xmmin     = w3d.xmmin
        self.f.xmmax     = w3d.xmmax
        self.f.ymmin     = w3d.ymmin
        self.f.ymmax     = w3d.ymmax
        self.f.deltaz    = self.deltaz
     
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

      self.f.set_verbosity(0)
      self.f.close()

  def saveplane(self):

    if self.f is None:
      self.f = 0
      self.initsaveplane()
      self.initsaveparticles()
      self.initsavephi()
      self.writeinitialdata()

    # --- Only save data if zplane is within the grid.
    if(self.zplane < w3d.zmmin+top.zbeam or
       self.zplane >= w3d.zmmax+top.zbeam): return

    # --- Only save data at specified frequency
    itt = max(1,nint(self.deltat/top.dt))
    if (top.it % itt) > 0: return

    if me == 0:
      # --- open file in append mode
      self.f = PW.PW(self.filename,'a')

    # --- save for each species
    for js in self.jslist:
      self.saveplanespecies(js)

    # --- Save phi at the plane
    self.saveplanephi()

    if me == 0:
      if(self.save_this_step):
        # --- Write start time and update end time
        if(self.it==1): self.f.tmin = top.time
        self.f.tmax = top.time

      # --- close file
      self.f.set_verbosity(0)
      self.f.close()

  def saveplanephi(self):
    if not self.lsavephi: return

    # --- phi is saved every time step whether or not there are particles saved
    # --- but only after saving has started.
    if(self.save_this_step):
      # --- get the two planes of phi to be saved
      iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)
      izz = max(1,nint(self.deltaz/w3d.dz))
      self.phi_save = transpose(array([transpose(getphi(iz=iz-izz)),
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

  def saveplanespecies(self,js):
    if not self.lsaveparticles: return

    j1 = top.pgroup.ins[js] - 1
    j2 = j1 + top.pgroup.nps[js]

    if j2 > j1:
      z = top.pgroup.zp[j1:j2]
      zold = top.pgroup.pid[j1:j2,self.zoldpid]

      # --- Find all of the particles which just crossed zplane.
      ii = compress(logical_and(zold < self.zplane,self.zplane <= z),
                    iota(j1,j2))

      # --- Get the data for those particles that crossed.
      xx = gatherarray(take(top.pgroup.xp,ii))
      yy = gatherarray(take(top.pgroup.yp,ii))
      zz = gatherarray(take(top.pgroup.zp,ii))
      ux = gatherarray(take(top.pgroup.uxp,ii))
      uy = gatherarray(take(top.pgroup.uyp,ii))
      uz = gatherarray(take(top.pgroup.uzp,ii))
      gi = gatherarray(take(top.pgroup.gaminv,ii))
      id = gatherarray(take(top.pgroup.pid,ii,axis=0))

      np_save = len(xx)
      if np_save > 0: self.save_this_step = true
    else:
      np_save = 0

    self.save_this_step = broadcast(self.save_this_step)

    if (self.save_this_step):
      self.it = self.it + 1

    if (np_save > 0 and self.save_this_step and me == 0):
      suffix = '%d_%d'%(self.it,js)
      try:
        self.f.write('xp'+suffix,    xx)
        self.f.write('yp'+suffix,    yy)
        self.f.write('zp'+suffix,    zz)
        self.f.write('uxp'+suffix,   ux)
        self.f.write('uyp'+suffix,   uy)
        self.f.write('uzp'+suffix,   uz)
        self.f.write('gaminv'+suffix,gi)
        self.f.write('pid'+suffix,   id)
      except:
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "ERROR: There was an error writing out the particle data"
        print "       This most likely means that an attempt was made"
        print "       to overwrite an existing file."
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    # --- The old z can now be reset
    if j2 > j1:
      zold[:] = z

