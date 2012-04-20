"""
Contains PlaneSave class that is used to save particle and field data at a
specified z plane. The data is used by PlaneRestore to continue the
simulation. The two simulations are linked together.
"""

__all__ = ['PlaneSave']

from warp import *
import cPickle
plane_save_version = "$Id: plane_save.py,v 1.24 2009/06/29 18:24:38 dave Exp $"

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
      self.started_saving_particles = True
    else:
      self.started_saving_particles = False

    # --- defines useful variables
    self.it = 0

    # --- Set so data is saved to file immediately after a field solve.
    installafterfs(self.saveplane)

    # --- Set the name of the file which will hold the data
    if filename is None:
      self.filename = arraytostr(top.runid)+".plane"
    else:
      self.filename = filename

    self.initted = False

  def write(self,name,val):
    # --- Check if _f is defined. If not, the create the file as needed.
    try:
      self._f
    except AttributeError:

      fileexists = os.access(self.filename,os.F_OK)
      if self.newfile or not fileexists:
        self._f = open(self.filename,'w')
        self.writeinitialdata()
      else:
        self._f = None

    if self._f is None or self._f.closed:
      self._f = open(self.filename,'a')

    # --- Write the data out as a named tuple
    cPickle.dump((name,val),self._f,-1)

  def flush(self):
    if me > 0: return
    if self._f is not None and not self._f.closed:
      self._f.flush()

  def fileclose(self):
    if me > 0: return
    try:
      if not self._f.closed:
        self._f.set_verbosity(0)
        self._f.close()
    except AttributeError:
      pass

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

    # --- save some initial data, and plane size and location and time step
    self.write('lsavephi',self.lsavephi)
    self.write('lsaveparticles',self.lsaveparticles)
    self.write('zplane',self.zplane)
    self.write('dt',top.dt)
    self.write('deltat',self.deltat)
    self.write('tmin',top.time)

    if self.lsaveparticles:
      self.write('npid',top.npid)
      # --- Write out particle quantities for each species
      for js in self.jslist:
        self.write('sq_%d'%js,top.pgroup.sq[js])
        self.write('sm_%d'%js,top.pgroup.sm[js])
        self.write('sw_%d'%js,top.pgroup.sw[js])

    if self.lsavephi:
      # --- Note that the file is already open
      self.write('deltaz',self.deltaz)
      self.write('nx_plane',w3d.nx)
      self.write('ny_plane',w3d.ny)
      self.write('ixa_plane',w3d.ix_axis)
      self.write('iya_plane',w3d.iy_axis)
      self.write('xmmin',w3d.xmmin)
      self.write('xmmax',w3d.xmmax)
      self.write('ymmin',w3d.ymmin)
      self.write('ymmax',w3d.ymmax)
      self.write('deltaz',self.deltaz)

      # --- set sym_plane and write it out
      if (w3d.l4symtry):
        sym_plane = 4
      elif (w3d.l2symtry):
        sym_plane = 2
      else:
        sym_plane = 1
      self.write('sym_plane',sym_plane)

      # --- Write out the solver geometry flag
      # --- This must be the last thing written out before the phi or particle
      # --- data, since it is used to flag the end of the initial data.
      self.write('solvergeom',w3d.solvergeom)

  def saveplane(self):

    if not self.initted:
      self.initsaveplane()
      self.initsaveparticles()
      self.initsavephi()
      self.initted = True

    # --- Only save data if zplane is within the grid.
    if(self.zplane < w3d.zmmin+top.zbeam or
       self.zplane >= w3d.zmmax+top.zbeam): return

    # --- Only save data at specified frequency
    itt = max(1,nint(self.deltat/top.dt))
    if (top.it % itt) > 0: return

    # --- save for each species
    for js in self.jslist:
      self.saveplanespecies(js)

    # --- Save phi at the plane
    self.saveplanephi()

    # --- Make sure the data is written out to the file
    if self.started_saving_particles:
      self.it += 1
      self.flush()

  def saveplanephi(self):
    if not self.lsavephi: return

    # --- phi is saved every time step whether or not there are particles saved
    # --- but only after saving has started.
    if self.started_saving_particles:
      # --- get the two planes of phi to be saved
      iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)
      izz = max(1,nint(self.deltaz/w3d.dz))
      self.phi_save = transpose(array([transpose(getphi(iz=iz-izz)),
                                       transpose(getphi(iz=iz))]))
      if me == 0:
        self.write('phiplane%09d'%self.it,self.phi_save)

  def saveplanespecies(self,js):
    if not self.lsaveparticles: return

    if top.pgroup.nps[js] > 0:

      j1 = top.pgroup.ins[js] - 1
      j2 = j1 + top.pgroup.nps[js]

      z = top.pgroup.zp[j1:j2]
      zold = top.pgroup.pid[j1:j2,self.zoldpid]

      # --- Find all of the particles which just crossed zplane.
      ii = logical_and(zold < self.zplane,self.zplane <= z)

      # --- Get the data for those particles that crossed the zplane.
      xx = top.pgroup.xp[j1:j2][ii]
      yy = top.pgroup.yp[j1:j2][ii]
      zz = top.pgroup.zp[j1:j2][ii]
      ux = top.pgroup.uxp[j1:j2][ii]
      uy = top.pgroup.uyp[j1:j2][ii]
      uz = top.pgroup.uzp[j1:j2][ii]
      gi = top.pgroup.gaminv[j1:j2][ii]
      id = top.pgroup.pid[j1:j2,:][ii,:]

      # --- The old z can now be reset
      zold[:] = z

    else:

      # --- If there are no particles, the particle arrays may be
      # --- unallocated, so just use zeros.
      xx = zeros(0,'d')
      yy = zeros(0,'d')
      zz = zeros(0,'d')
      ux = zeros(0,'d')
      uy = zeros(0,'d')
      uz = zeros(0,'d')
      gi = zeros(0,'d')
      id = zeros((0,top.npid),'d')

    # --- Gather the data from all of the processors (really only the one
    # --- or two where the saving plane is).
    xx = gatherarray(xx)
    yy = gatherarray(yy)
    zz = gatherarray(zz)
    ux = gatherarray(ux)
    uy = gatherarray(uy)
    uz = gatherarray(uz)
    gi = gatherarray(gi)
    id = gatherarray(id)

    np_save = len(xx)

    if not self.started_saving_particles:
      self.started_saving_particles = (globalmax(np_save) > 0)

    if not self.started_saving_particles:
      return

    if np_save > 0 and me == 0:
      suffix = '%09d_%d'%(self.it,js)
      self.write('xp'+suffix,xx)
      self.write('yp'+suffix,yy)
      self.write('zp'+suffix,zz)
      self.write('uxp'+suffix,ux)
      self.write('uyp'+suffix,uy)
      self.write('uzp'+suffix,uz)
      self.write('gaminv'+suffix,gi)
      self.write('pid'+suffix,id)


