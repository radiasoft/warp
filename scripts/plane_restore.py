"""
Contains PlaneRestore class that is used to restore particle and field data at a
specified z plane. The data have been saved by PlaneSave.
The two simulations are linked together.

###########################################################################
# Routines for restoring particle and phi data at a plane.
#
# beforefs restores the next plane.  It is automatically called every time
#          step before the field solve.
#
# afterfs resets the plane at iz-1 if iz=0 (that plane is changed in a field
# solve)
#
###########################################################################
# Note:                                                                   #
# The string plane_file must be set to the name of the file holding the   #
# plane data.                                                             #
# The variable new_zplane may be redefined.  It defaults to zmmin.  #
###########################################################################
# Assumptions:
#   new_zplane lies exactly on a grid cell.
###########################################################################
"""

from warp import *
plane_restore_version = "$Id: plane_restore.py,v 1.11 2007/12/12 23:39:01 dave Exp $"

class PlaneRestore:
  """
Saves the particle data and phi to a file just after the 
field solve is finished. The positions and velocities are staggered.
It is automatically called every time step after the field solve
Input:
  - filename=runid.plane: filename where data is stored
  - zplane: location where simulations are linked together. Units of meters
            relative to the lab frame. Defaults to w3d.zmmin.
  - js: species which are saved. Defaults to all species. Can be single
        integer or a list of integers.
  - l_restore_phi=1: flag for restoring phi or not.
  - lrestoreparticles=1: flag for restoring particles
  """

  def __init__(self,filename,zplane=None,js=None,
               l_restore_phi=1,lrestoreparticles=1):

    # --- Save some input values
    self.filename = filename
    self.zplane = zplane
    self.js = js
    self.l_restore_phi = l_restore_phi
    self.lrestoreparticles = lrestoreparticles

    # --- Install the routines that do the work.
    installbeforefs(self.restoreplane_bfs)
    installafterfs(self.restoreplane_afs)

    # --- Initialize self.f to None to flag that the initialization
    # --- has not yet happened.
    self.f = None

  def initrestoreplane(self):

    ##############################
    # --- Initialization stuff

    if self.zplane is None: self.zplane = w3d.zmmin

    # --- open the file which holds the data
    self.f = PR.PR(self.filename)
    self.zshift = self.zplane - self.f.zplane

    self.lsavephi = self.f.lsavephi
    self.lsaveparticles = self.f.lsaveparticles
    if not self.lsavephi: self.l_restore_phi = false
    if not self.lsaveparticles: self.lrestoreparticles = false

    # --- get time level of first plane and subtract 1
    self.it_restore = 0

    # --- get time step, tmin, tmax
    if top.dt != self.f.dt:
      print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      print "Warning: timestep size not the same as in saved run, top.dt will"
      print "         be adjusted to the value from the saved run"
      print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      top.dt = self.f.dt
    self.tmin = self.f.tmin
    self.tmax = self.f.tmax

    if self.lrestoreparticles:
      # --- initializes list of species
      if self.js is None:
        self.jslist = range(top.ns)
      else:
        try:
          list(self.js)
          self.jslist = self.js
        except TypeError:
          self.jslist= [self.js]

      # --- restore particle charge, mass, weight
      for js in self.jslist:
        top.pgroup.sq[js] = self.f.read('sq_%d'%js)
        top.pgroup.sm[js] = self.f.read('sm_%d'%js)
        top.pgroup.sw[js] = self.f.read('sw_%d'%js)

      # --- make sure that pid will be allocated
      top.npid = self.f.npid
      setuppgroup(top.pgroup)

    if self.l_restore_phi:
      # --- restore solver geometry of the saved data
      try:
        self.solvergeom = self.f.solvergeom
      except:
        self.solvergeom = w3d.XYZgeom

      # set up indices which specify transverse extent of saved and restored phi
      # _r for restored phi array, _s for saved phi array
      # '0' is minimum index, 'm' is maximum index

      #self.nx0_r = max(0,int(floor((self.f.xmmin - w3d.xmmin)/w3d.dx)))
      #self.ny0_r = max(0,int(floor((self.f.ymmin - w3d.ymmin)/w3d.dy)))
      #self.nxm_r = min(w3d.nx,int(floor((self.f.xmmax - w3d.xmmin)/w3d.dx)))
      #self.nym_r = min(w3d.ny,int(floor((self.f.ymmax - w3d.ymmin)/w3d.dy)))

      self.nx0_r = max(0, 0 - self.f.ixa_plane + w3d.ix_axis)
      self.ny0_r = max(0, 0 - self.f.iya_plane + w3d.iy_axis)
      self.nxm_r = min(w3d.nx, self.f.nx_plane - self.f.ixa_plane + w3d.ix_axis)
      self.nym_r = min(w3d.ny, self.f.ny_plane - self.f.iya_plane + w3d.iy_axis)
      self.nx0_s = self.nx0_r - w3d.ix_axis + self.f.ixa_plane
      self.ny0_s = self.ny0_r - w3d.iy_axis + self.f.iya_plane
      self.nxm_s = self.nxm_r - w3d.ix_axis + self.f.ixa_plane
      self.nym_s = self.nym_r - w3d.iy_axis + self.f.iya_plane

      # --- deal with symmetries
      # --- if saved is 2 or 4 fold symmetric and restored isn't,
      # --- lower half of restored is filled with inverted saved phi
      if ((self.f.sym_plane == 2 and (not w3d.l2symtry and not w3d.l4symtry)) or
          (self.f.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry))): 
        self.ny0_r2 = max(0, - self.f.ny_plane - self.f.iya_plane + w3d.iy_axis)
        self.nym_r2 = min(w3d.ny, 0 - self.f.iya_plane + w3d.iy_axis)
        self.ny0_s2 = - self.ny0_r + w3d.iy_axis + self.f.iya_plane
        self.nym_s2 =   self.nym_r - w3d.iy_axis + self.f.iya_plane
      if ((self.f.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry)) or
          (self.f.sym_plane == 4 and (    w3d.l2symtry and not w3d.l4symtry))):
        self.nx0_r2 = max(0, - self.f.nx_plane - self.f.ixa_plane + w3d.ix_axis)
        self.nxm_r2 = min(w3d.nx, 0 - self.f.ixa_plane + w3d.ix_axis)
        self.nx0_s2 = self.nxm_r - w3d.ix_axis + self.f.ixa_plane
        self.nxm_s2 = - self.nx0_r + w3d.ix_axis + self.f.ixa_plane

  ###########################################################################
  def disable_plane_restore(self):
    # for some reason, does not work!
    uninstallbeforefs(self.restoreplane_bfs)
    uninstallafterfs(self.restoreplane_afs)

  ###########################################################################
  # --- restore the next plane of data
  def restoreplane_bfs(self):

    # --- Do the initialization if it hasn't been done yet.
    if self.f is None:
      self.initrestoreplane()

    # --- restore only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- increment the timelevel of the plane
    self.it_restore = self.it_restore + 1
    it = self.it_restore

    # --- load particles for each species
    for js in self.jslist:
      self.restoreplaneparticles(js,it)

    # --- calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # --- load saved phi into the phi array
    try:
      self.restore_phi(iz,it)
    except:
      return

  def restoreplaneparticles(self,js=0,it=0):
    if not self.lrestoreparticles: return

    # --- put restored data into particle arrays, adjusting the z location
    # --- Check if data was written for this step.
    suffix = '%d_%d'%(it,js)
    if 'xp'+suffix in self.f.inquire_names():
      xx = self.f.read('xp'+suffix)
      yy = self.f.read('yp'+suffix)
      zz = self.f.read('zp'+suffix)+self.zshift
      ux = self.f.read('uxp'+suffix)
      uy = self.f.read('uyp'+suffix)
      uz = self.f.read('uzp'+suffix)
      gi = self.f.read('gaminv'+suffix)
      id = self.f.read('pid'+suffix)
      addparticles(xx,yy,zz,ux,uy,uz,gi,id,
                   js=js,
                   lallindomain=false,
                   lmomentum=true,
                   resetrho=false)

  def restoreplane_afs(self):
    # --- this routine resets the potential at the plane iz=-1 after the 
    # --- field solve if this is needed

    # --- restore only if between grid bounds
    if (self.zplane < w3d.zmmin+top.zbeam or
        self.zplane+top.zbeam >= w3d.zmmax): return

    # --- calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # --- reset phi at plane iz=-1 if zplane is at iz=0
    if (iz == 0):
      self.restore_phi(iz,self.it_restore)

  #######################################################################
  def restore_phi(self,iz,it):
    # --- return if flag indicates phi not to be restored
    if self.l_restore_phi is 0: return
    
    # --- Read in the phi data if it is available
    if 'phiplane%d'%it not in self.f.inquire_names(): return
    savedphi = self.f.read('phiplane%d'%it)

    if self.solvergeom == w3d.solvergeom:
      self.restore_phi_3d_to_3d(iz-top.izfsslave[me],it,savedphi,w3d.phi)
      if top.izpslave[me] != top.izfsslave[me]:
        self.restore_phi_3d_to_3d(iz-top.izpslave[me],it,savedphi,w3d.phip)
    elif self.solvergeom == w3d.RZgeom and w3d.solvergeom == w3d.XYZgeom:
      self.restore_phi_rz_to_3d(iz,it,savedphi,w3d.phi)
      if top.izpslave[me] != top.izfsslave[me]:
        self.restore_phi_rz_to_3d(iz,it,savedphi,w3d.phip)

  #######################################################################
  # This routine copies the saved phi plane into the current phi array
  # making use of different numbers of grid cells and differing symmetries.
  # Both saved and restored phi are 3-D.
  def restore_phi_3d_to_3d(self,iz,it,savedphi,phi):
    if iz < 0 or iz > w3d.nzlocal: return
    for i in range(2):
      grid2grid(phi[self.nx0_r:self.nxm_r+1,self.ny0_r:self.nym_r+1,iz+i],
                w3d.nx,w3d.ny,
                w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                savedphi[:,:,i],self.f.nx_plane,self.f.ny_plane,
                self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)
      
    if ((self.f.sym_plane == 2 and (not w3d.l2symtry and not w3d.l4symtry)) or
        (self.f.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry))):
    #     phi(self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz-1:iz)=
    #      self.f.phi_plane(nx0_s:nxm_s,nym_s2:ny0_s2:-1,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  w3d.nx,w3d.ny,
                  w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                  savedphi[:,::-1,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

    if ((self.f.sym_plane == 4 and ( w3d.l2symtry and not w3d.l4symtry))):
    #  phi(nx0_r2:nxm_r2,ny0_r:nym_r,iz-1:iz)=
    #    self.f.phi_plane(nx0_s2:nxm_s2:-1,ny0_s:nym_s,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  w3d.nx,w3d.ny,
                  w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                  savedphi[::-1,:,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

    if ((self.f.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry))):
    #  phi(nx0_r2:nxm_r2,ny0_r2:nym_r,iz-1:iz)=
    #    self.f.phi_plane(nx0_s2:nxm_s2:-1,ny0_s2:nym_s,)
      for i in range(2):
        grid2grid(phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz+i],
                  w3d.nx,w3d.ny,
                  w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                  savedphi[::-1,:,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

  #######################################################################
  # This routine copies the saved phi plane into the current phi array
  # making use of different numbers of grid cells and differing symmetries.
  # Saved phi is rz and restored phi is 3-D.
  def restore_phi_rz_to_3d(self,iz,it,savedphi,phi):
    if iz < 0 or iz > w3d.nzlocal: return
    try:
      self._rz_to_3d_inited
    except:
      self._rz_to_3d_inited = 1
      xmmin = w3d.xmmin + w3d.dx*self.nx0_r
      nx = self.nxm_r - self.nx0_r
      ymmin = w3d.ymmin + w3d.dy*self.ny0_r
      nx = self.nym_r - self.ny0_r
      xmesh,ymesh = getmesh2d(xmmin,w3d.dx,nx,ymmin,w3d.dy,ny)
      rmesh = sqrt(xmesh**2 + ymesh**2)
      dr = (self.f.xmmax - self.f.xmmin)/self.f.nx
      self.irmesh = int(rmesh/dr)
      self.wrmesh =     rmesh/dr  - self.irmesh
      self.wrmesh = where(self.irmesh >= self.f.nx,1,self.wrmesh)
      self.irmesh = where(self.irmesh >= self.f.nx,self.f.nx-1,self.irmesh)

    i1 = self.nx0_r
    i2 = self.nxm_r+1
    j1 = self.ny0_r
    j2 = self.nym_r+1
    for i in range(2):
      phi[i1:i2,j1:j2,iz+i] = (
          take(savephi[:,0,i],self.irmesh  )*(1.-self.wrmesh) + 
          take(savephi[:,0,i],self.irmesh+1)*self.wrmesh)

