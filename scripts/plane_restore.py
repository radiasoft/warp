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
# The variable new_zplane may be redefined.  It defaults to zmmin.        #
###########################################################################
# Assumptions:
#   Only one species of particles is saved.
#   new_zplane lies exactly on a grid cell.
###########################################################################

"""
Contains PlaneRestore class that is used to restore particle and field data at a
specified z plane. The data have been saved by PlaneSave.
The two simulations are linked together.
"""
from warp import *
plane_restore_version = "$Id: plane_restore.py,v 1.5 2003/08/07 23:35:56 dave Exp $"

class PlaneRestore:
  """
Saves the particle data and phi to a file just after the 
field solve is finished. The positions and velocities are staggered.
It is automatically called every time step after the field solve
Input:
  - filename=runid.plane: filename where data is stored
  - zplane: location where simulations are linked together. Units of meters
            relative to the lab frame. Note zplane must lie on a grid cell.
  - js: species which are saved. Defaults to all species. Can be single
        integer or a list of integers.
  - l_restore_phi=1: flag for restoring phi or not.
  """

  def __init__(self,filename,zplane=None,js=None,l_restore_phi=1):

    # --- Make sure that vbeamfrm is set to zero so that the grid doesn't move
    # --- out from under the restoration plane.
    top.vbeamfrm = 0.

    # set self.zplane
    if zplane is None:
      self.zplane = w3d.zmminglobal
    else:
      self.zplane = zplane

    # Only restore if zplane within the grid
    if(self.zplane<w3d.zmmin or self.zplane>=w3d.zmmax): return

    self.l_restore_phi = l_restore_phi

    ##############################
    # Initialization stuff

    # open the file which holds the data
    self.f = PR.PR(filename)
    self.zshift = self.zplane - self.f.zplane

    # get time level of first plane and subtract 1
    self.it_restore = 1

    # set record number to 0
    self.ir_plane = 0

    # get time step, tmin, tmax
    top.dt = self.f.dt
    self.tmin = self.f.tmin
    self.tmax = self.f.tmax

    # initializes list of species
    if type(js) == IntType:
      self.jslist= [js]
    elif js is None:
      self.jslist = range(top.ns)
    else:
      self.jslist = js

    # restore particle charge, mass, weight
    for js in self.jslist:
      top.sq[js] = self.f.read('sq_%08d'%js)
      top.sm[js] = self.f.read('sm_%08d'%js)
      top.sw[js] = self.f.read('sw_%08d'%js)

    # make sure that pid will be allocated
    if(top.npmaxi==1):
      top.npmax  = max(top.npmax,2)
      top.npmaxb = max(top.npmax,2)
      top.npmaxi = max(top.npmax,2)
      gchange('Particles')

    # restore solver geometry of the saved data
    try:
      self.solvergeom = self.f.solvergeom
    except:
      self.solvergeom = w3d.XYZgeom

    # restore only if between grid bounds
    if(self.zplane<w3d.zmmin or self.zplane>=w3d.zmmax): return
    
    # set up indices which specify transverse extent of saved and restored phi
    # _r for restored phi array, _s for saved phi array
    # '0' is minimum index, 'm' is maximum index
    self.nx0_r = max(0, 0 - self.f.ixa_plane + w3d.ix_axis)
    self.ny0_r = max(0, 0 - self.f.iya_plane + w3d.iy_axis)
    self.nxm_r = min(w3d.nx, self.f.nx_plane - self.f.ixa_plane + w3d.ix_axis)
    self.nym_r = min(w3d.ny, self.f.ny_plane - self.f.iya_plane + w3d.iy_axis)
    self.nx0_s = self.nx0_r - w3d.ix_axis + self.f.ixa_plane
    self.ny0_s = self.ny0_r - w3d.iy_axis + self.f.iya_plane
    self.nxm_s = self.nxm_r - w3d.ix_axis + self.f.ixa_plane
    self.nym_s = self.nym_r - w3d.iy_axis + self.f.iya_plane

    # deal with symmetries
    # if saved is 2 or 4 fold symmetric and restored isn't, lower half of restored
    # is filled with inverted saved phi
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

    installbeforefs(self.restoreplane_bfs)
    installafterfs(self.restoreplane_afs)

  ###########################################################################
  def disable_plane_restore(self):
    # for some reason, does not work!
    uninstallbeforefs(self.restoreplane_bfs)
    uninstallafterfs(self.restoreplane_afs)

  ###########################################################################
  # restore the next plane of data
  def restoreplane_bfs(self):
    # increment the timelevel of the plane
    self.it_restore = self.it_restore + 1
    it = self.it_restore

    # load particles for each species
    for js in self.jslist:
      self.restoreplane(js,it)

    # calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # load saved phi into the phi array
    try:
      self.restore_phi(iz,it)
    except:
      return

  def restoreplane(self,js=0,it=0):

    # get number of particles to be added and the starting point
    try:
      nn = self.f.read('np%08d_%08d'%(it,js))
    except:
      return
    
    # put restored data into particle arrays, adjusting the z location
    if (nn > 0):
      addparticles(self.f.read('xp%08d_%08d'%(it,js)),
                   self.f.read('yp%08d_%08d'%(it,js)),
                   self.f.read('zp%08d_%08d'%(it,js))+self.zshift,
                   self.f.read('uxp%08d_%08d'%(it,js)),
                   self.f.read('uyp%08d_%08d'%(it,js)),
                   self.f.read('uzp%08d_%08d'%(it,js)),
                   self.f.read('gaminv%08d_%08d'%(it,js)),
                   self.f.read('pid%08d_%08d'%(it,js)),
                   js,
                   lallindomain=true,
                   lmomentum=true,
                   resetrho=false)

  # this routine resets the potential at the plane iz=-1 after the field solve
  # if this is needed
  def restoreplane_afs(self):
    # calculate grid location of new_plane
    iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)

    # reset phi at plane iz=-1 if zplane is at iz=0
    if (iz == 0):
      try:
        self.restore_phi(iz,self.it_restore)
      except:
        return

  #######################################################################
  def restore_phi(self,iz,it):
    # return if flag indicates phi not to be restored
    if self.l_restore_phi is 0: return
    
    if self.solvergeom == w3d.solvergeom:
      restore_phi_3d_to_3d(iz,it)
    elif self.solvergeom == w3d.RZgeom and w3d.solvergeom == w3d.XYZgeom:
      restore_phi_rz_to_3d(iz,it)
    # Required since only phi is changed and not phip
    getphiforparticles

  #######################################################################
  # This routine copies the saved phi plane into the current phi array
  # making use of different numbers of grid cells and differing symmetries.
  # Both saved and restored phi are 3-D.
  def restore_phi_3d_to_3d(self,iz,it):
    savedphi = self.f.read('phiplane%08d'%it)
    for i in range(2):
      grid2grid(w3d.phi[self.nx0_r:self.nxm_r+1,self.ny0_r:self.nym_r+1,iz-1+i],
                w3d.nx,w3d.ny,
                w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                savedphi[:,:,i],self.f.nx_plane,self.f.ny_plane,
                self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)
      
    if ((self.f.sym_plane == 2 and (not w3d.l2symtry and not w3d.l4symtry)) or
        (self.f.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry))):
    #     phi(self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz-1:iz)=
    #      self.f.phi_plane(nx0_s:nxm_s,nym_s2:ny0_s2:-1,)
      for i in range(2):
        grid2grid(w3d.phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz-1+i],
                  w3d.nx,w3d.ny,
                  w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                  savedphi[:,::-1,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

    if ((self.f.sym_plane == 4 and ( w3d.l2symtry and not w3d.l4symtry))):
    #  phi(nx0_r2:nxm_r2,ny0_r:nym_r,iz-1:iz)=
    #    self.f.phi_plane(nx0_s2:nxm_s2:-1,ny0_s:nym_s,)
      for i in range(2):
        grid2grid(w3d.phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz-1+i],
                  w3d.nx,w3d.ny,
                  w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                  savedphi[::-1,:,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

    if ((self.f.sym_plane == 4 and (not w3d.l2symtry and not w3d.l4symtry))):
    #  phi(nx0_r2:nxm_r2,ny0_r2:nym_r,iz-1:iz)=
    #    self.f.phi_plane(nx0_s2:nxm_s2:-1,ny0_s2:nym_s,)
      for i in range(2):
        grid2grid(w3d.phi[self.nx0_r:self.nxm_r,self.ny0_r2:self.nym_r2,iz-1+i],
                  w3d.nx,w3d.ny,
                  w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
                  savedphi[::-1,:,i],self.f.nx_plane,self.f.ny_plane,
                  self.f.xmmin,self.f.xmmax,self.f.ymmin,self.f.ymmax)

  #######################################################################
  # This routine copies the saved phi plane into the current phi array
  # making use of different numbers of grid cells and differing symmetries.
  # Saved phi is rz and restored phi is 3-D.
  def restore_phi_rz_to_3d(self,iz,it):
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

    savedphi = self.f.read('phiplane%08d'%it)
    i1 = self.nx0_r
    i2 = self.nxm_r+1
    j1 = self.ny0_r
    j2 = self.nym_r+1
    for i in range(2):
      w3d.phi[i1:i2,j1:j2,iz-1+i] = (
          take(savephi[:,0,i],self.irmesh  )*(1.-self.wrmesh) + 
          take(savephi[:,0,i],self.irmesh+1)*self.wrmesh)

