"""
Contains PlaneSave class that is used to save particle and field data at a
specified z plane. The data is used by PlaneRestore to continue the
simulation. The two simulations are linked together.
"""
from warp import *
plane_save_version = "$Id: plane_save.py,v 1.1 2002/11/26 23:31:51 dave Exp $"

class PlaneSave:
  """
Saves the particle data and phi to a file just after the 
field solve is finished. The positions and velocities are staggered.
It is automatically called every time step after the field solve
Input:
  - zplane: location where simulations are linked together. Units of meters
            relative to the lab frame. Note zplane must lie on a grid cell.
  - filename=runid.plane: filename where data is stored
  - js: species which are saved. Defaults to all species. Can be single
        integer or a list of integers.
  """

  def __init__(self,zplane,filename=None,js=None):

    self.zplane = zplane

    # --- Arrays to find particles which cross the plane
    self.old_zp = zeros(top.zp.shape[0],'d')
    self.ii_zp  = arange(top.zp.shape[0])

    # --- Setup region of phi to be saved
    self.nx0 = 0
    self.nxm = w3d.nx
    self.ny0 = 0
    self.nym = max(w3d.ny,self.ny0+1)
    self.ixa_plane = w3d.ix_axis
    self.iya_plane = w3d.iy_axis
    self.nx_plane = self.nxm - self.nx0
    self.ny_plane = self.nym - self.ny0

    # --- Array used to save the potential
    self.phi_plane = zeros([self.nx_plane+1,self.ny_plane+1,2],'d')
    self.phi_plane = w3d.phi[self.nx0:self.nxm+1,self.ny0:self.nym+1,0:2]

    # create the file which will hold the data
    if filename is None:
      self.filename = arraytostr(top.runid)+".plane"
    else:
      self.filename = filename
    self.f = PW.PW(self.filename)

    # save plane size and location and time step
    self.f.zplane    = zplane
    self.f.nx_plane  = self.nx_plane
    self.f.ny_plane  = self.ny_plane
    self.f.ixa_plane = self.ixa_plane
    self.f.iya_plane = self.iya_plane
    self.f.dt        = top.dt

    # set sym_plane and write it out
    if (w3d.l4symtry):
      sym_plane = 4
    elif (w3d.l2symtry):
      sym_plane = 2
    else:
      sym_plane = 1
    self.f.sym_plane = sym_plane

    # initializes list of species
    if type(js) == IntType:
      self.jslist= [js]
    elif js is None:
      self.jslist = range(top.ns)
    else:
      self.jslist = js
    for js in self.jslist:
      self.f.write('sq_%08d'%js,top.sq[js])
      self.f.write('sm_%08d'%js,top.sm[js])
      self.f.write('sw_%08d'%js,top.sw[js])

    # defines useful variables
    self.it = 0
    self.np_save_tot = 0

    self.f.set_verbosity(0)
    self.f.close()

    installafterfs(self.saveplane)

  def saveplane(self):
    if(self.old_zp.shape <> top.zp.shape):
       self.old_zp = zeros(top.zp.shape[0],'d')
       self.ii_zp  = arange(top.zp.shape[0])

    # open file in append mode
    self.f = PW.PW(self.filename,'a')

    # save for each species
    for js in self.jslist:
      self.saveplanespecies(js)

    # close file
    self.f.set_verbosity(0)
    self.f.close()

  def saveplanespecies(self,js):

    # only save any data if there are particles
    if (top.nps[js] == 0):return

    # set range for particle arrays
    ipl = top.ins[js]-1
    ipu = top.ins[js]+top.nps[js]-1
  
    # Recalculate the previous axial location of the particles to check
    # if the particles crossed the saving plane.
    self.old_zp[ipl:ipu] = top.zp[ipl:ipu] - top.uzp[ipl:ipu]*top.gaminv[ipl:ipu]*top.dt

    # find indices of particles which crossed the plane in the last timestep
    ii = compress((self.old_zp[ipl:ipu]<=self.zplane) & (top.zp[ipl:ipu]>=self.zplane), self.ii_zp[ipl:ipu])

    np_save = shape(ii)[0]
    self.np_save_tot = self.np_save_tot + np_save

    # only save anything if there are particles to the right of zplane
    if(self.np_save_tot>0):

      # increment saving counter
      self.it = self.it + 1

      # save number of particles
      self.f.write('np%08d'%self.it,np_save)

      # if there are particles that crossed the plane, save the data
      if (np_save > 0):
        self.f.write('xp%08d_%08d'%(self.it,js),    take(top.xp,ii))
        self.f.write('yp%08d_%08d'%(self.it,js),    take(top.yp,ii))
        self.f.write('zp%08d_%08d'%(self.it,js),    take(top.zp,ii))
        self.f.write('uxp%08d_%08d'%(self.it,js),   take(top.uxp,ii))
        self.f.write('uyp%08d_%08d'%(self.it,js),   take(top.uyp,ii))
        self.f.write('uzp%08d_%08d'%(self.it,js),   take(top.uzp,ii))
        self.f.write('gaminv%08d_%08d'%(self.it,js),take(top.gaminv,ii))
        for id in range(top.npid):
          self.f.write('pid%08d_%08d_%08d'%(self.it,js,top.npid),take(top.pid,ii))

      # phi is saved every time step whether or not there are particles saved
      # get the two planes of phi to be saved
      iz = nint((self.zplane - top.zbeam - w3d.zmmin)/w3d.dz)
      self.f.write('phiplane%08d'%self.it,w3d.phi[self.nx0:self.nxm,self.ny0:self.nym,iz-1:iz])
