"""Plasma source model.
"""
from warp import *
from egun_like import gun
plasmasource_version = "$Id: plasmasource.py,v 1.1 2004/03/26 16:00:40 dave Exp $"

def plasmasourcedoc():
  import plasmasource
  print plasmasource.__doc__

class PlasmaSource:
  """
  - electrontemperature: in units of eV
  """
  def __init__(self,electrondensity,electrontemperature,
               plasmazmin,plasmazmax,plasmaxmin,plasmaxmax,
               sourcevolt):
    self.electrondensity = electrondensity
    self.electrontemperature = electrontemperature
    self.plasmazmin = plasmazmin
    self.plasmazmax = plasmazmax
    self.plasmaxmin = plasmaxmin
    self.plasmaxmax = plasmaxmax
    self.dx = w3d.dx
    self.dz = w3d.dz
    self.nz = int((self.plasmazmax - self.plasmazmin)/self.dz)
    self.nx = int((self.plasmaxmax - self.plasmaxmin)/self.dx)
    self.sourcevolt = sourcevolt
    self.gunstarted = false

  def go(self):
    # --- This version seems to work, starting with phi in plasma constant.

    if not self.gunstarted:
      frz.basegrid.phi[:,:int(-w3d.zmmin/w3d.dz)] = self.sourcevolt
      self.gunstarted = true

    gun(maxtime=1.,ipstep=3)
    self.rhoi = frz.basegrid.rho[0:self.nx+1,0:self.nz+1] + 0.
    frz.basegrid.rho[:,:] = where(frz.basegrid.phi[1:-1,1:-1]>self.sourcevolt,
                                  0.,frz.basegrid.rho)
    fieldsol(-1)
    self.dv = self.sourcevolt - frz.basegrid.phi[1:self.nx+2,1:self.nz+2]
    self.mb = exp(-minimum(100.,maximum(0.,self.dv)/self.electrontemperature))
    self.newrho = self.rhoi*(1. - self.mb)
    frz.basegrid.rho[0:self.nx+1,0:self.nz+1] = self.newrho
    fieldsol(-1)
    frz.basegrid.phi[:,:] = where(frz.basegrid.phi[:,:]>self.sourcevolt,
                                  self.sourcevolt,frz.basegrid.phi)
    fma()
    self.pplot()

  def pplot():
    w3d.phi[:,0,:] = frz.basegrid.phi[1:-1,:]
    ppgeneric(gridt=frz.basegrid.phi[1:self.nx+2,1:self.nz+2],
              xmin=self.plasmazmin,xmax=self.plasmazmin+self.nz*self.dz,
              ymin=self.plasmaxmin,ymax=self.plasmaxmin+self.nx*self.dx,
              cellarray=1)
    ppgeneric(gridt=frz.basegrid.phi[1:self.nx+2,1:self.nz+2],
              xmin=self.plasmazmin,xmax=self.plasmazmin+self.nz*self.dz,
              ymin=self.plasmaxmin,ymax=self.plasmaxmin-self.nx*self.dx,
              cellarray=1)
    pfzx(contours=100)
    ppzr()

