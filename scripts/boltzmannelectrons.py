from warp import *

class BoltzmanSolver:

  def __init__(self,rhoion,te,phi0,n=1):
    self.rhoion = rhoion
    self.te = te
    self.phi0 = phi0
    self.n = n
    installafterfs(self.solve)

  def solve(self):
    print "doing boltzmann solve"
    for i in range(self.n):
      f3d.sormaxit = 1
      fieldsol(-1)
      phi = w3d.phi[:,:,1:-1]
      rhoelectron = self.rhoion*exp(-abs(self.phi0-phi)/self.te)
      const = 0.5/(1./w3d.dx**2 + 1./w3d.dy**2 + 1./w3d.dz**2)
      phi[:,:,:] = phi - rhoelectron/eps0*const
      w3d.phi[:,:,0] = w3d.phi[:,:,1]
  

