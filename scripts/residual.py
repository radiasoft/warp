# Calculate the residual of Poisson's equation.  The residual has the units
# of Q/M**3 (same units as rho)
from warp import *
residual_version = "$Id: residual.py,v 1.3 2004/06/01 20:13:31 dave Exp $"

def getresidual(phi=None,rho=None,dx=None,dy=None,dz=None):
  if phi is None: phi = w3d.phi
  if rho is None: rho = w3d.rho
  if dx is None: dx = w3d.dx
  if dy is None: dy = w3d.dy
  if dz is None: dz = w3d.dz
  residual = zeros(shape(rho),"d")

  residual[1:-1,1:-1,1:-1] = (eps0*(
     (phi[0:-2,1:-1,2:-2] -
   2.*phi[1:-1,1:-1,2:-2] +
      phi[2:  ,1:-1,2:-2])/dx**2 +
     (phi[1:-1,0:-2,2:-2] -
   2.*phi[1:-1,1:-1,2:-2] +
      phi[1:-1,2:  ,2:-2])/dy**2 +
     (phi[1:-1,1:-1,1:-3] -
   2.*phi[1:-1,1:-1,2:-2] +
      phi[1:-1,1:-1,3:-1])/dz**2)+
      rho[1:-1,1:-1,1:-1])

  return residual

#print "Maximum residual = %e" % max(abs(ravel(residual)))
