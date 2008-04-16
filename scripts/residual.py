# Calculate the residual of Poisson's equation.  The residual has the units
# of volts, same as phi.
from warp import *
residual_version = "$Id: residual.py,v 1.6 2008/04/16 22:00:13 dave Exp $"

#print "Maximum residual = %e" % max(abs(ravel(residual)))

def getresidual(phi=None,rho=None,dx=None,dy=None,dz=None,
                iondensity=None,electrontemperature=None,plasmapotential=None):
  if phi is None: phi = w3d.phi
  if rho is None: rho = w3d.rho
  if dx is None: dx = w3d.dx
  if dy is None: dy = w3d.dy
  if dz is None: dz = w3d.dz
  if iondensity is None: iondensity = w3d.iondensity
  if electrontemperature is None: electrontemperature = w3d.electrontemperature
  if plasmapotential is None: plasmapotential = w3d.plasmapotential

  # --- Get an array that adds guard cells transversely.
  pshape = list(shape(phi))
  pshape[0] = pshape[0] + 2
  pshape[1] = pshape[1] + 2
  ppp = zeros(pshape,'d')
  ppp[1:-1,1:-1,:] = phi

  # --- Set the guard cells appropriately for the boundary conditions.
  # --- The guards cells in z are assumed already set.
  if w3d.boundxy == top.neumann:
    ppp[0,:,:] = ppp[2,:,:]
    ppp[-1,:,:] = ppp[-3,:,:]
    ppp[:,0,:] = ppp[:,2,:]
    ppp[:,-1,:] = ppp[:,-3,:]
  elif w3d.boundxy == top.periodic:
    ppp[0,:,:] = ppp[-3,:,:]
    ppp[-1,:,:] = ppp[2,:,:]
    ppp[:,0,:] = ppp[:,-3,:]
    ppp[:,-1,:] = ppp[:,2,:]

  # --- First calculate del squared phi plus rho/eps0
  eee = (
   (ppp[0:-2,1:-1,1:-1] + ppp[2:,1:-1,1:-1] - 2.*ppp[1:-1,1:-1,1:-1])/dx**2
  +(ppp[1:-1,0:-2,1:-1] + ppp[1:-1,2:,1:-1] - 2.*ppp[1:-1,1:-1,1:-1])/dy**2
  +(ppp[1:-1,1:-1,0:-2] + ppp[1:-1,1:-1,2:] - 2.*ppp[1:-1,1:-1,1:-1])/dz**2
  +rho/eps0)

  # --- If Boltzmann electron terms are given, then include the electrons
  # --- charge density.
  if iondensity != 0. and electrontemperature != 0.:
    eee = eee - (iondensity*exp(
                 (ppp[1:-1,1:-1,1:-1]-plasmapotential)/electrontemperature)
                 /eps0)


  # --- Rescale the error by the coefficient of ppp[1:-1,1:-1,1:-1].
  eee = eee/(2./dx**2+2./dy**2+2./dz**2)

  # --- Clear out the boundaries which are Dirichlet.
  if w3d.boundxy == top.dirichlet:
    eee[0,:,:] = 0.
    eee[-1,:,:] = 0.
    eee[:,0,:] = 0.
    eee[:,-1,:] = 0.
  if w3d.bound0 == top.dirichlet:
    eee[:,:,0] = 0.
  if w3d.boundnz == top.dirichlet:
    eee[:,:,-1] = 0.

  return eee

def getresidualXZ(phi=None,rho=None,dx=None,dz=None,
                  iondensity=None,electrontemperature=None,plasmapotential=None):
  if phi is None: phi = frz.basegrid.phi
  if rho is None: rho = frz.basegrid.rho
  if dx is None: dx = frz.basegrid.dr
  if dz is None: dz = frz.basegrid.dz
  if iondensity is None: iondensity = w3d.iondensity
  if electrontemperature is None: electrontemperature = w3d.electrontemperature
  if plasmapotential is None: plasmapotential = w3d.plasmapotential

  # --- First calculate del squared phi plus rho/eps0
  eee = (
   (phi[0:-2,1:-1] + phi[2:,1:-1] - 2.*phi[1:-1,1:-1])/dx**2
  +(phi[1:-1,0:-2] + phi[1:-1,2:] - 2.*phi[1:-1,1:-1])/dz**2
  +rho/eps0)

  # --- If Boltzmann electron terms are given, then include the electrons
  # --- charge density.
  if iondensity != 0. and electrontemperature != 0.:
    eee = eee - (iondensity*exp(
                 (ppp[1:-1,1:-1]-plasmapotential)/electrontemperature)
                 /eps0)


  # --- Rescale the error by the coefficient of phi[1:-1,1:-1].
  eee = eee/(2./dx**2+2./dz**2)

  # --- Clear out the boundaries which are Dirichlet.
  if w3d.boundxy == top.dirichlet:
    eee[0,:] = 0.
    eee[-1,:] = 0.
  if w3d.bound0 == top.dirichlet:
    eee[:,0] = 0.
  if w3d.boundnz == top.dirichlet:
    eee[:,-1] = 0.

  # --- Clear out the residual inside of conductors
  bnd = frz.basegrid.bndfirst
  if bnd.nb_conductors > 0:
    cnd = bnd.cndfirst
    while cnd is not None:
      for j,k in zip(cnd.jcond-1,cnd.kcond-1):
        eee[j,k] = 0.
      cnd = cnd.getpyobject('next')

  return eee

def getresidualRZ(phi=None,rho=None,dx=None,dz=None,xminodx=0.,solver=None):
  if solver is not None:
    phi = solver.phi[:,0,:]
    rho = solver.rho[:,0,:]
    dx = solver.dx
    dz = solver.dz
  if phi is None: phi = frz.basegrid.phi
  if rho is None: rho = frz.basegrid.rho
  if dx is None: dx = frz.basegrid.dr
  if dz is None: dz = frz.basegrid.dz

  dxsqi = 1./dx**2
  dzsqi = 1./dz**2
  const = 2.*(dxsqi+dzsqi)
  const0 = 2.*(2.*dxsqi + dzsqi)

  nx = phi.shape[0] - 3
  nz = phi.shape[-1] - 3
  rr,zz = getmesh2d(xminodx,1.,nx,0.,1.,nz)

  res = zeros(rho.shape,'d')

  ix1 = 0
  if xminodx == 0:
    ix1 = 1
    res[0,:] = (rho[0,:]/eps0
                  + (4*phi[2,1:-1])*dxsqi
                  + (phi[1,:-2] + phi[1,2:])*dzsqi
                  - phi[1,1:-1]*const0)

  res[ix1:,:] = (rho[ix1:,:]/eps0
        + ((rr[ix1:,:]-0.5)*phi[ix1:-2,1:-1] +
           (rr[ix1:,:]+0.5)*phi[ix1+2:,1:-1])*dxsqi/rr[ix1:,:]
        + (phi[ix1+1:-1,:-2] + phi[ix1+1:-1,2:])*dzsqi
        - phi[ix1+1:-1,1:-1]*const)

  # --- Clear out the boundaries which are Dirichlet.
  if w3d.boundxy == top.dirichlet:
    res[0,:] = 0.
    res[-1,:] = 0.
  if w3d.bound0 == top.dirichlet:
    res[:,0] = 0.
  if w3d.boundnz == top.dirichlet:
    res[:,-1] = 0.

  return res

