# Calculate the residual of Poisson's equation.  The residual has the units
# of Q/M**3 (same units as rho)
residual_version = "$Id: residual.py,v 1.2 2000/12/01 01:34:57 dave Exp $"

residual = zeros((w3d.nx+1,w3d.ny+1,w3d.nz+1),"d")

residual[1:w3d.nx,1:w3d.ny,1:w3d.nz] = (top.eps0*(
   (w3d.phi[0:w3d.nx-1,1:w3d.ny  ,2:w3d.nz+1] -
 2.*w3d.phi[1:w3d.nx  ,1:w3d.ny  ,2:w3d.nz+1] +
    w3d.phi[2:        ,1:w3d.ny  ,2:w3d.nz+1])/w3d.dx**2 +
   (w3d.phi[1:w3d.nx  ,0:w3d.ny-1,2:w3d.nz+1] -
 2.*w3d.phi[1:w3d.nx  ,1:w3d.ny  ,2:w3d.nz+1] +
    w3d.phi[1:w3d.nx  ,2:        ,2:w3d.nz+1])/w3d.dy**2 +
   (w3d.phi[1:w3d.nx  ,1:w3d.ny  ,1:w3d.nz  ] -
 2.*w3d.phi[1:w3d.nx  ,1:w3d.ny  ,2:w3d.nz+1] +
    w3d.phi[1:w3d.nx  ,1:w3d.ny  ,3:w3d.nz+2])/w3d.dz**2)+
    w3d.rho[1:w3d.nx  ,1:w3d.ny  ,1:w3d.nz  ])

#print "Maximum residual = %e" % max(abs(ravel(residual)))
