# Calculate the residual of Poisson's equation.  The residual has the units
# of Q/M**3 (same units as rho)
residual_version = "$Id: residual.py,v 1.1 2000/10/16 18:34:19 dave Exp $"

residual = zeros((w3d.nx+1,w3d.ny+1,w3d.nz+1),"d")

residual[1:w3d.nx,1:w3d.ny,1:w3d.nz] = (top.eps0*(
   (phi[0:w3d.nx-1,1:w3d.ny  ,1:w3d.nz  ] -
 2.*phi[1:w3d.nx  ,1:w3d.ny  ,1:w3d.nz  ] +
    phi[2:        ,1:w3d.ny  ,1:w3d.nz  ])/w3d.dx**2 +
   (phi[1:w3d.nx  ,0:w3d.ny-1,1:w3d.nz  ] -
 2.*phi[1:w3d.nx  ,1:w3d.ny  ,1:w3d.nz  ] +
    phi[1:w3d.nx  ,2:        ,1:w3d.nz  ])/w3d.dy**2 +
   (phi[1:w3d.nx  ,1:w3d.ny  ,0:w3d.nz-1] -
 2.*phi[1:w3d.nx  ,1:w3d.ny  ,1:w3d.nz  ] +
    phi[1:w3d.nx  ,1:w3d.ny  ,2:        ])/w3d.dz**2)+
    rho[1:w3d.nx  ,1:w3d.ny  ,1:w3d.nz  ])

print "Maximum residual = %e" % max(abs(ravel(residual)))
