__all__ = ['fieldsolvertesttools']
from warp import *

def gettestrho(solver,phi):

    # --- Setup the core slices
    cx = slice(1,-1)
    cy = slice(1,-1)
    cz = slice(1,-1)
    if solver.nx == 0: cx = 0
    if solver.ny == 0: cy = 0
    if solver.nz == 0: cz = 0

    phi = phi.copy()

    localbounds = solver.bounds.copy()
    fsdecomp = solver.fsdecomp
    ixproc = fsdecomp.ixproc
    iyproc = fsdecomp.iyproc
    izproc = fsdecomp.izproc
    if fsdecomp.ix[ixproc] > 0:                        localbounds[0] = -1
    if fsdecomp.ix[ixproc]+solver.nxlocal < solver.nx: localbounds[1] = -1
    if fsdecomp.iy[iyproc] > 0:                        localbounds[2] = -1
    if fsdecomp.iy[iyproc]+solver.nylocal < solver.ny: localbounds[3] = -1
    if fsdecomp.iz[izproc] > 0:                        localbounds[4] = -1
    if fsdecomp.iz[izproc]+solver.nzlocal < solver.nz: localbounds[5] = -1

    # --- phi in the guard cells needs to be modified for Neumann boundary
    # --- conditions. (Note that phi as defined is periodic.)
    if localbounds[0] == 1: phi[0,cy,cz] = phi[2,cy,cz]
    if localbounds[1] == 1: phi[-1,cy,cz] = phi[-3,cy,cz]
    if solver.ny > 0:
        if localbounds[2] == 1: phi[cx,0,cz] = phi[cx,2,cz]
        if localbounds[3] == 1: phi[cx,-1,cz] = phi[cx,-3,cz]
    if localbounds[4] == 1: phi[cx,cy,0] = phi[cx,cy,2]
    if localbounds[5] == 1: phi[cx,cy,-1] = phi[cx,cy,-3]

    # --- For Dirichlet boundaries, set phi in the boundaries consistent
    # --- with what the field solver should do.
    if localbounds[0] == 0: phi[0,cy,cz] = 2.*phi[1,cy,cz] - phi[2,cy,cz]
    if localbounds[1] == 0: phi[-1,cy,cz] = 2.*phi[-2,cy,cz] - phi[-3,cy,cz]
    if localbounds[2] == 0: phi[cx,0,cz] = 2.*phi[cx,1,cz] - phi[cx,2,cz]
    if localbounds[3] == 0: phi[cx,-1,cz] = 2.*phi[cx,-2,cz] - phi[cx,-3,cz]
    if localbounds[4] == 0: phi[cx,cy,0] = 2.*phi[cx,cy,1] - phi[cx,cy,2]
    if localbounds[5] == 0: phi[cx,cy,-1] = 2.*phi[cx,cy,-2] - phi[cx,cy,-3]

    dx,dy,dz = solver.dx,solver.dy,solver.dz
    if solver.ny > 0:
        rho = -eps0*((phi[:-2,cy,cz] - 2.*phi[cx,cy,cz] + phi[2:,cy,cz])/dx**2 +
                     (phi[cx,:-2,cz] - 2.*phi[cx,cy,cz] + phi[cx,2:,cz])/dy**2 +
                     (phi[cx,cy,:-2] - 2.*phi[cx,cy,cz] + phi[cx,cy,2:])/dz**2)
    elif solver.ny == 0:
        rho = -eps0*((phi[:-2,cy,cz] - 2.*phi[cx,cy,cz] + phi[2:,cy,cz])/dx**2 +
                     (phi[cx,cy,:-2] - 2.*phi[cx,cy,cz] + phi[cx,cy,2:])/dz**2)


    return rho,phi

def clearcorners(phi):
    # --- Clear out the corner guard cells.
    # --- These cells are ill defined and are not used so they are
    # --- ignored in these tests.
    phi[0,0,:] = 0.
    phi[0,-1,:] = 0.
    phi[0,:,0] = 0.
    phi[0,:,-1] = 0.
    phi[-1,0,:] = 0.
    phi[-1,-1,:] = 0.
    phi[-1,:,0] = 0.
    phi[-1,:,-1] = 0.
    phi[:,0,0] = 0.
    phi[:,0,-1] = 0.
    phi[:,-1,0] = 0.
    phi[:,-1,-1] = 0.

def fieldsolvertesttools(solver,testphi):
    rho,phi = gettestrho(solver,testphi)
    solver.rho[...] = rho
    solver.phi[...] = 0.

    # --- Set any Dirichlet boundary conditions.
    if solver.bounds[0] == 0: solver.phi[1,:,:] = phi[1,:,:]
    if solver.bounds[1] == 0: solver.phi[-2,:,:] = phi[-2,:,:]
    if solver.bounds[2] == 0: solver.phi[:,1,:] = phi[:,1,:]
    if solver.bounds[3] == 0: solver.phi[:,-2,:] = phi[:,-2,:]
    if solver.bounds[4] == 0: solver.phi[:,:,1] = phi[:,:,1]
    if solver.bounds[5] == 0: solver.phi[:,:,-2] = phi[:,:,-2]

    # --- Do the field solve
    fieldsolve()
    #solver.solve()

    if 0 not in solver.bounds:
        # --- Without Dirichlet boundaries, phi can change by an
        # --- arbitrary constant. solver.phi shifts by a small amount
        # --- (curiously, roughly by the error in single precision)
        # --- so subtract out an average value. For the parallel version,
        # --- the average is a global average.
        phiave = avend(getphi(bcast=1))
        solver.phi -= phiave
        cx = slice(1,-1)
        cy = slice(1,-1)
        cz = slice(1,-1)
        if solver.nx == 0: cx = 0
        if solver.ny == 0: cy = 0
        if solver.nz == 0: cz = 0
        phiave = avend(getdecomposedarray(phi[cx,cy,cz],solver=solver,bcast=1))
        phi -= phiave

    # --- Do some clean up, removing unused values.
    clearcorners(phi)
    clearcorners(solver.phi)

    # --- Find and return the error.
    phierror = maxnd(abs(solver.phi - phi))
    phierror = globalmax(phierror)

    return phierror

