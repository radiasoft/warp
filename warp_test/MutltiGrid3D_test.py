from warp import *
from fieldsolvertesttools import fieldsolvertesttools
import unittest

# --- Set four-character run id, comment lines, user's name.
top.pline2   = "Test of 3d Multigrid field solver"
top.pline1   = "for various boundary conditions"
top.runmaker = "David P. Grote"

# --- Turn off plot output
top.lprntpara = false
top.lpsplots = false

# --- Setup the grid size. It would be good to test a variety of sizes.
w3d.nx = 24*1
w3d.ny = 22*1
w3d.nz = 18*1

w3d.xmmin = -1.
w3d.xmmax =  1.
w3d.ymmin = -1.
w3d.ymmax =  1.
w3d.zmmin = -1.
w3d.zmmax =  1.

# --- Create and initialize the field solver
solver = MultiGrid()
registersolver(solver)
solver.mgverbose = -1
solver.mgtol = 1.e-12

# --- Initialize some things, including doing the decomposition for the
# --- the parallel version.
package('w3d')
generate()

# --- Fill rho with some charge. Note that these include the guard cells.
x,y,z = getmesh3d(solver.xmminlocal-solver.dx,solver.dx,solver.nxlocal+2,
                  solver.ymminlocal-solver.dy,solver.dy,solver.nylocal+2,
                  solver.zmminlocal-solver.dz,solver.dz,solver.nzlocal+2)

# --- The testphi is periodic and zero at the boundaries. This makes for
# --- easy comparison, but should still exersize all grid points.
testphi = sin(pi*x + 0.5)*sin(2*pi*y + 2.5)*sin(pi*z + 1.111)

class TestMultiGrid3DFunctions(unittest.TestCase):
    pass

# --- Dynamically generate the test methods, one for each set of boundary conditions.

# --- This dictionary will contain the boundary conditions already tested.
# --- This reduces repeats.
boundsdict = {}

# --- Loop over all possible types of boundary conditions.
# --- Note that if one side is periodic, both sides are made periodic.
for zp in [0,1,2]:
    for zm in [0,1,2]:
        if zp == 2: zm = 2
        if zm == 2: zp = 2

        for yp in [0,1,2]:
            for ym in [0,1,2]:
                if yp == 2: ym = 2
                if ym == 2: yp = 2

                for xp in [0,1,2]:
                    for xm in [0,1,2]:
                        if xp == 2: xm = 2
                        if xm == 2: xp = 2

                        boundsstr = '%d%d%d%d%d%d'%(xm,xp,ym,yp,zm,zp)
                        if boundsstr in boundsdict: continue
                        boundsdict[boundsstr] = 0

                        def doit(self,xm=xm,xp=xp,ym=ym,yp=yp,zm=zm,zp=zp):
                            solver.bounds[:] = [xm,xp,ym,yp,zm,zp]
                            phierror = fieldsolvertesttools(solver,testphi)
                            self.assertLess(phierror,1.e-10)
                            # --- Also check that the max number of iterations was not reached.
                            self.assertLess(solver.mgiters,solver.mgmaxiters)

                        setattr(TestMultiGrid3DFunctions,'test_MGwithbounds%s'%boundsstr,doit)

if __name__ == '__main__':
    unittest.main()

