from warp import *
from fieldsolvertesttools import fieldsolvertesttools
import unittest

# --- Set four-character run id, comment lines, user's name.
top.pline2   = "Test of 3d FFT field solvers"
top.pline1   = "for various conditions"
top.runmaker = "David P. Grote"

# --- Turn off plot output
top.lprntpara = false
top.lpsplots = false

# --- Setup the grid size. It would be good to test a variety of sizes.
w3d.nx = 32
w3d.ny = 32
w3d.nz = 32

w3d.xmmin = -1.
w3d.xmmax =  1.
w3d.ymmin = -1.
w3d.ymmax =  1.
w3d.zmmin = -1.
w3d.zmmax =  1.

# --- Initialize some things, including doing the decomposition for the
# --- the parallel version.
package('w3d')
generate()

# --- Fill rho with some charge. Note that these include the guard cells.
x,y,z = getmesh3d(w3d.xmminlocal-w3d.dx,w3d.dx,w3d.nxlocal+2,
                  w3d.ymminlocal-w3d.dy,w3d.dy,w3d.nylocal+2,
                  w3d.zmminlocal-w3d.dz,w3d.dz,w3d.nzlocal+2)

# --- The testphi is periodic and zero at the boundaries. This makes for
# --- easy comparison, but should still exersize all grid points.
testphi = sin(pi*x)*sin(2*pi*y)*sin(pi*z + 1.111)

solver = W3DFieldsolver()

class TestFFT3DFunctions(unittest.TestCase):
    pass

def createdoit(fstype,l4symtry,l2symtry,boundz):
    def doit(self,fstype=fstype,l4symtry=l4symtry,l2symtry=l2symtry,boundz=boundz):
        top.fstype = fstype
        w3d.bounds = [0,0,0,0,boundz,boundz]
        w3d.bound0 = boundz
        w3d.boundnz = boundz
        if l4symtry:
            w3d.bounds[0] = neumann
            w3d.bounds[2] = neumann
        if l2symtry:
            w3d.bounds[2] = neumann
        fieldsolve(1)
        phierror = fieldsolvertesttools(solver,testphi)
        self.assertLess(phierror,1.e-10)
    return doit

for l4symtry,l2symtry in [[False,False],[True,False],[False,True]]:
    for fstype,boundz in [[0,periodic],[0,dirichlet],[2,periodic],[4,dirichlet],[5,periodic],[6,periodic]]:
        setattr(TestFFT3DFunctions,'test_FFTfstype%dl4s%s2s%sboundz%s'%(fstype,l4symtry,l2symtry,boundz),createdoit(fstype,l4symtry,l2symtry,boundz))

if __name__ == '__main__':
    unittest.main()
