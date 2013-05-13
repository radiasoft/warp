from warp import *
from nose.tools import assert_almost_equal

# --- This is needed to turn off graphical output
top.lprntpara = false
top.lpsplots = false
def setup():
    """This dummy function replaces the setup from warp. Nose calls setup
    before running the test functions (and there is no way of preventing it)."""
    pass

# --- Create the beam species
beam = Species(type=Potassium,charge_state=+1,name="Beam species")

def ydipolefield(x,y,z):
    (ex,ey,ez,bx,by,bz) = zeros((6,len(x)))
    by[:] = 1.
    return (ex,ey,ez,bx,by,bz)

def xdipolefield(x,y,z):
    (ex,ey,ez,bx,by,bz) = zeros((6,len(x)))
    bx[:] = 1.
    return (ex,ey,ez,bx,by,bz)

addnewpyelem(zs=0.0,ze=0.2,fn=ydipolefield)
addnewpyelem(zs=0.1,ze=0.3,fn=xdipolefield)

w3d.xmmin = -1.
w3d.ymmin = -1.
w3d.zmmin = -1.
w3d.xmmax = +1.
w3d.ymmax = +1.
w3d.zmmax = +1.
w3d.nx = 4
w3d.ny = 4
w3d.nz = 10

top.fstype = -1

package('w3d')
generate()

beam.addparticles(x=[0.,0.,0.,0.],
                  y=[0.,0.,0.,0.],
                  z=[0.05,0.15,0.25,0.35],
                  lallindomain=true)

beam.getappliedfields()

def test_p0bx():
    assert_almost_equal(beam.bx[0],0.)

def test_p0by():
    assert_almost_equal(beam.by[0],1.)

def test_p1bx():
    assert_almost_equal(beam.bx[1],1.)

def test_p1by():
    assert_almost_equal(beam.by[1],1.)

def test_p2bx():
    assert_almost_equal(beam.bx[2],1.)

def test_p2by():
    assert_almost_equal(beam.by[2],0.)

def test_p3bx():
    assert_almost_equal(beam.bx[3],0.)

def test_p3by():
    assert_almost_equal(beam.by[3],0.)

