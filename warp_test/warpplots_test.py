from warp import *
import unittest

# --- Initialize the graphics
setup()

def setup():
    """This dummy function replaces the setup from warp and is only needed
    when using Nose. Nose calls setup before running the test functions
    (and there is no way of preventing it)."""
    pass

# --- Set four-character run id, comment lines, user's name.
top.pline2   = "Test of plotting routines"
top.pline1   = "S-G beam"
top.runmaker = "David P. Grote"

# --- Create the beam species
beam = Species(type=Potassium,charge_state=+1,name="Beam species")

# --- Set input parameters describing the beam, 72 to 17.
beam.a0       = 15.601528729428563*mm
beam.b0       =  8.759246730187062*mm
beam.emit     = 62.46985168627256e-06
beam.ap0      = 0.e0
beam.bp0      = 0.e0
beam.ibeam    = 2.e-03
beam.vbeam    = 0.e0
beam.ekin     = 80.e3
derivqty()
beam.vthz     = .5e0*beam.vbeam*beam.emit/sqrt(beam.a0*beam.b0) # Vthz ~ Vthperp

# +++ Set up arrays describing lattice.
# --- Set temp variables.
hlp     = 36.0e-2  # half lattice period length
piperad = 2.0e-2 # pipe radius
quadlen = 11.e-2   # quadrupole length
dbdx    = 8.4817660221033488

# --- Set general lattice variables.
top.tunelen   = 2.e0*hlp
env.zl        = -hlp
env.zu        = +hlp
env.dzenv     = top.tunelen/64

# --- Set up quadrupoles (matching, laps, and end quads are all the same).
top.zlatperi  = 2.e0*hlp
addnewquad(zs=    - quadlen/2.,
           ze=    + quadlen/2.,
           db=-dbdx)
addnewquad(zs=hlp - quadlen/2.,
           ze=hlp + quadlen/2.,
           db=+dbdx)
addnewquad(zs=2.*hlp - quadlen/2.,
           ze=2.*hlp + quadlen/2.,
           db=-dbdx)

# +++ Set input parameters describing the 3d simulation.
w3d.nx = 16
w3d.ny = 16
w3d.nz = 64
steps_p_perd = 50
top.dt = (top.tunelen/steps_p_perd)/beam.vbeam
top.prwall = piperad

# --- Set to finite beam.
top.pbound0  = top.pboundnz = periodic
top.pboundxy = absorb
w3d.xmmin = -piperad
w3d.xmmax =  piperad
w3d.ymmin = -piperad
w3d.ymmax =  piperad
w3d.zmmin = -hlp
w3d.zmmax = +hlp

# --- Set pulse length.
top.zimin = -hlp
top.zimax = +hlp

# --- Load Semi-Gaussian cigar beam.
top.npmax = 20000
w3d.distrbtn = "semigaus"
w3d.xrandom = "digitrev"
w3d.vtrandom = "digitrev"
w3d.vzrandom = "digitrev"
w3d.ldprfile = "polar"

# --- set up field solver
w3d.l4symtry = true
w3d.bound0 = periodic
w3d.boundnz = periodic
w3d.boundxy = dirichlet

# --- Run the envelope solver to provide data used to initialize particles.
package("env")
generate()
step()

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.).
package("w3d")
generate()

class TestPlotRoutines(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        os.remove('warpplots_test.000.cgm')
        os.remove('warpplots_test.000.cgmlog')

    def test_ppzx(self):
        ppzx()
        fma()
    def test_ppzy(self):
        ppzy()
        fma()
    def test_ppzvx(self):
        ppzvx()
        fma()
    def test_ppzvy(self):
        ppzvy()
        fma()
    def test_ppzvz(self):
        ppzvz()
        fma()
    def test_ppxxp(self):
        ppxxp()
        fma()
    def test_ppyyp(self):
        ppyyp()
        fma()

    def test_ppzxiy(self):
        ppzx(iy=1)
        fma()
    def test_ppzyix(self):
        ppzy(ix=2)
        fma()
    def test_ppzvxiy(self):
        ppzvx(iy=3)
        fma()
    def test_ppzvyix(self):
        ppzvy(ix=4)
        fma()
    def test_ppzvzix(self):
        ppzvz(ix=5)
        fma()
    def test_ppxxpiz(self):
        ppxxp(iz=6)
        fma()
    def test_ppyypiz(self):
        ppyyp(iz=7)
        fma()

    def test_pfzx(self):
        pfzx()
        fma()
    def test_pfzy(self):
        pfzy()
        fma()
    def test_pfxy(self):
        pfxy()
        fma()

    def test_pfzxiy(self):
        pfzx(iy=1)
        fma()
    def test_pfzyix(self):
        pfzy(ix=2)
        fma()
    def test_pfxyiz(self):
        pfxy(iz=3)
        fma()

    def test_ppzxwithbends(self):
        ppzx(withbends=True)
        fma()
    def test_ppzxcolordensitycontours(self):
        ppzx(color='density',contours=10)
        fma()
    def test_ppzxcolordensitycellarray(self):
        ppzx(color='density',cellarray=True)
        fma()
    def test_ppzxcolordensityfilledcontours(self):
        ppzx(color='density',contours=10,filled=1)
        fma()
    def test_ppzxcolordensityfilledcontoursleveloverlap(self):
        ppzx(color='density',contours=10,filled=1,leveloverlap=0.1)
        fma()

if __name__ == '__main__':
    unittest.main()

