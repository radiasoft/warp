from warp import *
from find_mgparam import *
from extpart import *

# --- Set four-character run id, comment lines, user's name.
top.runid          = "HE01"
top.pline2         = "HCX Injector:JUN.28.01:Vmarx=1.728MV,Vgate=72kV"
top.pline1         = "Child-Langmuir. Voltages:594:247:368:277:242"
top.runmaker       = "DPG and Enrique Henestroza"
# --- Invoke setup routine (THIS IS MANDATORY)
setup()
# --- Set input parameters describing the beam, 70 to 7
v_marx = 1728.0e3
v_gate = 72.e3
top.a0       =    5.08e-02
top.b0       =    5.08e-02
top.ap0      =    .0e0
top.bp0      =    .0e0
top.vbeam    =   1.0e0
top.emit     =    .0e0
top.ibeam    =    0.585
v_max        =    v_marx + v_gate
top.aion     =    39.1e0
top.zion     =    1.e0
top.lrelativ =    false
top.vthz     =    .0e0    *500.e0
top.vthperp  =    550.0
derivqty()
# --- set up arrays for lattice description
gaplen = 3.e0*2.54e-2  # gap between quad rod and plate, 3 inches
esq_platewid = .0254e0  # 1 inch
diode_len = .179e0
diode_volt = 593.573e3
extraction_volt = v_marx
quad_volt = array([247.114e3,368.424e3,277.067e3,241.822e3])
top.zlatperi  = 20.e0
top.tunelen   = (18.8e-2 + esq_platewid)*2.e0

# --- injector quadrupoles
top.quadzs[0] = - 12.22e0*2.54e-2/2.e0 # 12.22e0 inches
top.quadze[0] = + 12.22e0*2.54e-2/2.e0
top.quadap[0] = .12e0
top.quadvx[0] =  v_marx - diode_volt
top.quadvy[0] =  top.quadvx[0] - quad_volt[0]
top.quadrl[0] = 12.22e0*2.54e-2 - gaplen
top.quadgl[0] = gaplen
top.quadgp[0] = +1.e0
top.quadpa[0] = .075e0
top.qdelpar[0] = .025e0

top.quadzs[1] = top.quadze[0] + esq_platewid
top.quadze[1] = top.quadzs[1] + 18.02e0*2.54e-2
top.quadap[1] = .12e0
top.quadvy[1] = top.quadvy[0]
top.quadvx[1] = top.quadvy[1] - quad_volt[1]
top.quadrl[1] = 18.02e0*2.54e-2 - gaplen
top.quadgl[1] = gaplen
top.quadgp[1] = -1.e0
top.quadpa[1] = .10e0

top.quadzs[2] = top.quadze[1] + esq_platewid
top.quadze[2] = top.quadzs[2] + 18.8e0*2.54e-2
top.quadap[2] = .10e0
top.quadvx[2] = top.quadvx[1]
top.quadvy[2] = top.quadvx[2] - quad_volt[2]
top.quadrl[2] = 18.8e0*2.54e-2 - gaplen
top.quadgl[2] = gaplen
top.quadgp[2] = +1.e0
top.quadpa[2] = .10e0

top.quadzs[3] = top.quadze[2] + esq_platewid
top.quadze[3] = top.quadzs[3] + 18.8e0*2.54e-2
top.quadap[3] = .10e0
top.quadvy[3] = top.quadvy[2]
top.quadvx[3] = top.quadvy[3] - quad_volt[3]
top.quadrl[3] = 18.8e0*2.54e-2 - gaplen
top.quadgl[3] = gaplen
top.quadgp[3] = -1.e0
top.quadpa[3] = .10e0
top.qdelpar[3] = -.025e0

top.quadzs[4] = top.zlatperi + top.quadzs[0]
top.quadze[4] = top.zlatperi + top.quadze[0]
top.quadvx[4] = top.quadvx[0]
top.quadvy[4] = top.quadvy[0]
top.quadap[4] = top.quadap[0]
top.quadrl[4] = top.quadrl[0]
top.quadgl[4] = top.quadgl[0]
top.quadgp[4] = top.quadgp[0]
top.quadpa[4] = top.quadpa[0]
top.qdelpar[4] = top.qdelpar[0]

top.quadrr = 9./10.*top.quadap # This is to mimmic the egg-shaped rod profile. (esqround.10.8cmRod.aut)
#top.quadde = (top.quadvx - top.quadvy)#/top.quadap**2
top.quadpw[0:5] = esq_platewid
top.quadpr = 1.e0
f3d.rodfract = .5e0

# --- Set input parameters describing the 3d simulation
w3d.nx = 50;w3d.ny = 50;w3d.nz = 562
w3d.nx = 56;w3d.ny = 56;w3d.nz = 640
top.ibpush = 0                 # not mag quad focusing or bending
top.dt = 1.0e-9
top.allspecl = false
top.ifzmmnt = 2
top.prwall = 0.1
w3d.l4symtry = true
# --- Set to finite beam
top.periinz = false
top.stickyz = true
w3d.xmmin = -0.20e0
w3d.xmmax = +0.20e0
w3d.ymmin = -0.20e0
w3d.ymmax = +0.20e0
w3d.zmmin = 0.
w3d.zmmax=2.248 # Position of the slits
zmmax = w3d.zmmax

top.zlatstrt = w3d.zmmin + diode_len + esq_platewid - top.quadzs[0]
top.zimin = w3d.zmmin
top.zimax = w3d.zmmax
env.zl = w3d.zmmin
env.zu = w3d.zmmax
env.dzenv = 0.001
# --- Specify injection of the particles
top.npmax    = 0
top.inject   = 2
top.injctspc = 1
top.npinject = 2500
top.zinject  = w3d.zmmin
top.ainject  = top.a0
top.binject  = top.b0
top.rinject  = 20.32e-2 ##########CHANGE FOR DIFFERENT EMITTING SURFACE CURVATURE###########
top.apinject = 0.e0
top.bpinject = 0.e0
top.lvinject = false  # if false, source conductor input by user
top.jmaxinj = 1.0/top.pi/top.a0**2
# --- Set up some windows
top.xwindows[:,1] = [-0.001e0,.001e0]
top.ywindows[:,1] = [-0.001e0,.001e0]
top.zwindows[:,1] = [w3d.zmmin+0.015,w3d.zmmin+.017]
top.zwindows[:,2] = w3d.zmmin + diode_len/2.+array([0.000,0.002])
top.zwindows[:,3] = w3d.zmmin + diode_len+array([0.000,0.002])
top.zwindows[:,4] = top.zlatstrt + top.quadzs[0]+array([0.000,0.002])
top.zwindows[:,5] = top.zlatstrt + top.quadzs[1]+array([0.000,0.002])
top.zwindows[:,6] = top.zlatstrt + top.quadzs[2]+array([0.000,0.002])
top.zwindows[:,7] = top.zlatstrt + top.quadzs[3]+array([0.000,0.002])
top.zwindows[:,8] = top.zlatstrt + (top.quadzs[3]+top.quadze[3])/2.+array([0.000,0.002])
top.zwindows[:,9] = w3d.zmmax + array([-.010,-0.008])
# --- Select plot intervals, etc.
top.itplps[0:4]=0
top.itplfreq[0:4]=0
top.nhist=1
top.itmomnts[0:4]=[top.nhist,1000000,top.nhist,0]
top.lhlinechg = false
top.lhvzofz = false
top.lhcurrz = true
# --- set up multigrid
top.fstype = 7
w3d.bound0 = 1 #neumann boundary at iz=0
w3d.boundnz = 0 #dirichlet boundary at iz=w3d.nz
w3d.boundxy = 1 #neumann boundary
f3d.mgparam = 1.63125
f3d.downpasses = 4
f3d.uppasses = 4
f3d.mgtol = 1.e-3
f3d.mgmaxiters = 0

# --- Make some plots
def plotdiode(gridframe=0,axis='x'):
  plotquadoutline(gridframe=gridframe,axis=axis)
  if top.it > 0:
    pierce.draw()
    extractor.draw()
    plug.draw()
def myplots():
  #
  pfzx(contours=50,plotsg=0)
  ppzx(iy=w3d.iy_axis,wy=2)
  plotdiode(0,'x')
  limits(0.,zmmax,-.15,.15)
  fma()
  #
  pfzy(contours=50,plotsg=0)
  ppzy(ix=w3d.ix_axis,wx=2)
  plotdiode(0,'y')
  limits(0.,zmmax,-.15,.15)
  fma()
  #
  ppxy(iz=int(diode_len/w3d.dz),color="density",ncolor=200)
  limits(-.03,.03,-.03,.03)
  fma()
  #
  ppxy(iz=int(diode_len/w3d.dz),contours=200,filled=1,particles=0,nx=50,ny=50)
  limits(-.03,.03,-.03,.03)
  fma()
  #
  ppxxp(iz=int(diode_len/w3d.dz),color="density",ncolor=200,slope='a')
  limits(-.03,.03,-.02,.02)
  fma()
  #
  ppyyp(iz=int(diode_len/w3d.dz),color="density",ncolor=200,slope='a')
  limits(-.03,.03,-.02,.02)
  fma()
  #
  pprrp(iz=int(diode_len/w3d.dz),color="density",ncolor=200,slope='a')
  limits(.0,.03,-.02,.01)
  fma()
  #
  ppxy(iz=w3d.nz-1,color="density",ncolor=200)
  limits(-.06,.06,-.06,.06)
  fma()
  #
  ppxy(iz=w3d.nz-1,contours=200,filled=1,particles=0,nx=50,ny=50)
  limits(-.06,.06,-.06,.06)
  fma()
  #
  ppxxp(iz=w3d.nz-1,color="density",ncolor=200,slope='a')
  limits(-.06,.06,-.02,.02)
  fma()
  #
  ppyyp(iz=w3d.nz-1,color="density",ncolor=200,slope='a')
  limits(-.06,.06,-.02,.02)
  fma()
  #
  pprrp(iz=w3d.nz-1,color="density",ncolor=200,slope='a')
  limits(.0,.06,-.06,.06)
 
top.itplalways[0:3] = [0,1000000,100]
installplalways(myplots)




###########################################################################
package("w3d"); generate()

top.vbeamfrm = 0.e0
top.quads = false  # Do this because quadde is nonzero but we are not using "external" quads.
top.nplive = 1

# --- set right plane
w3d.phi[:,:,w3d.nz+1:] = top.quadvx[3]

top.vinject = v_max

###########################################################################
###########################################################################
# --- setup source conductor using srfrvinout
print "Setting up source"
xxxz1 = 0.00645245
xxxr1 = 0.0508
xxxz2 = 0.031496
xxxr2 = 0.080645
xxxz3 = 0.0359664
xxxr3 = 0.0980821
xxxz4 = -0.0002794
xxxr4 = 0.134328
xxxzc = -0.0002794
xxxrc = 0.0980821
xxxrnose = 0.0362458
xxxangle = 50.0

rmin = [0.,
        sqrt(top.rinject[0]**2 - (top.rinject[0] - xxxz1)**2),
        top.a0 + tan(xxxangle*top.pi/180.)*(xxxz2 - xxxz1),
        xxxrc]
zmin = w3d.zmmin+array([0.,xxxz1,xxxz2,xxxz3])
rcmin = [top.rinject[0],None,-xxxrnose]

rmax = [xxxr4,xxxr4,xxxrc]
zmax = w3d.zmmin+array([0.,xxxzc,xxxz3])
rcmax = [None,xxxrnose]

pierce = ZSrfrvInOut(zmin=w3d.zmmin,zmax=w3d.zmmin+xxxz3,voltage=top.vinject,
                     rminofzdata=rmin,rmaxofzdata=rmax,
                     zmindata=zmin,zmaxdata=zmax,
                     radmindata=rcmin,radmaxdata=rcmax,
                     condid=2)

installconductors(pierce)

###########################################################################

# --- setup extraction conductor
print "Setting up extraction ring"
rmin = [1.65227e-01,1.65227e-01,12.4993e-2,0.102286,0.091678526943706151,
        0.083591399873139335,
        0.0812038,0.0812038,0.0624865,0.0624865,0.0812038,0.0812038,0.0835914,0.09164809912811743,
        0.102286,0.11]
zmin = w3d.zmmin + array([0.,0.004318,4.45516e-02,4.45516e-02,0.050546,0.050546,0.0529336,0.0553212,
                          0.0553212,0.0584962,0.0584962,0.0608838,0.0632714,0.0632714,0.0693166,0.0693166],)
rcmin = [None,4.02336e-2,None,-1.23825e-2,None,-2.38760e-3,None,None,None,None,None,-2.38760e-3,None,-1.23825e-2,None]


rmax = [1.89992e-01,1.89992e-01,0.124993,0.11]
zmax = w3d.zmmin + array([0.,0.004318,0.0693166,0.0693166])
rcmax = [None,6.49986e-2,None]

extractor = ZSrfrvInOut(zmin=w3d.zmmin,zmax=w3d.zmmin+0.0693166,voltage=extraction_volt,
                        rminofzdata=rmin,rmaxofzdata=rmax,
                        zmindata=zmin,zmaxdata=zmax,
                        radmindata=rcmin,radmaxdata=rcmax,
                        condid=2)

###########################################################################

# --- setup extraction conductor PLUG
print "Setting up extraction ring PLUG"

rmin = [0.0599948,0.054991,0.054991,0.0599948,0.0674878]
zmin = w3d.zmmin + array([0.0459994,0.0510032,0.0534924,0.0584962,0.0584962])
rcmin = [-0.50038e-2,None,-0.50038e-2,None]

rmax = [0.062484,0.0674878,0.062484,0.0674878]
zmax = w3d.zmmin + array([0.0459994,0.0510032,0.0584962,0.0584962])
rcmax = [0.50038e-2,None,None]

plug = ZSrfrvInOut(zmin=0.0459994,zmax=0.0584962,voltage=extraction_volt,
                   rminofzdata=rmin,rmaxofzdata=rmax,
                   zmindata=zmin,zmaxdata=zmax,
                   radmindata=rcmin,radmaxdata=rcmax,
                   condid=3)


installconductors(pierce+extractor+plug)

# --- set so conductor is not recalculated
f3d.gridmode = 1

# --- Recalculate the fields
f3d.mgtol = 1.e-3
f3d.mgmaxiters = 100
fieldsol(-1)

