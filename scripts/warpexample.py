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
top.quadde = (top.quadvx - top.quadvy)#/top.quadap**2
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
# --- set up psor
top.fstype = 11
f3d.bound0 = 1 #neumann boundary at iz=0
f3d.boundnz = 0 #dirichlet boundary at iz=w3d.nz
f3d.boundxy = 1 #neumann boundary
f3d.sorparam =    1.957683
f3d.mgparam = 1.63125
f3d.downpasses = 4
f3d.uppasses = 4
f3d.lcndbndy = true
f3d.lplates  = true
f3d.sortol = 2.
f3d.mgtol = 1.e-3
f3d.mgmaxiters = 0
f3d.lchebshv= false

f3d.ncondmax = 1
f3d.ncndmax = 1


# --- Make some plots
def plotdiode(gridframe=0,axis='x'):
  plotquadoutline(gridframe=gridframe,axis=axis)
  if top.it > 0:
    plotsrfrv(sourcemin,0.,xxxz3,gridframe=gridframe)
    plotsrfrv(sourcemax,0.,xxxz3,gridframe=gridframe)
    plotsrfrv(sourcemin,0.,xxxz3,rscale=-1.,gridframe=gridframe)
    plotsrfrv(sourcemax,0.,xxxz3,rscale=-1.,gridframe=gridframe)
    plotsrfrv(extractmin,0.,0.0693166,gridframe=gridframe)
    plotsrfrv(extractmax,0.,0.0693166,gridframe=gridframe)
    plotsrfrv(extractmin,0.,0.0693166,rscale=-1.,gridframe=gridframe)
    plotsrfrv(extractmax,0.,0.0693166,rscale=-1.,gridframe=gridframe)
    plotsrfrv(plugmin,0.0459994,0.0584962,gridframe=gridframe)
    plotsrfrv(plugmax,0.0459994,0.0584962,gridframe=gridframe)
    plotsrfrv(plugmin,0.0459994,0.0584962,rscale=-1.,gridframe=gridframe)
    plotsrfrv(plugmax,0.0459994,0.0584962,rscale=-1.,gridframe=gridframe)
def myplots():
  #
  pfzx(contours=50,plotsg=0)
  pfzxi(contours=50,plotsg=0)
  ppzx(iy=w3d.iy_axis,wy=2)
  plotdiode(0,'x')
  limits(0.,zmmax,-.15,.15)
  fma()
  #
  pfzy(contours=50,plotsg=0)
  pfzyi(contours=50,plotsg=0)
  ppzy(ix=w3d.ix_axis,wx=2)
  plotdiode(0,'y')
  limits(0.,zmmax,-.15,.15)
  fma()
  #
  ppxy(iz=int(diode_len/w3d.dz),color="density",ncolor=240)
  limits(-.03,.03,-.03,.03)
  fma()
  #
  ppxy(iz=int(diode_len/w3d.dz),contours=240,filled=1,particles=0,nx=50,ny=50)
  limits(-.03,.03,-.03,.03)
  fma()
  #
  ppxxp(iz=int(diode_len/w3d.dz),color="density",ncolor=240,slope='a')
  limits(-.03,.03,-.02,.02)
  fma()
  #
  ppyyp(iz=int(diode_len/w3d.dz),color="density",ncolor=240,slope='a')
  limits(-.03,.03,-.02,.02)
  fma()
  #
  pprrp(iz=int(diode_len/w3d.dz),color="density",ncolor=240,slope='a')
  limits(.0,.03,-.02,.01)
  fma()
  #
  ppxy(iz=w3d.nz-1,color="density",ncolor=240)
  limits(-.06,.06,-.06,.06)
  fma()
  #
  ppxy(iz=w3d.nz-1,contours=240,filled=1,particles=0,nx=50,ny=50)
  limits(-.06,.06,-.06,.06)
  fma()
  #
  ppxxp(iz=w3d.nz-1,color="density",ncolor=240,slope='a')
  limits(-.06,.06,-.02,.02)
  fma()
  #
  ppyyp(iz=w3d.nz-1,color="density",ncolor=240,slope='a')
  limits(-.06,.06,-.02,.02)
  fma()
  #
  pprrp(iz=w3d.nz-1,color="density",ncolor=240,slope='a')
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
# --- All conductors generated from the lattice are numbered '1'.
f3d.condnumb[:f3d.ncond] = 1
f3d.ecnumb[:f3d.necndbdy] = 1
f3d.ocnumb[:f3d.nocndbdy] = 1
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

def sourcemin():
  zz = f3d.srfrv_z - w3d.zmmin
  if (zz <= 0.e0):
    f3d.srfrv_r = 0.e0
  elif (zz < xxxz1):
    f3d.srfrv_r = sqrt(max(0,top.rinject[0]**2 - (top.rinject[0] - zz)**2+1.e-6))
  elif (zz < xxxz2):
    f3d.srfrv_r = top.a0 + tan(xxxangle*top.pi/180.)*(zz - xxxz1)
  elif (zz < xxxz3):
    f3d.srfrv_r = xxxrc - sqrt(max(0,xxxrnose**2 - (zz-xxxzc)**2))
  else:
    f3d.srfrv_r = xxxrc

def sourcemax():
  zz = f3d.srfrv_z - w3d.zmmin
  if (zz < xxxzc):
    f3d.srfrv_r = xxxr4
  elif (zz < xxxz3):
    f3d.srfrv_r = xxxrc + sqrt(max(0,xxxrnose**2 - (zz-xxxzc)**2))
  else:
    f3d.srfrv_r = xxxrc

is1_start = f3d.ncond
io1_start = f3d.nocndbdy
ie1_start = f3d.necndbdy
# --- make sure there is a point at x=y=0 in the source
f3d.ncond = f3d.ncond + 1
f3d.ixcond[f3d.ncond] = 0
f3d.iycond[f3d.ncond] = 0
f3d.izcond[f3d.ncond] = 0
f3d.condvolt[f3d.ncond] = top.vinject[0]
srfrvinout("sourcemin","sourcemax",top.vinject,w3d.zmmin,w3d.zmmin+xxxz3,0.e0,0.e0,false,
         w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,true,
         w3d.zmmin,w3d.zmmax,top.zbeam,w3d.dx,w3d.dy,w3d.dz,w3d.nx,w3d.ny,w3d.nz,w3d.ix_axis,w3d.iy_axis,
         w3d.xmesh,w3d.ymesh,w3d.l2symtry,w3d.l4symtry)
is1_end = f3d.ncond
io1_end = f3d.nocndbdy
ie1_end = f3d.necndbdy
# --- Source and Pierce are numbered '2'
f3d.condnumb[is1_start:is1_end] = 2
f3d.ecnumb[ie1_start:ie1_end] = 2
f3d.ocnumb[io1_start:io1_end] = 2
###########################################################################

# --- setup extraction conductor
print "Setting up extraction ring"

def extractmin():
  zz = f3d.srfrv_z - w3d.zmmin
  if (zz < 0.004318):
    f3d.srfrv_r = 1.65227e-01
  elif (zz < 4.45516e-02):
    f3d.srfrv_r = 12.4993e-2 + sqrt(max(0,(4.02336e-2)**2 - (zz-0.004318)**2))
  elif (zz < 0.050546):
    f3d.srfrv_r = 10.2286e-2 - sqrt(max(0,1.23825e-2**2 - (0.0569341 - zz)**2))
  elif (zz < 0.0529336):
    f3d.srfrv_r = 0.0835914 - sqrt(max(0,2.38760e-3**2 - (0.0529336 - zz)**2))
  elif (zz < 0.0553212):
    f3d.srfrv_r = 0.0812038
  elif (zz < 0.0584962):
    f3d.srfrv_r = 0.0624865
  elif (zz < 0.0608838):
    f3d.srfrv_r = 0.0812038
  elif (zz < 0.0632714):
    f3d.srfrv_r = 0.0835914 - sqrt(max(0,2.38760e-3**2 - (0.0608838 - zz)**2))
  elif (zz < 0.0693166):
    f3d.srfrv_r = 0.102286 - sqrt(max(0,1.23825e-2**2 - (0.0569341 - zz)**2))
  else:
    f3d.srfrv_r = 0.11

def extractmax():
  zz = f3d.srfrv_z - w3d.zmmin
  if (zz < 0.004318):
    f3d.srfrv_r = 1.89992e-01
  elif (zz < 0.0693166):
    f3d.srfrv_r = 12.4993e-2 + sqrt(max(0,(6.49986e-2)**2 - (zz-0.004318)**2))
  else:
    f3d.srfrv_r = 0.11



is2_start = f3d.ncond
io2_start = f3d.nocndbdy
ie2_start = f3d.necndbdy
srfrvinout("extractmin","extractmax",extraction_volt,w3d.zmmin,w3d.zmmin+0.0693166,0.e0,0.e0,false,
         w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,true,
         w3d.zmmin,w3d.zmmax,top.zbeam,w3d.dx,w3d.dy,w3d.dz,w3d.nx,w3d.ny,w3d.nz,w3d.ix_axis,w3d.iy_axis,
         w3d.xmesh,w3d.ymesh,w3d.l2symtry,w3d.l4symtry)

###########################################################################

# --- setup extraction conductor PLUG
print "Setting up extraction ring PLUG"

def plugmin():
  zz = f3d.srfrv_z - w3d.zmmin
  if (zz == 0.0459994):
    f3d.srfrv_r = 0.0599948
  elif (zz < 0.0510032):
    f3d.srfrv_r = 0.0599948 - sqrt(max(0,0.50038e-2**2 - (0.0510032 - zz)**2+1.e-6))
  elif (zz < 0.0534924):
    f3d.srfrv_r = 0.054991
  elif (zz < 0.0584962):
    f3d.srfrv_r = 0.0599948 - sqrt(max(0,0.50038e-2**2 - (0.0534924 - zz)**2+1.e-6))
  else:
    f3d.srfrv_r = 0.0674878

def plugmax():
  zz = f3d.srfrv_z - w3d.zmmin
  if (zz == 0.0459994):
    f3d.srfrv_r = 0.0599948
  elif (zz < 0.0510032):
    f3d.srfrv_r = 0.062484 + sqrt(max(0,0.50038e-2**2 - (0.0510032 - zz)**2+1.e-6))
  elif (zz < 0.0584962):
    f3d.srfrv_r = 0.0674878

srfrvinout("plugmin","plugmax",extraction_volt,0.0459994,0.0584962,0.e0,0.e0,false,
         w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,true,
         w3d.zmmin,w3d.zmmax,top.zbeam,w3d.dx,w3d.dy,w3d.dz,w3d.nx,w3d.ny,w3d.nz,w3d.ix_axis,w3d.iy_axis,
         w3d.xmesh,w3d.ymesh,w3d.l2symtry,w3d.l4symtry)

is2_end = f3d.ncond
io2_end = f3d.nocndbdy
ie2_end = f3d.necndbdy

# --- Extraction electode is numbered '3'
f3d.condnumb[is2_start:is2_end] = 3
f3d.ecnumb[ie2_start:ie2_end] = 3
f3d.ocnumb[io2_start:io2_end] = 3

# --- set so conductor is not recalculated
f3d.gridmode = 1

if top.fstype == 7:
  subgrid_sor_to_mg(w3d.nx,w3d.ny,w3d.nz,w3d.dx,w3d.dy,w3d.dz,
                    w3d.l2symtry,w3d.l4symtry)
# --- Recalculate the fields
f3d.sortol = 1.e-3
f3d.mgmaxiters = 100
fieldsol(-1)

