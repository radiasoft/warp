#!/usr/local/python/bin/python
from warp import *
import namelist
import sys
import getopt
import hibeamlattice
import string
from hibeamdefaults import *
hibeam_version = "$Id: hibeam.py,v 1.3 2003/08/07 18:44:10 dave Exp $"


# --- Get the command line options.
runname = 'test'
latticefile = ''
wirefile = ''
inputfile = ''
for arg in sys.argv[1:]:
  mr = re.search('r=(\w*)',arg)
  ml = re.search('l=(\w*)',arg)
  mw = re.search('w=(\w*)',arg)
  mi = re.search('i=(\w*)',arg)
  if mr: runname = mr.group(1)
  if ml: latticefile = ml.group(1)
  if mw: wirefile = mw.group(1)
  if mi: inputfile = mi.group(1)
  mr = re.search('-r\s(\w*)',arg)
  ml = re.search('-l\s(\w*)',arg)
  mw = re.search('-w\s(\w*)',arg)
  mi = re.search('-i\s(\w*)',arg)
  if mr: runname = mr.group(1)
  if ml: latticefile = ml.group(1)
  if mw: wirefile = mw.group(1)
  if mi: inputfile = mi.group(1)

if not inputfile: inputfile = 'in'+runname

# --- Read in the namelist file.
namelist = namelist.NameList(inputfile)
try:
  hinit = namelist.namelists['hinit']
except KeyError:
  hinit = {}
try:
  ztime = namelist.namelists['ztime']
except KeyError:
  ztime = {}

for k,v in hinit.items():
  v = re.sub('\.f\w*','0',v) # Change .false. to 0
  v = re.sub('\.t\w*','1',v) # Change .true. to 1
  exec(k+'='+v)
for k,v in ztime.items():
  v = re.sub('\.f\w*','0',v) # Change .false. to 0
  v = re.sub('\.t\w*','1',v) # Change .true. to 1
  exec(k+'='+v)

# --- Parse the lattice file.
line = hibeamlattice.hibeamlattice(latticefile)

##########################################################################
##########################################################################
# --- Now the standard WARP input deck.

# --- Set four-character run id, comment lines, user's name.
top.runid          = runname[:4]
top.pline2         = runname
top.pline1         = 'Hibeam simulation'
top.runmaker       = ' '
# --- Invoke setup routine (THIS IS MANDATORY).
#numticks = 4
setup()
# --- Set input parameters describing the beam, 72 to 17.
# --- Parameters calculated with envelope code ignoring focusing effect
# --- of dipoles.
top.a0       = a
top.b0       = b
top.ap0      = ap
top.bp0      = bp
top.x0       = xoffset
top.y0       = yoffset
top.xp0      = xpoffset
top.yp0      = ypoffset
top.ibeam    = current0
top.emit     = enorm0
top.vbeam    = 0.e0
top.ekin     = emev0*1.e6
top.aion     = amu
top.zion     = charge
top.lrelativ = false
derivqty()
top.vthz     = dvzth

# +++ Set up arrays describing lattice.
zz = 0.
iq = -1
idrft = -1

for e in line:
  if e.type == 'quad':
    iq = iq + 1
    top.quadzs[iq] = zz
    top.quadze[iq] = zz + e.length
    top.quadde[iq] = e.gradient
    top.quadap[iq] = e.aperture
    top.quadrr[iq] = e.r_elem
    top.qoffx[iq] = e.offset_x*errordist(e.error_type)
    top.qoffy[iq] = e.offset_y*errordist(e.error_type)
    zz = zz + e.length
    if iq == top.nquad:
      top.nquad = top.nquad + 100
      gchange("Lattice")
  if e.type == 'hyperb':
    iq = iq + 1
    top.quadzs[iq] = zz
    top.quadze[iq] = zz + e.length
    top.quadde[iq] = e.gradient
    top.quadap[iq] = -e.aperture
    top.qoffx[iq] = e.offset_x*errordist(e.error_type)
    top.qoffy[iq] = e.offset_y*errordist(e.error_type)
    zz = zz + e.length
    if iq == top.nquad:
      top.nquad = top.nquad + 100
      gchange("Lattice")
  elif e.type == 'drift':
    idrft = idrft + 1
    top.drftzs[idrft] = zz
    top.drftze[idrft] = zz + e.length
    top.drftap[idrft] = e.aperture
    top.drftox[idrft] = e.offset_x*errordist(e.error_type)
    top.drftoy[idrft] = e.offset_y*errordist(e.error_type)
    zz = zz + e.length
    if idrft = top.ndrft:
      top.ndrft = top.ndrft + 100
      gchange("Lattice")
  elif e.type == 'box':
    idrft = idrft + 1
    top.drftzs[idrft] = zz
    top.drftze[idrft] = zz + e.length
    top.drftap[idrft] = -e.aperture
    top.drftox[idrft] = e.offset_x*errordist(e.error_type)
    top.drftoy[idrft] = e.offset_y*errordist(e.error_type)
    zz = zz + e.length
    if idrft = top.ndrft:
      top.ndrft = top.ndrft + 100
      gchange("Lattice")
  elif e.type == 'wire':
    idrft = idrft + 1
    top.drftzs[idrft] = zz
    top.drftze[idrft] = zz + e.length
    top.drftap[idrft] = e.aperture
    top.drftox[idrft] = e.offset_x*errordist(e.error_type)
    top.drftoy[idrft] = e.offset_y*errordist(e.error_type)
    zz = zz + e.length
    if idrft = top.ndrft:
      top.ndrft = top.ndrft + 100
      gchange("Lattice")

# --- Set general lattice variables.
top.zlatperi  = 2.*zz
top.zlatstrt  = 0.
if iq >= 2:
  top.tunelen = 0.5*(top.quadzs[2]+top.quadze[2]-top.quadzs[0]-top.quadze[0])
else:
  top.tunelen = top.zlatperi
env.zl        = zstart
env.zu        = max(zz,zmax)
env.dzenv     = zstep0

# +++ Set input parameters describing the 3d simulation.
w3d.nx = nx
w3d.ny = ny
top.ibpush = 0
wxy.ds = zstep0
# --- Set grid size
w3d.xmmin = -xdsize
w3d.xmmax =  xdsize
w3d.ymmin = -ydsize
w3d.ymmax =  ydsize
top.prwall = w3d.xmmax
# --- Load Semi-Gaussian cigar beam.
top.npmax = npart
w3d.distrbtn = "semigaus"
w3d.xrandom  = "digitrev"
w3d.vtrandom = "digitrev"
w3d.vzrandom = "digitrev"
w3d.ldprfile = "polar"
# --- Select plot intervals, etc.
top.nhist = max(1,nint(dzhist/wxy.ds))
top.izplfreq[0:3]=[0,100000,dzphaseplot]
top.zzplfreq[3:3+len(zphaseplot)]=zphaseplot
top.itmomnts[0:4]=[0,1000000,top.nhist,0]
# --- Select plots
top.ipxy[:,0] = always
top.iptrace[:,0] = always
top.icrhoxy[:,0] = always
top.icphixy[:,0] = always
# --- open a window
#winon()
# --- Run the envelope solver to provide data used to initialize particles.
#envgen();envexe()
#wxygen()
if l_env:
  package("env")
  generate()
  step()
package("wxy");generate()

step(nint(zmax/zstep0))
dump(runnume+'.dump')
