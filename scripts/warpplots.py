from warp import *
import RandomArray
import re
import os
warpplots_version = "$Id: warpplots.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

##########################################################################
# This setups the plot handling for warp.
##########################################################################

##########################################################################
warpplotsdocbasic = """
Basic graphics commands
winon(): creates X graphic windows
hcp(): send current plot to hard-copy file
fma(): do a frame advance
plg(): basic plotting routine
plp(): plots markers (dots) instead of lines
pla(): plots multi-dimensional array as a series of lines
plotc(): contour plots 2-D data
plotfc(): contour plots 2-D data with colored contour levels
limits(): sets plot limits in order left, right, bottom, top
mouse commmands: left button zoom in, middle shifts, right zoom out

These commands returns particle info based on selection criteria.
selectparticles(): return list of indices of particles selected
getn(): get number of particles selected
getx(), gety(), getz(), getr(), gettheta(): get particle position
getvx(), getvy(), getvz(): get particle velocity
getux(), getuy(), getuz(): get particle momentum/mass
getxp(), getyp(), getrp(): get tranverse normalized velocities

These return or set a slice out of the rho or phi array.
getrho(), getphi(), setrho(), setphi()

The following plot various particles projections.
ppzxy(), ppzx(), ppzy(), ppzxp(), ppzvx(), ppzyp(), ppzvy(), ppzvz(), ppxy()
ppxxp(), ppyyp(), ppxpyp(), ppxvx(), ppyvy(), ppxvz(), ppyvz(), pprrp()
pprvz(), pptrace()

The following plot various particles projections using color.
ppzxco(), ppzyco(), ppzxyco(), ppzvzco()

Plots arbitrary particle projections using color
ppco()

Plots various quantities versus z in beam frame
ppcurr(), ppegap(), pplchg(), ppvzofz(), ppezax(), ppphiax(), pprhoax()

A selection of history plots. Plot commands for any quantity for which the
history is saved can be obtained by importing histplots.py.
hpenvx(), hpenvy(), hpepsx(), hpepsy(), hpepsnx(), hpepsny()
hepsg(), hpepsh(), hpepsng(), hpepsnh(), hppnum()
hpxbar(), hpybar(), hpvzbar(), hpxprms(), hpyprms(), hpvzrms()

Plots solution from envelope code.
penv()

Plots contours of charge density (rho) or electrostatic potential (phi) in
various planes.
pcrhozy(), pcrhozx(), pcrhoxy()
pcphizy(), pcphizx(), pcphixy()
"""
##########################################################################
warpplotsdocmore = """
setup(): does the work needed to start writing plots to a file automatically
plotruninfo(): plots run info at bottom of plots (called by fma and hcp)
beforeplot[] list of functions called after each frame advance
afterplot[] list of functions called just before each frame advance
nf() = fma() Emulation of Basis command
sf() = redraw() Emulation of Basis command
warpplg(): for plotting data distributed across processors
warpplp(): for plotting particles distributed across processors
plothist(): convenience routine for plotting generic history data

Define variable names for various colors
fg, bg, white, black, red, green, blue, cyan, magenta, yellow 

Define variable names for plot markers
point, plus, star, circle

setup_subsets(): Create subsets for particle plots (negative window numbers)
settitles(): set plot titles
ptitles(): draw plot titles on the current frame

plps[]: list of interpreter defined plots controlled by itplps
plfreq[]: list of interpreter defined plots controlled by itplfreq
pltfld3d(): makes fields plots which have been turned on
onedplts(): makes 1-D plots which have been turned on
psplots(): makes particle phase space plots which have been turned on
"""
##########################################################################

def warpplotsdoc():
  print warpplotsdocbasic
  print warpplotsdocmore

##########################################################################

##########################################################################
# The setup routine does the work needed to start writing plots to a file
# automatically.
top.lpsplots = true
always = top.always
seldom = top.seldom
never = top.never
cgmlogfile = None
numframes = 0
if me == 0: pldefault(marks=0) # --- Set plot defaults, no line marks
def setup(makepsfile=0,prefix=None,cgmlog=1):
  """
Does the work needed to start writing plots to a file automatically
  - makepsfile=0 allows the specification of a ps file instead of cgm
  - prefix=None optional prefix to use for plotfile name instead of runid
  - cgmlog=1 Set to 0 to inhibit cgmlog file creation
"""
  # --- cgmlogfile is needed elsewhere
  global cgmlogfile
  # --- Only PE0 (or serial processor) should run this routine.
  if me > 0: return
  # --- Get next available plot file name.
  if not prefix: prefix = arraytostr(top.runid)
  if makepsfile:
    pname = getnextfilename(prefix,'ps')
  else:
    pname = getnextfilename(prefix,'cgm')
  # --- Create window(0), but have it only dump to the file pname for now.
  # --- Note that only plots made to window(0) are dumped to the file.
  window(0,display='',hcp=pname,dump=1)
  # --- Set so all fma's dump plot to file.
  hcpon()
  if cgmlog:
    # --- Create plot log file and write headint to it.
    plogname = getnextfilename(prefix,'cgmlog')
    cgmlogfile = open(plogname,"w")
    cgmlogfile.write("CGMLOG file for "+pname+"\n\n")

# --- Convenience function to open a window with default value specilized to
# --- WARP. By default, this opens up a window on the current display. If
# --- setup has been called, this just creates a window which is attached to
# --- the already created device. Otherwise, open a window attached to a
# --- new device.
def winon(winnum=0,dpi=100):
  """
Opens up an X window
  - winum=0 is the window number
  - dpi=100 is the dots per inch (either 100 or 75)
  """
  if winnum==0:
    window(winnum,dpi=dpi,display=os.environ['DISPLAY'])
  else:
    window(winnum,dpi=dpi)

##########################################################################
# Plot run info to the current plot and plot info to the log file.
framet=''
frameb=''
framel=''
framer=''
def plotruninfo():
  "Plot run info to the current plot and plot info to the log file"
  global numframes
  ss = (arraytostr(top.pline3)+'\n'+
        arraytostr(top.pline2)+'\n'+
        arraytostr(top.pline1))
  ss = re.sub(r'x10\|S2\|','e',ss)
  plt(ss,0.12,0.32)
  runmaker = arraytostr(top.runmaker)
  codeid = arraytostr(top.codeid)
  rundate = arraytostr(top.rundate)
  runtime = arraytostr(top.runtime)
  runid = arraytostr(top.runid)
  ss = '%-28s  %-8s  %-8s  %-9s  %-8s'%(runmaker,codeid,rundate,runtime,runid)
  plt(ss,0.12,0.28)
  if current_window()==0:
    # --- Only increment and print frame number and log if the active
    # --- device is window(0).
    numframes = numframes + 1
    plt(repr(numframes),0.68,0.9,justify='RA')
    if cgmlogfile:
      cgmlogfile.write('%d Step %d %s %s %s %s\n' %
                       (numframes,top.it,framet,frameb,framel,framer))

##########################################################################
# Frame advance and redraw routines. The fma routine from gist is replaced
# with one that prints informative text at the bottom of each frame just
# before the normal gist fma is called. Also created are alternate (Basis
# like) names for fma and redraw.
gistfma = fma
gisthcp = hcp
afterplot = []
beforeplot = []
def fma():
  """
Frame advance - plots run info on the bottom of the frame, gets graphics window
ready for next plot and sends image to hard copy file if one is opened. Checks
for before and after plot commands.
  """
  plotruninfo()
  for f in afterplot: f()
  gistfma()
  for f in beforeplot: f()
  oldlimits = limits()
def hcp():
  """
Hardcopy - plots run info on the bottom of the frame and sends image to hard
copy file.
  """
  for f in afterplot: f()
  plotruninfo()
  for f in beforeplot: f()
  gisthcp()

nf = fma
sf = redraw

##########################################################################
# Create the plotting routines. It is different in the serial and parallel
# versions.
if not lparallel:
  def warpplp(y,x,color="fg",type="none",marker="\1",msize=1.0):
    "Plots particles, same as plg but with different defaults"
    if len(y) > 0:
      plg(y,x,color=color,type=type,marker=marker,msize=msize)
  def warpplg(y,x,color="fg",type="solid",marks=0,marker=None,msize=1.0,
              width=1.0):
    "Same as plg but with different defaults"
    if len(y) > 0:
      plg(y,x,color=color,type=type,marks=marks,marker=marker,msize=msize,
          width=width)
else:
  warpplp = plotpart
  warpplg = plotarray

# --- Plot particles
circle = '\4'
star = '\3'
plus = '\2'
point = '\1'
def plp(y,x=None,color="fg",marker="\1",msize=1.0):
  "Plots particles, same as plg but with different defaults"
  if x:
    plg(y,x,type="none",marker=marker,color=color,msize=msize)
  else:
    plg(y,type="none",marker=marker,color=color,msize=msize)

# --- Plot history data. Convenience function that is only needed until
# --- the 'limited' capability is implemented.
def plothist(v,iw):
  """
Plots any history versus z
   - v is the history array
   - iw is the window number to plot
  """
  plg(v[iw,0:top.jhist+1],top.hzbeam[0:top.jhist+1])

# --- Simple interface to contour plotting. Only requires the 2-D array
# --- to be plotted.
def plotc(zz,xx=None,yy=None,ireg=None,color=None,levs=None,contours=8,
          filled=0):
  """
Simple interface to contour plotting, same arguments as plc
  - zz 2-D array to be plotted
  - xx, yy Optional axis. Can either be 1-D or 2-D.
  - ireg Optional region. Must be same shape as zz
  - color='fg'
  - contours=8 Optional number of levels or list of levels
  - filled=0 When 1, draws filled contours
  """
  s = shape(zz)
  if len(s) != 2:
    print 'First argument must be a 2-Dimensional array'
    return
  if not xx:
    xx = arange(s[0])[:,NewAxis]*ones(s[1],'d')
  elif len(shape(xx))==1:
    xx = xx[:,NewAxis]*ones(s[1],'d')
  if not yy:
    yy = arange(s[1])*ones(s[0],'d')[:,NewAxis]
  elif len(shape(yy))==1:
    yy = yy*ones(s[0],'d')[:,NewAxis]
  if not ireg:
    ireg = ones(s)
  if filled:
    plfc(zz,xx,yy,ireg,contours=contours)
  else:
    if not color: color = 'fg'
    if contours and not levs: levs=contours
    if levs:
      if type(levs) == type(1):
        mx = maxnd(zz)
        mn = minnd(zz)
        levs = iota(0,levs)*(mx-mn)/levs + mn
      plc(zz,xx,yy,color=color,levs=levs)
    else:
      plc(zz,xx,yy,color=color)
def plotfc(zz,xx=None,yy=None,ireg=None,contours=8):
  """
Simple interface to filled contour plotting, same arguments as plfc
  - zz 2-D array to be plotted
  - xx, yy Optional axis. Can either be 1-D or 2-D.
  - ireg Optional region. Must be same shape as zz
  - color='fg'
  - contours Optional number of levels or list of levels
  """
  plotc(zz,xx=xx,yy=yy,ireg=ireg,color=color,contours=contours,filled=1)

# --- Define variables names for the allowed colors
fg = 'fg'
bg = 'bg'
white = 'white'
black = 'black'
red = 'red'
green = 'green'
blue = 'blue'
cyan = 'cyan'
magenta = 'magenta'
yellow = 'yellow'

########################################################################
########################################################################
########################################################################
# The next part of this file contains Python routines which emulate compiled
# plot routines.
#
# Here are the plots available so far
#
# ppzx(), ppzy(), ppzxp(), ppzvx(), ppzyp(), ppzvy(), ppzvz()
# ppxy(), ppxxp(), ppyyp(), ppxpyp(), ppxvx(), ppyvy()
# ppxvz(), ppyvz(), pprvz(), ppzxy(), pptrace()
# ppco(y,x,z;uz,xmin,xmax,ymin,ymax,zmin,zmax)
# ppcurr, ppegap, pplchg, ppvzofz, ppezax, ppphiax, pprhoax
#
# The following only work properly serially
#
# ppzxco(), ppzyco(), ppzvzco(), ppzxyco()
#
##########################################################################
# List of available named colors.
color = ["red","green","blue","cyan","magenta","yellow"]
##########################################################################
# Setup the random subsets. This is only called when the first plot
# is made so that top.npsplt is known to be set.
# These routines can be called at any time by the user to add more subsets,
# for example subsets with more or fewer particles, or subsets based on
# the number of particles in a species other than 0.
psubset=[]
def setup_subsets(js=0):
  """
Adds plotting subset to the list
  - js=0 is the species to create a subset for
  """
  global psubset
  for i in xrange(0,len(top.npplot)):
    ntopick=min(top.nps[js],
                int(top.npplot[i]*float(top.nps[js])/float(top.np_s[js])+0.5))
    ii = arrayrange(top.nps[0])
    rr = top.nps[0]*RandomArray.random(top.nps[0])
    ii = compress(less(rr,ntopick),ii)
    psubset.append(ii.astype('i'))

# --- Old method which replicates the numbers selected by the fortran
# --- routine psubsets. Note that that method breaks down when
# --- nps/inclump[i] > npplot[i], when nps/inclump[i] particles will be plotted.
# --- This version is maintained in case a user wants a close comparison with
# --- the Basis version.
def setup_subsetsold(js=0):
  """Old subset calculator, do not use"""
  global psubset
  # --- Print warning if npsplt is zero, in which case the subsets won't work
  if (top.npsplt == 0):
    remark("WARNING: npsplt is zero, subsets not calculated")
    return
  d2 = len(top.npplot)
  # --- Create temp ii array and copy top.isubset into it. Then replicate
  # --- the data throughout ii.
  for i in xrange(0,d2):
    n = sum(top.isubset[:,i])
    nsets = int(min(top.np_s[js],top.npmax)/top.npsplt+1)
    ii = zeros(n*nsets,'i') + top.npmax
    ii[0:n] = nonzero(top.isubset[:,i])
    for j in xrange(1,nsets):
      ii[j*n:j*n+n] = ii[0:n] + j*top.npsplt
    ii = compress(less(ii,top.npmax),ii)
    psubset.append(ii)

########################################################################
# Note: Subtracted off 0.0337 from X position of titlel (10/21/99)
ptitle_placement = [
  [[0.3950, 0.8863], [0.3950, 0.3927], [0.1200, 0.6800], [0.6253, 0.6450]],
  [[0.3950, 0.8863], [0.3950, 0.3927], [0.1200, 0.6800], [0.6253, 0.6450]],
  [[0.2634, 0.8863], [0.2634, 0.6559], [0.1200, 0.8006], [0.3600, 0.7766]],
  [[0.5266, 0.8863], [0.5266, 0.6559], [0.3832, 0.8006], [0.6253, 0.7766]],
  [[0.2634, 0.6231], [0.2634, 0.3927], [0.1200, 0.5594], [0.3600, 0.5134]],
  [[0.5266, 0.6231], [0.5266, 0.3927], [0.3832, 0.5594], [0.6253, 0.5134]],
  [[0.2634, 0.8863], [0.2634, 0.3927], [0.1200, 0.6800], [0.3600, 0.6450]],
  [[0.5266, 0.8863], [0.5266, 0.3927], [0.3832, 0.6800], [0.6253, 0.6450]],
  [[0.3950, 0.8863], [0.3950, 0.6559], [0.1200, 0.8006], [0.6253, 0.7766]],
  [[0.3950, 0.6231], [0.3950, 0.3927], [0.1200, 0.5594], [0.6253, 0.5134]]]
default_titlet=""
default_titleb=""
default_titlel=""
default_titler=""
def settitles(titlet="",titleb="",titlel="",titler=""):
  "Sets titles which are plotted by ptitles"
  global default_titlet,default_titleb,default_titlel,default_titler
  default_titlet = titlet
  default_titleb = titleb
  default_titlel = titlel
  default_titler = titler
def ptitles(titlet="",titleb="",titlel="",titler="",v=1):
  "Plots titles, either uses input or titles set by settitles"
  global framet,frameb,framel,framer
  if titlet=="" and default_titlet: titlet = default_titlet
  if titleb=="" and default_titleb: titleb = default_titleb
  if titlel=="" and default_titlel: titlel = default_titlel
  if titler=="" and default_titler: titler = default_titler
  framet=titlet
  frameb=titleb
  framel=titlel
  framer=titler
  if titlet:
    plt(titlet,ptitle_placement[v-1][0][0],ptitle_placement[v-1][0][1],
        justify="CC",orient=0)
  if titleb:
    plt(titleb,ptitle_placement[v-1][1][0],ptitle_placement[v-1][1][1],
        justify="CC",orient=0)
  if titlel:
    plt(titlel,ptitle_placement[v-1][2][0],ptitle_placement[v-1][2][1],
        justify="CC",orient=1)
  if titler:
    plt(titler,ptitle_placement[v-1][3][0],ptitle_placement[v-1][3][1],
        justify="CC",orient=1)
  settitles()

##########################################################################
##########################   UTILITY ROUTINES  ###########################
##########################################################################
##########################################################################

#-------------------------------------------------------------------------
# This returns the indices of the particles selected.
def selectparticles(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None):
  """
Selects particles based on either subsets or windows. By default it selects
from window 0, getting all of the live partilces (whose uzp > 0).
  iw = 0     Window to chose from
  js = 1     Species to chose from
  win = top.zwindows+top.zbeam    Windows to use (in lab frame)
  z = top.zp Coordinate for range selection
  iz = -1    When 0 <= iz <= nz, picks particles within zmesh[iz]+-wz*dz
  wz = 1.    Width of window around zmesh[iz]
  zl = None  When specified, lower range of selection region
  zu = None  When specified, upper range of selection region
  """
  ir1 = top.ins[js-1]-1
  ir2 = top.ins[js-1]+top.nps[js-1]-1
  if ir2 <= ir1: return array([])
  if zl or zu:
    if not z: z = top.zp
    if not zl: zl = -top.largepos
    if not zu: zu = +top.largepos
    if zl > zu: print "Warning: zl > zu"
    ii=compress(logical_and(less(zl,z[ir1:ir2]),less(z[ir1:ir2],zu)),
                arrayrange(ir1,ir2))
  elif iz:
    if not z: z = top.zp
    if lparallel:
      zl = top.zmslmin[0] + iz*w3d.dz - wz*w3d.dz + top.zbeam
      zu = top.zmslmin[0] + iz*w3d.dz + wz*w3d.dz + top.zbeam
    else:
      zl = w3d.zmmin + iz*w3d.dz - wz*w3d.dz + top.zbeam
      zu = w3d.zmmin + iz*w3d.dz + wz*w3d.dz + top.zbeam
    ii=compress(logical_and(less(zl,z[ir1:ir2]),less(z[ir1:ir2],zu)),
                arrayrange(ir1,ir2))
  elif iw < 0:
    if psubset==[]: setup_subsets()
    ii = ir1 + compress(less(psubset[-iw-1],top.nps[js-1]),psubset[-iw-1])
  else:
    if not win: win = top.zwindows + top.zbeam
    if not z: z = top.zp
    ii=compress(logical_and(less(win[0,iw],z[ir1:ir2]),
                            less(z[ir1:ir2],win[1,iw])),
                arrayrange(ir1,ir2))
  ii = compress(not_equal(take(top.uzp,ii),0.),ii)
  return ii

#-------------------------------------------------------------------------
# The following return a specific coordinate of the selected particles
# More documetation added after they are declared.
#-------------------------------------------------------------------------
def getn(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns number of particles in selection."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  if lparallel and gather: return mpi_sum(len(ii))
  else: return len(ii)
#-------------------------------------------------------------------------
def getx(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the X positions."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.xp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def gety(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the Y positions."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.yp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getz(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the Z positions."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.zp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getr(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the R postions."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = sqrt(take(top.xp,ii)**2 + take(top.yp,ii)**2)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def gettheta(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the theta postions."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = arctan2(take(top.yp,ii),take(top.xp,ii))
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvx(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the X velocity."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.uxp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvy(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the Y velocity."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.uyp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvz(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the Z velocity."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.uzp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getux(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the X momentum over mass."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.uxp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getuy(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the Y momentum over mass."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.uyp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getuz(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the Z momentum over mass."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getxp(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the X velocity over the Z velocity (X')."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.uxp,ii)/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getyp(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the Y velocity over the Z velocity (Y')."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  result = take(top.uyp,ii)/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getrp(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,gather=1):
  "Returns the radial velocity over the Z velocity (R')."
  ii = selectparticles(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  tt = arctan2(take(top.yp,ii),take(top.xp,ii))
  result = (take(top.uxp,ii)*cos(tt)+take(top.uyp,ii)*sin(tt))/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
# Add the selectparticles documentation to each of the routines.
if sys.version[:5] != "1.5.1":
  if lparallel:
    _gatherdoc = "  gather=1 When 0, all data is gathered to PE0"
  else:
    _gatherdoc = ""
  getn.__doc__ = getn.__doc__ + selectparticles.__doc__ + _gatherdoc
  getx.__doc__ = getx.__doc__ + selectparticles.__doc__ + _gatherdoc
  gety.__doc__ = gety.__doc__ + selectparticles.__doc__ + _gatherdoc
  getz.__doc__ = getz.__doc__ + selectparticles.__doc__ + _gatherdoc
  getr.__doc__ = getr.__doc__ + selectparticles.__doc__ + _gatherdoc
  gettheta.__doc__ = gettheta.__doc__ + selectparticles.__doc__ + _gatherdoc
  getvx.__doc__ = getvx.__doc__ + selectparticles.__doc__ + _gatherdoc
  getvy.__doc__ = getvy.__doc__ + selectparticles.__doc__ + _gatherdoc
  getvz.__doc__ = getvz.__doc__ + selectparticles.__doc__ + _gatherdoc
  getux.__doc__ = getux.__doc__ + selectparticles.__doc__ + _gatherdoc
  getuy.__doc__ = getuy.__doc__ + selectparticles.__doc__ + _gatherdoc
  getuz.__doc__ = getuz.__doc__ + selectparticles.__doc__ + _gatherdoc
  getxp.__doc__ = getxp.__doc__ + selectparticles.__doc__ + _gatherdoc
  getyp.__doc__ = getyp.__doc__ + selectparticles.__doc__ + _gatherdoc
  getrp.__doc__ = getrp.__doc__ + selectparticles.__doc__ + _gatherdoc
#-------------------------------------------------------------------------

def fort_setplot(iw,js,ir1,ir2,l,b,w,win,z):
  # --- Attempt at general version of these routines.
  # --- Unused now though.
  if (iw < 0):
    if psubset==[]: setup_subsets()
    settitles(l+" vs "+b,b,l,
              "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles")
    ii = compress(less(psubset[-iw-1],top.nps[js]),psubset[-iw-1])
  else:
    settitles(l+" vs "+b,b,l,w+"window %d = %9.4e, %9.4e" % (iw,win[0],win[1]))
             #w+" window"+repr(iw)+" ="+repr(win[0])+","+repr(win[1]))
    ii=compress(logical_and(less(win[0],z[ir1:ir2]),less(z[ir1:ir2],win[1])),
                arrayrange(ir2-ir1))
  ii = compress(not_equal(take(top.uzp[ir1:ir2],ii),0.),ii)
  return ii

##########################################################################
def fort_setplotx(iw,ix,wx,js,l,b):
  ir1 = top.ins[js]-1
  ir2 = top.ins[js]+top.nps[js]-1
  if (0 <= ix <= w3d.nx):
    xl = w3d.xmesh[ix] - wx*w3d.dx
    xu = w3d.xmesh[ix] + wx*w3d.dx
    settitles(l+" vs "+b,b,l,"ix = %d, x range (%9.4e, %9.4e)" % (iw,xl,xu))
    ii=compress(logical_and(less(xl,top.xp[ir1:ir2]),less(top.xp[ir1:ir2],xu)),
                arrayrange(ir1,ir2))
  elif (iw < 0):
    if psubset==[]: setup_subsets()
    settitles(l+" vs "+b,b,l,
              "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles")
    ii = ir1 + compress(less(psubset[-iw-1],top.nps[js]),psubset[-iw-1])
  else:
    settitles(l+" vs "+b,b,l,
    "x window%d = %9.4e, %9.4e" % (iw,top.xwindows[0,iw],top.xwindows[1,iw]))
    ii=compress(logical_and(less(top.xwindows[0,iw],top.xp[ir1:ir2]),
                            less(top.xp[ir1:ir2],top.xwindows[1,iw])),
                arrayrange(ir1,ir2))
  ii = compress(not_equal(take(top.uzp,ii),0.),ii)
  return ii
def xwintitle(iw=0,ix=None,wx=None):
  if (0 <= ix <= w3d.nx):
    xl = w3d.xmesh[ix] - wx*w3d.dx
    xu = w3d.xmesh[ix] + wx*w3d.dx
    return "ix = %d, x range (%9.4e, %9.4e)" % (iw,xl,xu)
  elif iw == 0:
    return "All of beam"
  elif iw < 0:
    if psubset==[]: setup_subsets()
    return "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles"
  else:
    return "x window%d = %9.4e, %9.4e"% \
           (iw,top.xwindows[0,iw],top.xwindows[1,iw])

##########################################################################
def fort_setploty(iw,iy,wy,js,l,b):
  ir1 = top.ins[js]-1
  ir2 = top.ins[js]+top.nps[js]-1
  if (0 <= iy <= w3d.ny):
    yl = w3d.ymesh[iy] - wy*w3d.dy
    yu = w3d.ymesh[iy] + wy*w3d.dy
    settitles(l+" vs "+b,b,l,"iy = %d, y range (%9.4e, %9.4e)" % (iw,yl,yu))
    ii=compress(logical_and(less(yl,top.yp[ir1:ir2]),less(top.yp[ir1:ir2],yu)),
                arrayrange(ir1,ir2))
  elif (iw < 0):
    if psubset==[]: setup_subsets()
    settitles(l+" vs "+b,b,l,
              "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles")
    ii = ir1 + compress(less(psubset[-iw-1],top.nps[js]),psubset[-iw-1])
  else:
    settitles(l+" vs "+b,b,l,
    "y window%d = %9.4e, %9.4e" % (iw,top.ywindows[0,iw],top.ywindows[1,iw]))
    ii=compress(logical_and(less(top.ywindows[0,iw],top.yp[ir1:ir2]),
                            less(top.yp[ir1:ir2],top.ywindows[1,iw])),
                arrayrange(ir1,ir2))
  ii = compress(not_equal(take(top.uzp,ii),0.),ii)
  return ii
def ywintitle(iw=0,iy=None,wy=None):
  if (0 <= iy <= w3d.ny):
    yl = w3d.ymesh[iy] - wy*w3d.dy
    yu = w3d.ymesh[iy] + wy*w3d.dy
    return "y = %d, y range (%9.4e, %9.4e)" % (iw,yl,yu)
  elif iw == 0:
    return "All of beam"
  elif (iw < 0):
    if psubset==[]: setup_subsets()
    return "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles"
  else:
    return "y window%d = %9.4e, %9.4e"% \
           (iw,top.ywindows[0,iw],top.ywindows[1,iw])

##########################################################################
def zwintitle(iw=0,iz=None,wz=1.,zl=None,zu=None):
  if iz:
    # --- Assumes that this is done by left most processor
    zl = w3d.zmmin + iz*w3d.dz - wz*w3d.dz
    zu = w3d.zmmin + iz*w3d.dz + wz*w3d.dz
    return "iz = %d, z range (%9.4e, %9.4e)"%(iz,zl,zu)
  elif zl or zl:
    if not zl: zl = -top.largepos
    if not zu: zu = +top.largepos
    return "z range (%9.4e, %9.4e)"%(zl,zu)
  elif iw == 0:
    return "All of beam"
  elif iw < 0:
    return "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles"
  elif iw > 0:
    return "z window%d = %9.4e, %9.4e"% \
           (iw,top.zwindows[0,iw],top.zwindows[1,iw])
  else:
    return ""

##########################################################################
def fort_setplotr(iw,js,l,b):
  ir1 = top.ins[js]-1
  ir2 = top.ins[js]+top.nps[js]-1
  if (iw < 0):
    if psubset==[]: setup_subsets()
    settitles(l+" vs "+b,b,l,
              "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles")
    ii = ir1 + compress(less(psubset[-iw-1],top.nps[js]),psubset[-iw-1])
  else:
    settitles(l+" vs "+b,b,l,
    "r window%d = %9.4e, %9.4e" % (iw,top.rwindows[0,iw],top.rwindows[1,iw]))
    rr2 = (top.xp[ir1:ir2]**2+top.yp[ir1:ir2]**2)
    ii=compress(logical_and(less(top.rwindows[0,iw]**2,rr2),
                            less(rr2,top.rwindows[1,iw]**2)),
                arrayrange(ir1,ir2))
  ii = compress(not_equal(take(top.uzp,ii),0.),ii)
  return ii

########################################################################
########################################################################
########################################################################
########################################################################
def ppzxy(iw=0,js=0,color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-X and Z-Y in single page
     - iw=0 spatial window to use
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ytitler = ywintitle(iw)
  xtitler = xwintitle(iw)
  if titles and iw <= 0:
    ptitles(titler=ytitler)
    ytitler = ""
    xtitler = ""
  plsys(9)
  iiy = selectparticles(iw=iw,js=js,win=top.ywindows,z=top.yp)
  warpplp(take(top.xp,iiy),take(top.zp,iiy),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles("Y vs Z","Z","Y",ytitler,v=9)
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.xplmin,top.xplmax)
  plsys(10)
  iix = selectparticles(iw=iw,js=js,win=top.xwindows,z=top.xp)
  warpplp(take(top.yp,iix),take(top.zp,iix),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles("Y vs Z","Z","Y",xtitler,v=10)
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.yplmin,top.yplmax)

##########################################################################
def ppzx(iw=0,iy=-1,wy=1,js=0,
         color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-X
     - iw=0 spatial window to use
     - iy=-1 optional grid location to use
     - wy=1 range in dy around iy to get particles to plot
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = fort_setploty(iw,iy,wy,js,"X","Z")
  warpplp(take(top.xp,ii),take(top.zp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.xplmin,top.xplmax)

##########################################################################
def ppzy(iw=0,ix=-1,wx=1,js=0,color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-Y
     - iw=0 spatial window to use
     - ix=-1 optional grid location to use
     - wx=1 range in dx around ix to get particles to plot
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = fort_setplotx(iw,ix,wx,js,"Y","Z")
  warpplp(take(top.yp,ii),take(top.zp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.yplmin,top.yplmax)

##########################################################################
def ppzxp(iw=0,iy=-1,wy=1,js=0,
          color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-X'
     - iw=0 spatial window to use
     - iy=-1 optional grid location to use
     - wy=1 range in dy around iy to get particles to plot
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = fort_setploty(iw,iy,wy,js,"X'","Z")
  warpplp(take(top.uxp,ii)/take(top.uzp,ii),take(top.zp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.xpplmin,top.xpplmax)

##########################################################################
def ppzvx(iw=0,iy=-1,wy=1,js=0,
          color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-Vx
     - iw=0 spatial window to use
     - iy=-1 optional grid location to use
     - wy=1 range in dy around iy to get particles to plot
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = fort_setploty(iw,iy,wy,js,"Vx","Z")
  warpplp(take(top.uxp,ii)*take(top.gaminv,ii),take(top.zp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam)

##########################################################################
def ppzyp(iw=0,ix=-1,wx=1,js=0,
          color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-Y'
     - iw=0 spatial window to use
     - ix=-1 optional grid location to use
     - wx=1 range in dx around ix to get particles to plot
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = fort_setplotx(iw,ix,wx,js,"Y'","Z")
  warpplp(take(top.uyp,ii)/take(top.uzp,ii),take(top.zp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.ypplmin,top.ypplmax)

##########################################################################
def ppzvy(iw=0,ix=-1,wx=1,js=0,
          color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-Vy
     - iw=0 spatial window to use
     - ix=-1 optional grid location to use
     - wx=1 range in dx around ix to get particles to plot
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = fort_setplotx(iw,ix,wx,js,"Vy","Z")
  warpplp(take(top.uyp,ii)*take(top.gaminv,ii),take(top.zp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam)

##########################################################################
def ppzvz(iw=0,js=0,color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-Vz
     - iw=0 spatial window to use
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = fort_setplotr(iw,js,"Vz","Z")
  if (top.vzrng != 0.):
     vzmax = (1. + top.vtilt)*top.vbeam*(1.+top.vzrng) - top.vzshift
     vzmin = (1. - top.vtilt)*top.vbeam*(1.-top.vzrng) - top.vzshift
  else:
     vzmax = top.vzmaxp + 0.1*(top.vzmaxp-top.vzminp)
     vzmin = top.vzminp - 0.1*(top.vzmaxp-top.vzminp)
  warpplp(take(top.uzp,ii)*take(top.gaminv,ii),take(top.zp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,vzmin,vzmax)

##########################################################################
def ppxy(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
         color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots X-Y
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("Y vs X","X","Y",titler)
  warpplp(take(top.yp,ii),take(top.xp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  if (lframe): limits(top.xplmin,top.xplmax,top.yplmin,top.yplmax)

##########################################################################
def ppxxp(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
          color="fg",marker="\1",msize=1.0,lframe=0,slope=0.,titles=1):
  """Plots X-X'
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - slope=0 slope to subtract (plots x'-slope*x,x)
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("X' vs X","X","X'",titler)
  xpo = 0.
  if type(slope) == type(''):
    if 0 <= iz <= w3d.nz:
      slope = (top.xxpbarz[iz]-top.xbarz[iz]*top.xpbarz[iz])/top.xrmsz[iz]**2
      xpo = top.xpbarz[iz]-slope*top.xbarz[iz]
    else:
      iiw = max(0,iw)
      slope = (top.xxpbar[iiw]-top.xbar[iiw]*top.xpbar[iiw])/top.xrms[iiw]**2
      xpo = top.xpbar[iiw]-slope*top.xbar[iiw]
  warpplp(take(top.uxp,ii)/take(top.uzp,ii)-
               slope*take(top.xp,ii)-xpo,take(top.xp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)

##########################################################################
def ppyyp(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
          color="fg",marker="\1",msize=1.0,lframe=0,slope=0.,titles=1):
  """Plots Y-Y'
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - slope=0 slope to subtract (plots y'-slope*y,y)
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("Y' vs Y","Y","Y'",titler)
  ypo = 0.
  if type(slope) == type(''):
    if 0 <= iz <= w3d.nz:
      slope = (top.yypbarz[iz]-top.ybarz[iz]*top.ypbarz[iz])/top.yrmsz[iz]**2
      ypo = top.ypbarz[iz]-slope*top.ybarz[iz]
    else:
      iiw = max(0,iw)
      slope = (top.yypbar[iiw]-top.ybar[iiw]*top.ypbar[iiw])/top.yrms[iiw]**2
      ypo = top.ypbar[iiw]-slope*top.ybar[iiw]
  warpplp(take(top.uyp,ii)/take(top.uzp,ii)-
               slope*take(top.yp,ii)-ypo,take(top.yp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.yplmin,top.yplmax,top.ypplmin,top.ypplmax)

##########################################################################
def ppxpyp(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
           color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots X'-Y'
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("Y' vs X'","X'","Y'",titler)
  warpplp(take(top.uyp,ii)/take(top.uzp,ii),take(top.uxp,ii)/take(top.uzp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)

##########################################################################
def ppxvx(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
          color="fg",marker="\1",msize=1.0,lframe=0,slope=0.,titles=1):
  """Plots X-Vx
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - slope=0 slope to subtract (plots vx-slope*x,x)
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("Vx vs X","X","Vx",titler)
  vxo = 0.
  if type(slope) == type(''):
    if 0 <= iz <= w3d.nz:
      slope = ((top.xxpbarz[iz]-top.xbarz[iz]*top.xpbarz[iz])/
               top.xrmsz[iz]**2*top.vzbarz[iz])
      vxo = top.vxbarz[iz]-slope*top.xbarz[iz]
    else:
      iiw = max(0,iw)
      slope = ((top.xxpbar[iiw]-top.xbar[iiw]*top.xpbar[iiw])/
               top.xrms[iiw]**2*top.vzbar[iiw])
      vxo = top.vxbar[iiw]-slope*top.xbar[iiw]
  warpplp(take(top.uxp*top.gaminv,ii)-
               slope*take(top.xp,ii)-vxo,take(top.xp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  #if (lframe): limits(top.xplmin,top.xplmax)

##########################################################################
def ppyvy(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
          color="fg",marker="\1",msize=1.0,lframe=0,slope=0.,titles=1):
  """Plots Y-Vy
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - slope=0 slope to subtract (plots vy-slope*y,y)
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("Vy vs Y","Y","Vy",titler)
  vyo = 0.
  if type(slope) == type(''):
    if 0 <= iz <= w3d.nz:
      slope = ((top.yypbarz[iz]-top.ybarz[iz]*top.ypbarz[iz])/
               top.yrmsz[iz]**2*top.vzbarz[iz])
      vyo = top.vybarz[iz]-slope*top.ybarz[iz]
    else:
      iiw = max(0,iw)
      slope = ((top.yypbar[iiw]-top.ybar[iiw]*top.ypbar[iiw])/
               top.yrms[iiw]**2*top.vzbar[iiw])
      vyo = top.vybar[iiw]-slope*top.ybar[iiw]
  warpplp(take(top.uyp*top.gaminv,ii)-
               slope*take(top.yp,ii)-vyo, take(top.yp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  #if (lframe): limits(top.yplmin,top.yplmax)

##########################################################################
def ppxvz(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
          color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots X-Vz
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("Vz vs X","X","Vz",titler)
  if (top.vzrng != 0.):
     vzmax = (1. + top.vtilt)*top.vbeam*(1.+top.vzrng) - top.vzshift
     vzmin = (1. - top.vtilt)*top.vbeam*(1.-top.vzrng) - top.vzshift
  else:
     vzmax = top.vzmaxp + 0.1*(top.vzmaxp-top.vzminp)
     vzmin = top.vzminp - 0.1*(top.vzmaxp-top.vzminp)
  warpplp(take(top.uzp,ii)*take(top.gaminv,ii),take(top.xp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.xplmin,top.xplmax,vzmin,vzmax)

##########################################################################
def ppyvz(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
          color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Y-Vz
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("Vz vs Y","Y","Vz",titler)
  if (top.vzrng != 0.):
     vzmax = (1. + top.vtilt)*top.vbeam*(1.+top.vzrng) - top.vzshift
     vzmin = (1. - top.vtilt)*top.vbeam*(1.-top.vzrng) - top.vzshift
  else:
     vzmax = top.vzmaxp + 0.1*(top.vzmaxp-top.vzminp)
     vzmin = top.vzminp - 0.1*(top.vzmaxp-top.vzminp)
  warpplp(take(top.uzp,ii)*tak(top.gaminv,ii),take(top.yp,ii),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.yplmin,top.yplmax,vzmin,vzmax)

##########################################################################
def pprrp(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
          color="fg",marker="\1",msize=1.0,lframe=0,slope=0.,scale=0,titles=1):
  """Plots R-R'
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - slope=0 slope to subtract (plots r'-slope*r,r)
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("R' vs R","R","R'",titler)
  xscale = 1.
  yscale = 1.
  xpscale = 1.
  ypscale = 1.
  if scale:
    iiw = max(0,iw)
    xscale = 2.*top.xrms[iiw]
    yscale = 2.*top.yrms[iiw]
    xpscale = 2.*top.vxrms[iiw]/top.vzbar[iiw]
    ypscale = 2.*top.vyrms[iiw]/top.vzbar[iiw]
  rr = sqrt((take(top.xp,ii)/xscale)**2 + (take(top.yp,ii)/yscale)**2)
  tt = arctan2(take(top.yp,ii)/ypscale,take(top.xp,ii)/xpscale)
  rp = ((take(top.uxp,ii)/xpscale*cos(tt)+
         take(top.uyp,ii)/ypscale*sin(tt))/take(top.uzp,ii))
  if type(slope) == type(''):
    slope = ave(rr*rp)/ave(rr**2)
  warpplp(rp-slope*rr,rr,color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)

##########################################################################
def pprvz(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
          color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots R-Vz
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - zl=None when specified, lower range of selection region
     - zu=None when specified, upper range of selection region
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
  settitles("Vz vs R","R","Vz",titler)
  if (top.vzrng != 0.):
     vzmax = (1. + top.vtilt)*top.vbeam*(1.+top.vzrng) - top.vzshift
     vzmin = (1. - top.vtilt)*top.vbeam*(1.-top.vzrng) - top.vzshift
  else:
     vzmax = top.vzmaxp + 0.1*(top.vzmaxp-top.vzminp)
     vzmin = top.vzminp - 0.1*(top.vzmaxp-top.vzminp)
  warpplp(take(top.uzp,ii)*take(top.gaminv,ii),
          sqrt(take(top.xp,ii)**2 + take(top.yp,ii)**2),
          color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles()
  settitles() 
  if (lframe): limits(0.,sup(top.xplmax,top.yplmax),vzmin,vzmax)

##########################################################################
def pptrace(iw=0,iz=None,wz=1,js=0,zl=None,zu=None,
            color="fg",marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots X-Y, X-X', Y'-Y, Y'-X' in single page
     - iw=0 spatial window to use
     - iz=-1 optional grid location to use
     - wz=1 range in dz around iz to get particles to plot
     - js=0 species to plot
     - color='fg' particle color
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  ii = selectparticles(iw=iw,js=js,iz=iz,wz=wz,zl=zl,zu=zu)
  x = take(top.xp,ii)
  y = take(top.yp,ii)
  u = take(top.uxp,ii)/take(top.uzp,ii)
  v = take(top.uyp,ii)/take(top.uzp,ii)
  if titles: ptitles(titler=zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu))
 
  plsys(3)
  warpplp(y,x,color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles("Y vs X","X","Y","",3)
  if (lframe): limits(top.xplmin,top.xplmax,top.yplmin,top.yplmax)
 
  plsys(4)
  warpplp(y,v,color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles("Y vs Y'","Y'","Y","",4)
  if (lframe): limits(top.ypplmin,top.ypplmax,top.yplmin,top.yplmax)
 
  plsys(5)
  warpplp(u,x,color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles("X' vs X","X","X'","",5)
  if (lframe): limits(top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)
 
  plsys(6)
  warpplp(u,v,color=color,type="none",marker=marker,msize=msize)
  if titles: ptitles("X' vs Y'","Y'","X'","",6)
  if (lframe): limits(top.ypplmin,top.ypplmax,top.xpplmin,top.xpplmax)

##########################################################################
def ppzxco(js=0,marker="\1",msize=1.0,lframe=0,sys=1,titles=1):
  """Plots Z-X with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - sys=1 specifies section of plot frame to use
     - titles=1 specifies whether or not to plot titles"""
  inp=top.nps[js]/top.ncolor
  istep=top.nskipcol*top.nstepcol
  #if (lframadv) nf
  istart = top.ins[js]-1
  if (inp < istep): istep = 1
  for ij in xrange(1,istep+1,top.nskipcol*2):
    for ic in xrange(1,top.ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.xp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
            color=color[ic%len(color)],type="none",marker=marker,msize=msize)
    for ic in xrange(top.ncolor,0,-1):
      irs1 = istart+ij+top.nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-top.nskipcol-inp*(ic-1))/istep
      warpplp(take(top.xp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],type="none",marker=marker,msize=msize)
  if titles: ptitles(" ","Z","X"," ",sys)
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.xplmin,top.xplmax)

##########################################################################
def ppzyco(js=0,marker="\1",msize=1.0,lframe=0,sys=1,titles=1):
  """Plots Z-Y with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - sys=1 specifies section of plot frame to use
     - titles=1 specifies whether or not to plot titles"""
  inp=top.nps[js]/top.ncolor
  istep=top.nskipcol*top.nstepcol
  #if (lframadv) nf
  istart = top.ins[js]-1
  if (inp < istep): istep = 1
  for ij in xrange(1,istep+1,top.nskipcol*2):
    for ic in xrange(1,top.ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.yp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
            color=color[ic%len(color)],type="none",marker=marker,msize=msize)
    for ic in xrange(top.ncolor,0,-1):
      irs1 = istart+ij+top.nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-top.nskipcol-inp*(ic-1))/istep
      warpplp(take(top.yp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
           color=color[ic%len(color)],type="none",marker=marker,msize=msize)
  if titles: ptitles(" ","Z","Y"," ",sys)
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.yplmin,top.yplmax)

##########################################################################
def ppzxyco(js=0,marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-X and Z-Y in single frame with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  plsys(9)
  ppzxco(js,sys=9)
  plsys(10)
  ppzyco(js,sys=10)

##########################################################################
def ppzvzco(js=0,marker="\1",msize=1.0,lframe=0,titles=1):
  """Plots Z-Vz with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  inp=top.nps[js]/top.ncolor
  istep=top.nskipcol*top.nstepcol
  #if (lframadv) nf
  if (top.vzrng != 0.):
     vzmax = (1. + top.vtilt)*top.vbeam*(1.+top.vzrng) - top.vzshift
     vzmin = (1. - top.vtilt)*top.vbeam*(1.-top.vzrng) - top.vzshift
  else:
     vzmax = top.vzmaxp + 0.1*(top.vzmaxp-top.vzminp)
     vzmin = top.vzminp - 0.1*(top.vzmaxp-top.vzminp)
  istart = top.ins[js]-1
  for ij in xrange(1,istep+1,top.nskipcol*2):
    for ic in xrange(1,top.ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.uzp[irs1:irs2:irs3]*top.gaminv[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
            color=color[ic%len(color)],type="none",marker=marker,msize=msize)
    for ic in xrange(top.ncolor,0,-1):
      irs1 = istart+ij+top.nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-top.nskipcol-inp*(ic-1))/istep
      warpplp(take(top.uzp[irs1:irs2:irs3]*top.gaminv[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
            color=color[ic%len(color)],type="none",marker=marker,msize=msize)
  if titles: ptitles("Vz vs Z","Z","Vz"," ")
  if (lframe): limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,vzmin,vzmax)

##########################################################################
def ppco(y,x,z,uz=1.,xmino=None,xmaxo=None,ymino=None,ymaxo=None,
         zmino=None,zmaxo=None,
         marker="\1",msize=1.0,lframe=0):
  """Plots y versus x with color based in z
     - y is y coordinate
     - x is x coordinate
     - z is used to calculate the color
     - xmino, xmaxo, ymino, ymaxo, zmino, zmaxo optional bounds
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  rx = ravel(x)
  ry = ravel(y)
  rz = ravel(z)
  if not lparallel:
    if not xmino: xmino = min(rx)
    if not xmaxo: xmaxo = max(rx)
    if not ymino: ymino = min(ry)
    if not ymaxo: ymaxo = max(ry)
    if not zmino: zmino = min(rz)
    if not zmaxo: zmaxo = max(rz)
  else:
    if not xmino: xmino = globalmin(rx)
    if not xmaxo: xmaxo = globalmax(rx)
    if not ymino: ymino = globalmin(ry)
    if not ymaxo: ymaxo = globalmax(ry)
    if not zmino: zmino = globalmin(rz)
    if not zmaxo: zmaxo = globalmax(rz)
  dd = (zmaxo - zmino)/top.ncolor
  #if (lframadv) nf
  for ic in xrange(1,top.ncolor+1):
    ii = compress(logical_and(logical_and(not_equal(uz,0.),
           less(zmino+(ic-1)*dd,rz)),less(rz,zmino+ic*dd)), iota(0,len(rx)))
    warpplp(take(y,ii),take(x,ii),
            color=color[ic%len(color)],type="none",marker=marker,msize=msize)
  if (lframe): limits(xmino,xmaxo,ymino,ymaxo)

##########################################################################
def ppcurr(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """Plots current along z-axis
     - color='fg' particle color
     - marks=0 turns on identifying marks on the curve
     - marker=None marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  warpplg(top.curr,top.zplmesh,color=color,
          marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Beam Current","Z")
##########################################################################
def ppegap(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """Plots smeared Ez along z-axis
     - color='fg' particle color
     - marks=0 turns on identifying marks on the curve
     - marker=None marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  warpplg(top.egap,top.zplmesh,color=color,
          marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Gap Electric Field","Z")
##########################################################################
def pplchg(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """Plots linecharge along the z-axis
     - color='fg' particle color
     - marks=0 turns on identifying marks on the curve
     - marker=None marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  warpplg(top.linechg,top.zplmesh,color=color,
          marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Line Charge","Z")
##########################################################################
def ppvzofz(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """Plots Vz along the z-axis
     - color='fg' particle color
     - marks=0 turns on identifying marks on the curve
     - marker=None marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  warpplg(top.vzofz,top.zplmesh,color=color,
          marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Mean Axial Velocity","Z")
##########################################################################
def ppezax(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """Plots Self Ez along the z-axis
     - color='fg' particle color
     - marks=0 turns on identifying marks on the curve
     - marker=None marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  warpplg(top.ezax,top.zplmesh,color=color,
          marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Z Electric Field on Axis","Z")
##########################################################################
def ppphiax(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """Plots electrostatic potential along the z-axis
     - color='fg' particle color
     - marks=0 turns on identifying marks on the curve
     - marker=None marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  warpplg(top.phiax,top.zplmesh,color=color,
          marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Electrostatic Potential on Axis","Z")
  if ((top.phiplmin != 0.0)&(top.phiplmax == 0.0)):
    if (lframe): limits(top.zzmin,top.zzmax,top.phiplmin)
  elif ((top.phiplmin == 0.0)&(top.phiplmax != 0.0)):
    if (lframe): limits(top.zzmin,top.zzmax,max(top.phiax),top.phiplmax)
  elif ((top.phiplmin != 0.0)&(top.phiplmax != 0.0)):
    if (lframe): limits(top.zzmin,top.zzmax,top.phiplmin,top.phiplmax)
##########################################################################
def pprhoax(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """Plots space-charge density along the z-axis
     - color='fg' particle color
     - marks=0 turns on identifying marks on the curve
     - marker=None marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  warpplg(top.rhoax,top.zplmesh,color=color,
          marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Charge Density on Axis","Z")
##########################################################################

# To be implemented
#ppzx4
#ppzy4
#ppzxp4
#ppzyp4
#ppzvz4
#ppxy4
#ppxxp4
#ppyyp4
#ppxpyp4
#ppxxpco
#ppyypco


##########################################################################
# This routine allows plotting of multi-dimensioned arrays.
def pla(y,x=None,color="fg",type="solid",marks=0,marker=None,msize=1.0,
        width=1.):
  "Same as plg but can plot multidimensional array"
  if x:
    if shape(x)!=shape(y) and (len(shape(x))==1 and shape(x)[0]!=shape(y)[0]):
      raise TypeError,"pla: x must either be the same shape as y, or it must be 1-D with the same len as the 1st dimension of y."
  else:
    x = arange(0,y.shape[0],1,'d')
  if len(shape(x)) > 2:
    # This next 'if' is needed since non-contiguous arrays cannot be reshaped.
    if x.iscontiguous():
      xx = x
    else:
      xx = x + 0.
    # Change xx into a 2-D array, with all of the upper dims lumped into
    # the second dimension.
    xx.shape = (xx.shape[0],product(array(xx.shape[1:])))
  elif len(shape(x)) == 2:
    # The input x is usable as is.
    xx = x
  else:
    # Extend xx into a 2-D array, with a second dimension of length 1.
    xx=x[:,NewAxis]
  if len(shape(y)) > 2:
    # This next 'if' is needed since non-contiguous arrays cannot be reshaped.
    if y.iscontiguous():
      yy = y
    else:
      yy = y + 0.
    # Change yy into a 2-D array, with all of the upper dims lumped into
    # the second dimension.
    yy.shape = (yy.shape[0],product(array(yy.shape[1:])))
  elif len(shape(y)) == 2:
    # The input y is usable as is.
    yy = y
  else:
    # Extend yy into a 2-D array, with a second dimension of length 1.
    yy=y[:,NewAxis]
  n = shape(xx)[1]
  for i in xrange(yy.shape[1]):
    warpplg(yy[:,i],xx[:,i%n],color=color,type=type,
            marks=marks,marker=marker,msize=msize,width=width)

##########################################################################
##########################################################################
# Some history plotting routines
##########################################################################
##########################################################################
def hzaxis(lzshift=0,istep=1,lvsz=1,iw=0):
  """
This is a convenience function returning an array for the z axis.
  - lzshift specifies whether the z-axis is shifted by the window locaton
  - istep=1
  - lvsz=true when true, plot versus z, otherwise versus time
  """
  if lvsz:
    if lzshift:
      return top.hzbeam[:top.jhist+1:istep] +  \
                             0.5*(top.zwindows[0,iw]+top.zwindows[1,iw])
    else:
      return top.hzbeam[:top.jhist+1:istep]
  else:
      return top.thist[:top.jhist+1:istep]
##########################################################################
def hpenvx(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1):
  """
Plots history of X envelope
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time"""
  if not me==0: return
  plg(2.*top.hxrms[iw,:top.jhist+1:istep],hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("X Beam Edge","Z","X (m)",zwintitle(iw))
##########################################################################
def hpenvy(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1):
  """
Plots history of Y envelope
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time"""
  if not me==0: return
  plg(2.*top.hyrms[iw,:top.jhist+1:istep],hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Y Beam Edge","Z","Y (m)",zwintitle(iw))
##########################################################################
def hpepsx(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of X emittance
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hepsx[iw,0]
  else: s = 1.
  plg(top.hepsx[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Unnormalized X emittance","Z","(pi-m-rad)",zwintitle(iw))
##########################################################################
def hpepsy(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of Y emittance
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hepsy[iw,0]
  else: s = 1.
  plg(top.hepsy[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Unnormalized Y emittance","Z","(pi-m-rad)",zwintitle(iw))
##########################################################################
def hpepsnx(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
            titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of X normalized emittance
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hepsnx[iw,0]
  else: s = 1.
  plg(top.hepsnx[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Normalized X emittance","Z","(pi-mm-mrad)",zwintitle(iw))
##########################################################################
def hpepsny(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
            titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of Y normalized emittance
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hepsny[iw,0]
  else: s = 1.
  plg(top.hepsny[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Normalized Y emittance","Z","(pi-mm-mrad)",zwintitle(iw))
##########################################################################
def hpepsg(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of 1st generalized emittance
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hepsg[iw,0]
  else: s = 1.
  plg(top.hepsg[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles:
    ptitles("Unnormalized generalized emittance","Z","(pi-m-rad)",zwintitle(iw))
##########################################################################
def hpepsh(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of 2nd generalized emittance
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hepsh[iw,0]
  else: s = 1.
  plg(top.hepsh[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles:
    ptitles("Unnormalized generalized emittance","Z","(pi-m-rad)",zwintitle(iw))
##########################################################################
def hpepsng(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
            titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of 1st generalized normalized emittance
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hepsng[iw,0]
  else: s = 1.
  plg(top.hepsng[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles:
    ptitles("Normalized generalized emittance","Z","(pi-mm-mrad)",zwintitle(iw))
##########################################################################
def hpepsnh(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
            titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of 2nd generalized normalized emittance
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hepsnh[iw,0]
  else: s = 1.
  plg(top.hepsnh[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles:
    ptitles("Normalized generalized emittance","Z","(pi-mm-mrad)",zwintitle(iw))
##########################################################################
def hppnum(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of the number of particles
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hpnum[iw,0]
  else: s = 1.
  plg(top.hpnum[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Number of particles","Z","#",zwintitle(iw))
##########################################################################
def hpxbar(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of X bar
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hxbar[iw,0]
  else: s = 1.
  plg(top.hxbar[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("X bar","Z","X (m)",zwintitle(iw))
##########################################################################
def hpybar(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
           titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of Y bar
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hybar[iw,0]
  else: s = 1.
  plg(top.hybar[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Y bar","Z","Y (m)",zwintitle(iw))
##########################################################################
def hpvzbar(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
            titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of average Vz
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hvzbar[iw,0]
  else: s = 1.
  plg(top.hvzbar[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Average Vz","Z","Vz (m/s)",zwintitle(iw))
##########################################################################
def hpxprms(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
            titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of X' RMS
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hxprms[iw,0]/top.hvzbar[iw,0]
  else: s = 1.
  plg(top.hvxrms[iw,:top.jhist+1:istep]/top.hvzbar[iw,:top.jhist+1:istep]/s,
      hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("X' RMS","Z","X' (rad)",zwintitle(iw))
##########################################################################
def hpyprms(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
            titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of Y' RMS
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hyprms[iw,0]/top.hvzbar[iw,0]
  else: s = 1.
  plg(top.hvyrms[iw,:top.jhist+1:istep]/top.hvzbar[iw,:top.jhist+1:istep]/s,
      hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Y' RMS","Z","Y' (rad)",zwintitle(iw))
##########################################################################
def hpvzrms(iw=0,istep=1,color="fg",marks=0,marker=None,msize=1.0,lframe=0,
            titles=1,lzshift=0,lvsz=1,lnormalized=0):
  """
Plots history of Vz RMS
  - iw=0 spatial window to use
  - istep=1
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  - lzshift=0 specifies whether the z-axis is shifted by the window
              location
  - lvsz=true when true, plot versus z, otherwise versus time
  - lnormalized=false when true, divide by initial value"""
  if not me==0: return
  if lnormalized: s = top.hvzrms[iw,0]
  else: s = 1.
  plg(top.hvzrms[iw,:top.jhist+1:istep]/s,hzaxis(lzshift,istep,lvsz,iw),
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Vz RMS","Z","Vz (m/s)",zwintitle(iw))
##########################################################################
def penv(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots a and b envelope
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles"""
  if not me==0: return
  plg(env.aenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(env.benv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Envelope","Z")
##########################################################################
##########################################################################
# --- These functions returns or sets slices of phi and rho.
##########################################################################
def getrho(ix=None,iy=None,iz=None,bcast=0):
  """Returns slices of rho, the charge density array. The shape of the object
returned depends on the number of arguments specified, which can be from none
to all three.
  - ix = None
  - iy = None
  - iz = None
  - bcast=0 When 1, the result is broadcast to all of the processors
  """
  if not lparallel:
    if ix == None and iy == None and iz == None:
      return w3d.rho
    if ix != None and iy == None and iz == None:
      return w3d.rho[ix,:,:]
    if ix == None and iy != None and iz == None:
      return w3d.rho[:,iy,:]
    if ix == None and iy == None and iz != None:
      return w3d.rho[:,:,iz]
    if ix != None and iy != None and iz == None:
      return w3d.rho[ix,iy,:]
    if ix != None and iy == None and iz != None:
      return w3d.rho[ix,:,iz]
    if ix == None and iy != None and iz != None:
      return w3d.rho[:,iy,iz]
    if ix != None and iy != None and iz != None:
      return w3d.rho[ix,iy,iz]
  else:
    if me < npes-1:
      ppp = w3d.rho[:,:,:-top.grid_overlap]
    else:
      ppp = w3d.rho[:,:,:]
    if ix != None and iy == None:
      ppp = ppp[ix,:,:]
    elif ix == None and iy != None:
      ppp = ppp[:,iy,:]
    elif ix != None and iy != None:
      ppp = ppp[ix,iy,:]
    if iz == None:
      ppp = transpose(gatherarray(transpose(ppp)))
    else:
      pe = convertiztope(iz)
      if pe == None: return None
      if me == pe: ppp = ppp[...,iz-top.izslave[me+1]]
      if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
    if bcast: ppp = mpi.bcast(ppp)
    return ppp
# --------------------------------------------------------------------------
def setrho(val,ix=None,iy=None,iz=None):
  """Sets slices of rho, the charge density array. The shape of the input
object depends on the number of arguments specified, which can be from none
to all three.
  - val input array (must be supplied)
  - ix = None
  - iy = None
  - iz = None
  """
  if not lparallel:
    if ix == None and iy == None and iz == None:
      w3d.rho[:,:,:] = val
    if ix != None and iy == None and iz == None:
      w3d.rho[ix,:,:] = val
    if ix == None and iy != None and iz == None:
      w3d.rho[:,iy,:] = val
    if ix == None and iy == None and iz != None:
      w3d.rho[:,:,iz] = val
    if ix != None and iy != None and iz == None:
      w3d.rho[ix,iy,:] = val
    if ix != None and iy == None and iz != None:
      w3d.rho[ix,:,iz] = val
    if ix == None and iy != None and iz != None:
      w3d.rho[:,iy,iz] = val
    if ix != None and iy != None and iz != None:
      w3d.rho[ix,iy,iz] = val
  else:
    print "Warning, setrho this is not yet implemented in parallel"
   #if me < npes-1:
   #  ppp = w3d.rho[:,:,:-top.grid_overlap]
   #else:
   #  ppp = w3d.rho[:,:,:]
   #if ix != None and iy == None:
   #  ppp = ppp[ix,:,:]
   #elif ix == None and iy != None:
   #  ppp = ppp[:,iy,:]
   #elif ix != None and iy != None:
   #  ppp = ppp[ix,iy,:]
   #if iz == None:
   #  ppp = transpose(gatherarray(transpose(ppp)))
   #else:
   #  pe = convertiztope(iz)
   #  if pe == None: return None
   #  if me == pe: ppp = ppp[...,iz-top.izslave[me+1]+1]
   #  if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
   #if bcast: ppp = mpi.bcast(ppp)
   #return ppp
# --------------------------------------------------------------------------
def getphi(ix=None,iy=None,iz=None,bcast=0):
  """Returns slices of phi, the electrostatic potential array. The shape of
the object returned depends on the number of arguments specified, which can
be from none to all three.
  - ix = None
  - iy = None
  - iz = None Value is relative to the fortran indexing, so iz ranges
              from -1 to nz+1
  - bcast=0 When 1, the result is broadcast to all of the processors
  """
  if not lparallel:
    if ix == None and iy == None and iz == None:
      return w3d.phi[:,:,1:-1]
    if ix != None and iy == None and iz == None:
      return w3d.phi[ix,:,1:-1]
    if ix == None and iy != None and iz == None:
      return w3d.phi[:,iy,1:-1]
    if ix == None and iy == None and iz != None:
      return w3d.phi[:,:,iz+1]
    if ix != None and iy != None and iz == None:
      return w3d.phi[ix,iy,1:-1]
    if ix != None and iy == None and iz != None:
      return w3d.phi[ix,:,iz+1]
    if ix == None and iy != None and iz != None:
      return w3d.phi[:,iy,iz+1]
    if ix != None and iy != None and iz != None:
      return w3d.phi[ix,iy,iz+1]
  else:
    if me < npes-1:
      ppp = w3d.phi[:,:,1:w3d.nz-top.grid_overlap+2]
    else:
      ppp = w3d.phi[:,:,1:-1]
    if ix != None and iy == None:
      ppp = ppp[ix,:,:]
    elif ix == None and iy != None:
      ppp = ppp[:,iy,:]
    elif ix != None and iy != None:
      ppp = ppp[ix,iy,:]
    if iz == None:
      ppp = transpose(gatherarray(transpose(ppp)))
    else:
      pe = convertiztope(iz)
      if pe == None: return None
      if me == pe: ppp = ppp[...,iz-top.izslave[me+1]]
      if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
    if bcast: ppp = mpi.bcast(ppp)
    return ppp
# --------------------------------------------------------------------------
def setphi(val,ix=None,iy=None,iz=None):
  """Sets slices of phi, the electrostatic potential array. The shape of
the input object depends on the number of arguments specified, which can
be from none to all three.
  - val input array (must be supplied)
  - ix = None
  - iy = None
  - iz = None Value is relative to the fortran indexing, so iz ranges
              from -1 to nz+1
  """
  if not lparallel:
    if ix == None and iy == None and iz == None:
      w3d.phi[:,:,1:-1] = val
    if ix != None and iy == None and iz == None:
      w3d.phi[ix,:,1:-1] = val
    if ix == None and iy != None and iz == None:
      w3d.phi[:,iy,1:-1] = val
    if ix == None and iy == None and iz != None:
      w3d.phi[:,:,iz+1] = val
    if ix != None and iy != None and iz == None:
      w3d.phi[ix,iy,1:-1] = val
    if ix != None and iy == None and iz != None:
      w3d.phi[ix,:,iz+1] = val
    if ix == None and iy != None and iz != None:
      w3d.phi[:,iy,iz+1] = val
    if ix != None and iy != None and iz != None:
      w3d.phi[ix,iy,iz+1] = val
  else:
    print "Warning, setphi this is not yet implemented in parallel"
   #if me < npes-1:
   #  ppp = w3d.phi[:,:,1:w3d.nz-top.grid_overlap+2]
   #else:
   #  ppp = w3d.phi[:,:,1:-1]
   #if ix != None and iy == None:
   #  ppp = ppp[ix,:,:]
   #elif ix == None and iy != None:
   #  ppp = ppp[:,iy,:]
   #elif ix != None and iy != None:
   #  ppp = ppp[ix,iy,:]
   #if iz == None:
   #  ppp = transpose(gatherarray(transpose(ppp)))
   #else:
   #  pe = convertiztope(iz)
   #  if pe == None: return None
   #  if me == pe: ppp = ppp[...,iz-top.izslave[me+1]+1]
   #  if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
   #if bcast: ppp = mpi.bcast(ppp)
   #return ppp
##########################################################################
##########################################################################
def pcrhozy(ix=None,contours=20,titles=1,filled=1,color=None):
  """
Plots contours of charge density in the Z-Y plane
  - ix=w3d.ix_axis X index of plane
  - contours=20 number or list of contours
  - titles=1 specifies whether or not to plot titles"""
  if not ix: ix = w3d.ix_axis
  plotc(getrho(ix=ix),w3d.ymesh,w3d.zmesh,
        contours=contours,filled=filled,color=color)
  if titles: ptitles("Charge density in z-y plane","Z","Y","ix = "+repr(ix))
##########################################################################
def pcrhozx(iy=None,contours=20,titles=1,filled=1,color=None):
  """Plots contours of charge density in the Z-X plane
     - iy=w3d.iy_axis Y index of plane
     - contours=20 number or list of contours
     - titles=1 specifies whether or not to plot titles"""
  if not iy: iy = w3d.iy_axis
  zz=w3d.zmesh*ones(w3d.nz+1,'d')[:,NewAxis]
  xx=w3d.xmesh[:,NewAxis]*ones(w3d.nx+1,'d')
  ireg=ones((w3d.nz+1,w3d.nx+1))
  plotc(getrho(iy=iy),w3d.xmesh,w3d.zmesh,
        contours=contours,filled=filled,color=color)
  if titles: ptitles("Charge density in z-x plane","Z","X","iy = "+repr(iy))
##########################################################################
def pcrhoxy(iz=None,contours=20,titles=1,filled=1,color=None):
  """Plots contours of charge density in the X-Y plane
     - iz=w3d.iz_axis Z index of plane
     - contours=20 number or list of contours
     - titles=1 specifies whether or not to plot titles"""
  if not iz: iz = w3d.iz_axis
  xx=w3d.xmesh*ones(w3d.nx+1,'d')[:,NewAxis]
  yy=w3d.ymesh[:,NewAxis]*ones(w3d.ny+1,'d')
  ireg=ones((w3d.nx+1,w3d.ny+1))
  plotc(getrho(iz=iz),w3d.ymesh,w3d.xmesh,
        contours=contours,filled=filled,color=color)
  if titles: ptitles("Charge density in x-y plane","X","Y","iz = "+repr(iz))
##########################################################################
def pcphizy(ix=None,contours=20,titles=1,filled=1,color=None):
  """Plots contours of electrostatic potential in the Z-Y plane
     - ix=w3d.ix_axis X index of plane
     - contours=20 number or list of contours
     - titles=1 specifies whether or not to plot titles"""
  if not ix: ix = w3d.ix_axis
  zz=w3d.zmesh*ones(w3d.nz+1,'d')[:,NewAxis]
  yy=w3d.ymesh[:,NewAxis]*ones(w3d.ny+1,'d')
  ireg=ones((w3d.nz+1,w3d.ny+1))
  plotc(getphi(ix=ix),w3d.ymesh,w3d.zmesh,
        contours=contours,filled=filled,color=color)
  if titles: ptitles("Potential in z-y plane","Z","Y","ix = "+repr(ix))
##########################################################################
def pcphizx(iy=None,contours=20,titles=1,filled=1,color=None):
  """Plots contours of electrostatic potential in the Z-X plane
     - iy=w3d.iy_axis Y index of plane
     - contours=20 number or list of contours
     - titles=1 specifies whether or not to plot titles"""
  if not iy: iy = w3d.iy_axis
  zz=w3d.zmesh*ones(w3d.nz+1,'d')[:,NewAxis]
  xx=w3d.xmesh[:,NewAxis]*ones(w3d.nx+1,'d')
  ireg=ones((w3d.nz+1,w3d.nx+1))
  plotc(getphi(iy=iy),w3d.xmesh,w3d.zmesh,
        contours=contours,filled=filled,color=color)
  if titles: ptitles("Potential in z-x plane","Z","X","iy = "+repr(iy))
##########################################################################
def pcphixy(iz=None,contours=20,titles=1,filled=1,color=None):
  """Plots contours of electrostatic potential in the X-Y plane
     - iz=w3d.iz_axis Z index of plane
     - contours=20 number or list of contours
     - titles=1 specifies whether or not to plot titles"""
  if not iz: iz = w3d.iz_axis
  xx=w3d.xmesh*ones(w3d.nx+1,'d')[:,NewAxis]
  yy=w3d.ymesh[:,NewAxis]*ones(w3d.ny+1,'d')
  ireg=ones((w3d.nx+1,w3d.ny+1))
  plotc(getphi(iz=iz),w3d.ymesh,w3d.xmesh,
        contours=contours,filled=filled,color=color)
  if titles: ptitles("Potential in x-y plane","X","Y","iz = "+repr(iz))
##########################################################################
##########################################################################
plps = []
plfreq = []
def pltfld3d(fld='phi',freqflag=always):
  """Makes fields plots which have been turned on
     - fld='phi' quantity to plot, either 'phi' or 'rho'
     - freqflag=always frequency flag, either always, seldom, or never"""
  currentwindow = current_window()
  window(0)
  nwindows = 9
  for i in xrange(nwindows):
    if (top.icrhoxy[i] == freqflag and fld == "rho"): pcrhoxy[i]
    if (top.icrhozx[i] == freqflag and fld == "rho"): pcrhozx[i]
    if (top.icrhozy[i] == freqflag and fld == "rho"): pcrhozy[i]
    if (top.icphixy[i] == freqflag and fld == "phi"): pcphixy[i]
    if (top.icphizx[i] == freqflag and fld == "phi"): pcphizx[i]
    if (top.icphizy[i] == freqflag and fld == "phi"): pcphizy[i]
  #if (top.icrhoxy4 == freqflag and fld == "rho"): pcrhoxy4
  #if (top.icrhozx4 == freqflag and fld == "rho"): pcrhozx4
  #if (top.icrhozy4 == freqflag and fld == "rho"): pcrhozy4
  #if (top.icphixy4 == freqflag and fld == "phi"): pcphixy4
  #if (top.icphizx4 == freqflag and fld == "phi"): pcphizx4
  #if (top.icphizy4 == freqflag and fld == "phi"): pcphizy4
  oldlimits = limits()
  window(currentwindow)

def onedplts(freqflag=always):
  """Makes 1-D plots which have been turned on
     - freqflag=always frequency flag, either always, seldom, or never"""
  currentwindow = current_window()
  window(0)
  if freqflag == top.ipcurr: ppcurr()
  if freqflag == top.ipegap: ppegap()
  if freqflag == top.iplchg: pplchg()
  if freqflag == top.ipvzofz: ppvzofz()
  if freqflag == top.iprhoax: pprhoax()
  if freqflag == top.ipphiax: ppphiax()
  if freqflag == top.ipezax: ppezax()
  oldlimits = limits()
  window(currentwindow)

def psplots(freqflag=always,js=0):
  """Makes particle phase space plots which have been turned on
     - freqflag=always frequency flag, either always, seldom, or never
     - js=0 specifies the species of particles to plot"""
# --- Phase space plots, both "frequent" ones and others
# --- Do z-x,y 2-to-a-page subset and all-particle plots
  # --- Save current device and set active device to window(0). This
  # --- ensures that plots created by this routine will be dumped to
  # --- the appropriate plot file.
  currentwindow = current_window()
  window(0)

  nsubsets = 3
  nwindows = 9

  for i in xrange(-nsubsets,1):
    if (top.ipzxy[i] == freqflag):
      ppzxy(i,lframe=true)
      fma()

# --- Do z-x,y 2-to-a-page in color, skipping NSKIPCOL particles
  if (top.ipzxyco == freqflag):
    ppzxyco(js,lframe=true)
    fma()

# --- Do z-vz in color, skipping NSKIPCOL particles
  if (top.ipzvzco == freqflag):
    ppzvzco(js,lframe=true)
    fma()

# --- Do x-xp in color, skipping NSKIPCOL particles
  for i in xrange(nwindows+1):
   if (top.ipxxpco[i] == freqflag):
     ppxxpco(i,lframe=true)
     fma()

# --- Do y-yp in color, skipping NSKIPCOL particles
  for i in xrange(nwindows+1):
   if (top.ipyypco[i] == freqflag):
     ppyypco(i,lframe=true)
     fma()

# --- Do z-x and z-xp subset and y-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzx[i] == freqflag):
      ppzx(i,lframe=true)
      fma()
  #if (top.ipzx4 == freqflag):
    #ppzx4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzxp[i] == freqflag):
      ppzxp(i,lframe=true)
      fma()
  #if (top.ipzxp4 == freqflag):
    #ppzxp4

# --- Do z-y and z-yp subset and x-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzy[i] == freqflag):
      ppzy(i,lframe=true)
      fma()
  #if (top.ipzy4 == freqflag):
    #ppzy4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzyp[i] == freqflag):
      ppzyp(i,lframe=true)
      fma()
  #if (top.ipzyp4 == freqflag):
    #ppzyp4
    #fma()

# --- Do z-vz subset and r-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzvz[i] == freqflag):
      ppzvz(i,lframe=true)
      fma()
  #if (top.ipzvz4 == freqflag):
    #ppzvz4
    #fma()

# --- Do transverse phase-space subset and z-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxy[i] == freqflag):
      ppxy(i,lframe=true)
      fma()
  #if (top.ipxy4 == freqflag):
    #ppxy4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxxp[i] == freqflag):
      ppxxp(i,lframe=true)
      fma()
  #if (top.ipxxp4 == freqflag):
    #ppxxp4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipyyp[i] == freqflag):
      ppyyp(i,lframe=true)
      fma()
  #if (top.ipyyp4 == freqflag):
    #ppyyp4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxpyp[i] == freqflag):
      ppxpyp(i,lframe=true)
      fma()
  #if (top.ipxpyp4 == freqflag):
    #ppxpyp4
    #fma()

# --- Do trace-space z-window plots
  for i in xrange(nwindows+1):
    if (top.iptrace[i] == freqflag and i >= 0):
      pptrace(i,lframe=true)
      fma()

# --- Do the user defined plots
  oldlimits = limits()
  if freqflag == always:
    for p in plfreq:
      p()
      fma()
      oldlimits = limits()
  if freqflag == seldom:
    for p in plps:
      p()
      fma()
      oldlimits = limits()

# --- Reset the current window to it previous value.
  oldlimits = limits()
  window(currentwindow)
