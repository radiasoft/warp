from warp import *
import RandomArray
import re
import os
if me == 0:
  try:
    import pl3d
    import plwf
  except ImportError:
    pass
warpplots_version = "$Id: warpplots.py,v 1.62 2001/12/21 18:52:26 dave Exp $"

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
getgaminv(): get gamma inverse

These return or set a slice out of the rho or phi array.
getrho(), getphi(), setrho(), setphi()

The following plot various particles projections.
ppzxy(), ppzx(), ppzy(), ppzr(), ppzxp(), ppzvx(), ppzyp(), ppzvy(), ppzvz(), ppxy()
ppxxp(), ppyyp(), ppxpyp(), ppxvx(), ppyvy(), ppxvz(), ppyvz(), pprrp()
pprvz(), pptrace()

The following plot various particles projections using color.
ppzxco(), ppzyco(), ppzxyco(), ppzvzco()

Plots arbitrary particle projections using color
ppco()

Plots various quantities versus z in beam frame (see pzplotsdoc())

Run histplotsdoc() for a list of history plots.

Plots solution from envelope code.
penv()

Plots contours of charge density (rho) or electrostatic potential (phi) in
various planes.
pcrhozy(), pcrhozx(), pcrhoxy()
pcphizy(), pcphizx(), pcphixy()

Dynamically view any 3-D surface plot
viewsurface
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
clear_subsets(): Clears the subsets for particle plots (negative window numbers)
settitles(): set plot titles
ptitles(): draw plot titles on the current frame

plps[]: list of interpreter defined plots controlled by itplps
plfreq[]: list of interpreter defined plots controlled by itplfreq
plseldom[]: list of interpreter defined plots controlled by itplseldom
plalways[]: list of interpreter defined plots controlled by itplalways
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
def setup(makepsfile=0,prefix=None,cgmlog=1,runcomments=''):
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
  print "Plot file name",pname
  # --- Set so all fma's dump plot to file.
  hcpon()
  if cgmlog:
    # --- Create plot log file and write headint to it.
    plogname = getnextfilename(prefix,'cgmlog')
    cgmlogfile = open(plogname,"w")
    cgmlogfile.write("CGMLOG file for "+pname+"\n\n")
  # --- Print the versions to the plot file.
  plt(versionstext()+'\n'+runcomments,0.15,0.88,justify="LT")
  fma()

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
def fma(legend=1):
  """
Frame advance - plots run info on the bottom of the frame, gets graphics window
ready for next plot and sends image to hard copy file if one is opened. Checks
for before and after plot commands.
  - legend=1: when set to 0, the text at the frame bottom is omitted
  """
  if legend: plotruninfo()
  for f in afterplot: f()
  gistfma()
  for f in beforeplot: f()
  oldlimits = limits()
def hcp(legend=1):
  """
Hardcopy - plots run info on the bottom of the frame and sends image to hard
copy file.
  - legend=1: when set to 0, the text at the frame bottom is omitted
  """
  for f in afterplot: f()
  if legend: plotruninfo()
  for f in beforeplot: f()
  gisthcp()

nf = fma
sf = redraw

##########################################################################
# Create the plotting routines. It is different in the serial and parallel
# versions.
if not lparallel:
  def warpplp(y,x,color="fg",linetype="none",marker="\1",msize=1.0):
    "Plots particles, same as plg but with different defaults"
    if len(y) > 0:
      plg(y,x,color=color,type=linetype,marker=marker,msize=msize)
  def warpplg(y,x,color="fg",linetype="solid",marks=0,marker=None,msize=1.0,
              width=1.0):
    "Same as plg but with different defaults"
    if len(y) > 0:
      plg(y,x,color=color,type=linetype,marks=marks,marker=marker,msize=msize,
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
  """Plots particles, same as plg but with different defaults so it plots
markers instead of lines"""
  if len(y) == 0: return
  if x is not None:
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
          filled=0,width=1.,linetype='solid'):
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
  if type(contours) == ListType: contours = array(contours)
  if type(contours) == TupleType: contours = array(contours)
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
      plc(zz,xx,yy,color=color,levs=levs,width=width,type=linetype)
    else:
      plc(zz,xx,yy,color=color,width=width,type=linetype)
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
# ppzx(), ppzy(), ppzr(), ppzxp(), ppzvx(), ppzyp(), ppzvy(), ppzvz()
# ppxy(), ppxxp(), ppyyp(), ppxpyp(), ppxvx(), ppyvy()
# ppxvz(), ppyvz(), pprvz(), ppzxy(), pptrace()
# ppco(y,x,z;uz,xmin,xmax,ymin,ymax,zmin,zmax)
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
  if lparallel:
    totalnp = parallelsum(top.nps[js])
    fracnp = float(top.nps[js])/float(totalnp)
  else:
    fracnp = 1.
  for i in xrange(0,len(top.npplot)):
    ntopick=min(top.nps[js],int(top.npplot[i]*fracnp+0.5))
    ii = arrayrange(top.nps[0])
    rr = top.nps[0]*RandomArray.random(top.nps[0])
    ii = compress(less(rr,ntopick),ii)
    psubset.append(ii.astype('i'))
#----------------------------------------------------------------------------
def clear_subsets():
  "Clears the particle subsets so that they can be updated."
  global psubset
  psubset = []

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
def checkarguments(input,arglist):
  "Compare inputs against and argument list and return list of bad arguments"
  inputcopy = input.copy()
  for i in inputcopy.keys():
    if i in arglist.keys(): del inputcopy[i]
  return inputcopy


#-------------------------------------------------------------------------
# This returns the indices of the particles selected.
def selectparticles(iw=0,kwdict={},**kw):
  """
Selects particles based on either subsets or windows. By default it selects
from window 0, getting all of the live partilces (whose uzp > 0).
  - iw=0: Window to chose from
  - js=0: Species to chose from
  - jslist=None: List of Species to choose from, e.g. [0,3,4]; -1 for all specs
  - win=top.zwindows+top.zbeam: Windows to use (in lab frame)
  - z=top.zp: Coordinate for range selection
  - ix=-1: When 0 <= ix <= nx, picks particles within xmesh[ix]+-wx*dx
  - wx=1.: Width of window around xmesh[ix]
  - iy=-1: When 0 <= iy <= ny, picks particles within ymesh[iy]+-wy*dy
  - wy=1.: Width of window around ymesh[iy]
  - iz=-1: When 0 <= iz <= nz, picks particles within zmesh[iz]+-wz*dz
  - wz=1.: Width of window around zmesh[iz]
  - zl=None: When specified, lower range of selection region
  - zu=None: When specified, upper range of selection region
  """
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {"js":0,"jslist":None,"win":None,"z":None,
                "ix":None,"wx":1.,"iy":None,"wy":1.,"iz":None,"wz":1.,
                "zl":None,"zu":None,'checkargs':0,'allowbadargs':0}

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in kwdefaults.
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.
  badargs = checkarguments(kwvalues,kwdefaults)
  if checkargs: return badargs
  if badargs and not allowbadargs:
    raise "bad argument ",string.join(badargs.keys())

#  If jslist defined, call selectparticles repeatedly for each species on the list

  del kwvalues['jslist']    # Remove list so selectparticles is subsequently
                            # called with one species at a time
  if jslist is not None:
    if jslist == -1:    jslist = range(0,top.ns)
    partlist = array([])
    for js in jslist:
        kwvalues['js'] = js
        newparts = selectparticles(iw, kwvalues)
        partlist = array(list(partlist)+list(newparts))
    return partlist

  ir1 = top.ins[js]-1
  ir2 = top.ins[js]+top.nps[js]-1
  if ir2 <= ir1: return array([])
  if zl is not None or zu is not None:
    if z is None: z = top.zp
    if zl is None: zl = -top.largepos
    if zu is None: zu = +top.largepos
    if zl > zu: print "Warning: zl > zu"
    ii=compress(logical_and(less(zl,z[ir1:ir2]),less(z[ir1:ir2],zu)),
                arrayrange(ir1,ir2))
  elif ix is not None:
    xl = w3d.xmmin + ix*w3d.dx - wx*w3d.dx
    xu = w3d.xmmin + ix*w3d.dx + wx*w3d.dx
    ii=compress(logical_and(less(xl,top.xp[ir1:ir2]),less(top.xp[ir1:ir2],xu)),
                arrayrange(ir1,ir2))
  elif iy is not None:
    yl = w3d.ymmin + iy*w3d.dy - wy*w3d.dy
    yu = w3d.ymmin + iy*w3d.dy + wy*w3d.dy
    ii=compress(logical_and(less(yl,top.yp[ir1:ir2]),less(top.yp[ir1:ir2],yu)),
                arrayrange(ir1,ir2))
  elif iz is not None:
    z = top.zp
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
    if -iw > len(psubset): raise "Bad window number"
    ii = ir1 + compress(less(psubset[-iw-1],top.nps[js]),psubset[-iw-1])
  elif iw == 0:
    ii = xrange(ir1,ir2)
  else:
    if win is None: win = top.zwindows[:,iw] + top.zbeam
    if len(shape(win)) == 2: win = win[:,iw]
    if z is None: z = top.zp
    ii=compress(logical_and(less(win[0],z[ir1:ir2]),less(z[ir1:ir2],win[1])),
                arrayrange(ir1,ir2))
  ii = compress(not_equal(take(top.uzp,ii),0.),ii)
  return ii

#-------------------------------------------------------------------------
# The following return a specific coordinate of the selected particles
# More documetation added after they are declared.
#-------------------------------------------------------------------------
def getn(iw=0,gather=1,**kw):
  "Returns number of particles in selection."
  ii = selectparticles(iw=iw,kwdict=kw)
  if lparallel and gather: return globalsum(len(ii))
  else: return len(ii)
#-------------------------------------------------------------------------
def getx(iw=0,gather=1,**kw):
  "Returns the X positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.xp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def gety(iw=0,gather=1,**kw):
  "Returns the Y positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.yp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getz(iw=0,gather=1,**kw):
  "Returns the Z positions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.zp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getr(iw=0,gather=1,**kw):
  "Returns the R postions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = sqrt(take(top.xp,ii)**2 + take(top.yp,ii)**2)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def gettheta(iw=0,gather=1,**kw):
  "Returns the theta postions."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = arctan2(take(top.yp,ii),take(top.xp,ii))
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvx(iw=0,gather=1,**kw):
  "Returns the X velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uxp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvy(iw=0,gather=1,**kw):
  "Returns the Y velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uyp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getvz(iw=0,gather=1,**kw):
  "Returns the Z velocity."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uzp*top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getux(iw=0,gather=1,**kw):
  "Returns the X momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uxp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getuy(iw=0,gather=1,**kw):
  "Returns the Y momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uyp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getuz(iw=0,gather=1,**kw):
  "Returns the Z momentum over mass."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getxp(iw=0,gather=1,**kw):
  "Returns the X velocity over the Z velocity (X')."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uxp,ii)/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getyp(iw=0,gather=1,**kw):
  "Returns the Y velocity over the Z velocity (Y')."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.uyp,ii)/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getrp(iw=0,gather=1,**kw):
  "Returns the radial velocity over the Z velocity (R')."
  ii = selectparticles(iw=iw,kwdict=kw)
  tt = arctan2(take(top.yp,ii),take(top.xp,ii))
  result = (take(top.uxp,ii)*cos(tt)+take(top.uyp,ii)*sin(tt))/take(top.uzp,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
def getgaminv(iw=0,gather=1,**kw):
  "Returns the gamma inverse."
  ii = selectparticles(iw=iw,kwdict=kw)
  result = take(top.gaminv,ii)
  if lparallel and gather: return gatherarray(result)
  else: return result
#-------------------------------------------------------------------------
# Add the selectparticles documentation to each of the routines.
if sys.version[:5] != "1.5.1":
  if lparallel:
    _gatherdoc = "  gather=1 When 1, all data is gathered to PE0"
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

##########################################################################
def getxxpslope(iw=0,iz=-1):
  """
Calculates the x-x' slope based on either the window moments in window iw
or the zmoments at iz. This returns a tuple containing (slope,offset,vz).
The product slope*vz gives the slope for x-vx.
  """
  if not lparallel:
    if 0 <= iz <= w3d.nz:
      slope = (top.xxpbarz[iz]-top.xbarz[iz]*top.xpbarz[iz])/top.xrmsz[iz]**2
      offset = top.xpbarz[iz]-slope*top.xbarz[iz]
      vz = top.vzbarz[iz]
    else:
      iiw = max(0,iw)
      slope = (top.xxpbar[iiw]-top.xbar[iiw]*top.xpbar[iiw])/top.xrms[iiw]**2
      offset = top.xpbar[iiw]-slope*top.xbar[iiw]
      vz = top.vzbar[iiw]
  else:
    if 0 <= iz <= w3d.nzfull:
      pe = convertizptope(iz)
      if me == pe:
        iz = iz - top.izpslave[me]
        slope = (top.xxpbarz[iz]-top.xbarz[iz]*top.xpbarz[iz])/top.xrmsz[iz]**2
        offset = top.xpbarz[iz]-slope*top.xbarz[iz]
        vz = top.vzbarz[iz]
      else:
        (slope,offset,vz) = (0.,0.,0.)
      (slope,offset,vz) = tuple(broadcast(array([slope,offset,vz]),pe))
    else:
      iiw = max(0,iw)
      pe = convertiwtope(iiw)
      if me == pe:
        slope = (top.xxpbar[iiw]-top.xbar[iiw]*top.xpbar[iiw])/top.xrms[iiw]**2
        offset = top.xpbar[iiw]-slope*top.xbar[iiw]
        vz = top.vzbar[iiw]
      else:
        (slope,offset,vz) = (0.,0.,0.)
      (slope,offset,vz) = tuple(broadcast(array([slope,offset,vz]),pe))
  return (slope,offset,vz)
#-------------------------------------------------------------------------
def getyypslope(iw=0,iz=-1):
  """
Calculates the y-y' slope based on either the window moments in window iw
or the zmoments at iz. This returns a tuple containing (slope,offset,vz).
The product slope*vz gives the slope for y-vy.
  """
  if not lparallel:
    if 0 <= iz <= w3d.nz:
      slope = (top.yypbarz[iz]-top.ybarz[iz]*top.ypbarz[iz])/top.yrmsz[iz]**2
      offset = top.ypbarz[iz]-slope*top.ybarz[iz]
      vz = top.vzbarz[iz]
    else:
      iiw = max(0,iw)
      slope = (top.yypbar[iiw]-top.ybar[iiw]*top.ypbar[iiw])/top.yrms[iiw]**2
      offset = top.ypbar[iiw]-slope*top.ybar[iiw]
      vz = top.vzbar[iiw]
  else:
    if 0 <= iz <= w3d.nzfull:
      pe = convertizptope(iz)
      if me == pe:
        iz = iz - top.izpslave[me]
        slope = (top.yypbarz[iz]-top.ybarz[iz]*top.ypbarz[iz])/top.yrmsz[iz]**2
        offset = top.ypbarz[iz]-slope*top.ybarz[iz]
        vz = top.vzbarz[iz]
      else:
        (slope,offset,vz) = (0.,0.,0.)
      (slope,offset,vz) = tuple(broadcast(array([slope,offset,vz]),pe))
    else:
      iiw = max(0,iw)
      pe = convertiwtope(iiw)
      if me == pe:
        slope = (top.yypbar[iiw]-top.ybar[iiw]*top.ypbar[iiw])/top.yrms[iiw]**2
        offset = top.ypbar[iiw]-slope*top.ybar[iiw]
        vz = top.vzbar[iiw]
      else:
        (slope,offset,vz) = (0.,0.,0.)
      (slope,offset,vz) = tuple(broadcast(array([slope,offset,vz]),pe))
  return (slope,offset,vz)
#-------------------------------------------------------------------------
def getvzrange():
  "Returns a tuple containg the Vz range for plots"
  if (top.vzrng != 0.):
     vzmax = (1. + top.vtilt)*top.vbeam*(1.+top.vzrng) - top.vzshift
     vzmin = (1. - top.vtilt)*top.vbeam*(1.-top.vzrng) - top.vzshift
  else:
     vzmax = top.vzmaxp + 0.1*(top.vzmaxp-top.vzminp)
     vzmin = top.vzminp - 0.1*(top.vzmaxp-top.vzminp)
  return (vzmin,vzmax)

##########################################################################
def pptitleright(iw=0,kwdict={},**kw):
  "Returns right plot title. Takes same arguments as selectparticles"
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {"js":0,"win":None,"z":None,
                "ix":None,"wx":1.,"iy":None,"wy":1.,"iz":None,"wz":1.,
                "zl":None,"zu":None,"slope":0,
                'checkargs':0,'allowbadargs':0}

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in kwdefaults.
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.
  badargs = checkarguments(kwvalues,kwdefaults)
  if checkargs: return badargs
  if badargs and not allowbadargs:
    raise "bad argument ",string.join(badargs.keys())

  # --- Return appropriate right title
  if zl is not None or zu is not None:
    if z is None: prefix = ""
    else: prefix = "z "
    if zl is None: zl = -top.largepos
    if zu is None: zu = +top.largepos
    result = prefix+"range (%9.4e, %9.4e)"%(zl,zu)
  elif ix is not None:
    xl = w3d.xmmin + ix*w3d.dx - wx*w3d.dx
    xu = w3d.xmmin + ix*w3d.dx + wx*w3d.dx
    result = "ix = %d, x range (%9.4e, %9.4e)"%(ix,xl,xu)
  elif iy is not None:
    yl = w3d.ymmin + iy*w3d.dy - wy*w3d.dy
    yu = w3d.ymmin + iy*w3d.dy + wy*w3d.dy
    result = "iy = %d, y range (%9.4e, %9.4e)"%(iy,yl,yu)
  elif iz is not None:
    if lparallel:
      zl = top.zmslmin[0] + iz*w3d.dz - wz*w3d.dz + top.zbeam
      zu = top.zmslmin[0] + iz*w3d.dz + wz*w3d.dz + top.zbeam
    else:
      zl = w3d.zmmin + iz*w3d.dz - wz*w3d.dz + top.zbeam
      zu = w3d.zmmin + iz*w3d.dz + wz*w3d.dz + top.zbeam
    result = "iz = %d, z range (%9.4e, %9.4e)"%(iz,zl,zu)
  elif iw < 0:
    if psubset==[]: setup_subsets()
    result = "subset "+repr(-iw)+": "+repr(len(psubset[-iw-1]))+" particles"
  else:
    if win is None:
      win = top.zwindows[:,iw] + top.zbeam
      prefix = "z "
    else:
      prefix = ""
    if len(shape(win)) == 2: win = win[:,iw]
    result = prefix+"window%d = %9.4e, %9.4e"%(iw,win[0],win[1])
  if slope != 0:
    result = result + ", slope=%7.4f"%slope
  return result

#-------------------------------------------------------------------------
def ppmoments(text):
  "Plots text in upper right hand corner of the plot"
  plt(text,0.61,.855,justify="RT",height=12,font="courierB")

#############################################################################
#############################################################################
#############################################################################
def ppgeneric_doc(x,y):
  doc = selectparticles.__doc__ + """
  - zz: optional third particle data quantity - when supplied, it is deposited
       on a grid and that is used for contour levels.
  - grid: optional grid to plot (instead of deriving grid from particle data)
  - nx, ny: grid size, defaults to 20x20
  - slope=0.: slope to subtract from %(y)s coordinate (%(y)s-slope*%(x)s)
  - offset=0.: %(y)s-offset of particles
  - titles=1: when true, plot the titles
  - lframe=0: when true, the plot limits are set to the plmin and plmax input
              arguments, which default to the plmin and plmax variables from
              the group InDiag
  - pplimits=None: a tuple of (xplmin, xplmax, yplmin, yplmax), limits of plot
                   range (used when lframe=1)
  - xmin, xmax, ymin, ymax: extrema of density grid, defaults to particle
                            extrema (x for %(x)s and y for %(y)s)
  - particles=0: when true, plot particles
  - uselog=0: when true, logarithmic levels of the number density are used
  - color='fg': color of particles, when=='density', color by number density
  - ncolor=None: when plotting color by number density, number of colors to use,
                 defaults to top.ncolor
  - denmin, denmax: set extrema for coloring particles by number density
  - chopped=None: only particles where r < chopped*maxdensity/density
                  are plotted, where r is a random number between 0 and 1
                  and density is the density at the particle location
  - marker=dot: particle marker to plot
  - msize=1.: scaled size of marker
  - hash=0: flag to turn on or off the hash plot
  - line_scale=.9: scaling factor on line length
  - hcolor='fg': color of hash marks
  - width=1.0: width of hash marks
  - contours=None: number of countours to plot
  - filled=0: when true, plot filled contours (assumes contours is set)
  - ccolor='fg': contour color (when not filled)
  - cellarray=0: when true, plot grid as cell array
  - ctop=240, cmin=minnd(grid), cmax=maxnd(grid): control palette of cellarray
  - ldensityscale=0: when true, scale the density by its max.
  - view=1: view window to use (experts only)
  - colbarunitless=0: when true, color-bar scale is unitless
  - colbarlinear=1: when true, the color-bar is laid out linearly in density,
                    otherwise each contour level gets an equal sized area.
                    Only in effect when a list of colorbars is specified.
  - surface=0: when true, a 3-d surface plot is made of the gridded data
               The view can be changed dynamically with the function viewsurface
  """
  return doc%vars()
#-------------------------------------------------------------------------
def ppgeneric(y=None,x=None,kwdict={},**kw):
  """
Generic particle plotting routine. Allows plotting of particle points, density
contours, and/or density hash marks.
Note that either the x and y coordinates or the grid must be passed in.
  - y, x: optional particle data (instead of using inputted grid)
  """
  # --- Complete dictionary of possible keywords and their default values
  kwdefaults = {'zz':None,'grid':None,'nx':20,'ny':20,'slope':0.,
                'offset':0.,'titles':1,'lframe':0,
                'xmin':None,'xmax':None,'ymin':None,'ymax':None,
                'pplimits':('e','e','e','e'),
                'particles':0,'uselog':0,'color':'fg','ncolor':top.ncolor,
                'usepalette':1,'marker':'\1','msize':1.0,
                'denmin':None,'denmax':None,'chopped':None,
                'hash':0,'line_scale':.9,'hcolor':'fg','width':1.0,
                'contours':None,'filled':0,'ccolor':'fg',
                'cellarray':0,'ctop':199,
                'cmin':-top.largepos,'cmax':top.largepos,
                'ldensityscale':0,
                'view':1,'colbarunitless':0,'colbarlinear':1,'surface':0,
                'checkargs':0,'allowbadargs':0}

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in kwdefaults.
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.
  badargs = checkarguments(kwvalues,kwdefaults)
  if checkargs: return badargs
  assert (not badargs or allowbadargs), \
         "bad argument %s"%string.join(badargs.keys())

  # --- Do some error checking on the consistency of the input
  assert (type(grid) == ArrayType or \
          (type(x) == ArrayType and type(y) == ArrayType)), \
         "either the grid and/or both x and y must be specified"
  assert (not particles or (type(x) == ArrayType and type(y) == ArrayType)), \
         "both x and y must be specified if particles are to be plotted"
  assert ((type(x) != ArrayType and type(y) != ArrayType) or len(x) == len(y)),\
         "both x and y must be of the same length"
  assert (zz is None) or (type(zz) == ArrayType and len(zz) == len(x)),\
         "zz must be the same length as x"
  assert (type(slope) != StringType),"slope must be a number"
  assert (zz is None) or (grid is None),\
         "only one of zz and grid can be specified"

  # --- If there are no particles and no grid to plot, just return
  if type(x) == ArrayType and type(y) == ArrayType: np = globalsum(len(x))
  else: np = 0
  if np == 0 and grid is None: return

  # --- Make sure that nothing is not plotted over a surface plot
  if surface:
    particles = 0
    contours = 0
    hash = 0
    lframe = 0
    titles = 0

  # -- Set the plotting view window
  plsys(view)

  # --- Make sure that the grid size nx and ny are consistent with grid
  # --- is one is input
  if type(grid) == ArrayType:
    nx = shape(grid)[0] - 1
    ny = shape(grid)[1] - 1

  # --- Calculate extrema of the particles
  if type(x) == ArrayType and type(y) == ArrayType:
    # --- Get slope subtracted value of y
    yms = y - x*slope - offset
    # --- Get mins and maxs of particles that were not supplied by the user.
    if lparallel:
      if xmin is None: xmintemp = globalmin(x)
      if xmax is None: xmaxtemp = globalmax(x)
      if ymin is None: ymintemp = globalmin(yms)
      if ymax is None: ymaxtemp = globalmax(yms)
    else:
      xmintemp = 0.
      xmaxtemp = 0.
      ymintemp = 0.
      ymaxtemp = 0.
      if xmin is None and len(x) > 0: xmintemp = min(x)
      if xmax is None and len(x) > 0: xmaxtemp = max(x)
      if ymin is None and len(yms) > 0: ymintemp = min(yms)
      if ymax is None and len(yms) > 0: ymaxtemp = max(yms)
    # --- When neither the min or max are supplied by the user, extend
    # --- extrema by one grid cell so that all particles are within the
    # --- limits of the grid. This is the most common case.
    if xmin is None and xmax is None:
      xmintemp = xmintemp - (xmaxtemp-xmintemp)/(nx-2)
      xmaxtemp = xmaxtemp + (xmaxtemp-xmintemp)/(nx-2)
    if ymin is None and ymax is None:
      ymintemp = ymintemp - (ymaxtemp-ymintemp)/(ny-2)
      ymaxtemp = ymaxtemp + (ymaxtemp-ymintemp)/(ny-2)
    # --- Now set main versions of min and max
    if xmin is None: xmin = xmintemp
    if xmax is None: xmax = xmaxtemp
    if ymin is None: ymin = ymintemp
    if ymax is None: ymax = ymaxtemp
  else:
    # --- If no particles are inputted and the extrema are not set, then
    # --- can only make a guess.
    if xmin is None: xmin = 0
    if xmax is None: xmax = nx
    if ymin is None: ymin = 0
    if ymax is None: ymax = ny

  # --- Get grid cell sizes
  dx = (xmax-xmin)/nx
  dy = (ymax-ymin)/ny

  # --- If the grid is needed for the plot and it was not passed in, generate
  # --- it from the inputted particle data (if there was any)
  if type(grid) != ArrayType and \
     (hash or contours or color=='density' or chopped or surface or cellarray):
    if zz is None:
      densitygrid = 1

      # --- Create space for data
      grid = fzeros((1+nx,1+ny),'d')

      # --- Deposit the density onto the grid.
      setgrid2d(len(x),x,yms,nx,ny,grid,xmin,xmax,ymin,ymax)

      # --- If parallel, do a reduction on the grid
      if lparallel: grid = parallelsum(grid)

    else:
      densitygrid = 0

      # --- Create space for data
      grid = fzeros((1+nx,1+ny),'d')
      gridcount = fzeros((1+nx,1+ny),'d')

      # --- Deposit the data onto the grid. itask is 1 so that the parallel
      # --- version can be done properly.
      deposgrid2d(1,len(x),x,yms,zz,nx,ny,grid,gridcount,xmin,xmax,ymin,ymax)

      # --- If parallel, do a reduction on the grid
      if lparallel:
        grid = parallelsum(grid)
        gridcount = parallelsum(gridcount)

      # --- Divide out the particle counts by hand.
      grid = grid/where(greater(gridcount,0.),gridcount,1.)

  elif (hash or contours or color=='density' or chopped or surface):
    densitygrid = 0
 
  # --- Scale the grid by its maximum if requested.
  if ldensityscale and grid is not None:
    grid[:,:] = grid/maxnd(grid)

  # --- If using logarithmic number density levels, take the log of the grid
  # --- data. The original grid is left unchanged since that is still needed
  # --- by some operations below.
  if uselog:
    if densitygrid:
      # --- Take the log, raising all values below 0.1 to 0.1. The
      # --- threshold is used so that none of the elements are zero.
      # --- That value 0.1 is used since values any smaller do not have
      # --- much meaning since a value of 1.0 means that there is already
      # --- only one particle in that cell.
      grid1 = log10(where(less(grid,0.1),0.1,grid))
    else:
      # --- Before taking the log of the user supplied grid data, make sure
      # --- that there are no negative values. Zero is ok since they will
      # --- be replaced with a minimum value.
      dmax = maxnd(grid)
      dmin = minnd(where(equal(grid,0.),dmax,grid))
      if dmin <= 0.:
        raise "Can't take log since the grid has negative values"
      grid1 = log(where(less(grid,dmin/10.),dmin/10.,grid))
  else:
    grid1 = grid

  # --- Get grid mesh if it is needed
  if contours or hash or surface or cellarray:
    xmesh = xmin + dx*arange(nx+1)[:,NewAxis]*ones(ny+1,'d')
    ymesh = ymin + dy*arange(ny+1)*ones(nx+1,'d')[:,NewAxis]

  # --- Make filled contour plot of grid first since it covers everything
  # --- plotted before it.
  if contours and filled:
    if maxnd(grid1) != minnd(grid1):
      plotc(transpose(grid1),transpose(ymesh),transpose(xmesh),
            color=ccolor,contours=contours,filled=filled)

  # --- Make cell-array plot. This also is done early since it covers anything
  # --- done before it.
  if cellarray:
    cmin = max(cmin,minnd(grid1))
    cmax = min(cmax,maxnd(grid1))
    pli(transpose(grid1),xmin,ymin,xmax,ymax,top=ctop,cmin=cmin,cmax=cmax)

  # --- Plot particles
  if particles:
    if color == 'density':
      z1 = zeros(len(x),'d')
      getgrid2d(len(x),x,yms,z1,nx,ny,grid1,xmin,xmax,ymin,ymax)
    if chopped:
      dd = zeros(len(x),'d')
      getgrid2d(len(x),x,yms,dd,nx,ny,grid,xmin,xmax,ymin,ymax)
      maxdensity = maxnd(grid)
      npart = len(x)
      ipick = less(RandomArray.random(shape(x)),maxdensity*chopped/dd)
      x = compress(ipick,x)
      yms = compress(ipick,yms)
      if color == 'density':
        z1 = compress(ipick,z1)
    if color == 'density':
      # --- Plot particles with color based on the density from the grid.
      if denmin is None: demmin = minnd(grid)
      if denmax is None: demmax = maxnd(grid)
      ppco(yms,x,z1,uz=1.,marker=marker,msize=msize,lframe=0,
           xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=denmin,zmax=denmax,
           ncolor=ncolor,usepalette=usepalette)
    else:
      # --- Plot particles as a solid color.
      warpplp(yms,x,color=color,marker=marker,msize=msize)

  # --- Now plot unfilled contours, which are easier to see on top of the
  # --- particles
  if contours and not filled:
    if maxnd(grid1) != minnd(grid1):
      plotc(transpose(grid1),transpose(ymesh),transpose(xmesh),
            color=ccolor,contours=contours,filled=filled)

  # --- Plot hash last since it easiest seen on top of everything else.
  if hash:
    # --- Set line length
    sss = line_scale*(xmax-xmin)/nx/maxnd(grid1)
    # --- Make plot of tick marks
    for ix in range(nx+1):
      for iy in range(ny+1):
        plg(ymesh[ix,iy]+zeros(2),xmesh[ix,iy]+array([0.,sss*grid1[ix,iy]]),
            color=hcolor,width=width)

  # --- Add colorbar if needed
  if (contours and filled==1) or (color == 'density' and len(x) > 0) or \
     (cellarray):
    if (contours and filled==1):
      try:
        nc = len(contours) + 1
        levs = contours
      except TypeError:
        nc = contours + 1
        levs = None
    elif (color == 'density' and len(x) > 0):
      nc = ncolor + 1
      levs = None
    elif (cellarray):
      nc = ctop
      levs = None
    if colbarunitless:
      dmax = 1.0
      dmin = 0.
    elif cellarray:
      dmin = cmin
      dmax = cmax
    else:
      dmax = maxnd(grid1)
      dmin = minnd(grid1)
    colorbar(dmin,dmax,uselog=uselog,ncolor=nc,view=view,levs=levs,
             colbarlinear=colbarlinear,ctop=ctop)

  # --- Make surface plot
  if surface and me == 0:
    pl3d.orient3()
    pl3d.light3()
    plwf.plwf(grid1,xmesh,ymesh,fill=grid1,edges=0)
    [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
    #limits(xmin3,xmax3,ymin3,ymax3)

  # --- Finish off the plot, adding titles and setting the frame limits.
  if titles: ptitles(v=view)
  settitles() 
  if (lframe):
    limits(pplimits[0],pplimits[1],pplimits[2],pplimits[3])
if sys.version[:5] != "1.5.1":
  ppgeneric.__doc__ = ppgeneric.__doc__ + ppgeneric_doc('x','y')


#############################################################################
#############################################################################
# ed williams' colorbar stuff / modified for Warp by J.-L. Vay on 01/22/2001
def nicelevels(z,n=8) :
  """nicelevels(z,n=8) finds approximately n "nice values"
between min(z) and max(z) for axis labels. n defaults to eight.
  """
  zmax = max(ravel(z))
  zmin = min(ravel(z))
  if zmin == zmax: return array([zmin,zmax])
  finest = abs(zmax - zmin)/float (n)
  # blows up on zmin=zmax
  unit = 10.**floor (log10 (finest))
  finest = finest/unit
  if   finest > 5.0: finest = 10.
  elif finest > 2.:  finest = 5.
  elif finest > 1.:  finest = 2.
  unit = unit*finest
  cmin = unit*ceil(zmin/unit)
  if (abs(cmin - zmin) < 0.01*unit) :
     cmin = cmin
  cmax = unit*floor(zmax/unit)
  if (abs(cmax - zmax) < 0.01*unit) :
     cmax = cmax
  n = int(((cmax - cmin)/unit + 0.5) + 1)
  levs = cmin + arange(n)*unit
  llist = nonzero(less(abs(levs),0.1*unit))
  if len(llist) > 0:
     array_set(levs,llist,0.0)
  return levs

#-----------------------------------------------------------------------
def colorbar(zmin,zmax,uselog=0,ncolor=100,view=1,levs=None,colbarlinear=1,
             ctop=199):
  """
Plots a color bar to the right of the plot square labelled by the z
values from zmin to zmax.
  - zmin, zmax: lower and upper range for color bar
  - uselog=0: when true, labels are printed in the form 10^x
  - ncolor=100: default number of colors to include
  - view=1: specifies the view that is associated with the color bar
  - levs: an optional list of color levels
  - ctop=199: number of colors from palette to use
  """
  plsys(0)
  xmin = 0.66
  xmax = 0.68
  ymax = 0.85
  ymin = 0.44
  # --- draw the bar
  if colbarlinear and levs is not None:
    # --- Use the contour plotting routine plfc for this case. The data
    # --- plotted is uniformly spaced between zmin and zmax. The contour
    # --- levels are those specified. The result is that the colorbar
    # --- shows the contours levels by their values relative to zmin and zmax.
    plotval = span(zmin,zmax,255)[:,NewAxis]*ones(2)
    xx = array([xmin,xmax])*ones(255)[:,NewAxis]
    yy = span(ymin,ymax,255)[:,NewAxis]*ones(2)
    ireg = ones((255,2))
    plfc(plotval,yy,xx,ireg,contours=array(levs))
  else:
    # --- Use cell array plotting for this case. All of the colors get a block
    # --- of the same size. If levs is not specified, the uniform spacing 
    # --- matches the uniform spacing of the contours. If levs is specified,
    # --- each equal sized block represents one contour level, independent of
    # --- the range of the level relative to other levels.
    if type(zmin) == type(zmax) == type(1):
       plotval = arange(zmin,zmax+1,typecode='b')[:,NewAxis]*ones(2)
    else:
       plotval = (arange(ncolor)/(ncolor-1.))[:,NewAxis]*ones(2)
    pli(plotval,xmin,ymin,xmax,ymax,top=ctop)
  # --- Draw a black box around it
  pldj([xmin,xmin,xmin,xmax],[ymin,ymax,ymin,ymin],
       [xmax,xmax,xmin,xmax],[ymin,ymax,ymax,ymax])
  # --- Generate nice levels for the labels and tick marks.
  if levs is None:
    # --- Use the nicelevels routine to get evenly spaced labels.
    nicelevs = nicelevels(array([zmin,zmax]))
  else:
    # --- If there are less than 15 specified contour levels, put a label
    # --- at each of the labels. If there are more, pick out roughly 15
    # --- evenly spaced values. Also, if the levels do not extend beyond
    # --- zmin and zmax, add labels at those points too.
    nicelevs = levs
    if zmin < levs[0]:  nicelevs = array([zmin] + list(nicelevs))
    if zmax > levs[-1]: nicelevs = array(list(nicelevs) + [zmax])
  llev = len(nicelevs)
  # --- Create the labels
  labels = []
  # --- Calculate the location of the labels.
  if not colbarlinear and levs is not None:
    # --- The ys are evenly spaced
    ys = ymin + arange(len(nicelevs))/(len(levs)+1.)*(ymax - ymin)
    # --- If the lowest level is less than zmin, then bump up the y's
    # --- by one block size.
    if levs[0] < zmin: ys = ys + 1./(len(levs)+1.)*(ymax - ymin)
  elif llev==2 and (nicelevs[0] == nicelevs[1]):
    ys = array([ymin,ymax])
  else:
    ys = ymin + (ymax - ymin)*(nicelevs - zmin)/(zmax - zmin)
  # --- Plot the labels, skipping ones that are too close together.
  if uselog: ss = " 10^%.5g"
  else:      ss = " %.5g"
  ylast = 0.
  for i in xrange(llev):
    if ys[i] - ylast > (ymax-ymin)/30:
      plt(ss%nicelevs[i],xmax+0.005,ys[i]-0.005)
      ylast = ys[i]
  # --- Plot the tick marks
  pldj(llev*[xmax],ys,llev*[xmax+0.005],ys)
  # --- Return to plot system 1.
  plsys(view)

#############################################################################
#############################################################################
def changepalette(returnpalette=0,help=0):
  """
Dynamically change the color palette.
  - returnpalette=0: when true, returns tuple of (red, green, blue)
  - help=0: when true, prints this message
Mouse actions:
  Button 1: shifts a point, compressing and stretching the rest of the colors
  Button 2: reset palette to original
  Button 3: shifts a point, sliding the colors up and down
  Shift Button 1: reverse the palette
  Control Button 1: add black point
  Control Button 3: add white point
  """
  # --- Print out help if wanted
  if help: print changepalette.__doc__
  # --- min's and max's are the same as in the colorbar routine
  xmin = 0.66
  xmax = 0.68
  ymax = 0.85
  ymin = 0.44
  # --- Create storate arrays
  # --- rr, gg, bb hold the original palette
  rr = zeros(200,'b')
  gg = zeros(200,'b')
  bb = zeros(200,'b')
  palette(rr,gg,bb,query=1)
  # --- newrr, newgg, newbb hold the new palette
  newrr = zeros(200,'b')
  newgg = zeros(200,'b')
  newbb = zeros(200,'b')
  # --- position relative to the original palette
  cc = arange(0,200)*1.
  newcc = arange(0,200)*1.
  # --- List of black and white points
  blacklist = []
  whitelist = []
  while 1:
    mm = mouse(0,0,"")
    if mm == None: break
    # --- Get mouse positions. Skip if outside the colorbar
    (x1, y1, x2, y2) = tuple(mm[:4])
    if x1 < xmin or x1 > xmax or x2 < xmin or x2 > xmax: continue
    if y1 < ymin or y1 > ymax or y2 < ymin or y2 > ymax: continue

    if mm[9] == 1 and mm[10] == 0:
      # --- Button 1, no keys
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      i2 = nint((y2 - ymin)/(ymax - ymin)*200)
      up = (ymax - y1)/(ymax - y2)
      down = (y1 - ymin)/(y2 - ymin)
      for i in xrange(1,i2):
        iold = int(i*down)
        wold =     i*down - iold
        newcc[i] = cc[iold]*(1.-wold) + cc[iold+1]*wold
      for i in xrange(i2,199):
        iold = 199 - int((199-i)*up)
        wold = iold - (199 -    ((199-i)*up))
        newcc[i] = cc[iold]*(1.-wold) + cc[iold-1]*wold

    if mm[9] == 2:
      # --- Button 2, no keys
      # --- Restore original palette
      newcc = arange(0,200)*1.
      blacklist = []
      whitelist = []

    if mm[9] == 3:
      # --- Button 3, no keys
      # --- slide whole palette
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      i2 = nint((y2 - ymin)/(ymax - ymin)*200)
      for i in xrange(0,200):
        iold = i - (i2 - i1)
        if iold < 0: newcc[i] = cc[0]
        elif iold > 199: newcc[i] = cc[-1]
        else: newcc[i] = cc[iold]

    if mm[9] == 1 and mm[10] == 1:
      # --- Button 1, shift
      # --- Reverse the palette
      newcc[:] = cc[::-1]

    if mm[9] == 1 and mm[10] == 4:
      # --- button 1, control
      # --- Add black point
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      blacklist.append(i1)

    if mm[9] == 3 and mm[10] == 4:
      # --- button 3, control
      # --- Add white point
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      whitelist.append(i1)

    # --- Calculate the new palette based on the position relative to the
    # --- original palette.
    for i in xrange(0,200):
      ii = int(newcc[i])
      ww =     newcc[i]  - ii
      iip1 = min(ii+1,199)
      newrr[i] = (nint(rr[ii]*(1.-ww) + rr[iip1]*ww))
      newgg[i] = (nint(gg[ii]*(1.-ww) + gg[iip1]*ww))
      newbb[i] = (nint(bb[ii]*(1.-ww) + bb[iip1]*ww))
    for ii in blacklist: (newrr[ii], newgg[ii], newbb[ii]) = 0,0,0
    for ii in whitelist: (newrr[ii], newgg[ii], newbb[ii]) = 255,255,255
    cc[:] = newcc
    palette(newrr,newgg,newbb)
  if returnpalette: return (newrr,newgg,newbb)

#############################################################################
#############################################################################
def viewsurface(scale=4.,gnomon=1):
  """
Dynamically view a surface plot. The mouse is used to change to view angle.
With button 1 pushed, the horizontal movement changes the z angle, and
vertical the y angle. With button 2 pressed, horizontal changes the x angle.
When finished, press return in the python window.
  - scale=4.: multiplicative factor to convert mouse movement to angle change
  """
  pl3d.gnomon(gnomon)
  [xmin3min,xmax3max,ymin3min,ymax3max,sys] = limits()
  while 1:
    mm = mouse(0,0,"")
    if mm == None: break
    (xa, ya, za) = (0.,0.,0.)
    if mm[9] == 1:
      ya = - (mm[3] - mm[1])*scale
      za = - (mm[2] - mm[0])*scale
    if mm[9] == 3:
      xa = (mm[2] - mm[0])*scale
    pl3d.rot3(xa,ya,za)
    [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
    xmin3min = min(xmin3min,xmin3)
    xmax3max = max(xmax3max,xmax3)
    ymin3min = min(ymin3min,ymin3)
    ymax3max = max(ymax3max,ymax3)
    limits(xmin3min,xmax3max,ymin3min,ymax3max)
  pl3d.gnomon(gnomon)

def _viewsurfacetest(scale=4.,gnomon=1):
  """
Dynamically view a surface plot. The mouse is used to change to view angle.
With button 1 pushed, the horizontal movement changes the z angle, and
vertical the y angle. With button 2 pressed, horizontal changes the x angle.
When finished, press return in the python window.
  - scale=4.: multiplicative factor to convert mouse movement to angle change
  """
  pl3d.gnomon(gnomon)
  pl3d.orient3(phi=0.,theta=0.)
  [xmin3min,xmax3max,ymin3min,ymax3max] = pl3d.draw3(1)
  phi = 0.
  theta = 0.
  (xa, ya, za) = (0.,0.,0.)
  while 1:
    mm = mouse(0,0,"")
    if mm == None: break
    dphi   = (mm[3] - mm[1])*scale
    dtheta = (mm[2] - mm[0])*scale
    print theta,phi
    newxa = xa + dtheta*sin(phi)*cos(theta) + dphi*cos(phi)*cos(theta)
    newya = ya + dtheta*sin(phi)*sin(theta) + dphi*cos(phi)*sin(theta)
    newza = za + dtheta*cos(phi)*cos(theta) + dphi*sin(phi)*sin(theta)
    phi = xa*cos(za) + ya*sin(za)
    theta = za
    pl3d.rot3(newxa-xa,newya-ya,newza-za)
    xa = newxa
    ya = newya
    za = newza
    [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
    xmin3min = min(xmin3min,xmin3)
    xmax3max = max(xmax3max,xmax3)
    ymin3min = min(ymin3min,ymin3)
    ymax3max = max(ymax3max,ymax3)
    limits(xmin3min,xmax3max,ymin3min,ymax3max)
  pl3d.gnomon(gnomon)

#############################################################################
#############################################################################
def ppmultispecies(pp,args,kw):
  """checks if js defined and assign it to a list if plotting multispecies.
  Also assign colors accordingly
  """
  if 'js' in kw.keys():
    js = kw['js']
    if js != -1 and type(js) != ListType:
      return false
    else:
      if js == -1: js = range(top.ns)
      if 'color' in kw.keys(): color = kw['color']
      else: color = range(0,240,240/len(js))
      for i in xrange(len(js)):
        kw['js'] = js[i]
        kw['color'] = color[i]
        apply(pp,args,kw)
      return true
  else:
    return false

########################################################################
########################################################################
########################################################################
########################################################################
def checkparticleplotarguments(kw):
  """Convenience routine to check arguments of particle plot routines.
Warning: this has the side affect of adding the arguement allowbadargs to
the kw dictionary. This is done since the calls to these functions here to
make the plots may have unused arguements since the entire kw list passed
into each of the pp plotting routines is passed into each of these
functions.
  """
  badargs = selectparticles(checkargs=1,kwdict=kw)
  badargs = pptitleright(checkargs=1,kwdict=badargs)
  badargs = ppgeneric(checkargs=1,kwdict=badargs)
  kw['allowbadargs'] = 1
  if badargs: raise "bad arguments ",string.join(badargs.keys())
########################################################################
def ppzxy(iw=0,particles=1,**kw):
  "Plots Z-X and Z-Y in single page"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzxy,(iw,particles),kw): return
  kw['particles'] = particles
  kw['view'] = 9
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  ii = selectparticles(iw=iw,win=top.ywindows,z=top.yp,kwdict=kw)
  settitles("X vs Z","Z","X",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.xp,ii),take(top.zp,ii),kwdict=kw)

  kw['view'] = 10
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.yplmin,top.yplmax)
  ii = selectparticles(iw=iw,win=top.xwindows,z=top.xp,kwdict=kw)
  settitles("Y vs Z","Z","Y",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.yp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzxy.__doc__ = ppzxy.__doc__ + ppgeneric_doc('z','x')

##########################################################################
def ppzx(iw=0,particles=1,**kw):
  "Plots Z-X"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzx,(iw,particles),kw): return
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("X vs Z","Z","X",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.xp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzx.__doc__ = ppzx.__doc__ + ppgeneric_doc('z','x')

##########################################################################
def ppzy(iw=0,particles=1,**kw):
  "Plots Z-Y"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzy,(iw,particles),kw): return
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.yplmin,top.yplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Y vs Z","Z","Y",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.yp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzy.__doc__ = ppzy.__doc__ + ppgeneric_doc('z','y')

##########################################################################
def ppzr(iw=0,particles=1,**kw):
  "Plots Z-R"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzr,(iw,particles),kw): return
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("R vs Z","Z","R",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(sqrt(top.xp**2+top.yp**2),ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzr.__doc__ = ppzr.__doc__ + ppgeneric_doc('z','r')

##########################################################################
def ppzxp(iw=0,particles=1,**kw):
  "Plots Z-X'"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzxp,(iw,particles),kw): return
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin,top.xpplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("X' vs Z","Z","X'",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uxp,ii)/take(top.uzp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzxp.__doc__ = ppzxp.__doc__ + ppgeneric_doc('z',"x'")

##########################################################################
def ppzvx(iw=0,particles=1,**kw):
  "Plots Z-Vx"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvx,(iw,particles),kw): return
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Vx vs Z","Z","Vx",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uxp,ii)*take(top.gaminv,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzvx.__doc__ = ppzvx.__doc__ + ppgeneric_doc('z',"vx")

##########################################################################
def ppzyp(iw=0,particles=1,**kw):
  "Plots Z-Y'"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzyp,(iw,particles),kw): return
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.ypplmin,top.ypplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Y' vs Z","Z","Y'",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uyp,ii)/take(top.uzp,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzyp.__doc__ = ppzyp.__doc__ + ppgeneric_doc('z',"y'")

##########################################################################
def ppzvy(iw=0,particles=1,**kw):
  "Plots Z-Vy"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvy,(iw,particles),kw): return
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.ypplmin*top.vbeam,top.ypplmax*top.vbeam)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Vy vs Z","Z","Vy",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uyp,ii)*take(top.gaminv,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzvy.__doc__ = ppzvy.__doc__ + ppgeneric_doc('z',"vy")

##########################################################################
def ppzvz(iw=0,particles=1,**kw):
  "Plots Z-Vz"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvz,(iw,particles),kw): return
  kw['particles'] = particles
  (vzmin,vzmax) = getvzrange()
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,vzmin,vzmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Vz vs Z","Z","Vz",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uzp,ii)*take(top.gaminv,ii),take(top.zp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppzvz.__doc__ = ppzvz.__doc__ + ppgeneric_doc('z',"vz")

##########################################################################
def ppxy(iw=0,particles=1,**kw):
  "Plots X-Y"
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxy,(iw,particles),kw): return
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,top.yplmin,top.yplmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Y vs X","X","Y",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.yp,ii),take(top.xp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxy.__doc__ = ppxy.__doc__ + ppgeneric_doc('x','y')

##########################################################################
def ppxxp(iw=0,iz=None,slope=0.,offset=0.,particles=1,lmoments=0,**kw):
  "Plots X-X'. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxxp,(iw,iz,slope,offset,particles,lmoments),kw): return
  if type(slope) == type(''): (slope,offset,vz) = getxxpslope(iw=iw,iz=iz)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)
  kw['iz'] = iz
  kw['slope'] = slope
  kw['offset'] = offset
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("X' vs X","X","X'",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uxp,ii)/take(top.uzp,ii),take(top.xp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxxp.__doc__ = ppxxp.__doc__ + ppgeneric_doc("x","x'")

##########################################################################
def ppyyp(iw=0,iz=None,slope=0.,offset=0.,particles=1,**kw):
  "Plots Y-Y'. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyyp,(iw,iz,slope,offset,particles),kw): return
  if type(slope) == type(''): (slope,offset,vz) = getyypslope(iw=iw,iz=iz)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax)
  kw['iz'] = iz
  kw['slope'] = slope
  kw['offset'] = offset
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Y' vs Y","Y","Y'",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uyp,ii)/take(top.uzp,ii),take(top.yp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppyyp.__doc__ = ppyyp.__doc__ + ppgeneric_doc("y","y'")

##########################################################################
def ppxpyp(iw=0,iz=None,slope=0.,offset=0.,particles=1,**kw):
  "Plots X'-Y'. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxpyp,(iw,iz,slope,offset,particles),kw): return
  if type(slope) == type(''):
    (xslope,xoffset,vz) = getxxpslope(iw=iw,iz=iz)
    (yslope,yoffset,vz) = getyypslope(iw=iw,iz=iz)
  else:
    (xslope,xoffset) = (slope,0.)
    (yslope,yoffset) = (slope,0.)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)
  kw['iz'] = iz
  settitles("Y' vs X'","X'","Y'",pptitleright(iw=iw,kwdict=kw))
  ii = selectparticles(iw=iw,kwdict=kw)
  xp = take(top.uxp,ii)/take(top.uzp,ii) - xslope*take(top.xp,ii) - xoffset
  yp = take(top.uyp,ii)/take(top.uzp,ii) - yslope*take(top.yp,ii) - yoffset
  ppgeneric(yp,xp,kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxpyp.__doc__ = ppxpyp.__doc__ + ppgeneric_doc("x'","y'")

##########################################################################
def ppxvx(iw=0,iz=None,slope=0.,offset=0.,particles=1,**kw):
  "Plots X-Vx. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxvx,(iw,iz,slope,offset,particles),kw): return
  if type(slope) == type(''):
    (slope,offset,vz) = getxxpslope(iw=iw,iz=iz)
    kw['slope'] = slope*vz
    kw['offset'] = offset*vz
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw['iz'] = iz
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Vx vs X","X","Vx",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uxp,ii),take(top.xp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxvx.__doc__ = ppxvx.__doc__ + ppgeneric_doc("x","Vx")

##########################################################################
def ppyvy(iw=0,iz=None,slope=0.,offset=0.,particles=1,**kw):
  "Plots Y-Vy. If slope='auto', it is calculated from the moments."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyvy,(iw,iz,slope,offset,particles),kw): return
  if type(slope) == type(''):
    (slope,offset,vz) = getyypslope(iw=iw,iz=iz)
    kw['slope'] = slope*vz
    kw['offset'] = offset*vz
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,
                      top.ypplmin*top.vbeam,top.ypplmax*top.vbeam)
  kw['iz'] = iz
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Vy vs Y","Y","Vy",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uyp,ii)*take(top.gaminv,ii),take(top.yp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppyvy.__doc__ = ppyvy.__doc__ + ppgeneric_doc("y","Vy")

##########################################################################
def ppxvz(iw=0,particles=1,**kw):
  "Plots X-Vz."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxvz,(iw,particles),kw): return
  (vzmin,vzmax) = getvzrange()
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,vzmin,vzmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Vz vs X","X","Vz",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uzp,ii)*take(top.gaminv,ii),take(top.xp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppxvz.__doc__ = ppxvz.__doc__ + ppgeneric_doc("x","Vz")

##########################################################################
def ppyvz(iw=0,particles=1,**kw):
  "Plots Y-Vz."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyvz,(iw,particles),kw): return
  (vzmin,vzmax) = getvzrange()
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,vzmin,vzmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Vz vs Y","Y","Vz",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(take(top.uzp,ii)*take(top.gaminv,ii),take(top.yp,ii),kwdict=kw)
if sys.version[:5] != "1.5.1":
  ppyvz.__doc__ = ppyvz.__doc__ + ppgeneric_doc("y","Vz")

##########################################################################
def pprrp(iw=0,scale=0,slope=0.,particles=1,**kw):
  """Plots R-R', If slope='auto', it is calculated from the moments.
  - scale=0: when true, scale particle by 2*rms"""
  checkparticleplotarguments(kw)
  if ppmultispecies(pprrp,(iw,scale,slope,particles),kw): return
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
  ii = selectparticles(iw=iw,kwdict=kw)
  xx = take(top.xp,ii)/xscale
  yy = take(top.yp,ii)/yscale
  xp = take(top.uxp,ii)/take(top.uzp,ii)/xpscale
  yp = take(top.uyp,ii)/take(top.uzp,ii)/ypscale
  rr = sqrt(xx**2 + yy**2)
  tt = arctan2(yy,xx)
  rp = xp*cos(tt) + yp*sin(tt)
  if type(slope) == type(''):
    if lparallel:
      aversq = globalave(rr**2)
      averrp = globalave(rr*rp)
    else:
      aversq = ave(rr**2)
      averrp = ave(rr*rp)
    if aversq > 0.:
      slope = averrp/aversq
    else:
      slope = 0.
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                      top.xpplmin/xpscale,top.xpplmax/ypscale)
  kw['slope'] = slope
  settitles("R' vs R","R","R'",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(rp,rr,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pprrp.__doc__ = pprrp.__doc__ + ppgeneric_doc("r","r'")

##########################################################################
def pprvz(iw=0,particles=1,**kw):
  "Plots R-Vz"
  checkparticleplotarguments(kw)
  if ppmultispecies(pprvz,(iw,particles),kw): return
  (vzmin,vzmax) = getvzrange()
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (0.,max(top.xplmax,top.yplmax),vzmin,vzmax)
  ii = selectparticles(iw=iw,kwdict=kw)
  settitles("Vz vs R","R","Vz",pptitleright(iw=iw,kwdict=kw))
  rr = sqrt(take(top.xp,ii)**2 + take(top.yp,ii)**2)
  vz = take(top.uzp,ii)*take(top.gaminv,ii)
  ppgeneric(vz,rr,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pprvz.__doc__ = pprvz.__doc__ + ppgeneric_doc("r","vz")

##########################################################################
def pptrace(iw=0,slope=0.,iz=-1,particles=1,titles=1,pplimits=None,**kw):
  """
Plots X-Y, X-X', Y-Y', Y'-X' in single page
If slope='auto', it is calculated from the moments for X-X' and Y-Y' plots.
pplimits can be a list of up to four tuples, one for each phase space plot.
If any of the tuples are empty, the limits used will be the usual ones for
that plot.
  """
  checkparticleplotarguments(kw)
  if ppmultispecies(pptrace,(iw,slope,iz,particles,titles,pplimits),kw): return
  ii = selectparticles(iw=iw,kwdict=kw)
  x = take(top.xp,ii)
  y = take(top.yp,ii)
  xp = take(top.uxp,ii)/take(top.uzp,ii)
  yp = take(top.uyp,ii)/take(top.uzp,ii)
  titler=pptitleright(iw=iw,kwdict=kw)
  kw['iz'] = iz
  kw['particles'] = particles
  kw['titles'] = titles
  if titles: ptitles(titler=titler)
  defaultpplimits = [(top.xplmin,top.xplmax,top.yplmin,top.yplmax),
                     (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax),
                     (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax),
                     (top.ypplmin,top.ypplmax,top.xpplmin,top.xpplmax)]
  if pplimits is None:
    pplimits = defaultpplimits
  else:
    kw['lframe'] = 1
    if type(pplimits[0]) != type(()):
      pplimits = 4*[pplimits]
    else:
      for i in xrange(4):
        if i == len(pplimits): pplimits.append(defaultpplimits[i])
        if not pplimits[i]: pplimits[i] = defaultpplimits[i]

  kw['view'] = 3
  kw['pplimits'] = pplimits[0]
  settitles("Y vs X","X","Y")
  ppgeneric(y,x,kwdict=kw)
 
  kw['view'] = 4
  kw['pplimits'] = pplimits[1]
  if type(slope)==type(''): kw['slope'] = getyypslope(iw=iw,iz=iz)[0]
  settitles("Y' vs Y","Y","Y'")
  ppgeneric(yp,y,kwdict=kw)
 
  kw['view'] = 5
  kw['pplimits'] = pplimits[2]
  if type(slope)==type(''): kw['slope'] = getxxpslope(iw=iw,iz=iz)[0]
  settitles("X' vs X","X","X'")
  ppgeneric(xp,x,kwdict=kw)
 
  kw['view'] = 6
  kw['pplimits'] = pplimits[3]
  if type(slope)==type(''): kw['slope'] = 0.
  settitles("X' vs Y'","Y'","X'")
  ppgeneric(xp,yp,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pptrace.__doc__ = pptrace.__doc__ + ppgeneric_doc("x","x'")

##########################################################################
##########################################################################
##########################################################################
def ppzxco(js=0,marker="\1",msize=1.0,lframe=0,sys=1,titles=1,
           ncolor=None,nskipcol=None,nstepcol=None):
  """Plots Z-X with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - sys=1 specifies section of plot frame to use
     - titles=1 specifies whether or not to plot titles"""
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: ncolor = top.nskipcol
  if nstepcol is None: ncolor = top.nstepcol
  inp=top.nps[js]/ncolor
  istep=nskipcol*nstepcol
  #if (lframadv) nf
  istart = top.ins[js]-1
  if (inp < istep): istep = 1
  for ij in xrange(1,istep+1,nskipcol*2):
    for ic in xrange(1,ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.xp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
    for ic in xrange(ncolor,0,-1):
      irs1 = istart+ij+nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-nskipcol-inp*(ic-1))/istep
      warpplp(take(top.xp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
  if titles: ptitles(" ","Z","X"," ",sys)
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.xplmin,top.xplmax)

##########################################################################
def ppzyco(js=0,marker="\1",msize=1.0,lframe=0,sys=1,titles=1,
           ncolor=None,nskipcol=None,nstepcol=None):
  """Plots Z-Y with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - sys=1 specifies section of plot frame to use
     - titles=1 specifies whether or not to plot titles"""
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: ncolor = top.nskipcol
  if nstepcol is None: ncolor = top.nstepcol
  inp=top.nps[js]/ncolor
  istep=nskipcol*nstepcol
  #if (lframadv) nf
  istart = top.ins[js]-1
  if (inp < istep): istep = 1
  for ij in xrange(1,istep+1,nskipcol*2):
    for ic in xrange(1,ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.yp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
         color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
    for ic in xrange(ncolor,0,-1):
      irs1 = istart+ij+nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-nskipcol-inp*(ic-1))/istep
      warpplp(take(top.yp[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
           color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
  if titles: ptitles(" ","Z","Y"," ",sys)
  if (lframe):
    limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,top.yplmin,top.yplmax)

##########################################################################
def ppzxyco(js=0,marker="\1",msize=1.0,lframe=0,titles=1,
            ncolor=None,nskipcol=None,nstepcol=None):
  """Plots Z-X and Z-Y in single frame with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  plsys(9)
  ppzxco(js,sys=9,marker=marker,msize=msize,lframe=lframe,titles=titles,
         ncolor=ncolor,nskipcol=nskipcol,nstepcol=nstepcol)
  plsys(10)
  ppzyco(js,sys=10,marker=marker,msize=msize,lframe=lframe,titles=titles,
         ncolor=ncolor,nskipcol=nskipcol,nstepcol=nstepcol)

##########################################################################
def ppzvzco(js=0,marker="\1",msize=1.0,lframe=0,titles=1,
            ncolor=None,nskipcol=None,nstepcol=None):
  """Plots Z-Vz with color based in paricle index
     - js=0 species to plot
     - marker='\1' marker type (see gist manual for the list)
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles"""
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: ncolor = top.nskipcol
  if nstepcol is None: ncolor = top.nstepcol
  inp=top.nps[js]/ncolor
  istep=nskipcol*nstepcol
  #if (lframadv) nf
  (vzmin,vzmax) = getvzrange()
  istart = top.ins[js]-1
  for ij in xrange(1,istep+1,nskipcol*2):
    for ic in xrange(1,ncolor+1):
      irs1 = istart+ij+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii=compress(not_equal(top.uzp[irs1:irs2:irs3],0.),iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-inp*(ic-1))/istep
      warpplp(take(top.uzp[irs1:irs2:irs3]*top.gaminv[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
         color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
    for ic in xrange(ncolor,0,-1):
      irs1 = istart+ij+nskipcol+inp*(ic-1)
      irs2 = istart+inp*ic
      irs3 = istep
      ii = compress(not_equal(top.uzp[irs1:irs2:irs3],0.),
                    iota(irs1,irs2,irs3))
      ii = (ii-istart-ij-nskipcol-inp*(ic-1))/istep
      warpplp(take(top.uzp[irs1:irs2:irs3]*top.gaminv[irs1:irs2:irs3],ii),
              take(top.zp[irs1:irs2:irs3],ii),
         color=color[ic%len(color)],linetype="none",marker=marker,msize=msize)
  if titles: ptitles("Vz vs Z","Z","Vz"," ")
  if (lframe): limits(top.zplmin+top.zbeam,top.zplmax+top.zbeam,vzmin,vzmax)

##########################################################################
def ppco(y,x,z,uz=1.,xmin=None,xmax=None,ymin=None,ymax=None,
         zmin=None,zmax=None,ncolor=None,usepalette=1,
         marker="\1",msize=1.0,lframe=0):
  """Plots y versus x with color based in z
     - y is y coordinate
     - x is x coordinate
     - z is used to calculate the color
     - xmin, xmax, ymin, ymax, zmin, zmax optional bounds
     - ncolor is number of colors to use, defaults to top.ncolor
     - usepalette=1 when true, uses palette, otherwise uses colors in array
                    color
     - msize=1.0 marker size
     - lframe=0 specifies whether or not to set plot limits
     - titles=1 specifies whether or not to plot titles
  """

  # --- If there are not particles to plot, just return
  np = globalsum(len(x))
  if np == 0: return

  # --- Make sure the lengths of the input are the same
  assert (len(y) == len(x) == len(z)),"x, y, and z must all be the same length"

  # --- This routine can be expensive in parallel when there are many
  # --- colors since synchronization is needed for each color.
  # --- So, if there arn't too many particles, transfer everything to PE0
  # --- and let it do the work.
  if np < 1000000:
    alllocal = 1
    y = gatherarray(y)
    x = gatherarray(x)
    z = gatherarray(z)
    if type(uz) == ArrayType: uz = gatherarray(uz)
  else:
    alllocal = 0

  # --- Make sure arrays are 1-D
  rx = ravel(x)
  ry = ravel(y)
  rz = ravel(z)

  # --- Find extrema
  if xmin is None: xmin = globalmin(rx)
  if xmax is None: xmax = globalmax(rx)
  if ymin is None: ymin = globalmin(ry)
  if ymax is None: ymax = globalmax(ry)
  if zmin is None: zmin = globalmin(rz)
  if zmax is None: zmax = globalmax(rz)

  if ncolor is None: ncolor = top.ncolor
  dd = (zmax - zmin)/ncolor
  #if (lframadv) nf
  for ic in xrange(1,ncolor+1):
    ii = compress(logical_and(logical_and(not_equal(uz,0.),
           less(zmin+(ic-1)*dd,rz)),less(rz,zmin+ic*dd)), iota(0,len(rx)))
    if usepalette:
      c = int(240*ic/ncolor)
    else:
      c = color[ic%len(color)]
    if alllocal:
      plp(take(y,ii),take(x,ii),color=c,marker=marker,msize=msize)
    else:
      warpplp(take(y,ii),take(x,ii),
              color=c,linetype="none",marker=marker,msize=msize)
  if (lframe): limits(xmin,xmax,ymin,ymax)


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
def pla(y,x=None,color="fg",linetype="solid",marks=0,marker=None,msize=1.0,
        width=1.,decomposed=0):
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
    if decomposed:
      warpplg(yy[:,i],xx[:,i%n],color=color,linetype=linetype,
              marks=marks,marker=marker,msize=msize,width=width)
    else:
      if len(yy[:,i]) > 0:
        plg(yy[:,i],xx[:,i%n],color=color,type=linetype,
            marks=marks,marker=marker,msize=msize,width=width)

##########################################################################
##########################################################################
# history plotting routines have been replaced by those in histplots.py
##########################################################################
##########################################################################

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
  - bcast=0: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  """
  if not lparallel:
    if ix is None     and iy is None     and iz is None    :
      return w3d.rho
    if ix is not None and iy is None     and iz is None    :
      return w3d.rho[ix,:,:]
    if ix is None     and iy is not None and iz is None    :
      return w3d.rho[:,iy,:]
    if ix is None     and iy is None     and iz is not None:
      return w3d.rho[:,:,iz]
    if ix is not None and iy is not None and iz is None    :
      return w3d.rho[ix,iy,:]
    if ix is not None and iy is None     and iz is not None:
      return w3d.rho[ix,:,iz]
    if ix is None     and iy is not None and iz is not None:
      return w3d.rho[:,iy,iz]
    if ix is not None and iy is not None and iz is not None:
      return w3d.rho[ix,iy,iz]
  else:
    iz1 = top.izfsslave[me] - top.izslave[me]
    if me < npes-1:
      iz2 = top.izfsslave[me+1] - top.izslave[me]
    else:
      iz2 = iz1 + top.nzfsslave[me] + 1
    ppp = w3d.rho[:,:,iz1:iz2]
    if ix is not None and iy is None:
      ppp = ppp[ix,:,:]
    elif ix is None and iy is not None:
      ppp = ppp[:,iy,:]
    elif ix is not None and iy is not None:
      ppp = ppp[ix,iy,:]
    if iz is None:
      ppp = transpose(gatherarray(transpose(ppp)))
    else:
      pe = convertizfstope(iz)
      if pe is None: return None
      if me == pe: ppp = ppp[...,iz-top.izfsslave[me]]
      else:        ppp = zeros(shape(ppp[...,0]),'d')
      if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
    if bcast: ppp = broadcast(ppp)
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
    if ix is None     and iy is None     and iz is None    :
      w3d.rho[:,:,:] = val
    if ix is not None and iy is None     and iz is None    :
      w3d.rho[ix,:,:] = val
    if ix is None     and iy is not None and iz is None    :
      w3d.rho[:,iy,:] = val
    if ix is None     and iy is None     and iz is not None:
      w3d.rho[:,:,iz] = val
    if ix is not None and iy is not None and iz is None    :
      w3d.rho[ix,iy,:] = val
    if ix is not None and iy is None     and iz is not None:
      w3d.rho[ix,:,iz] = val
    if ix is None     and iy is not None and iz is not None:
      w3d.rho[:,iy,iz] = val
    if ix is not None and iy is not None and iz is not None:
      w3d.rho[ix,iy,iz] = val
  else:
    print "Warning, setrho this is not yet implemented in parallel"
   #if me < npes-1:
   #  ppp = w3d.rho[:,:,:-top.grid_overlap]
   #else:
   #  ppp = w3d.rho[:,:,:]
   #if ix is not None and iy is None    :
   #  ppp = ppp[ix,:,:]
   #elif ix is None and iy is not None:
   #  ppp = ppp[:,iy,:]
   #elif ix is not None and iy is not None:
   #  ppp = ppp[ix,iy,:]
   #if iz is None:
   #  ppp = transpose(gatherarray(transpose(ppp)))
   #else:
   #  pe = convertiztope(iz)
   #  if pe is None: return None
   #  if me == pe: ppp = ppp[...,iz-top.izslave[me+1]+1]
   #  if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
   #if bcast: ppp = broadcast(ppp)
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
  - bcast=0: When 1, the result is broadcast to all of the processors
             (otherwise returns None to all but PE0
  """
  if not lparallel:
    if ix is None     and iy is None     and iz is None    :
      return w3d.phi[:,:,1:-1]
    if ix is not None and iy is None     and iz is None    :
      return w3d.phi[ix,:,1:-1]
    if ix is None     and iy is not None and iz is None    :
      return w3d.phi[:,iy,1:-1]
    if ix is None     and iy is None     and iz is not None:
      return w3d.phi[:,:,iz+1]
    if ix is not None and iy is not None and iz is None    :
      return w3d.phi[ix,iy,1:-1]
    if ix is not None and iy is None     and iz is not None:
      return w3d.phi[ix,:,iz+1]
    if ix is None     and iy is not None and iz is not None:
      return w3d.phi[:,iy,iz+1]
    if ix is not None and iy is not None and iz is not None:
      return w3d.phi[ix,iy,iz+1]
  else:
    iz1 = top.izfsslave[me] - top.izslave[me]
    if me < npes-1:
      iz2 = top.izfsslave[me+1] - top.izslave[me]
    else:
      iz2 = iz1 + top.nzfsslave[me] + 1
    ppp = w3d.phi[:,:,iz1+1:iz2+1]
    if ix is not None and iy is None:
      ppp = ppp[ix,:,:]
    elif ix is None and iy is not None:
      ppp = ppp[:,iy,:]
    elif ix is not None and iy is not None:
      ppp = ppp[ix,iy,:]
    if iz is None:
      ppp = transpose(gatherarray(transpose(ppp)))
    else:
      pe = convertizfstope(iz)
      if pe is None: return None
      if me == pe: ppp = ppp[...,iz-top.izfsslave[me]]
      else:        ppp = zeros(shape(ppp[...,0]),'d')
      if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
    if bcast: ppp = broadcast(ppp)
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
    if ix is None     and iy is None     and iz is None    :
      w3d.phi[:,:,1:-1] = val
    if ix is not None and iy is None     and iz is None    :
      w3d.phi[ix,:,1:-1] = val
    if ix is None     and iy is not None and iz is None    :
      w3d.phi[:,iy,1:-1] = val
    if ix is None     and iy is None     and iz is not None:
      w3d.phi[:,:,iz+1] = val
    if ix is not None and iy is not None and iz is None    :
      w3d.phi[ix,iy,1:-1] = val
    if ix is not None and iy is None     and iz is not None:
      w3d.phi[ix,:,iz+1] = val
    if ix is None     and iy is not None and iz is not None:
      w3d.phi[:,iy,iz+1] = val
    if ix is not None and iy is not None and iz is not None:
      w3d.phi[ix,iy,iz+1] = val
  else:
    print "Warning, setphi this is not yet implemented in parallel"
   #if me < npes-1:
   #  ppp = w3d.phi[:,:,1:w3d.nz-top.grid_overlap+2]
   #else:
   #  ppp = w3d.phi[:,:,1:-1]
   #if ix is not None and iy is None:
   #  ppp = ppp[ix,:,:]
   #elif ix is None and iy is not None:
   #  ppp = ppp[:,iy,:]
   #elif ix is not None and iy is not None:
   #  ppp = ppp[ix,iy,:]
   #if iz is None:
   #  ppp = transpose(gatherarray(transpose(ppp)))
   #else:
   #  pe = convertiztope(iz)
   #  if pe is None: return None
   #  if me == pe: ppp = ppp[...,iz-top.izslave[me+1]+1]
   #  if (me == pe or me == 0) and (pe != 0): ppp = getarray(pe,ppp,0)
   #if bcast: ppp = broadcast(ppp)
   #return ppp
##########################################################################
##########################################################################
def pcrhozy(ix=None,**kw):
  """Plots contours of charge density in the Z-Y plane
     - ix=w3d.ix_axis X index of plane"""
  if ix is None: ix = w3d.ix_axis
  if not kw.has_key('xmin'): kw['xmin'] = w3d.zmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.zmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if not kw.has_key('contours'): kw['contours'] = 20
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.zmmin,w3d.zmmax,w3d.ymmin,w3d.ymmax)
  settitles("Charge density in z-y plane","Z","Y","ix = "+repr(ix))
  rrr = getrho(ix=ix)
  if rrr is None: rrr = zeros((w3d.ny,w3d.nzfull),'d')
  ppgeneric(grid=transpose(rrr),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcrhozy.__doc__ = pcrhozy.__doc__ + ppgeneric_doc("z","y")
##########################################################################
def pcrhozx(iy=None,**kw):
  """Plots contours of charge density in the Z-X plane
     - iy=w3d.iy_axis Y index of plane"""
  if iy is None: iy = w3d.iy_axis
  if not kw.has_key('xmin'): kw['xmin'] = w3d.zmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.zmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.xmmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.xmmax
  if not kw.has_key('contours'): kw['contours'] = 20
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.zmmin,w3d.zmmax,w3d.xmmin,w3d.xmmax)
  settitles("Charge density in z-x plane","Z","X","iy = "+repr(iy))
  rrr = getrho(iy=iy)
  if rrr is None: rrr = zeros((w3d.nx,w3d.nzfull),'d')
  ppgeneric(grid=transpose(rrr),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcrhozx.__doc__ = pcrhozx.__doc__ + ppgeneric_doc("z","x")
##########################################################################
def pcrhoxy(iz=None,**kw):
  """Plots contours of charge density in the X-Y plane
     - iz=w3d.iz_axis Z index of plane"""
  if iz is None: iz = w3d.iz_axis + top.izslave[me]
  if not kw.has_key('xmin'): kw['xmin'] = w3d.xmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.xmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if not kw.has_key('contours'): kw['contours'] = 20
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax)
  settitles("Charge density in x-y plane","X","Y","iz = "+repr(iz))
  rrr = getrho(iz=iz)
  if rrr is None: rrr = zeros((w3d.nx,w3d.ny),'d')
  ppgeneric(grid=rrr,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcrhoxy.__doc__ = pcrhoxy.__doc__ + ppgeneric_doc("x","y")
##########################################################################
def pcphizy(ix=None,**kw):
  """Plots contours of electrostatic potential in the Z-Y plane
     - ix=w3d.ix_axis X index of plane"""
  if ix is None: ix = w3d.ix_axis
  if not kw.has_key('xmin'): kw['xmin'] = w3d.zmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.zmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if not kw.has_key('contours'): kw['contours'] = 20
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.zmmin,w3d.zmmax,w3d.ymmin,w3d.ymmax)
  settitles("Charge density in z-y plane","Z","Y","ix = "+repr(ix))
  ppp = getphi(ix=ix)
  if ppp is None: ppp = zeros((w3d.ny,w3d.nzfull),'d')
  ppgeneric(grid=transpose(ppp),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcphizy.__doc__ = pcphizy.__doc__ + ppgeneric_doc("z","y")
##########################################################################
def pcphizx(iy=None,**kw):
  """Plots contours of electrostatic potential in the Z-X plane
     - iy=w3d.iy_axis Y index of plane"""
  if iy is None: iy = w3d.iy_axis
  if not kw.has_key('xmin'): kw['xmin'] = w3d.zmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.zmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.xmmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.xmmax
  if not kw.has_key('contours'): kw['contours'] = 20
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.zmmin,w3d.zmmax,w3d.xmmin,w3d.xmmax)
  settitles("Charge density in z-x plane","Z","X","iy = "+repr(iy))
  ppp = getphi(iy=iy)
  if ppp is None: ppp = zeros((w3d.nx,w3d.nzfull),'d')
  ppgeneric(grid=transpose(ppp),kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcphizx.__doc__ = pcphizx.__doc__ + ppgeneric_doc("z","x")
##########################################################################
def pcphixy(iz=None,**kw):
  """Plots contours of electrostatic potential in the X-Y plane
     - iz=w3d.iz_axis Z index of plane"""
  if iz is None: iz = w3d.iz_axis + top.izslave[me]
  if not kw.has_key('xmin'): kw['xmin'] = w3d.xmmin
  if not kw.has_key('xmax'): kw['xmax'] = w3d.xmmax
  if not kw.has_key('ymin'): kw['ymin'] = w3d.ymmin
  if not kw.has_key('ymax'): kw['ymax'] = w3d.ymmax
  if not kw.has_key('contours'): kw['contours'] = 20
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (w3d.zmmin,w3d.zmmax,w3d.xmmin,w3d.xmmax)
  settitles("Charge density in x-y plane","X","Y","iz = "+repr(iz))
  ppp = getphi(iz=iz)
  if ppp is None: ppp = zeros((w3d.nx,w3d.ny),'d')
  ppgeneric(grid=ppp,kwdict=kw)
if sys.version[:5] != "1.5.1":
  pcphixy.__doc__ = pcphixy.__doc__ + ppgeneric_doc("x","y")
##########################################################################
##########################################################################
def ppdecomposition(scale=1.,minscale=0.,gap=0.2):
  """Shows the domain decomposition in a graphical way. For each
processor, the total mesh extent is plotted as a filled rectangle
covering the z-length and with height determined by 'scale' and the
number of processors. Another filled rectangle is plotted in the top
half showing the particle domains, and one on the lower half shows the
field domain.
  - scale=1.: the maximum vertical extent of the graph
  - minscale=0.: the minimum vertical extent of the graph
  - gap=0.2: fractional vertical gap between rectangles
  """
  z = []
  x = []
  y = []
  dd = 1.*scale/top.maxslaves
  mm = 1. - gap
  for i in xrange(top.maxslaves):
    z = z + [1.]
    zmin = top.zmslmin[i]
    zmax = top.zmslmax[i]
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([-mm,-mm,mm,mm,-mm]))
  for i in xrange(top.maxslaves):
    z = z + [2.]
    zmin = top.zpslmin[i]
    zmax = top.zpslmax[i]
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([0,0,mm,mm,0]))
  for i in xrange(top.maxslaves):
    z = z + [3.]
    zmin = top.izfsslave[i]*w3d.dz
    zmax = top.izfsslave[i]*w3d.dz + top.nzfsslave[i]*w3d.dz
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([-mm,-mm,0,0,-mm]))
  plfp(array(z),y,x,5*ones(len(z)),cmin=0,cmax=4)
  for i in xrange(len(z)):
    pldj(x[i*5:i*5+4],y[i*5:i*5+4],x[i*5+1:i*5+5],y[i*5+1:i*5+5])
      

##########################################################################
##########################################################################
plps = []
plfreq = []
plseldom = []
plalways = []
#--------------------------------------------------------------------------
def installplseldom(f):
  "Adds a function to the list of functions called when seldom plots are made"
  plseldom.append(f)
#--------------------------------------------------------------------------
def uninstallplseldom(f):
  """Removes the function from the list of functions called when seldom plots
     are made"""
  if f in plseldom:
    plseldom.remove(f)
  else:
    raise 'Warning: uninstallplseldom: no such function had been installed'
#--------------------------------------------------------------------------
def installplalways(f):
  "Adds a function to the list of functions called when always plots are made"
  plalways.append(f)
#--------------------------------------------------------------------------
def uninstallplalways(f):
  """Removes the function from the list of functions called when always plots
     are made"""
  if f in plalways:
    plalways.remove(f)
  else:
    raise 'Warning: uninstallplalways: no such function had been installed'

##########################################################################
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

##########################################################################
def onedplts(freqflag=always):
  """Makes 1-D plots which have been turned on
     - freqflag=always frequency flag, either always, seldom, or never"""
  currentwindow = current_window()
  window(0)
  if freqflag == top.ipcurr: pzcurr()
  if freqflag == top.ipegap: pzegap()
  if freqflag == top.iplchg: pzlchg()
  if freqflag == top.ipvzofz: pzvzofz()
  if freqflag == top.iprhoax: pzrhoax()
  if freqflag == top.ipphiax: pzphiax()
  if freqflag == top.ipezax: pzezax()
  oldlimits = limits()
  window(currentwindow)

##########################################################################
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
    for p in plalways:
      p()
      fma()
      oldlimits = limits()
    for p in plfreq:
      p()
      fma()
      oldlimits = limits()
  if freqflag == seldom:
    for p in plseldom:
      p()
      fma()
      oldlimits = limits()
    for p in plps:
      p()
      fma()
      oldlimits = limits()

# --- Reset the current window to it previous value.
  oldlimits = limits()
  window(currentwindow)
