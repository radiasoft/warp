from warp import *
from appendablearray import *
import cPickle
extpart_version = "$Id: extpart.py,v 1.5 2001/07/19 21:41:58 dave Exp $"

def extpartdoc():
  print """
Creates a class for handling extrapolated particle windows, ExtPart. Type
doc(ExtPart) for more help.
Two functions are available for saving the object in a file.
 - dumpExtPart(object,filename)
 - restoreExtPart(object,filename)
  """

############################################################################
class ExtPart:
  """This class defines a container to setup and keep track of extropolated
particle data. It can optionally accumulate the data over multiple time steps.
The creator options are:
 - iz: grid location where the extrapolated data is saved.
 - nepmax: max size of the arrays. Defaults to 3*top.pnumz[iz] if non-zero,
           otherwise 10000.
 - laccumulate=0: when true, particles are accumulated over multiple steps.

Available methods:
 - setaccumulate(v=1): Turns accumulation on or off. If turned on, all
                       currently saved data is deleted.
 - clear(): Will clear the existing data.
 - disable(): Turns off collecting of data
 - enable(): Turns on collecting of data (only needed after disable)

The follow all take an optional argument to specify species number.
 - nep: Get number of particles
 - t: Get time at which particle was saved
 - x, y, ux, uy, uz, vx, vy, vz, xp, yp, r, theta, rp: Get the various
     coordinates or velocity of particles

The following are available plot routines. All take an optional argument to
specify species number. Additional arguments are the same as the 'pp' plotting
routines (such as ppxxp).
 - pxy, pxxp, pyyp, pxpyp, prrp, ptx, pty, ptxp, ptyp, ptux, ptuy, ptuz, ptvx
 - ptvy, ptvz, ptrace
  """

  def __init__(self,iz,nepmax=None,laccumulate=0):
    # --- Save input values, getting default values when needed
    self.iz = iz
    self.izlocal = iz - top.izslave[me]
    self.laccumulate = laccumulate
    if nepmax is None:
      self.nepmax = 10000
      if top.allocated("pnumz") and 0 <= self.izlocal <= top.nzmmnt:
        if top.pnumz[self.izlocal] > 0: self.nepmax = top.pnumz[self.izlocal]*3
    else:
      self.nepmax = nepmax
    # --- Add this new window to the ExtPart group in top
    self.enabled = 0
    self.enable()
    # --- Setup empty arrays for accumulation if laccumulate if true.
    # --- Otherwise, the arrays will just point to the data in ExtPart.
    self.setuparrays()

  def setuparrays(self):
    if self.laccumulate:
      self.tep = []
      self.xep = []
      self.yep = []
      self.uxep = []
      self.uyep = []
      self.uzep = []
      for js in xrange(top.ns):
        self.tep.append(AppendableArray(self.nepmax,'d',self.nepmax))
        self.xep.append(AppendableArray(self.nepmax,'d',self.nepmax))
        self.yep.append(AppendableArray(self.nepmax,'d',self.nepmax))
        self.uxep.append(AppendableArray(self.nepmax,'d',self.nepmax))
        self.uyep.append(AppendableArray(self.nepmax,'d',self.nepmax))
        self.uzep.append(AppendableArray(self.nepmax,'d',self.nepmax))
    else:
      self.tep = top.ns*[None]
      self.xep = top.ns*[None]
      self.yep = top.ns*[None]
      self.uxep = top.ns*[None]
      self.uyep = top.ns*[None]
      self.uzep = top.ns*[None]

  def clear(self):
    self.setuparrays()

  def getid(self):
    assert self.enabled,"This window is disabled and there is no associated id"
    for i in xrange(top.nepwin):
      if top.izepwin[i] == self.izlocal: return i
    raise "Uh Ooh! Somehow the window was deleted! I can't continue!"

  def enable(self):
    # --- Add this window to the list
    # --- Only add this iz to the list if it is not already there.
    # --- Note that it is not an error to have more than one instance
    # --- have the same value of iz. For example one could be accumulating
    # --- while another isn't.
    if self.enabled: return
    top.nepwin = top.nepwin + 1
    if top.nepmax < self.nepmax: top.nepmax = self.nepmax
    err = gchange("ExtPart")
    top.izepwin[-1] = self.izlocal
    # --- Set so accumulate method is called after time steps
    installafterstep(self.accumulate)
    self.enabled = 1

  def disable(self):
    if not self.enabled: return
    # --- Set so accumulate method is not called after time steps
    uninstallafterstep(self.accumulate)
    # --- Remove this window from the list.
    top.nepwin = top.nepwin - 1
    for i in xrange(self.getid(),top.nepwin):
      top.izepwin[i] = top.izepwin[i+1]
    gchange("ExtPart")
    self.enabled = 0

  def accumulate(self):
    if globalmax(maxnd(top.nep) == top.nepmax):
      print "************* WARNING *************"
      print "**** Not enough space was allocated for the ExtPart arrays."
      print "**** The data will not be correct and will not be saved."
      print "**** A guess will be made as to how much to increase the size"
      print "**** of the arrays. Please run another timestep to accumulate new"
      print "**** data"
      if top.allocated("pnumz") and 0 <= self.izlocal <= top.nzmmnt:
        guess = 3*top.pnumz[self.izlocal]
      else:
        guess = 0
      # --- Only do this on if the relation is true. This avoids unnecessarily
      # --- increasing the size of the arrays on processors where no data
      # --- is gathered.
      if maxnd(top.nep) == top.nepmax:
        top.nepmax = max(2*top.nepmax,guess)
        err = gchange("ExtPart")
      return
    for js in xrange(top.ns):
      id = self.getid()
      # --- Gather the data.
      # --- In parallel, the data is gathered in PE0, return empty arrays
      # --- on other processors. In serial, the arrays are just returned as is.
      nn = top.nep[id,js]
      t = gatherarray(top.tep[:nn,id,js]+0.,othersempty=1)
      x = gatherarray(top.xep[:nn,id,js]+0.,othersempty=1)
      y = gatherarray(top.yep[:nn,id,js]+0.,othersempty=1)
      ux = gatherarray(top.uxep[:nn,id,js]+0.,othersempty=1)
      uy = gatherarray(top.uyep[:nn,id,js]+0.,othersempty=1)
      uz = gatherarray(top.uzep[:nn,id,js]+0.,othersempty=1)
      if self.laccumulate:
        self.tep[js].append(t)
        self.xep[js].append(x)
        self.yep[js].append(y)
        self.uxep[js].append(ux)
        self.uyep[js].append(uy)
        self.uzep[js].append(uz)
      else:
        self.tep[js] = t
        self.xep[js] = x
        self.yep[js] = y
        self.uxep[js] = ux
        self.uyep[js] = uy
        self.uzep[js] = uz

  def setaccumulate(self,v=1):
    self.laccumulate = v
    if self.laccumulate: self.setuparrays()

  def t(self,js=0): return self.tep[js][:]
  def x(self,js=0): return self.xep[js][:]
  def y(self,js=0): return self.yep[js][:]
  def ux(self,js=0): return self.uxep[js][:]
  def uy(self,js=0): return self.uyep[js][:]
  def uz(self,js=0): return self.uzep[js][:]
  def vx(self,js=0): return self.uxep[js][:]
  def vy(self,js=0): return self.uyep[js][:]
  def vz(self,js=0): return self.uzep[js][:]
  def xp(self,js=0): return self.ux(js)/self.uz(js)
  def yp(self,js=0): return self.uy(js)/self.uz(js)
  def r(self,js=0): return sqrt(self.x(js)**2 + self.y(js)**2)
  def theta(self,js=0): return arctan2(self.y(js),self.x(js))
  def rp(self,js=0): return self.xp(js)*cos(self.theta(js)) + \
                            self.yp(js)*sin(self.theta(js))
  def nep(self,js=0): return len(self.tep[js][:])

  ############################################################################
  ############################################################################
  # --- Define plotting routines for the extrapolated particles.

  def checkplotargs(self,kw):
    """Convenience routine to check arguments of particle plot routines.
Warning: this has the side affect of adding the arguement allowbadargs to
the kw dictionary. This is done since the calls to these functions here to
make the plots may have unused arguements since the entire kw list passed
into each of the pp plotting routines is passed into each of these
functions.
    """
    badargs = ppgeneric(checkargs=1,kwdict=kw)
    kw['allowbadargs'] = 1
    if badargs: raise "bad arguments ",string.join(badargs.keys())

  def titleright(self):
    if lparallel:
      return "iz = %d (z = %f m)"%(self.iz,top.zmslmin[0]+self.iz*w3d.dz)
    else:
      return "iz = %d (z = %f m)"%(self.iz,w3d.zmesh[self.iz])

  ############################################################################
  def pxy(self,js=0,particles=1,**kw):
    """Plots X-Y for extraploated particles"""
    self.checkplotargs(kw)
    x = self.x(js)
    y = self.y(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (top.xplmin,top.xplmax,top.yplmin,top.yplmax)
    settitles("Y vs X","X","Y",self.titleright())
    ppgeneric(y,x,kwdict=kw)

  ############################################################################
  def pxxp(self,js=0,slope=0.,offset=0.,particles=1,**kw):
    """Plots X-X' for extraploated particles"""
    self.checkplotargs(kw)
    x = self.x(js)
    xp = self.xp(js)
    if type(slope) == type(''):
      slope = (ave(x*xp)-ave(x)*ave(xp))/(ave(x*x) - ave(x)**2)
      offset = ave(xp)-slope*ave(x)
    kw['slope'] = slope
    kw['offset'] = offset
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)
    settitles("X' vs X","X","X'",self.titleright())
    ppgeneric(xp,x,kwdict=kw)

  ############################################################################
  def pyyp(self,js=0,slope=0.,offset=0.,particles=1,**kw):
    """Plots Y-Y' for extraploated particles"""
    self.checkplotargs(kw)
    y = self.y(js)
    yp = self.yp(js)
    if type(slope) == type(''):
      slope = (ave(y*yp)-ave(y)*ave(yp))/(ave(y*y) - ave(y)**2)
      offset = ave(yp)-slope*ave(y)
    kw['slope'] = slope
    kw['offset'] = offset
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax)
    settitles("Y' vs Y","Y","Y'",self.titleright())
    ppgeneric(yp,y,kwdict=kw)

  ############################################################################
  def pxpyp(self,js=0,particles=1,**kw):
    """Plots X'-Y' for extraploated particles"""
    self.checkplotargs(kw)
    xp = self.xp(js)
    yp = self.yp(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)
    settitles("Y' vs X'","X'","Y'",self.titleright())
    ppgeneric(yp,xp,kwdict=kw)

  ############################################################################
  def prrp(self,js=0,scale=0.,slope=0.,offset=0.,particles=1,**kw):
    """Plots R-R' for extraploated particles"""
    self.checkplotargs(kw)
    x = self.x(js)
    y = self.y(js)
    xp = self.xp(js)
    yp = self.yp(js)
    xscale = 1.
    yscale = 1.
    xpscale = 1.
    ypscale = 1.
    if scale:
      xscale = 2.*sqrt(ave(x*x) - ave(x)**2)
      yscale = 2.*sqrt(ave(y*y) - ave(y)**2)
      xpscale = 2.*sqrt(ave(xp*xp) - ave(xp)**2)
      ypscale = 2.*sqrt(ave(yp*yp) - ave(yp)**2)
    x = x/xscale
    y = y/yscale
    xp = xp/xpscale
    yp = yp/ypscale
    r = sqrt(x**2 + y**2)
    t = arctan2(y,x)
    rp = xp*cos(t) + yp*sin(t)
    if type(slope) == type(''): slope = ave(r*rp)/ave(r*r)
    kw['slope'] = slope
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                        top.xpplmin/xpscale,top.xpplmax/ypscale)
    settitles("R' vs R","R","R'",self.titleright())
    ppgeneric(rp,r,kwdict=kw)

  ############################################################################
  def ptx(self,js=0,particles=1,**kw):
    """Plots time-X for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    x = self.x(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = ('e','e',top.xplmin,top.xplmax)
    settitles("X vs time","time","X",self.titleright())
    ppgeneric(x,t,kwdict=kw)

  ############################################################################
  def pty(self,js=0,particles=1,**kw):
    """Plots time-Y for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    y = self.y(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = ('e','e',top.yplmin,top.yplmax)
    settitles("Y vs time","time","Y",self.titleright())
    ppgeneric(y,t,kwdict=kw)

  ############################################################################
  def ptxp(self,js=0,particles=1,**kw):
    """Plots time-X' for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    xp = self.xp(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = ('e','e',top.xpplmin,top.xpplmax)
    settitles("X' vs time","time","X'",self.titleright())
    ppgeneric(xp,t,kwdict=kw)

  ############################################################################
  def ptyp(self,js=0,particles=1,**kw):
    """Plots time-Y' for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    yp = self.yp(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    else:
      kw['pplimits'] = ('e','e',top.ypplmin,top.ypplmax)
    settitles("Y' vs time","time","Y'",self.titleright())
    ppgeneric(yp,t,kwdict=kw)

  ############################################################################
  def ptux(self,js=0,particles=1,**kw):
    """Plots time-ux for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    ux = self.ux(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("ux vs time","time","ux",self.titleright())
    ppgeneric(ux,t,kwdict=kw)

  ############################################################################
  def ptuy(self,js=0,particles=1,**kw):
    """Plots time-uy for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    uy = self.uy(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("uy vs time","time","uy",self.titleright())
    ppgeneric(uy,t,kwdict=kw)

  ############################################################################
  def ptuz(self,js=0,particles=1,**kw):
    """Plots time-uz for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    uz = self.uz(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("uz vs time","time","uz",self.titleright())
    ppgeneric(uz,t,kwdict=kw)

  ############################################################################
  def ptvx(self,js=0,particles=1,**kw):
    """Plots time-Vx for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    vx = self.vx(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("Vx vs time","time","Vx",self.titleright())
    ppgeneric(vx,t,kwdict=kw)

  ############################################################################
  def ptvy(self,js=0,particles=1,**kw):
    """Plots time-Vy for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    vy = self.vy(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("Vy vs time","time","Vy",self.titleright())
    ppgeneric(vy,t,kwdict=kw)

  ############################################################################
  def ptvz(self,js=0,particles=1,**kw):
    """Plots time-Vz for extraploated particles"""
    self.checkplotargs(kw)
    t = self.t(js)
    vz = self.vz(js)
    kw['particles'] = particles
    if 'pplimits' in kw.keys():
      kw['lframe'] = 1
    settitles("Vz vs time","time","Vz",self.titleright())
    ppgeneric(vz,t,kwdict=kw)

  ############################################################################
  def ptrace(self,js=0,slope=0.,particles=1,pplimits=None,**kw):
    """
Plots X-Y, X-X', Y-Y', Y'-X' in single page
If slope='auto', it is calculated from the moments for X-X' and Y-Y' plots.
pplimits can be a list of up to four tuples, one for each phase space plot.
If any of the tuples are empty, the limits used will be the usual ones for
that plot.
    """
    self.checkplotargs(kw)
    x = self.x(js)
    y = self.y(js)
    xp = self.xp(js)
    yp = self.yp(js)
    titler = self.titleright()
    kw['particles'] = particles
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
    if type(slope)==type(''): kw['slope'] = 0.
    settitles("Y vs X","X","Y",titler)
    ppgeneric(y,x,kwdict=kw)
 
    kw['view'] = 4
    kw['pplimits'] = pplimits[1]
    if type(slope)==type(''):
      kw['slope'] = (ave(y*yp)-ave(y)*ave(yp))/(ave(y*y) - ave(y)**2)
    settitles("Y' vs Y","Y","Y'",titler)
    ppgeneric(yp,y,kwdict=kw)

    kw['view'] = 5
    kw['pplimits'] = pplimits[2]
    if type(slope)==type(''):
      kw['slope'] = (ave(x*xp)-ave(x)*ave(xp))/(ave(x*x) - ave(x)**2)
    settitles("X' vs X","X","X'",titler)
    ppgeneric(xp,x,kwdict=kw)
 
    kw['view'] = 6
    kw['pplimits'] = pplimits[3]
    if type(slope)==type(''): kw['slope'] = 0.
    settitles("X' vs Y'","Y'","X'",titler)
    ppgeneric(xp,yp,kwdict=kw)

##############################################################################
def dumpExtPart(object,filename):
  """Dump the saved extrapolated data to a file
 - filename: The name of the file to save the data in"""
  if me == 0:
    # --- Only PE0 writes the object to the file since it is the processor
    # --- where the data is gathered.
    ff = open(filename,'w')
    cPickle.dump(object,ff,1)
    ff.close()

def restoreExtPart(object,filename):
  """Restore extrapolated data from the given file"""
  if me == 0:
    # --- Only PE0 wrote the object to the file since it is the processor
    # --- where the data was gathered.
    ff = open(filename,'r')
    result = cPickle.load(ff)
    ff.close()
    result.enable()
    # --- Get the value of iz
    iz = result.iz
  else:
    # --- Create temp iz
    iz = 0
  # --- PE0 broadcasts its value of iz to all of the other processors
  # --- whcih create new instances of the ExtPart class.
  iz = broadcast(iz)
  if me > 0: result = ExtPart(iz)
  return result






