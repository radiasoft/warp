from warp import *
from appendablearray import *
import cPickle
extpart_version = "$Id: extpart.py,v 1.2 2001/06/26 00:11:34 dave Exp $"

def extpartdoc():
  print """
Sets up windows in which extrapolated particle data is saved. The data is
extrapolated during the moments calculation. The following routines are
available:
 - addextpart: Add a window to save extrapolated particles in
 - deleteextpart: Delete an extrapolated particle window
 - disableextpart: Disables an extrapolated particle window
 - enableextpart: Enables an extrapolated particle window
 - setaccumulateextpart: Turn accumulation in a window on or off
 - clearextpart: Clear the data saved in a window
 - getnep: Get the number of particles saved in a window
 - dumpextpart: Dump the saved extrapolated data to a file
 - restoreextpart: Restore extrapolated data from the given file
Routines which return the particle coordinates, all prefixed with 'extpart':
 - extpartt, x, y, xp, yp, vx, vy, vz, r, theta, rp
As well as plots of any two coordinates, all prefixed with 'ep'. For example
 - epxy, eprrp, ...
  """

############################################################################
class _ExtPart:
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
 - setuparrays(): Will clear the existing data.
  """
  def __init__(self,iz,nepmax=None,laccumulate=0):
    # --- Save input values, getting default values when needed
    self.iz = iz
    self.laccumulate = laccumulate
    if nepmax is None:
      self.nepmax = 10000
      if top.allocated("pnumz"):
        if top.pnumz[iz] > 0: self.nepmax = top.pnumz[iz]*3
    else:
      self.nepmax = nepmax
    # --- Add this new window to the ExtPart group in top
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

  def getid(self):
    assert self.enabled,"This window is disabled and there is no associated id"
    for i in xrange(top.nepwin):
      if top.izepwin[i] == self.iz: return i
    raise "Oh Ooh! Somehow the window was deleted! I can't continue!"

  def enable(self):
    # --- Add this window to the list
    # --- Only add this iz to the list if it is not already there.
    # --- Note that it is not an error to have more than one instance
    # --- have the same value of iz. For example one could be accumulating
    # --- while another isn't.
    if top.nepwin == 0 or self.iz not in top.izepwin:
      top.nepwin = top.nepwin + 1
      if top.nepmax < self.nepmax: top.nepmax = self.nepmax
      err = gchange("ExtPart")
      top.izepwin[-1] = self.iz
    # --- Set so accumulate method is called after time steps
    installafterstep(self.accumulate)
    self.enabled = 1

  def disable(self):
    self.enabled = 0
    # --- Set so accumulate method is not called after time steps
    uninstallafterstep(self.accumulate)
    # --- Remove this window from the list.
    top.nepwin = top.nepwin - 1
    for i in xrange(self.getid(),top.nepwin):
      top.izepwin[i] = top.izepwin[i+1]
    gchange("ExtPart")

  def accumulate(self):
    if maxnd(top.nep) == top.nepmax:
      print "************* WARNING *************"
      print "**** Not enough space was allocated for the ExtPart arrays."
      print "**** The data will not be correct and will not be saved."
      print "**** A guess will be made as to how much to increase the size"
      print "**** of the arrays. Please run another timestep to accumulate new"
      print "**** data"
      if top.allocated('pnumz'): guess = 3*top.pnumz[self.iz]
      else:                       guess = 0
      top.nepmax = max(2*top.nepmax,guess)
      err = gchange("ExtPart")
      return
    id = self.getid()
    if self.laccumulate:
      for js in xrange(top.ns):
        nn = top.nep[id,js]
        if nn > 0:
          self.tep[js].append(top.tep[:nn,id,js])
          self.xep[js].append(top.xep[:nn,id,js])
          self.yep[js].append(top.yep[:nn,id,js])
          self.uxep[js].append(top.uxep[:nn,id,js])
          self.uyep[js].append(top.uyep[:nn,id,js])
          self.uzep[js].append(top.uzep[:nn,id,js])
    else:
      for js in xrange(top.ns):
        nn = top.nep[id,js]
        if nn > 0:
          self.tep[js] = top.tep[:nn,id,js] + 0.
          self.xep[js] = top.xep[:nn,id,js] + 0.
          self.yep[js] = top.yep[:nn,id,js] + 0.
          self.uxep[js] = top.uxep[:nn,id,js] + 0.
          self.uyep[js] = top.uyep[:nn,id,js] + 0.
          self.uzep[js] = top.uzep[:nn,id,js] + 0.

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
  def rp(self,js=0):
    return self.xp(js)*cos(self.theta(js)) + self.yp(js)*sin(self.theta(js))
  def nep(self,js=0): return len(self.tep[js][:])


##############################################################################
##############################################################################
##############################################################################
# --- Create a list of instances of ExtPart.
_extpartlist = []

def _convertidtoiz(id):
  assert (id >= 0) and (id < len(_extpartlist)), \
    "id must be within the range 0 to less then the number of windows"
  return _extpartlist[id].iz

def _convertiztoid(iz):
  for i in xrange(len(_extpartlist)):
    if _extpartlist[i].iz == iz:
      return i
  raise "No windows found at the value of iz = %d"%iz

##############################################################################
##############################################################################
##############################################################################
# --- Define routines which are callable outside the module.

def addextpart(iz,nepmax=None,laccumulate=0):
  """Add a window to save extrapolated particles.
 - iz: grid location where the extrapolated data is saved.
 - nepmax: max size of the arrays. Defaults to 3*top.pnumz[iz] if non-zero,
   otherwise 10000.
 - laccumulate=0: when true, particles are accumulated over multiple steps.
  """
  _extpartlist.append(_ExtPart(iz,nepmax,laccumulate))

def deleteextpart(iz=None,id=None):
  """Delete an extrapolated particle window. One of iz or id must be
specified."""
  assert (iz is not None) or (id is not None),"Either iz or id must be input"
  if iz is None: iz = _convertidtoiz(id) # Check validity of input id
  if id is None: id = _convertiztoid(iz)
  _extpartlist[id].disable()
  del _extpartlist[id]

def disableextpart(iz=None,id=None):
  """Disables an extrapolated particle window. No more data will be collected,
but the current data is saved. One of iz or id must be specified."""
  assert (iz is not None) or (id is not None),"Either iz or id must be input"
  if iz is None: iz = _convertidtoiz(id) # Check validity of input id
  if id is None: id = _convertiztoid(iz)
  _extpartlist[id].disable()

def enableextpart(iz=None,id=None):
  """Enables an extrapolated particle window. New data will now be collected.
One of iz or id must be specified."""
  assert (iz is not None) or (id is not None),"Either iz or id must be input"
  if iz is None: iz = _convertidtoiz(id) # Check validity of input id
  if id is None: id = _convertiztoid(iz)
  _extpartlist[id].enable()

def setaccumulateextpart(iz=None,id=0,v=1):
  """Turn accumulation in a window on or off"""
  if iz is not None: id = _convertiztoid(iz)
  _extpartlist[id].setaccumulate(v)

def clearextpart(iz=None,id=0):
  """Clear the data saved in a window"""
  if iz is not None: id = _convertiztoid(iz)
  _extpartlist[id].setuparrays()

def getnep(iz=None,id=0,js=0):
  """Get the number of particles saved in a windows"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].nep(js)

def izlistextpart():
  """Returns a list of iz's for the windows"""
  # --- This list should always be identical to top.izepwin
  result = []
  for e in _extpartlist: result.append(e.iz)
  return result

def dumpextpart(filename):
  """Dump the saved extrapolated data to a file
 - filename: The name of the file to save the data in"""
  if me == 0:
    ff = open(filename,'w')
    cPickle.dump(_extpartlist,ff,1)
    ff.close()

def restoreextpart(filename):
  """Restore extrapolated data from the given file"""
  global _extpartlist
  if me == 0:
    ff = open(filename,'r')
    _extpartlist = cPickle.load(ff)
    ff.close()
  for e in _extpartlist: e.enable()


############################################################################
# --- Define functions which return the particle data
def extpartiz(iz=None,id=0,js=0):
  if iz is not None: id = _convertiztoid(iz) # --- or "return iz"
  return _extpartlist[id].iz
def extpartt(iz=None,id=0,js=0):
  """Returns x given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].t(js)
def extpartx(iz=None,id=0,js=0):
  """Returns x given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].x(js)
def extparty(iz=None,id=0,js=0):
  """Returns y given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].y(js)
def extpartux(iz=None,id=0,js=0):
  """Returns ux given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].ux(js)
def extpartuy(iz=None,id=0,js=0):
  """Returns uy given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].uy(js)
def extpartuz(iz=None,id=0,js=0):
  """Returns uz given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].uz(js)
def extpartvx(iz=None,id=0,js=0):
  """Returns vx given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].vx(js)
def extpartvy(iz=None,id=0,js=0):
  """Returns vy given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].vy(js)
def extpartvz(iz=None,id=0,js=0):
  """Returns vz given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].vz(js)
def extpartxp(iz=None,id=0,js=0):
  """Returns xp given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].xp(js)
def extpartyp(iz=None,id=0,js=0):
  """Returns yp given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].yp(js)
def extpartr(iz=None,id=0,js=0):
  """Returns r given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].r(js)
def extparttheta(iz=None,id=0,js=0):
  """Returns theta given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].theta(js)
def extpartrp(iz=None,id=0,js=0):
  """Returns rp given iz or id and optionally js (defaults to 0)"""
  if iz is not None: id = _convertiztoid(iz)
  return _extpartlist[id].rp(js)

############################################################################
############################################################################
############################################################################
# --- Define plotting routines for the extrapolated particles.

def _checkepplotarguments(kw):
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

############################################################################
def epxy(iz=None,id=0,js=0,particles=1,**kw):
  """Plots X-Y for extraploated particles"""
  _checkepplotarguments(kw)
  x = extpartx(iz=iz,id=id,js=js)
  y = extparty(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,top.yplmin,top.yplmax)
  settitles("Y vs X","X","Y","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(y,x,kwdict=kw)

############################################################################
def epxxp(iz=None,id=0,js=0,slope=0.,offset=0.,particles=1,**kw):
  """Plots X-X' for extraploated particles"""
  _checkepplotarguments(kw)
  x = extpartx(iz=iz,id=id,js=js)
  xp = extpartxp(iz=iz,id=id,js=js)
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
  settitles("X' vs X","X","X'","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(xp,x,kwdict=kw)

############################################################################
def epyyp(iz=None,id=0,js=0,slope=0.,offset=0.,particles=1,**kw):
  """Plots Y-Y' for extraploated particles"""
  _checkepplotarguments(kw)
  y = extparty(iz=iz,id=id,js=js)
  yp = extpartyp(iz=iz,id=id,js=js)
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
  settitles("Y' vs Y","Y","Y'","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(yp,y,kwdict=kw)

############################################################################
def epxpyp(iz=None,id=0,js=0,particles=1,**kw):
  """Plots X'-Y' for extraploated particles"""
  _checkepplotarguments(kw)
  xp = extpartxp(iz=iz,id=id,js=js)
  yp = extpartyp(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)
  settitles("Y' vs X'","X'","Y'","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(yp,xp,kwdict=kw)

############################################################################
def eprrp(iz=None,id=0,js=0,scale=0.,slope=0.,offset=0.,particles=1,**kw):
  """Plots R-R' for extraploated particles"""
  _checkepplotarguments(kw)
  x = extpartx(iz=iz,id=id,js=js)
  y = extparty(iz=iz,id=id,js=js)
  xp = extpartxp(iz=iz,id=id,js=js)
  yp = extpartyp(iz=iz,id=id,js=js)
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
  settitles("R' vs R","R","R'","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(rp,r,kwdict=kw)

############################################################################
def eptx(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-X for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  x = extpartx(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = ('e','e',top.xplmin,top.xplmax)
  settitles("X vs time","time","X","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(x,t,kwdict=kw)

############################################################################
def epty(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-Y for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  y = extparty(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = ('e','e',top.yplmin,top.yplmax)
  settitles("Y vs time","time","Y","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(y,t,kwdict=kw)

############################################################################
def eptxp(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-X' for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  xp = extpartxp(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = ('e','e',top.xpplmin,top.xpplmax)
  settitles("X' vs time","time","X'","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(xp,t,kwdict=kw)

############################################################################
def eptyp(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-Y' for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  yp = extpartyp(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  else:
    kw['pplimits'] = ('e','e',top.ypplmin,top.ypplmax)
  settitles("Y' vs time","time","Y'","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(yp,t,kwdict=kw)

############################################################################
def eptux(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-ux for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  ux = extpartux(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  settitles("ux vs time","time","ux","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(ux,t,kwdict=kw)

############################################################################
def eptuy(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-uy for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  uy = extpartuy(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  settitles("uy vs time","time","uy","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(uy,t,kwdict=kw)

############################################################################
def eptuz(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-uz for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  uz = extpartuz(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  settitles("uz vs time","time","uz","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(uz,t,kwdict=kw)

############################################################################
def eptvx(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-Vx for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  vx = extpartvx(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  settitles("Vx vs time","time","Vx","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(vx,t,kwdict=kw)

############################################################################
def eptvy(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-Vy for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  vy = extpartvy(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  settitles("Vy vs time","time","Vy","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(vy,t,kwdict=kw)

############################################################################
def eptvz(iz=None,id=0,js=0,particles=1,**kw):
  """Plots time-Vz for extraploated particles"""
  _checkepplotarguments(kw)
  t = extpartt(iz=iz,id=id,js=js)
  vz = extpartvz(iz=iz,id=id,js=js)
  kw['particles'] = particles
  if 'pplimits' in kw.keys():
    kw['lframe'] = 1
  settitles("Vz vs time","time","Vz","iz = %d"%extpartiz(iz=iz,id=id,js=js))
  ppgeneric(vz,t,kwdict=kw)

############################################################################
def eptrace(iz=None,id=0,js=0,slope=0.,particles=1,pplimits=None,**kw):
  """
Plots X-Y, X-X', Y-Y', Y'-X' in single page
If slope='auto', it is calculated from the moments for X-X' and Y-Y' plots.
pplimits can be a list of up to four tuples, one for each phase space plot.
If any of the tuples are empty, the limits used will be the usual ones for
that plot.
  """
  _checkepplotarguments(kw)
  x = extpartx(iz=iz,id=id,js=js)
  y = extparty(iz=iz,id=id,js=js)
  xp = extpartxp(iz=iz,id=id,js=js)
  yp = extpartyp(iz=iz,id=id,js=js)
  titler = "iz = %d"%extpartiz(iz=iz,id=id,js=js)
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

