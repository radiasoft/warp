from warp import *
from drifts import *
import __main__
realboundaries_version = "$Id: realboundaries.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

# --- Make plot command to plot particles and conductors.
def ppxycond(iw=0):
  "X-Y plot of particles including conductors"
  ppxy(iw)
  plp(fxy.ycond,fxy.xcond,color='red',msize=2.)

############################################################################
# --- Create classes to keep track of the element statuses:
# --- Save the id of the element that the beam is near, and whether or not the
# --- beam is in the element.
class ElemStatus:
  def isin(self):
    if not self.elems: return 0
    if self.zs[self.id] <= top.zbeam and top.zbeam <= self.ze[self.id]:
      return 1
    return 0
  def getid(self):
    offset = ((int(10.+(top.zbeam-top.zlatstrt)/top.zlatperi) - 10 )*
              top.zlatperi + top.zlatstrt)
    id = 0
    for ielt in xrange(1,len(self.ze)):
      zm = offset + 0.5*(self.ze[ielt-1] + self.zs[ielt])
      if zm < top.zbeam: id = ielt
    return id
  def update(self):
    self.elems = eval(self.sstring)
    if not self.elems: return 0
    id = self.getid()
    if self.id == id and self.inelem == self.isin():
      changed = 0
    else:
      changed = 1
    self.id = id
    self.inelem = self.isin()
    return changed
  def __init__(self,zs,ze,ap,ox,oy,sstring,rods):
    self.zs = zs
    self.ze = ze
    self.ap = ap
    self.ox = ox
    self.oy = oy
    self.sstring = sstring
    self.elems = eval(self.sstring)
    self.rods = rods*ones(len(self.zs))
    self.id = self.getid()
    self.inelem = self.isin()
  def getparams(self):
    if (not self.elems) or (not self.isin()): return ()
    self.offsetx = self.ox[self.id]
    self.offsety = self.oy[self.id]
    self.aperture = self.ap[self.id]
    return (self.aperture,self.offsetx,self.offsety)
  def arerods(self):
    return self.rods[self.id]

# --- Create a class to keep track of all of the element statusis
elementstati = []
class ElementStati:
  def isin(self):
    for s in self.stati:
      if s.isin(): return 1
    return 0
  def update(self):
    changed = 0
    for s in self.stati:
      if s.update(): changed = 1
    return changed
  def getparams(self):
    for s in self.stati:
      p = s.getparams()
      if p:
        self.aperture = p[0]
        self.offsetx  = p[1]
        self.offsety  = p[2]
        return (self.aperture,self.offsetx,self.offsety)
    return ()
  def arerods(self):
    for s in self.stati:
      if s.isin() and s.arerods(): return 1
    return 0
  def __init__(self):
    self.stati = []
    if top.quads:
      self.stati.append(ElemStatus(top.quadzs,top.quadze,top.quadap,
                                   top.qoffx,top.qoffy,'top.quads',
                                   greater(abs(top.quadde),0.)))
    if top.heles:
      self.stati.append(ElemStatus(top.helezs,top.heleze,top.heleap,
                                   top.heleox,top.heleoy,'top.heles',
                                   greater(abs(top.heleae[0,:]),0.)))
    if top.emlts:
      self.stati.append(ElemStatus(top.emltzs,top.emltze,top.emltap,
                                   top.emltox,top.emltoy,'top.emlts',1))
    if top.mmlts:
      self.stati.append(ElemStatus(top.mmltzs,top.mmltze,top.mmltap,
                                   top.mmltox,top.mmltoy,'top.mmlts',0))
    if top.pgrds:
      self.stati.append(ElemStatus(top.pgrdzs,top.pgrdze,top.pgrdap,
                                   top.pgrdox,top.pgrdoy,'top.pgrds',1))
    if top.bgrds:
      self.stati.append(ElemStatus(top.bgrdzs,top.bgrdze,top.bgrdap,
                                   top.bgrdox,top.bgrdoy,'top.bgrds',0))
    if top.drfts:
      self.stati.append(ElemStatus(top.drftzs,top.drftze,top.drftap,
				   top.drftox,top.drftoy,'top.drfts',0))

############################################################################
# --- Setup capacity matrix to include round pipe.
def makeroundpipe(ap,ox,oy):
  "Creates capacity matrix for a round pipe"
  niq = 2*nint((15./7.*ap)/w3d.dx)
  fxy.ncxymax = niq
  fxy.ncxy = fxy.ncxymax
  gchange("CapMatxy",0)
  symm_fact = 1.
  if (w3d.l2symtry): symm_fact = 0.5
  if (w3d.l4symtry): symm_fact = 0.25
  fxy.vcond[:] = 0.0
  fxy.xcond[:] = 15./7.*ap*cos(symm_fact*2.*pi*iota(0,niq-1)/niq)
  fxy.ycond[:] = 15./7.*ap*sin(symm_fact*2.*pi*iota(0,niq-1)/niq)
  fxy.xcond = fxy.xcond + ox
  fxy.ycond = fxy.ycond + oy
  fieldsol(0)
  ncxymaxsave.append(fxy.ncxymax)
  ncxysave.append(fxy.ncxy)
  xcondsave.append((fxy.xcond - ox)/(15./7.*ap))
  ycondsave.append((fxy.ycond - oy)/(15./7.*ap))
  vcondsave.append(fxy.vcond*1.+0.)
  cmatxysave.append(fxy.cmatxy*(w3d.xmmax-w3d.xmmin)*(w3d.ymmax-w3d.ymmin))
  kpvtxysave.append(fxy.kpvtxy*1+0)

# --- Setup capacity matrix to include quadrupole rods.
def makeroundrods(ap,ox,oy):
  "Creates capacity matrix for a set of four quadrupole rods"
  niq = nint((16./7.*ap)/w3d.dx)
  fxy.ncxymax = 4*niq
  fxy.ncxy = fxy.ncxymax
  gchange("CapMatxy",0)
  #symm_fact = 1.
  #if (w3d.l2symtry): symm_fact = 0.5
  #if (w3d.l4symtry): symm_fact = 0.25
  fxy.vcond[:] = 0.0
  xx = 8./7.*ap*cos(1.*pi*iota(0,niq-1)/niq)
  yy = 8./7.*ap*sin(1.*pi*iota(0,niq-1)/niq)
  fxy.xcond[0    :  niq] = ap + 8./7.*ap - yy
  fxy.ycond[0    :  niq] = xx
  fxy.xcond[  niq:2*niq] = -ap - 8./7.*ap + yy
  fxy.ycond[  niq:2*niq] = xx
  fxy.xcond[2*niq:3*niq] = xx
  fxy.ycond[2*niq:3*niq] = ap + 8./7.*ap - yy
  fxy.xcond[3*niq:4*niq] = xx
  fxy.ycond[3*niq:4*niq] = -ap - 8./7.*ap + yy
  fxy.xcond = fxy.xcond + ox
  fxy.ycond = fxy.ycond + oy
  fieldsol(0)
  ncxymaxsave.append(fxy.ncxymax)
  ncxysave.append(fxy.ncxy)
  xcondsave.append((fxy.xcond - ox)/(15./7.*ap))
  ycondsave.append((fxy.ycond - oy)/(15./7.*ap))
  vcondsave.append(fxy.vcond*1.+0.)
  cmatxysave.append(fxy.cmatxy*(15./7.*ap)**2)
  kpvtxysave.append(fxy.kpvtxy*1+0)

############################################################################
# --- This routine redefines the mesh to match the size of the element.
# --- The default is to extend far enough to cover half of quadrupole rods
# --- plus a little extra space to keep the conducting points away from
# --- the mesh edge.
def setmesh():
  "Set mesh and particle scraping parameters and resets FFT data"
  (elementap, elementox, elementoy) = elementstati[0].getparams()
  w3d.xmmin = elementox - 15./7.*elementap*(1. + 1./(w3d.nx/2-1))
  w3d.xmmax = elementox + 15./7.*elementap*(1. + 1./(w3d.nx/2-1))
  w3d.ymmin = elementoy - 15./7.*elementap*(1. + 1./(w3d.ny/2-1))
  w3d.ymmax = elementoy + 15./7.*elementap*(1. + 1./(w3d.ny/2-1))
  w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
  w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
  w3d.xmesh = iota(0,w3d.nx)*w3d.dx + w3d.xmmin
  w3d.ymesh = iota(0,w3d.ny)*w3d.dy + w3d.ymmin
  top.prwall = elementap
  top.prwallx = elementox
  top.prwally = elementoy
  top.prwallz = top.prwall
  top.prwallxz = top.prwallx
  top.prwallyz = top.prwally
  # --- Recalculate ksqx and ksqy (notice that for this part, the capacity
  # --- matrix is turned off since the matrix does not need to be calculated)
  if top.fstype >= 0:
    ff = top.fstype
    top.fstype = 0
    fieldsol(1)
    top.fstype = ff

def resetgrid():
  "Resets grid to match current element aperture"
  # --- Update the element status, which returns whether or not it actually
  # --- changed. If it did not change, the return, doing nothing.
  if not elementstati[0].update(): return
  # --- Setup the mesh.
  setmesh()
  if top.fstype >= 0:
    # --- Get the current element parameters
    (elementap, elementox, elementoy) = elementstati[0].getparams()
    # --- Check if in an electric quadrupole with rods.
    if elementstati[0].arerods():
      ic = 1
      if len(ncxymaxsave) == 1:
        makeroundrods(elementap,elementox,elementoy)
    else:
      ic = 0
    # --- Reset the capacity matrix data.
    fxy.ncxymax = ncxymaxsave[ic]
    fxy.ncxy    = ncxysave[ic]
    gchange("CapMatxy",0)
    # --- Retreive the data and scale it appropriately.
    fxy.xcond[:]   = xcondsave[ic]*15./7.*elementap + elementox
    fxy.ycond[:]   = ycondsave[ic]*15./7.*elementap + elementoy
    fxy.vcond[:]   = vcondsave[ic]
    fxy.cmatxy[:,:]=cmatxysave[ic]/((w3d.xmmax-w3d.xmmin)*(w3d.ymmax-w3d.ymmin))
    fxy.kpvtxy[:]  = kpvtxysave[ic]
    # --- Reload the charge density.
    w3d.rho = 0.
    loadrho()
    # --- Now recalculate the correct self field plus conductors.
    fieldsol(-1)

############################################################################
# --- Now, do initialization work.
# --- Create list of capacity matrices and associated data.
ncxymaxsave = []
ncxysave    = []
xcondsave   = []
ycondsave   = []
vcondsave   = []
cmatxysave  = []
kpvtxysave  = []

def initrealboundaries():
  "Initializes the things needed for using the real-boundaries module"
  #global elementstati
  global _mainstep

  # --- Get the drift elements if needed.
  initdrifts()

  # --- Turn on the capacity matrix field solver.
  top.fstype = 1

  # --- Create instance of element status class
  elementstati.append(ElementStati())
  (elementap, elementox, elementoy) = elementstati[0].getparams()

  # --- Initialize round pipe
  setmesh()
  makeroundpipe(elementap, elementox, elementoy)

  # --- Set grid to initial values
  resetgrid()

  # --- Add resetgrid to beforestep
  beforestep.append(resetgrid)
  # --- Replace the step command with one that always calls resetgrid.
  #if not __main__.__dict__.has_key('ctlstep'):
    #__main__.__dict__['ctlstep'] = __main__.__dict__['step']
  #_mainstep = __main__.__dict__['step']
  #def rqstep(n=1):
    #for i in xrange(n):
      #resetgrid()
      #_mainstep(1)
  #__main__.__dict__['step'] = rqstep
