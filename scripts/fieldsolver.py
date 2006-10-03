"""Generic class describing the interface needed for a field solver.
"""
from warp import *

class FieldSolver(object):

  def setrhopforparticles(self,irhopndtscopies,indts,iselfb):
    if irhopndtscopies is not None:
      rhopndts = self.getrhop(lndts=1)
      self.rhop = rhopndts[...,irhopndtscopies,indts]
    else:
      rhopselfb = self.getrhop(lselfb=1)
      self.rhop = rhopselfb[...,iselfb]

  def aftersetrhop(self,lzero):
    "Anything that needs to be done to rhop after the deposition"
    pass

  def loadrho(self,lzero=true,**kw):
    'Charge deposition, uses particles from top directly'
    self.allocatedataarrays()
    if lzero: self.zerorhop()

    for js,i,n,q,w in zip(arange(top.pgroup.ns),top.pgroup.ins-1,
                          top.pgroup.nps,top.pgroup.sq,top.pgroup.sw):
      if n > 0:
        args = [top.pgroup.xp[i:i+n],top.pgroup.yp[i:i+n],
                top.pgroup.zp[i:i+n],top.pgroup.uzp[i:i+n],
                q,w*top.pgroup.dtscale[js]]
        if top.pgroup.ldts[js]:
          indts = top.ndtstorho[top.pgroup.ndts[js]-1]
          self.setrhopforparticles(0,indts,None)
          self.setrhop(*args)
          self.aftersetrhop(lzero)
        if top.pgroup.iselfb[js] > -1:
          iselfb = top.pgroup.iselfb[js]
          self.setrhopforparticles(None,None,iselfb)
          self.setrhop(*args)
          self.aftersetrhop(lzero)

    self.averagerhopwithsubcycling()

  def getphip(self,lndts=None,lselfb=None):
    if lndts is None and lselfb is None:
      return self.phiparray
    elif lndts is not None:
      return self.phiparray[...,:top.nsndtsphi]
    elif lselfb is not None:
      return self.phiparray[...,top.nsndtsphi:]

  def getrhop(self,lndts=None,lselfb=None):
    if lndts is None and lselfb is None:
      return self.rhoparray
    elif lndts is not None:
      rhop = self.rhoparray[...,:top.nrhopndtscopies*top.nsndts]
      trhop = transpose(rhop)
      sss = list(shape(rhop)[:-1])+[top.nrhopndtscopies,top.nsndts]
      sss.reverse()
      trhop.shape = sss
      return transpose(trhop)
    elif lselfb is not None:
      return self.rhoparray[...,top.nrhopndtscopies*top.nsndts:]

  def loadj(self,lzero=true,**kw):
    'Charge deposition, uses particles from top directly'
    pass

  def setphipforparticles(self,indts,iselfb):
    if indts is not None:
      phipndts = self.getphip(lndts=1)
      self.phip = phipndts[...,min(indts,w3d.nsndtsphi3d-1)]
    else:
      phipselfb = self.getphip(lselfb=1)
      self.phip = phipselfb[...,iselfb]

  def fetche(self,**kw):
    'Fetches the E field, uses arrays from w3d module FieldSolveAPI'
    if w3d.api_xlf2:
      w3d.xfsapi = top.pgroup.xp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.yfsapi = top.pgroup.yp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
      w3d.zfsapi = top.pgroup.zp[w3d.ipminapi-1:w3d.ipminapi-1+w3d.ipapi]
    args = [w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
            w3d.exfsapi,w3d.eyfsapi,w3d.ezfsapi]

    js = w3d.jsapi

    tmpnsndts = getnsndtsforsubcycling()
    indts = min(tmpnsndts-1,top.ndtstorho[top.pgroup.ndts[js]-1])
    self.setphipforparticles(indts,None)
    self.fetchefrompositions(*args,**kw)

    # add self B correction as needed
    if top.pgroup.iselfb[js] > -1:
      exfsapispecies = zeros(shape(w3d.exfsapi),'d')
      eyfsapispecies = zeros(shape(w3d.eyfsapi),'d')
      ezfsapispecies = zeros(shape(w3d.ezfsapi),'d')
      self.setphipforparticles(None,top.pgroup.iselfb[js])
      args = [w3d.xfsapi,w3d.yfsapi,w3d.zfsapi,
              exfsapispecies,eyfsapispecies,ezfsapispecies]
      self.fetchefrompositions(*args,**kw)
      w3d.exfsapi[...] += exfsapispecies
      w3d.eyfsapi[...] += eyfsapispecies

  def fetchb(self):
    'Fetches the B field, uses arrays from w3d module FieldSolveAPI'
    pass
  def fetchphi(self):
    'Fetches the potential, uses arrays from w3d module FieldSolveAPI'
    pass

  def setrhoandphiforfieldsolve(self,irhopndtscopies,indts,iselfb):
    if irhopndtscopies is not None:
      rhopndts = self.getrhop(lndts=1)
      phipndts = self.getphip(lndts=1)
      self.rho = rhopndts[...,irhopndtscopies,indts]
      self.phi = phipndts[...,min(indts,w3d.nsndtsphi3d-1)]
    else:
      rhopselfb = self.getrhop(lselfb=1)
      phipselfb = self.getphip(lselfb=1)
      self.rho = rhopselfb[...,iselfb]
      self.phi = phipselfb[...,iselfb]

  def getphipforparticles(self,indts,iselfb):
    if indts is not None:
      phipndts = self.getphip(lndts=1)
      phipndts[...,min(indts,w3d.nsndtsphi3d-1)] = self.phi[...]
    else:
      phipselfb = self.getphip(lselfb=1)
      phipselfb[...,iselfb] = self.phi[...]

  def selfbcorrection(self,js):
    "scale phispecies by -(1-1/gamma*2) store into top.fselfb"
    phipselfb = self.getphip(lselfb=1)
    phipselfb[...,js] *= top.fselfb[js]

  def dosolveonphi(self,iwhich,irhopndtscopies,indts,iselfb):
    "points rho and phi appropriately and call the solving routine"
    self.setrhoandphiforfieldsolve(irhopndtscopies,indts,iselfb)
    self.dosolve(iwhich,irhopndtscopies,indts,iselfb)
    if iselfb is not None: self.selfbcorrection(iselfb)
    self.getphipforparticles(indts,iselfb)

  def solve(self,iwhich=0):
    self.allocatedataarrays()
    # --- Loop over the subcyling groups and do any field solves that
    # --- are necessary.
    # --- Do loop in reverse order so that rho and phi end up with the arrays
    # --- for the speices with the smallest timestep.
    tmpnsndts = getnsndtsforsubcycling()
    for indts in range(tmpnsndts-1,-1,-1):
      if (not top.ldts[indts] and
          (top.ndtsaveraging == 0 and not sum(top.ldts))): continue
      self.dosolveonphi(iwhich,top.nrhopndtscopies-1,indts,None)

    # --- Solve for phi for groups which require the self B correction
    for js in range(top.nsselfb):
      self.dosolveonphi(iwhich,None,None,js)

  def installconductor(self,conductor):
    'Installs conductor. Does nothing if not defined'
    pass

  def getlinecharge(self):
    'Gets line-charge'
    raise 'getlinecharge must be defined'
  def getese(self):
    'Gets electrostatic energy, rho*phi'
    raise 'getese must be defined'

  def pfzx(self,**kw):
    'Plots potential in z-x plane. Does nothing if not defined.'
    pass
  def pfzy(self,**kw):
    'Plots potential in z-y plane. Does nothing if not defined.'
    pass
  def pfxy(self,**kw):
    'Plots potential in x-y plane. Does nothing if not defined.'
    pass

  def getpdims(self):
    raise """getpdims must be supplied - it should return a list of the dimensions
of the rho and field arrays"""

  def allocatedataarrays(self):
    # --- Get base dimension of the arrays
    rhopdims,phipdims = self.getpdims()

    # --- Setup rho and phi arrays, including extra copies for subcycling
    # --- and self B corrections.
    setupSubcycling(top.pgroup)
    setupSelfB(top.pgroup)

    extrarhopdim = top.nrhopndtscopies*top.nsndts + top.nsselfb
    dims = list(rhopdims) + [extrarhopdim]
    if 'rhoparray' not in self.__dict__ or shape(self.rhoparray) != tuple(dims):
      self.rhoparray = fzeros(dims,'d')

    extraphidim = top.nsndtsphi + top.nsselfb
    dims = list(phipdims) + [extraphidim]
    if 'phiparray' not in self.__dict__ or shape(self.phiparray) != tuple(dims):
      self.phiparray = fzeros(dims,'d')

  def zerorhop(self):
    if top.ndtsaveraging == 0:
      self.zerorhopwithsampledsubcycling()
    elif top.ndtsaveraging == 1:
      self.zerorhopwithfullvsubcycling()
    elif top.ndtsaveraging == 2:
      self.zerorhopwithhalfvsubcycling()
    self.zerorhopselfb()

  def zerorhopwithsampledsubcycling(self):
    # --- Zero the rhop copy for species when the positions
    # --- are advanced.
    # --- rhopndts(...,2,in1) holds the old rhop which is still needed
    # --- for the faster advanced groups.
    rhopndts = self.getrhop(lndts=1)
    for in1 in range(top.nsndts):
      if top.ldts[in1]:
        if top.nrhopndtscopies == 2:
          rhopndts[...,1,in1] = rhopndts[...,0,in1]
        rhopndts[...,0,in1] = 0.

  def zerorhopwithfullvsubcycling(self):
    raise "fullv subcycling not yet implemented"

  def zerorhopwithhalfvsubcycling(self):
    raise "halfv subcycling not yet implemented"

  def averagerhopwithsubcycling(self):
    if top.ndtsaveraging == 0:
      self.averagerhopwithsampledsubcycling()
    elif top.ndtsaveraging == 1:
      self.averagerhopwithfullvsubcycling()
    elif top.ndtsaveraging == 2:
      self.averagerhopwithhalfvsubcycling()

  def averagerhopwithsampledsubcycling(self):
    if top.ndtsmax == 1: return

    rhopndts = self.getrhop(lndts=1)

    # --- During the generate, do the copy of the new rhop to the old rhop
    # --- for group 0, which is normally done during the zerorhop call.
    if top.it == 0: rhopndts[...,1,0] = rhopndts[...,0,0]

    # --- Save the rhop where the fastest particle's rhop is. For now,
    # --- assume that this is in1=0
    for in1 in range(1,top.nsndts):
      if top.it == 0:
        # --- At top.it==0, before the first step, always add the new rhop.
        rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,0,in1])
      elif top.rhotondts[in1]%2 == 1:
        # --- Use the rhop that is closest in time to the current time.
        if ((top.it-1)%top.rhotondts[in1] > top.rhotondts[in1]/2. or
            (top.it-1)%top.rhotondts[in1] == 0):
          rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,0,in1])
        else:
          rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,1,in1])
      else:
        # --- When ndts is even, at the mid point of the step, take the
        # --- average of the old and the new
        # --- Otherwise, use the rhop that is closest in time to the currentc
        # --- time.
        if (top.it-1)%top.rhotondts[in1] == top.rhotondts[in1]/2:
          rhopndts[...,1,0] = (rhopndts[...,1,0] +
               0.5*(rhopndts[...,0,in1] + rhopndts[...,1,in1]))
        elif ((top.it-1)%top.rhotondts[in1] > top.rhotondts[in1]/2. or
              (top.it-1)%top.rhotondts[in1] == 0):
          rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,0,in1])
        else:
          rhopndts[...,1,0] = (rhopndts[...,1,0] + rhopndts[...,1,in1])

  def zerorhopselfb(self):
    rhopselfb = self.getrhop(lselfb=1)
    rhopselfb[...] = 0.

