"""Generates the applied field from the helix accelerating structure
"""
from warp import *
from plot_conductor import setconductorvoltage
from appendablearray import AppendableArray

class LinearHelix:
  """
Creates helix with constant Ez in the pulse.
 - zs,ze: extent over which the full helix extends
 - ap: aperture of helix
 - length: length of pulse - length of the region where Ez is.
 - V0: voltage change along pulse
 - vz: velocity of pulse
 - ts: time at which pulse end is at zs
 - nz=10000: number of grid points
 - lengthfinal: final length of the pulse
 - zoffset=0.: offset of z start of pulse head relative to zs
 - startgaplen=0.2: length of gap at start of helix
 - finalgaplen=0.2: length of gap at end of helix
 - ezoft,ezduration=None: lists of ez and the durations of each time block
 - vvdata=None: Specifies voltage as a function of time.
                The change of voltage in each time slice is given by
                V(t) = Vt*(v1*t/dt + v2*(t/dt)**2 + v3*(t/dt)**3 + ... )
                where Vt and dt are either directly from vvdata
                (vvdatatype='delta')
                or Vt=(V(t+dt)-V(t)) and dt=(T(i+1)-T(i)
                (vvdatatype='absolute').
                vvdata is a list of parameters for the time slices. Each set
                of data is V, dt, v1, v2, ..., where
                V is either the Vt, or V(t), and dt is either dt or T
                according to the above.
                dt is the length of the time slice
                v1, v2, ... are the fractions assigned to each polynomial term
 - vvdatatype='delta': When 'delta', vvdata is assumed to be the changes in
                       voltage during each time slice, or when 'absolute',
                       the voltage at the start at each time slice.
 - vvfunc=None: A function which returns the voltage and electric field given
                a single argument, the time.
 - startgaptype=0: specified what to do about starting gap
                   0 means inductive coils, which gives an inverse voltage
                     change in the firing region
                   1 means induction cell, which gives the same Ez in the
                     firing range.
 - model=0: Specifies the model to use to apply the waveforms
            0 - uses purely applied fields, via the accl and emlt lattice
                elements
            1 - the voltage is used as a boundary condition on the field
                solve - a condid must be specified
            2 - same as 1, but curcuit equations are used to advance the
                voltage rather than simple advection
 - condid=None: id conductor where the boundary condition is applied for
                models 1 and 2
 - inductance=None: required for model=2
 - capacitance=None: required for model=2
 - resistance=None: required for model=2
 - terminductance=None: required for model=2
 - termcapacitance=None: required for model=2
 - termresistance=None: required for model=2
 - pitch=None: pitch of helix, now only used in calculation of the B field.
               Default to dzwave.
 - lbeamloading=1: turn on or off beam loading for model=2
 - ncircuit=10: number of circuit steps to take per particle step
 - nspirals=None: number of spirals in helix, normally (ze-zs)/(vz*dt)
 - lsavehist=0: when true, saves history of voltages
 - lcalcunloaded=0: when true, in a side calculation, calculates the voltages
                    without beam loading (only applies when lbeamloading=1).
 - lcalcnowave=0: when true, in a side calculation, calculates the voltages
                  with only beam loading
 - lcalcbfield=0: when true, calculates the B field by solving the
                  magnetostatic Poisson equation using the helix current.
                  This only works with model==2.
  """
  def __init__(self,zs,ze,ap=None,length=None,V0=None,vz=None,ts=0.,nz=10000,
                    lengthfinal=None,
                    zoffset=0.,startgaplen=0.2,finalgaplen=0.2,
                    ezoft=None,ezduration=None,vvdata=None,vvdatatype='delta',
                    vvfunc=None,
                    startgaptype=0,model=0,condid=None,
                    inductance=None,capacitance=None,resistance=None,
                    terminductance=None,termcapacitance=None,
                    termresistance=None,pitch=None,
                    lbeamloading=1,ncircuit=10,nspirals=None,lsavehist=0,
                    lcalcunloaded=0,lcalcnowave=0,lcalcbfield=0):
    # --- Check input values
    if model in [1,2]:
      assert condid is not None,\
             "If model 1 or 2 is used, a conid must be specified"
      assert startgaptype == 0,\
             "Model 1 can only be used when startgaptype is 0"

    # --- Save input quantities
    self.zs = zs
    self.ze = ze
    self.ap = ap
    self.length = length
    self.V0 = V0
    if lengthfinal is None: self.lengthfinal = length
    else:                   self.lengthfinal = lengthfinal
    self.vz = vz
    self.ts = ts # --- Maybe reset below based on other input data
    self.zoffset = zoffset
    self.startgaplen = startgaplen
    self.finalgaplen = finalgaplen
    self.startgaptype = startgaptype
    self.model = model
    self.condid = condid
    self.inductance = inductance
    self.capacitance = capacitance
    self.resistance = resistance
    self.termcapacitance = where(termcapacitance is None,capacitance,
                                 termcapacitance)
    self.terminductance = where(terminductance is None,inductance,
                                terminductance)
    self.termresistance = termresistance
    self.pitch = pitch
    self.lbeamloading = lbeamloading
    self.ncircuit = ncircuit
    self.nspirals = nspirals
    self.lsavehist = lsavehist
    self.lcalcunloaded = lcalcunloaded
    self.lcalcnowave = lcalcnowave
    self.lcalcbfield = lcalcbfield

    self.ezoft = ezoft
    self.ezduration = ezduration

    if V0 is not None and length is not None:
      vvdata = [[self.V0,length/vz]]
    elif ezoft is not None and ezduration is not None:
      voft = array(ezoft)*self.vz*array(ezduration)
      vvdata = zip(voft,ezduration)

    self.vvdata = vvdata
    self.vvdatatype = vvdatatype
    self.vvfunc = vvfunc

    if self.vvdata is not None:
      self.ioft = 0
      self.vslice = [self.vvdata[i][0] for i in range(len(self.vvdata))]
      self.tslice = [self.vvdata[i][1] for i in range(len(self.vvdata))]
      if self.vvdatatype == 'delta':
        self.ts = ts
        toft = cumsum(self.tslice) + self.ts
        self.tlist = array([self.ts] + list(toft))
        voft = cumsum(self.vslice)
        self.vlist = array([0.] + list(voft))
      else:
        self.tlist = self.tslice
        self.vlist = self.vslice
        self.ts = self.tlist[0]

    self.vstart = 0.
    self.vend = 0.
    self.vvprevious = 0.

    self.setuparrays()
    if self.model == 0: self.setuplattice()
    #elif self.model == 1: self.setupboundarycondition()

    if self.pitch is None: self.pitch = self.dzwave

    self.hv = [0.]
    self.hez = [0.]
    self.ht = [self.ts]


    if self.model == 0:
      self.setlatticefield()
      installbeforefs(self.setlatticefield)
    elif self.model == 1:
      self.setboundaries()
      installbeforefs(self.setboundaries)
    elif self.model == 2:
      self.advanceboundaries()
      installbeforefs(self.advanceboundaries)
      if self.lcalcbfield: top.bfstype = 7

  def setuparrays(self):

    self.zswave = self.zs + self.startgaplen
    self.zewave = self.ze - self.finalgaplen
    if self.nspirals is None:
      self.dzwave = self.vz*top.dt
    else:
      self.dzwave = (self.zewave - self.zswave)/self.nspirals
    self.nzwave = int((self.zewave - self.zswave)/self.dzwave+1)
    self.vwave = zeros(1+self.nzwave,'d')
    self.ezwave = zeros(1+self.nzwave,'d')
    self.zwave = iota(0,self.nzwave)*self.dzwave + self.zswave

    if self.lsavehist:
      self.hvwave = AppendableArray(initunit=self.vwave)

    if self.model == 2:
      self.iwave = zeros(self.nzwave,'d')
      self.dzwave = (self.zewave - self.zswave)/self.nzwave
      self.zsterm = self.ze - self.finalgaplen
      self.zeterm = self.ze
      # --- Setup so dzterm ~ dzwave
      self.nzterm = nint((self.zeterm - self.zsterm)/self.dzwave)
      self.dzterm = (self.zeterm - self.zsterm)/self.nzterm
      self.zterm = iota(0,self.nzterm)*self.dzterm + self.zsterm
      self.vterm = zeros(1+self.nzterm,'d')
      self.iterm = zeros(self.nzterm,'d')
      if self.lbeamloading:
        self.phiawave = zeros(1+self.nzwave,'d')
        self.phiaterm = zeros(1+self.nzterm,'d')
        r = int(self.ap/w3d.dx)*w3d.dx
        if (self.ap - r) < w3d.dx/10.: r = r - w3d.dx
        self.xwave = r + zeros(1+self.nzwave,'d')
        self.xterm = r + zeros(1+self.nzterm,'d')
        self.ywave =     zeros(1+self.nzwave,'d')
        self.yterm =     zeros(1+self.nzterm,'d')
        self.lambdawaveprev = zeros(1+self.nzwave,'d')
        self.lambdatermprev = zeros(1+self.nzterm,'d')
      if self.lcalcunloaded:
        self.vwavenoload = zeros(1+self.nzwave,'d')
        self.iwavenoload = zeros(self.nzwave,'d')
        self.vtermnoload = zeros(1+self.nzterm,'d')
        self.itermnoload = zeros(self.nzterm,'d')
      if self.lcalcnowave:
        self.vwavenowave = zeros(1+self.nzwave,'d')
        self.iwavenowave = zeros(self.nzwave,'d')
        self.vtermnowave = zeros(1+self.nzterm,'d')
        self.itermnowave = zeros(self.nzterm,'d')

      self.delibeamwave = zeros(1+self.nzwave,'d')
      self.delibeamterm = zeros(1+self.nzterm,'d')
      if self.lsavehist:
        self.hiwave = AppendableArray(initunit=self.iwave)
        self.hvterm = AppendableArray(initunit=self.vterm)
        self.hiterm = AppendableArray(initunit=self.iterm)
        if self.lbeamloading:
          self.hlambdawave = AppendableArray(initunit=self.lambdawaveprev)
          self.hlambdaterm = AppendableArray(initunit=self.lambdatermprev)
          if self.lcalcunloaded:
            self.hvwavenoload = AppendableArray(initunit=self.vwavenoload)
            self.hiwavenoload = AppendableArray(initunit=self.iwavenoload)
            self.hvtermnoload = AppendableArray(initunit=self.vtermnoload)
            self.hitermnoload = AppendableArray(initunit=self.itermnoload)
          if self.lcalcnowave:
            self.hvwavenowave = AppendableArray(initunit=self.vwavenowave)
            self.hiwavenowave = AppendableArray(initunit=self.iwavenowave)
            self.hvtermnowave = AppendableArray(initunit=self.vtermnowave)
            self.hitermnowave = AppendableArray(initunit=self.itermnowave)


  def setuplattice(self):
    # --- Setup emlt element which is used to apply the moving pulse
    top.nemlt = top.nemlt + 1
    self.eid = top.nemlt
    gchange("Lattice")
    top.emltzs[self.eid] = self.zs + self.startgaplen
    top.emltze[self.eid] = self.ze - self.finalgaplen
    top.emltid[self.eid] = self.eid + 1
    top.emltap[self.eid] = self.ap
    top.nemltsets = top.nemltsets + 1
    top.nesmult = 1
    top.nzemltmax = max(self.nzwave,top.nzemltmax)
    gchange('Mult_data')
    top.nzemlt[self.eid] = self.nzwave
    top.dzemlt[self.eid] = self.dzwave
    top.emlt_n[0] = 0
    top.emlt_v[0] = 0
    top.esemlt[:,:,self.eid] = 0.

    # --- Setup the accl element which gives the acceleration at the start
    top.naccl = top.naccl + 1
    gchange("Lattice")
    self.asid = top.naccl
    top.acclzs[self.asid] = self.zs
    top.acclze[self.asid] = self.zs + self.startgaplen
    top.acclez[self.asid] = 0.
    top.acclap[self.asid] = self.ap

    # --- Setup the accl element which gives the acceleration at the end
    top.naccl = top.naccl + 1
    gchange("Lattice")
    self.aeid = top.naccl
    top.acclzs[self.aeid] = self.ze - self.finalgaplen
    top.acclze[self.aeid] = self.ze
    top.acclez[self.aeid] = 0.
    top.acclap[self.aeid] = self.ap

    resetlat()

  def getvoltagefromvvdata(self,time):

    # --- Check if the next time slice is starting
    if time < self.tlist[self.ioft]: self.ioft = 0
    while (self.ioft < len(self.vlist)-1 and
           time > self.tlist[self.ioft+1]):
      self.ioft += 1

    # --- Get current value of voltage at helix start.
    if self.ioft < len(self.vlist)-1:
      # --- Voltage at the start of the current time slice.
      vbase = self.vlist[self.ioft]
      # --- Time relative to start of time slice
      tt = time - self.tlist[self.ioft]
      # --- Get polynomial cofficients
      if len(self.vvdata[self.ioft]) == 2:
        vp = [1.]
      else:
        vp = self.vvdata[self.ioft][2:]

      dv = self.vlist[self.ioft+1] - self.vlist[self.ioft]
      dt = self.tlist[self.ioft+1] - self.tlist[self.ioft]

      # --- Now the voltage can be calculated
      vv = 0.
      for i in range(len(vp)):
        vv = vv + vp[i]*(tt/dt)**(i+1)
      vv = vbase + vv*dv

      # --- Calculate dV/dt, to get Ez = dV/dt / vc
      vdot = 0.
      for i in range(len(vp)):
        vdot = vdot + (i+1)*vp[i]*tt**i/dt**(i+1)
      ez = vdot*dv/self.vz

      self.vvprevious = vv

    else:
      vv = self.vvprevious
      ez = 0.

    return vv,ez

  def getvoltagefromvvfunc(self,time):
    return self.vvfunc(time)

  def getvoltage(self,time,advect=1):

    if time < self.ts: return 0.,0.

    if self.vvdata is not None: vv,ez = self.getvoltagefromvvdata(time)
    if self.vvfunc is not None: vv,ez = self.getvoltagefromvvfunc(time)

    # --- Advect the vwave and ezwave meshes
    if advect:
      self.vwave[1:] = self.vwave[:-1]
      self.vwave[0] = vv
      self.ezwave[1:] = self.ezwave[:-1]
      self.ezwave[0] = ez

      # --- Save the time histories
      self.hv.append(vv)
      self.hez.append(ez)
      self.ht.append(time)

    return vv,ez

  def setlatticefield(self,time=None):
    "Used for model == 0"

    if time is None: time = top.time
    if time < self.ts: return

    vv,ez = self.getvoltage(time)

    # --- Adjust voltage at start
    top.accls = 1
    self.vstart = vv
    if self.startgaptype == 0:
      top.acclez[self.asid] = -self.vstart/self.startgaplen
    elif self.startgaptype == 1:
      top.acclez[self.asid] = ez

    # --- Adjust voltage at end
    iz = top.nzemlt[self.eid]
    self.vend = self.vend + top.esemlt[iz,0,self.eid]*self.vz*top.dt
    top.acclez[self.aeid] = self.vend/self.finalgaplen

    # --- Advect top.esemlt
    top.esemlt[1:iz+1,:,self.eid] = top.esemlt[:iz,:,self.eid].copy()
    top.esemlt[0,0,self.eid] = ez

    # --- locations where the ez changes
    zlist = self.zswave + (top.time - self.tlist)*self.vz
    self.zlist = maximum(self.zswave,minimum(self.zewave,zlist))

  def setboundaries(self,time=None):
    "Used for model == 1"

    if time is None: time = top.time
    if time < self.ts: return

    vv,ez = self.getvoltage(time)

    # --- Interpolate from helix mesh to the field grid
    iz  = int((w3d.zmeshlocal + top.zgrid - self.zswave)/self.dzwave)
    wz0 =     (w3d.zmeshlocal + top.zgrid - self.zswave)/self.dzwave - iz
    wz1 = 1. - wz0
    
    # --- Zero out any points beyond the helix mesh
    wz0 = where((0 <= iz) & (iz < self.nzwave),wz0,0.)
    wz1 = where((0 <= iz) & (iz < self.nzwave),wz1,0.)
    iz  = where((0 <= iz) & (iz < self.nzwave),iz,0)

    vgrid = take(self.vwave,iz)*wz1 + take(self.vwave,iz+1)*wz0

    # --- Add in voltages of start and end regions.
    iz1 = max(0,int((self.zs - top.zgrid - w3d.zmmin)/w3d.dz) + 1)
    iz2 = min(w3d.nz,
              int((self.zs + self.startgaplen - top.zgrid - w3d.zmmin)/w3d.dz))
    if iz1 <= iz2:
      zz = iota(iz1,iz2)*w3d.dz + top.zgrid + w3d.zmmin
      vgrid[iz1:iz2+1] = (zz - self.zs)/self.startgaplen*self.vwave[0]
    iz1 = max(0,
              int((self.ze - self.finalgaplen - top.zgrid - w3d.zmmin)/w3d.dz))
    iz2 = min(w3d.nz,int((self.ze - top.zgrid - w3d.zmmin)/w3d.dz) + 1)
    if iz1 <= iz2:
      zz = iota(iz1,iz2)*w3d.dz + top.zgrid + w3d.zmmin
      vgrid[iz1:iz2+1] = (self.ze - zz)/self.finalgaplen*self.vwave[-1]

    setconductorvoltage(vgrid,condid=self.condid)

  def getbeamloading(self):
    # --- Save phi which includes helix boundary conditions
    fullphi = frz.basegrid.phi.copy()

    # --- Set helix conductor to ground and redo the field solve
    setconductorvoltage(0.,condid=self.condid)
    fieldsol(-1)

    # --- Get phi at grid point just below pipe radius
    getgrid3d(1+self.nzwave,self.xwave,self.ywave,self.zwave,
              self.phiawave,w3d.nx,w3d.ny,w3d.nz,w3d.phi[:,:,1:-1],
              w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
              w3d.zmminlocal+top.zgrid,w3d.zmmaxlocal+top.zgrid,
              w3d.l2symtry,w3d.l4symtry)
    getgrid3d(1+self.nzterm,self.xterm,self.yterm,self.zterm,
              self.phiaterm,w3d.nx,w3d.ny,w3d.nz,w3d.phi[:,:,1:-1],
              w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,
              w3d.zmminlocal+top.zgrid,w3d.zmmaxlocal+top.zgrid,
              w3d.l2symtry,w3d.l4symtry)

    self.lambdawave = 2.*pi*self.ap*self.phiawave/(self.ap - self.xwave)*eps0
    self.lambdaterm = 2.*pi*self.ap*self.phiaterm/(self.ap - self.xterm)*eps0
    self.delibeamwave = (self.lambdawave - self.lambdawaveprev)/top.dt
    self.delibeamterm = (self.lambdaterm - self.lambdatermprev)/top.dt
    self.lambdawaveprev = self.lambdawave 
    self.lambdatermprev = self.lambdaterm 

    self.hlambdawave.append(self.lambdawave)
    self.hlambdaterm.append(self.lambdaterm)

    frz.basegrid.phi[:,:] = fullphi
    fieldsol(-1)

  def advancewaveonestep(self,updatev=None,
                         iwave=None,vwave=None,iterm=None,vterm=None,
                         lbeamloading=None):
    if iwave is None: iwave = self.iwave
    if vwave is None: vwave = self.vwave
    if iterm is None: iterm = self.iterm
    if vterm is None: vterm = self.vterm
    if lbeamloading is None: lbeamloading = self.lbeamloading

    dt = top.dt/self.ncircuit

    tfacw = self.dzwave/sqrt(self.inductance*self.capacitance)/self.vz
    Lw = self.inductance*tfacw
    Cw = self.capacitance*tfacw
    Rw = self.resistance #*tfacw

    tfact = self.dzterm/sqrt(self.terminductance*self.termcapacitance)/self.vz
    Lt = self.terminductance*tfact
    Ct = self.termcapacitance*tfact
    #Rt = self.termresistance*(arange(self.nzterm)/(self.nzterm-1.))**1
    Rt = self.termresistance #*tfact

    # --- Advance of the voltage and current one time step

    # --- Half advance of the current.
    iwave[:] = (iwave + 0.5*dt/Lw*(vwave[:-1] - vwave[1:]))/(1. + 0.5*dt*Rw/Lw)
    iterm[:] = (iterm + 0.5*dt/Lt*(vterm[:-1] - vterm[1:]))/(1. + 0.5*dt*Rt/Lt)

    if updatev is not None: vwave[0] = updatev
    
    vwave[1:-1] = vwave[1:-1] + dt/Cw*(iwave[:-1] - iwave[1:])
    if lbeamloading: vwave[1:-1] = vwave[1:-1] + dt/Cw*self.delibeamwave[1:-1]
    vwave[-1] = vwave[-1] + dt/Cw*(iwave[-1] - iterm[0])
    if lbeamloading: vwave[-1] = vwave[-1] + dt/Cw*self.delibeamwave[-1]
   #vwave[-1] = vwave[-1] + dt/Cw*(iwave[-1] - vwave[-1]/Rw)
    vterm[0] = vwave[-1]
    vterm[1:-1] = vterm[1:-1] + dt/Ct*(iterm[:-1] - iterm[1:])
    if lbeamloading: vterm[1:-1] = vterm[1:-1] + dt/Ct*self.delibeamterm[1:-1]
    # --- Note that vterm[-1] is always zero, at ground.

    # --- Half advance of the current.
    iwave[:] = iwave + 0.5*dt/Lw*(vwave[:-1] - vwave[1:] - Rw*iwave)
    iterm[:] = iterm + 0.5*dt/Lt*(vterm[:-1] - vterm[1:] - Rt*iterm)

  def advancewave(self,n=None,updatev=1,time=None):
    if n is None: n = self.ncircuit

    if self.lbeamloading: self.getbeamloading()
    dt = top.dt/self.ncircuit
    for i in range(n):
      if updatev: vv,ez = self.getvoltage(time+(i+1)*dt,advect=0)
      else: vv = None
      self.advancewaveonestep(updatev=vv)
      if self.lcalcunloaded:
        self.advancewaveonestep(updatev=vv,
                                iwave=self.iwavenoload,vwave=self.vwavenoload,
                                iterm=self.itermnoload,vterm=self.vtermnoload,
                                lbeamloading=0)
      if self.lcalcnowave:
        self.advancewaveonestep(updatev=None,
                                iwave=self.iwavenowave,vwave=self.vwavenowave,
                                iterm=self.itermnowave,vterm=self.vtermnowave,
                                lbeamloading=1)


  def advanceboundaries(self,time=None):
    "Used for model == 2"

    if time is None: time = top.time
    if time < self.ts: return

    self.advancewave(time=time)
    if self.lsavehist:
      self.hv.append(self.vwave[0])
      self.ht.append(time)
      self.hvwave.append(self.vwave)
      self.hvterm.append(self.vterm)
      self.hiwave.append(self.iwave)
      self.hiterm.append(self.iterm)
      if self.lbeamloading:
        if self.lcalcunloaded:
          self.hvwavenoload.append(self.vwavenoload)
          self.hvtermnoload.append(self.vtermnoload)
          self.hiwavenoload.append(self.iwavenoload)
          self.hitermnoload.append(self.itermnoload)
        if self.lcalcnowave:
          self.hvwavenowave.append(self.vwavenowave)
          self.hvtermnowave.append(self.vtermnowave)
          self.hiwavenowave.append(self.iwavenowave)
          self.hitermnowave.append(self.itermnowave)

    vgrid = self.interptogrid(self.vwave,self.vterm)
    setconductorvoltage(vgrid,condid=self.condid)
    self.depositcurrentdensity()

  def interptogrid(self,wave,waveterm,ngp=0):

    # --- Add in wave part of helix
    # --- Interpolate from helix mesh to the field grid
    iz  = (w3d.zmeshlocal + top.zgrid - self.zswave)/self.dzwave
    if ngp:
      wz1 = ones(1+w3d.nz,'d')
      wz0 = zeros(1+w3d.nz,'d')
    else:
      wz0 = (w3d.zmeshlocal + top.zgrid - self.zswave)/self.dzwave - int(iz)
      wz1 = 1. - wz0
    
    # --- Zero out any points beyond the helix mesh
    wz0 = where((0. <= iz) & (iz <= self.nzwave),wz0,0.)
    wz1 = where((0. <= iz) & (iz <= self.nzwave),wz1,0.)
    iz  = where((0. <= iz) & (iz <= self.nzwave),iz,0.)
    iz = int(iz)
    iz = where(iz == len(wave),iz-1,iz) # --- needed for ngp
    izp1 = iz + 1
    izp1 = where(izp1 > len(wave)-1,iz,izp1)

    wgrid = take(wave,iz)*wz1 + take(wave,izp1)*wz0

    # --- Add in termination part of helix
    # --- Interpolate from helix mesh to the field grid
    iz  = (w3d.zmeshlocal + top.zgrid - self.zsterm)/self.dzterm
    if ngp:
      wz1 = ones(1+w3d.nz,'d')
      wz0 = zeros(1+w3d.nz,'d')
    else:
      wz0 = (w3d.zmeshlocal + top.zgrid - self.zsterm)/self.dzterm - int(iz)
      wz1 = 1. - wz0
    
    # --- Zero out any points beyond the helix mesh
    wz0 = where((0. <= iz) & (iz <= self.nzterm),wz0,0.)
    wz1 = where((0. <= iz) & (iz <= self.nzterm),wz1,0.)
    iz  = where((0. <= iz) & (iz <= self.nzterm),iz,0.)
    iz = int(iz)
    iz = where(iz == len(waveterm),iz-1,iz) # --- needed for ngp
    izp1 = iz + 1
    izp1 = where(izp1 > len(waveterm)-1,iz,izp1)

    wgrid = wgrid + take(waveterm,iz)*wz1 + take(waveterm,izp1)*wz0

    # --- Add in wave in start region.
    iz1 = int(max(0.,(self.zs - top.zgrid - w3d.zmmin)/w3d.dz))
    iz2 = int(min(w3d.nz,
              (self.zs + self.startgaplen - top.zgrid - w3d.zmmin)/w3d.dz))
    if iz1 <= iz2:
      if ngp:
        wgrid[iz1:iz2+1] = wave[0]
      else:
        zz = iota(iz1,iz2)*w3d.dz + top.zgrid + w3d.zmmin
        wgrid[iz1:iz2+1] = (zz - self.zs)/self.startgaplen*wave[0]

    return wgrid

  def depositcurrentdensity(self):
    if not self.lcalcbfield: return

    # --- Get current on field solve grid
    igrid = self.interptogrid(self.iwave,self.iterm,ngp=1)
    jzgrid = igrid/(2.*pi*self.ap*w3d.dx)
    jthetagrid = igrid/(w3d.dx*w3d.dz)

    # --- Put current on the grid
    # --- Put it in both bfield and bfieldp just to make sure...
    rr = self.ap/w3d.dx
    ir = int(rr)
    wr = rr - ir
    f3d.bfield.j[1,ir,0,:]    = jthetagrid*(1. - wr)
    f3d.bfieldp.j[1,ir,0,:]   = jthetagrid*(1. - wr)
    f3d.bfield.j[1,ir+1,0,:]  = jthetagrid*wr
    f3d.bfieldp.j[1,ir+1,0,:] = jthetagrid*wr
    f3d.bfield.j[2,ir,0,:]    = jzgrid*(1. - wr)
    f3d.bfieldp.j[2,ir,0,:]   = jzgrid*(1. - wr)
    f3d.bfield.j[2,ir+1,0,:]  = jzgrid*wr
    f3d.bfieldp.j[2,ir+1,0,:] = jzgrid*wr

  def plotvoltage(self,color='fg',scale=1.):
    if self.startgaptype == 0:
      # --- In mode 0, the helix starts at ground
      zz = [self.zs]
      vv = [0.]
    elif self.startgaptype == 1:
      # --- In mode 1, the Ez in the start gap is set
      zz = [self.zs]
      vv = [self.vstart + self.startgaplen*self.ezwave[0]]

    # --- Add interior of helix
    zz = zz + list(self.zwave[:-1])
    vv = vv + list(self.vwave[:-1])

    if self.model == 2:
      zz = zz + list(self.zterm[:-1])
      vv = vv + list(self.vterm[:-1])
    else:
      # --- Add in end of helix interior
      zz.append(self.ze-self.finalgaplen)
      vv.append(self.vwave[-1])

    # --- Add end of helix
    zz.append(self.ze)
    vv.append(0.)

    zz = array(zz)
    vv = array(vv)
    plg(vv*scale,zz,color=color)

  def plotez(self,color='fg',scale=1.):
    zz = [self.zs]
    if self.startgaptype == 0: ezstart = -self.vstart/self.startgaplen
    else:                      ezstart = self.ezwave[0]
    ez = [ezstart]

    # --- Add interior of helix
    zz = zz + list(self.zwave[:-1])
    ez = ez + list(self.ezwave[:-1])

    # --- Add in end of helix interior
    zz.append(self.ze-self.finalgaplen)
    ez.append(self.ezwave[-1])

    # --- Add end of helix
    zz.append(self.ze)
    ez.append(self.vwave[-1]/self.finalgaplen)

    zz = array(zz)
    ez = array(ez)
    plg(ez*scale,zz,color=color)

  def plotvoft(self,ts=0.,te=None,dt=None,color='fg'):
    if te is None: te = self.tlist[-1]
    if dt is None: dt = top.dt
    if self.vvdata is not None: ioftsave = self.ioft
    vvprevioussave = self.vvprevious
    
    self.ioft = 0

    hv = []
    hez = []
    ht = []

    for time in arrayrange(ts,te,dt):
      vv,ez = self.getvoltage(time,advect=0)
      hv.append(vv)
      hez.append(ez)
      ht.append(time)

    plg(hv,ht,color=color)

    if self.vvdata is not None: self.ioft = ioftsave
    self.vvprevious = vvprevioussave

  def plotezoft(self,ts=0.,te=None,dt=None,color='fg'):
    if te is None: te = self.tlist[-1]
    if dt is None: dt = top.dt
    ioftsave = self.ioft
    vvprevioussave = self.vvprevious
    
    self.ioft = 0

    hv = []
    hez = []
    ht = []

    for time in arrayrange(ts,te,dt):
      vv,ez = self.getvoltage(time,advect=0)
      hv.append(vv)
      hez.append(ez)
      ht.append(time)

    plg(hez,ht,color=color)

    self.ioft = ioftsave
    self.vvprevious = vvprevioussave

  def hpcell(self,wave,term,cmin=None,cmax=None):
    ppgeneric(wave[...],xmin=self.ts,xmax=top.time,
              ymin=self.zswave,ymax=self.zewave,
              cmin=cmin,cmax=cmax)
    ppgeneric(term[...],xmin=self.ts,xmax=top.time,
              ymin=self.zsterm,ymax=self.zeterm,
              cmin=cmin,cmax=cmax)

  def hpmplot(self,title,wave,term,nlines=100,navg=0,offset=0.,ltranspose=0):
    if ltranspose:
      wave = transpose(wave)
      term = transpose(term)
    mountainplot1(title,wave,dz=self.dzwave,zmmin=self.zswave,ifvst=0,
                  nlines=nlines,navg=navg,offset=offset)
    mountainplot1(title,term,dz=self.dzterm,zmmin=self.zsterm,ifvst=0,
                  nlines=nlines*self.nzterm/self.nzwave,navg=navg,offset=offset)

  def hpvwave(self,cmin=None,cmax=None,noload=0,nowave=0):
    if not noload and not nowave:
      self.hpcell(self.hvwave,self.hvterm,cmin=cmin,cmax=cmax)
    elif noload:
      self.hpcell(self.hvwavenoload,self.hvtermnoload,cmin=cmin,cmax=cmax)
    elif nowave:
      self.hpcell(self.hvwavenowave,self.hvtermnowave,cmin=cmin,cmax=cmax)

  def ifnone(self,s,a):
    if a is None: return ''
    else:         return s%a

  def printdata(self,name=None):

    text = (
      self.ifnone('name = %s\n',name) +
      self.ifnone('zs = %f\n',self.zs) +
      self.ifnone('ze = %f\n',self.ze) +
      self.ifnone('ap = %f\n',self.ap) +
      self.ifnone('length = %f\n',self.length) +
      self.ifnone('V0 = %e\n',self.V0) +
      self.ifnone('lengthfinal = %f\n',self.lengthfinal) +
      self.ifnone('vz = %e\n',self.vz) +
      self.ifnone('ts = %e\n',self.ts) +
      self.ifnone('zoffset = %f\n',self.zoffset) +
      self.ifnone('startgaplen = %f\n',self.startgaplen) +
      self.ifnone('finalgaplen = %f\n',self.finalgaplen) +
      self.ifnone('startgaptype = %f\n',self.startgaptype) +
      self.ifnone('model = %d\n',self.model) +
      self.ifnone('condid = %d\n',self.condid) +
      self.ifnone('inductance = %e\n',self.inductance) +
      self.ifnone('capacitance = %e\n',self.capacitance) +
      self.ifnone('resistance = %e\n',self.resistance) +
      self.ifnone('terminductance = %e\n',self.terminductance) +
      self.ifnone('termcapacitance = %e\n',self.termcapacitance) +
      self.ifnone('termresistance = %e\n',self.termresistance) +
      self.ifnone('lbeamloading = %d\n',self.lbeamloading) +
      self.ifnone('ncircuit = %d\n',self.ncircuit) +
      self.ifnone('nspirals = %e\n',self.nspirals) +
      self.ifnone('lsavehist = %d\n',self.lsavehist) +
      self.ifnone('lcalcunloaded = %d\n',self.lcalcunloaded) +
      self.ifnone('lcalcnowave = %d\n',self.lcalcnowave)
      )

    plt(text,0.12,0.88,justify="LT")
    fma()

    text = ''
#   if self.ezoft is not None and self.ezduration is not None:
#     text = text + (
#   'ezoft = %s\n'%array2string(array(self.ezoft),max_line_width=60,separator=', ') +
#   'ezduration = %s\n'%array2string(array(self.ezduration),max_line_width=60,separator=', '))

    if self.vvdata is not None:
      ii = min(10,len(self.vvdata))
      text = text + 'vvdata = %s\n'%array2string(array(self.vvdata[:ii]),max_line_width=60,separator=', ')

    if text:
      plt(text,0.12,0.88,justify="LT")
      fma()

