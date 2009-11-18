"""Sets voltage of source so that a constant current is emitted, emulating Lampel-Tieffenback
procedure, but using the full simulation instead of 1-D approximation.
"""
from warp import *
from timedependentvoltage import TimeVoltage

constantcurrentinjection_version = "$Id: constantcurrentinjection.py,v 1.8 2009/11/18 22:09:32 jlvay Exp $"
def constantcurrentinjectiondoc():
  import constantcurrentinjection
  print constantcurrentinjection.__doc__

###########################################################################################
###########################################################################################
class ConstantCurrentRiseTime:
  """
Sets voltage of source so that a constant current is emitted, emulating
the Lampel-Tieffenback procedure, but using the full simulation instead
of 1-D approximation.  

Here is an example if its use...

ccrt = ConstantCurrentRiseTime(sourceid=1,
                               currentdensity=top.ibeam/(pi*top.a0*top.b0),
                               sourcevolt=top.ekin,
                               otherids=[2,3,4,5],
                               othervolts=[100.,200.,150.,0.])

NOTE: This code assumes that RZ geometry is being used and
      that w3d.l_inj_rz = true.

Input arguements:
 - sourceid: conductor id of the source
 - currentdensity: the current density to be injected
 - sourcevolt: the steady-state voltage on the source
 - otherids: a list of all other conductor ids in the system
 - othervolts: a list of the voltages on all other conductors in the system
               (matching otherids)
 - endplatevolt=0.: voltage on plate at end of diode. In a triode
                    configuration, or a diode with a guard plate, this should
                    be the voltage on the guard plate.

This does the equivalent of setting the following parameters
frz.l_find_rise_time = true
frz.inj_phi_eq = ekininit - vplate[0]
frz.v_max = ekininit
frz.calc_a = 3
  """

  def __init__(self,sourceid,currentdensity,sourcevolt,otherids=[],othervolts=[],
                    othercontrolledids=[],othercontrolledvolts=[],endplatevolt=0.,
                    l_setvinject=1):
    self.currentdensity = currentdensity
    self.sourceid = sourceid
    self.sourcevolt = sourcevolt
    self.otherids = otherids
    self.othervolts = othervolts
    self.othercontrolledids = othercontrolledids
    self.othercontrolledvolts = othercontrolledvolts
    self.endplatevolt = endplatevolt
    self.l_setvinject = l_setvinject
    self.hphiref = []

    # --- Initialize parameters
    self.setphiref(currentdensity)
    print 'phiref = ',self.phiref
    if w3d.l_inj_rz:
      self.ww = 2.*pi*iota(0,w3d.inj_nx)*w3d.inj_dx**2*w3d.inj_area[:,0,0]
      self.ww[0] = 0.25*pi*w3d.inj_dx**2
    else:
      self.ww = w3d.inj_dx*w3d.inj_dy*w3d.inj_area[:,:,0]
    self.wwsum = sum(self.ww)

    # --- First set all but source to ground
    setconductorvoltage(sourcevolt-self.endplatevolt,condid=self.sourceid)
    for i,id in enumerate(othercontrolledids):
      setconductorvoltage(self.othercontrolledvolts[i],condid=id)
    for id in otherids:
      setconductorvoltage(0.,condid=id)

    # --- Calculate and save a copy of the fields
    fieldsol(-1)
    solver = getregisteredsolver()
    if solver is None:
      solver = frz.basegrid
    self.phisave = solver.phi.copy()
    try:
      self.blocklists = solver.blocklists
      for i in range(1,len(self.blocklists)):
          for c in self.blocklists[i]:
            if c.isactive:
              c.phisave = c.phi.copy()
    except:
      self.blocklists = []

    # --- Calculate phiv
    if self.l_setvinject:top.vinject = sourcevolt-self.endplatevolt
    getinj_phi()
    if w3d.l_inj_rz:
      self.phiv = sum(self.ww*w3d.inj_phi[:,0,0])/self.wwsum
    else:
      self.phiv = sum(self.ww*w3d.inj_phi[:,:,0])/self.wwsum

    # --- Now reset voltages on plates and set source to voltage of endplatevolt
    setconductorvoltage(self.endplatevolt,condid=self.sourceid)
    for id,v in map(None,otherids,othervolts):
      setconductorvoltage(v,condid=id)
    for id in othercontrolledids:
      setconductorvoltage(self.endplatevolt,condid=id)

    # --- Setup histories
    self.hsourcevolt = AppendableArray(typecode='d')
    self.hafact = AppendableArray(typecode='d')
    self.hphirho = AppendableArray(typecode='d')
    self.hnp = AppendableArray(typecode='d')
    self.htime = AppendableArray(typecode='d')

    # --- Make initial call
    fieldsol(-1)
    self.setsourcevolt()

    # --- Install setsourcevolt so it is called after a fieldsolve
    installafterfs(self.setsourcevolt)

  def setsourcevolt(self):
    if top.inject==0:return
    # --- Calculate phirho
    #fieldsol(-1) # done afterfs
    if self.l_setvinject:top.vinject = self.endplatevolt
    getinj_phi()
    if w3d.l_inj_rz:
      phirho = -sum(self.ww*w3d.inj_phi[:,0,0])/self.wwsum
    else:
      phirho = -sum(self.ww*w3d.inj_phi[:,:,0])/self.wwsum

    self.afact = (self.phiref + phirho)/self.phiv
    solver = getregisteredsolver()
    if solver is None:
      solver = frz.basegrid
    solver.phi += self.afact*self.phisave
    for i in range(1,len(self.blocklists)):
      for c in self.blocklists[i]:
        if c.isactive:
          c.phi += self.afact*c.phisave
    voltage = (self.endplatevolt +
                   self.afact*(self.sourcevolt-self.endplatevolt))
    if self.l_setvinject:
      top.vinject = voltage
    setconductorvoltage(voltage,condid=self.sourceid)
    for id,v in map(None,self.othercontrolledids,self.othercontrolledvolts):
      voltage = (self.endplatevolt +
                     self.afact*(v-self.endplatevolt))
      setconductorvoltage(voltage,condid=id)
    fieldsol(-1)
    getinj_phi()

    # --- Leave the conductor voltages at the endpatevolt so that
    # --- the automatic field solve will not include the diode voltage which
    # --- is added in in this routine.
    setconductorvoltage(self.endplatevolt,condid=self.sourceid)
    for id in self.othercontrolledids:
      setconductorvoltage(self.endplatevolt,condid=id)

    self.hsourcevolt.append(voltage)
    self.hafact.append(self.afact)
    self.hphirho.append(phirho)
    self.hnp.append(getn())
    self.htime.append(top.time)
    import os


  def disable(self):
    if isinstalledafterfs(self.setsourcevolt):
      uninstallafterfs(self.setsourcevolt)
      setconductorvoltage(self.sourcevolt,condid=self.sourceid)
      for id,v in map(None,self.othercontrolledids,self.othercontrolledvolts):
        setconductorvoltage(v,condid=id)

  def plothist(self,tscale=1.,vscale=1.,tunits='s',vunits='V',title=1):
    pla(self.hsourcevolt[:]*vscale,self.htime[:]*tscale)
    if title:ptitles('','Time ('+tunits+')','Voltage ('+vunits+')')
    
  def setphiref(self,currentdensity=None):
    if top.inject==0:return
    if currentdensity is None:
      self.currentdensity = self.currentdensityfunc(top.time)
    else:
      self.currentdensity = currentdensity

    chi = 4*eps0/9*sqrt(2*top.pgroup.sq[0]/top.pgroup.sm[0])
    d = top.inj_d[0]*w3d.inj_dz
    if top.inject in [2,3]:
      # --- space-charge limited injection
      self.phiref = (self.currentdensity/chi)**(2./3.)*d**(4./3.)
    elif top.inject in [4,5]:
      # --- source limited injection
      kT = top.boltzmann*top.tempinject[0]
      A0 = ((4.*pi*top.pgroup.sm[0]*top.boltzmann**2*top.pgroup.sq[0])/(top.planck**3)) * 0.5
      def currentdensfunc(V=0.,J0=0.): 
        te_const = A0 * top.tempinject[0]**2 
        dw = sqrt((echarge**3*V/d)/(4.*pi*eps0))
        J_s = te_const * exp(-(top.workinject-dw)/kT)
        if top.inject==4:
          return J_s
        elif top.inject==5:
          J_cl = chi/d**2 * V**1.5
          return J_cl * (1.-exp(-J_s/J_cl))
      self.phiref = bisection(currentdensfunc,1.e-10,1.e8,f0=self.currentdensity)
      self.currentdensfunc = currentdensfunc
    else:
      raise("Error in ConstantCurrentRiseTime, top.inject<1 or top.inject>5.")

    self.hphiref.append(self.phiref)
    
class SpecifiedCurrentRiseTime(ConstantCurrentRiseTime):
  """
Same as constantcurrentinjection but allows to specify time dependent current profile.
  """
  def __init__(self,sourceid,currentdensityfunc,sourcevolt,otherids=[],othervolts=[],
                    othercontrolledids=[],othercontrolledvolts=[],endplatevolt=0.,
                    l_setvinject=1):
    ConstantCurrentRiseTime.__init__(self,sourceid,currentdensityfunc(top.time),sourcevolt,
                                          otherids,othervolts,
                                          othercontrolledids,othercontrolledvolts,endplatevolt,
                                          l_setvinject)
    self.currentdensityfunc = currentdensityfunc
    installbeforefs(self.setphiref)
    self.hphiref = AppendableArray(typecode='d')
    

