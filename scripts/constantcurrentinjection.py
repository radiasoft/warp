"""Sets voltage of source so that a constant current is emitted, emulating Lampel-Tieffenback
procedure, but using the full simulation instead of 1-D approximation.
"""
from warp import *
from timedependentvoltage import TimeVoltage

constantcurrentinjection_version = "$Id: constantcurrentinjection.py,v 1.1 2005/02/16 19:51:22 dave Exp $"
def constantcurrentinjectiondoc():
  import constantcurrentinjection
  print constantcurrentinjection.__doc__

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

  def __init__(self,sourceid,currentdensity,sourcevolt,otherids,othervolts,
                    endplatevolt=0.):
    self.sourceid = sourceid
    self.otherids = otherids
    self.currentdensity = currentdensity
    self.sourcevolt = sourcevolt
    self.othervolts = othervolts
    self.endplatevolt = endplatevolt

    # --- Initialize parameters
    chi = 4*eps0/9*sqrt(2*top.sq[0]/top.sm[0])
    self.phiref = (currentdensity/chi)**(2./3.)*(top.inj_d[0]*w3d.inj_dz)**(4./3.)
    self.ww = 2.*pi*iota(0,w3d.inj_nx)*w3d.inj_dx**2*w3d.inj_area[:,0,0]
    self.ww[0] = 0.25*pi*w3d.inj_dx**2
    self.wwsum = sum(self.ww)

    # --- First set all but source to ground
    setconductorvoltage(sourcevolt-self.endplatevolt,condid=self.sourceid)
    for id in otherids:
      setconductorvoltage(0.,condid=id)

    # --- Calculate and save a copy of the fields
    fieldsol(-1)
    self.phisave = frz.basegrid.phi + 0.

    # --- Calculate phiv
    top.vinject = sourcevolt-self.endplatevolt
    getinj_phi()
    self.phiv = sum(self.ww*w3d.inj_phi[:,0,0])/self.wwsum

    # --- Now reset voltages on plates and set source to ground
    setconductorvoltage(self.endplatevolt,condid=self.sourceid)
    for id,v in map(None,otherids,othervolts):
      setconductorvoltage(v,condid=id)

    # --- Setup histories
    self.hsourcevolt = []
    self.hafact = []
    self.hphirho = []

    # --- Make initial call
    fieldsol(-1)
    self.setsourcevolt()

    # --- Install setsourcevolt so it is called after a fieldsolve
    installafterfs(self.setsourcevolt)

  def setsourcevolt(self):
    # --- Calculate phirho
    #fieldsol(-1) # done afterfs
    top.vinject = self.endplatevolt
    getinj_phi()
    phirho = -sum(self.ww*w3d.inj_phi[:,0,0])/self.wwsum

    self.afact = (self.phiref + phirho)/self.phiv
    frz.basegrid.phi[:,:] = frz.basegrid.phi + self.afact*self.phisave
    top.vinject = (self.endplatevolt +
                   self.afact*(self.sourcevolt-self.endplatevolt))
    setconductorvoltage(top.vinject[0],condid=self.sourceid)
    fieldsol(-1)
    getinj_phi()

    # --- Leave the conductor voltages at the endpatevolt so that
    # --- the automatic field solve will not include the diode voltage which
    # --- is added in in this routine.
    setconductorvoltage(self.endplatevolt,condid=self.sourceid)

    self.hsourcevolt.append(top.vinject[0])
    self.hafact.append(self.afact)
    self.hphirho.append(phirho)

  def disable(self):
    uninstallafterfs(self.setsourcevolt)
    setconductorvoltage(self.sourcevolt,condid=self.sourceid)
