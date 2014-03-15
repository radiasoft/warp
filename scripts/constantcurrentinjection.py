"""Sets voltage of source so that a constant current is emitted, emulating Lampel-Tieffenback
procedure, but using the full simulation instead of 1-D approximation.
"""
from warp import *
from timedependentvoltage import TimeVoltage


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
        that w3d.l_inj_rz = true or w3d.l_inj_rz_grid = true.

  Input arguements:
   - sourceid: conductor object or id of the source
   - currentdensity: the current density to be injected
   - sourcevolt: the steady-state voltage on the source or extactor plate,
                 the primary conductor with a time dependent voltage
   - otherids: a list of other conductor objects or ids in the system.
               Conductors that will have a fixed voltage that is not ground
               should be included in this list.
   - othervolts: the steady-state voltages of the other conductors
                 (matching otherids)
   - othercontrolledids: other conductors objects or ids that will have a time
                         dependent voltage that is proporational to the source
                         voltage
   - othercontrolledvolts: the steady-state voltages of the other controlled
                           conductors
   - endplatevolt=0.: voltage on plate at end of diode. In a triode
                      configuration, or a diode with a guard plate, this should
                      be the voltage on the guard plate.
   - l_setvinject=1: when true, top.vinject will be set to the same voltage
                     that the source is set to. This should be true when
                     the primary conductor with a time dependent voltage is
                     the source.
                     If the primary conductor is an extractor plate, this
                     should be false, and top.vinject must be set to the source
                     voltage before this is called.
   - maxvoltagechangerate=None: optional maximum rate at which the source voltage
                                can be varied. If this is not specified, the
                                change in the voltage in unconstrained. This
                                will mostly affect the initial voltage. Without
                                the constraint, the initial voltage will just be
                                whatever it comes out to be. With the
                                constraint, the initial voltage will be relative
                                to the given value, sourcevoltinit.
   - sourcevoltinit=endplatevolt: Initial value of the source voltage. This is
                                  only used when maxvoltagechangerate is
                                  specified.

  This does the equivalent of setting the following parameters
  frz.l_find_rise_time = true
  frz.inj_phi_eq = ekininit - vplate[0]
  frz.v_max = ekininit
  frz.calc_a = 3
    """

    def __init__(self,sourceid,currentdensity,sourcevolt,
                 otherids=[],othervolts=[],
                 othercontrolledids=[],othercontrolledvolts=[],endplatevolt=0.,
                 l_setvinject=1,
                 maxvoltagechangerate=None,sourcevoltinit=None):
        self.currentdensity = currentdensity
        self.sourceid = sourceid
        self.sourcevolt = sourcevolt
        self.otherids = otherids
        self.othervolts = othervolts
        self.othercontrolledids = othercontrolledids
        self.othercontrolledvolts = othercontrolledvolts
        self.endplatevolt = endplatevolt
        self.l_setvinject = l_setvinject
        self.maxvoltagechangerate = maxvoltagechangerate
        if sourcevoltinit is None:
            self.sourcevoltinit = endplatevolt
        else:
            self.sourcevoltinit = sourcevoltinit

        # --- Initialize parameters
        self.hphiref = []
        self.setphiref(currentdensity)
        print 'phiref = ',self.phiref
        if w3d.l_inj_rz or w3d.l_inj_rz_grid:
            self.ww = 2.*pi*iota(0,w3d.inj_nx)*w3d.inj_dx**2*w3d.inj_area[:,0,0]
            self.ww[0] = 0.25*pi*w3d.inj_dx**2
        else:
            self.ww = w3d.inj_dx*w3d.inj_dy*w3d.inj_area[:,:,0]
        self.wwsum = sum(self.ww)

        # --- Calculate the potentials with only the source and controlled
        # --- conductors charged to their steady-state value and everything
        # --- else at endplatevolt. Though everything is shifted down by
        # --- endplatevolt for convenience, since it is assumed that there
        # --- will be more conductors around, but all at ground.
        # --- This is needed so that the other conductors won't influence
        # --- the calculation here, which is supposed to be only the affect
        # --- of the source (and other controlled conductors).
        setconductorvoltage(sourcevolt-self.endplatevolt,condid=self.sourceid)
        for id,v in zip(othercontrolledids,othercontrolledvolts):
            setconductorvoltage(v-self.endplatevolt,condid=id)
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
        if self.l_setvinject: top.vinject = sourcevolt
        top.vinject = top.vinject - self.endplatevolt
        getinj_phi()
        top.vinject = top.vinject + self.endplatevolt
        if w3d.l_inj_rz or w3d.l_inj_rz_grid:
            self.phiv = sum(self.ww*w3d.inj_phi[:,0,0])/self.wwsum
        else:
            self.phiv = sum(self.ww*w3d.inj_phi[:,:,0])/self.wwsum

        # --- Restore the voltage on the other conductors.
        for id,v in zip(otherids,othervolts):
            setconductorvoltage(v,condid=id)

        # --- Setup histories
        self.hsourcevolt = AppendableArray(typecode='d')
        self.hafact = AppendableArray(typecode='d')
        self.hphirho = AppendableArray(typecode='d')
        self.hnp = AppendableArray(typecode='d')
        self.htime = AppendableArray(typecode='d')
        self.hafactunconstrained = AppendableArray(typecode='d')

        # --- Make initial call
        self.setsourcevolt()

        # --- Install setsourcevolt so it is called before a fieldsolve
        # --- Calling setsourcevolt beforefs leaves the voltages on the
        # --- conductors at their correct values.
        installbeforefs(self.setsourcevolt)

    def setsourcevolt(self):
        if top.inject==0: return

        # --- Set the conductor voltages at the endpatevolt so that
        # --- the field solve will not include the diode voltage.
        # --- It will include the potential from the space charge, plus
        # --- any voltage that leaks through the diode end plate.
        setconductorvoltage(self.endplatevolt,condid=self.sourceid)
        for id in self.othercontrolledids:
            setconductorvoltage(self.endplatevolt,condid=id)

        # --- Calculate phirho
        fieldsol(-1)

        if self.l_setvinject:
            top.vinject = self.endplatevolt
        getinj_phi()

        if w3d.l_inj_rz or w3d.l_inj_rz_grid:
            phirho = -sum(self.ww*w3d.inj_phi[:,0,0])/self.wwsum
        else:
            phirho = -sum(self.ww*w3d.inj_phi[:,:,0])/self.wwsum

        self.afact = (self.phiref + phirho)/self.phiv
        self.hafactunconstrained.append(self.afact)

        voltage = (self.endplatevolt +
                       self.afact*(self.sourcevolt-self.endplatevolt))

        # --- The voltage may be constrained by maxvoltagechangerate.
        # --- Note that this modifies self.afact.
        voltage = self.constrainvoltage(voltage)

        solver = getregisteredsolver()
        if solver is None:
            solver = frz.basegrid
        solver.phi += self.afact*self.phisave
        for i in range(1,len(self.blocklists)):
            for c in self.blocklists[i]:
                if c.isactive:
                    c.phi += self.afact*c.phisave

        if self.l_setvinject:
            top.vinject = voltage
        setconductorvoltage(voltage,condid=self.sourceid)
        for id,v in zip(self.othercontrolledids,self.othercontrolledvolts):
            ocvoltage = (self.endplatevolt + voltage*(v-self.endplatevolt)/
                                                (self.sourcevolt-self.endplatevolt))
                           #self.afact*(v-self.endplatevolt))
            setconductorvoltage(ocvoltage,condid=id)

        self.hsourcevolt.append(voltage)
        self.hafact.append(self.afact)
        self.hphirho.append(phirho)
        self.hnp.append(getn())
        self.htime.append(top.time)

    def disable(self):
        if isinstalledbeforefs(self.setsourcevolt):
            uninstallbeforefs(self.setsourcevolt)
            setconductorvoltage(self.sourcevolt,condid=self.sourceid)
            for id,v in zip(self.othercontrolledids,self.othercontrolledvolts):
                setconductorvoltage(v,condid=id)

    def plothist(self,tscale=1.,vscale=1.,tunits='s',vunits='V',title=1):
        pla(self.hsourcevolt[:]*vscale,self.htime[:]*tscale)
        if title:ptitles('','Time ('+tunits+')','Voltage ('+vunits+')')

    def setphiref(self,currentdensity=None):
        if top.inject==0:return
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
            raise Exception("Error in ConstantCurrentRiseTime, top.inject<1 or top.inject>5.")

        self.hphiref.append(self.phiref)

    def constrainvoltage(self,voltage):
        # --- Do nothing if there is no constraint.
        if self.maxvoltagechangerate is None: return voltage

        if len(self.hsourcevolt) == 0:
            # --- Use the initial voltage if this is the first call.
            voltage = self.sourcevoltinit
        else:
            # --- Constrain the voltage to be within the rate of change of the
            # --- previous voltage.
            vprevious = self.hsourcevolt[-1]
            dvmax = self.maxvoltagechangerate*top.dt
            voltage = min(voltage,vprevious+dvmax)
            voltage = max(voltage,vprevious-dvmax)
            #if len(self.hsourcevolt) > 1:
            #  dv = 2*self.hsourcevolt[-1] - self.hsourcevolt[-2]
            #  voltage = max(-2.+dv,min(+2.+dv,voltage))

        # --- Adjust afact so that it would have calculated the constrained
        # --- voltage. This is only needed if there are other controlled voltages.
        self.afact = ((voltage - self.endplatevolt)/
                      (self.sourcevolt-self.endplatevolt))

        return voltage

class SpecifiedCurrentRiseTime(ConstantCurrentRiseTime):
    """
  Same as constantcurrentinjection but allows to specify time dependent current profile.
    """
    def __init__(self,sourceid,currentdensityfunc,sourcevolt,
                 otherids=[],othervolts=[],
                 othercontrolledids=[],othercontrolledvolts=[],endplatevolt=0.,
                 l_setvinject=1,
                 maxvoltagechangerate=None,sourcevoltinit=None):
        self.currentdensityfunc = currentdensityfunc
        ConstantCurrentRiseTime.__init__(self,sourceid,
                                         currentdensityfunc(top.time),sourcevolt,
                                         otherids,othervolts,
                                         othercontrolledids,othercontrolledvolts,
                                         endplatevolt,
                                         l_setvinject,
                                         maxvoltagechangerate,
                                         sourcevoltinit)

    def setsourcevolt(self):
        self.setphiref(self.currentdensityfunc(top.time))
        ConstantCurrentRiseTime.setsourcevolt(self)
