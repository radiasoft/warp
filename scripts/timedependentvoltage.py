from warp import *

class TimeVoltage:
  """
Makes the voltage on a conductor time dependent. One of several methods can be
specified to calculate the voltage as a function of time.
Input for constructor:
 - condid=0: Id (or list of id's) of the conductor which is to be varied
 - discrete=true: z locations for plus/minus z subgrid points are round up/down.
 - setvinject=false: when true, top.vinject is set to the same voltage
 - doitnow=false: when true, applies the voltage and recalculates the fields
                  when object created

 - tieffenback=false: when true, use the Tieffenback profile, with the
                      following parameters. Both maxvoltage and risetime
                      must be specified.
 - minvoltage=0: starting voltage
 - maxvoltage: ending voltage (voltage of flattop)
 - risetime: rise time of the pulse
 - flattime=largepos: flat time, at maximum voltage
 - falltime=largepos: at end of pulse

 - voltdata: sequence of voltage values
 - timedata: times at which the voltages are specified
   Linear interpolation is done between values.

 - voltfunc: user supplied function which takes a single argument, the current
             time. It returns the voltage at that time.
  """

  def __init__(self,condid=0,discrete=true,
                setvinject=false,doitnow=false,
                tieffenback=false,minvoltage=0.,maxvoltage=None,risetime=None,
                flattime=top.largepos,falltime=top.largepos,
                voltdata=None,timedata=None,
                voltfunc=None,
                aftervoltset=None):
    assert ((not tieffenback) or
           (tieffenback and maxvoltage is not None and risetime is not None)),\
           "If Tieffenback profile is used, both a maxvoltage and a risetime must be specified"

    if type(condid) == IntType: self.condid = [condid]
    else:                       self.condid = condid
    self.tieffenback = tieffenback
    self.minvoltage = minvoltage
    self.maxvoltage = maxvoltage
    self.risetime = risetime
    self.flattime = flattime
    self.falltime = falltime
    self.voltdata = voltdata
    self.timedata = timedata
    self.voltfunc = voltfunc
    self.discrete = discrete
    self.setvinject = setvinject
    self.aftervoltset = aftervoltset

    # --- Select method of calculating the voltage
    self.getvolt = None
    if self.tieffenback: self.getvolt = self.tieffenbackvoltage
    elif self.voltdata is not None: self.getvolt = self.voltagefromdata
    elif self.voltfunc is not None: self.getvolt = self.voltfunc
    assert self.getvolt is not None, \
           "At least one method of calculating the voltage must be specified"

    # --- Save history of voltage
    self.hvolt = []
    self.htime = []

    # --- Now, set so the function is called before each field solve
    # --- to apply the voltage
    self.applyvoltage()
    installbeforefs(self.applyvoltage)

    # --- Do it now if requested, both applying the voltage and calculating
    # --- the fields.
    if doitnow: fieldsol(-1)

  def disable(self):
    uninstallbeforefs(self.applyvoltage)

  def applyvoltage(self,time=None):
    # --- AMR is being used, this routine will install itself there to
    # --- be called before the field solve. This is so that the voltages
    # --- on the conductors get set properly for the updated set of
    # --- meshrefined blocks.
    import __main__
    if isinstalledbeforefs(self.applyvoltage):
      if 'AMRtree' in __main__.__dict__:
        __main__.__dict__['AMRtree'].installbeforefs(self.applyvoltage)
        uninstallbeforefs(self.applyvoltage)
        return
    if time is None: time = top.time
    volt = self.getvolt(time)
    # --- Get the appropriate conductors object to install into.
    if 'AMRtree' in __main__.__dict__:
      cond = __main__.__dict__['AMRtree'].getconductors()
    else:
      cond = f3d.conductors
    for c in self.condid:
      setconductorvoltage(volt,c,discrete=self.discrete,
                          setvinject=self.setvinject,
                          conductors=cond)
    self.hvolt.append(volt)
    self.htime.append(time)
    if self.aftervoltset is not None: self.aftervoltset(volt)

  def tieffenbackvoltage(self,time):
    risetime = self.risetime
    flattime = self.flattime
    falltime = self.falltime
    minvoltage = self.minvoltage
    maxvoltage = self.maxvoltage
    if time < risetime:
      tau = time/risetime
      volt = ((maxvoltage - minvoltage)*
              (4.e0/3.e0*tau - 1.e0/3.e0*tau**4) + minvoltage)
    elif risetime <= time and time <= risetime+flattime:
      volt = maxvoltage
    elif risetime+flattime < time < risetime+flattime+falltime:
      tau = (risetime+flattime+falltime-time)/falltime
      volt = ((maxvoltage - minvoltage)*
              (4.e0/3.e0*tau - 1.e0/3.e0*tau**4) + minvoltage)
    elif time >= risetime+flattime+falltime:
      volt = minvoltage
    return volt
    
  def voltagefromdata(self,time):
    if time <= self.timedata[0]: return self.voltdata[0]
    if time >= self.timedata[-1]: return self.voltdata[-1]
    try:
      test = self.index
    except AttributeError:
      self.index = 0
    while self.timedata[self.index+1] < time: self.index = self.index + 1
    wt = ((time - self.timedata[self.index])/
          (self.timedata[self.index+1] - self.timedata[self.index]))
    volt = (self.voltdata[self.index  ]*(1. - wt) +
            self.voltdata[self.index+1]*wt)
    return volt




