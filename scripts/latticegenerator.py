from warp import *
from lattice import *
import cPickle
latticegenerator_version = "$Id: latticegenerator.py,v 1.9 2009/02/05 18:18:40 dave Exp $"
######################################################################
# Lattice builder
# 
######################################################################

class LatticeGenerator:
  """Generates a lattice
  ion_mass:
  beam_duration:
  ekinmid:
  charge_per_beam:
  gap_len=.01:
  ngappoints=1000: Number of slices in the beam
  nendpoints=ngappoints/10: Additional points added to end of acclet array
  nhlpswithoutgaps=0:
  loadandfire=1:
  firetime=0.:
  risetime=None:
  unfireandunload=1:
  stoptime=1.e36:
  falltime=None:
  accel_gradient=0.:
  firstquadsign=+1:
  lattice=None: An existing lattice can be passed in
  nhist=1:
  zhist=1:
  lverbose=1:
  maxgapgradient=1.e6:
  luservgap=None:
  straight=0.8: Fraction of the beam which is not ends (used in Eears calc)
  icharge=4: CIRCE charge model to use in calculation of Ez for Eears
  lfixed=false: CIRCE flag for using newer version of code in getezbeam
  """
  def __init__(s,ion_mass,beam_duration,ekinmid,charge_per_beam,gap_len=.01,
               ngappoints=1000,nendpoints=None,
               nhlpswithoutgaps=0,
               loadandfire=1,firetime=0.,risetime=None,
               unfireandunload=1,stoptime=1.e36,falltime=None,
               accel_gradient=0.,firstquadsign=+1,
               lattice=None,uselattice=0,
               nhist=1,zhist=1,lverbose=1,
               maxgapgradient=1.e6,luservgap=None,
               straight=0.8,icharge=4,lfixed=false,setears=1):
    s.ion_mass = ion_mass
    s.beam_duration = beam_duration
    s.ekinmid = ekinmid
    s.charge_per_beam = charge_per_beam
    s.gap_len = gap_len
    s.ngappoints = ngappoints
    if nendpoints is None: s.nendpoints = ngappoints/10
    else: s.nendpoints = nendpoints
    s.ntaccl = s.nendpoints + s.ngappoints + s.nendpoints
    s.nhlpswithoutgaps = nhlpswithoutgaps
    s.loadandfire = loadandfire
    s.firetime = firetime
    s.risetime = risetime
    s.unfireandunload = unfireandunload
    s.stoptime = stoptime
    s.falltime = falltime
    s.accel_gradient = accel_gradient
    s.firstquadsign = firstquadsign
    s.lattice = lattice
    s.uselattice = uselattice
    if s.lattice is not None: s.uselattice = 1
    s.nhist = nhist
    s.zhist = zhist
    s.lverbose = lverbose
    s.maxgapgradient = maxgapgradient
    s.luservgap = None
    s.straight = straight
    s.icharge = icharge
    s.lfixed = lfixed
    s.setears = setears
    # --- element counters
    s.iquad = 0
    s.iaccl = 0
    s.idrft = 0
    # --- Setup element attributes
    if not s.uselattice:
      s.latlast = 0.
      s.drftzs = []
      s.drftze = []
      s.drftap = []
      s.helezs = []
      s.heleze = []
      s.heleae = []
      s.heleam = []
      s.heleap = []
      s.acclzs = []
      s.acclze = []
      s.acclts = []
      s.accldt = []
      s.acclet = []
      s.acclap = []
    # --- Setup beam
    s.nn = s.ngappoints
    s.nno2 = s.ngappoints/2
    s.itime = 1.*iota(0,s.nn)/s.nn
    s.beamtime = s.itime*s.beam_duration
    s.vzbeam = s.ekintov(s.ekinmid)*ones(s.nn+1,'d')
    s.vzmid = s.vzbeam[s.nno2]
    s.ekin = s.ekinmid*ones(s.nn+1,'d')
    s.chargemid = 1.1*s.charge_per_beam/(s.nn+1)
    s.earsinit()

  #-----------------------------------------------------------------------
  # --- More initialization that requires the user defined routines
  def moreinit(s):

    # --- Make sure that this routine is only called once.
    try:
      s.moreinitcalled = s.moreinitcalled + 1
    except AttributeError:
      s.moreinitcalled = 0
    if s.moreinitcalled: return

    # --- Setup first quad
    s.linechgmid = s.chargemid/ \
                  (0.5*(s.beamtime[s.nno2+1]-s.beamtime[s.nno2-1])*s.vzmid)
    s.perveancemid = s.linechgmid/(4.*pi*eps0*s.ekinmid)
    s.hlp = s.amean(s)*sqrt((1.-cos(s.sigma(s)*pi/180.))/(2.*s.perveancemid))
    s.zlast = 0.
    s.ihlp = 0

    #for iq in range(s.nhlpswithoutgaps):
    s.adddrft((1.-s.occupancy(s))*s.hlp/2.)
    s.addquad(s.zlast,s.hlp,s.ekinmid)
    s.adddrft((1.-s.occupancy(s))*s.hlp/2.-s.gap_len/2.)

    # --- Load and fire risetime default is transit time of one lattice
    # --- period
    if s.loadandfire and not s.risetime:
      s.risetime = s.hlp/s.vzmid
    if s.unfireandunload and not s.falltime:
      s.falltime = s.risetime

    # --- Setup history arrays
    s.savehist(linit=1)

    # --- Get beam_duration before first gap
    s.advancebeam(len=1.) #+s.nhlpswithoutgaps)
    s.zlast = s.zlast + s.hlp*(1.) #+s.nhlpswithoutgaps)
    s.savehist()

  #-----------------------------------------------------------------------
  def smoother(s):
    if s <= 0:
      return 0.
    elif s < 1.:
  #    return 0.5*(1. - tanh(cos(pi*s)/sin(pi*s)))
       return s
    else:
      return 1.

  #---------------------------------------------------------------------------
  # --- Write data out to a pickle file
  def saveself(s,filename):
    ff = open(filename,'w')
    cPickle.dump(s,ff,1)
    ff.close()

  #---------------------------------------------------------------------------
  # --- Print out a line of info
  def oneliner(s):
    print "hlp %d %6.4f beta %6.4f V %6.1f tau %6.4f tilt %4.2f"%(
           s.ihlp,s.hlp,s.vzmid/clight,
           s.ekinmid/1.e6,s.beam_duration*1.e6,
           (s.vzbeam[-1]-s.vzbeam[0])/s.vzbeam[s.nno2])

  #---------------------------------------------------------------------------
  # --- Convenienct functions
  def ekintov(s,ekin):
    ke = jperev*ekin/(s.ion_mass*amu*clight**2)
    gammabar = 1. + ke
    return clight*sqrt((2*ke+ke**2)/gammabar**2)
  def vtoekin(s,v):
    u = (v/clight)**2
    gammabar = 1./sqrt(1. - u)
    return (s.ion_mass*amu*clight**2)*(u/(sqrt(1.-u)+1.-u))/jperev
  def vtogamma(s,v):
    u = (v/clight)**2
    return 1./sqrt(1. - u)
  def ekintogamma(s,ekin):
    ke = jperev*ekin/(s.ion_mass*amu*clight**2)
    return 1. + ke

  #-----------------------------------------------------------------------
  def advancebeam(s,gapez=None,len=1.):
    if gapez:
      accldt = (s.beamtime[-1] - s.beamtime[0])/s.ngappoints
      acclts = s.beamtime[0] - s.nendpoints*accldt
      ia = ((s.beamtime - acclts)/accldt).astype(long)
      wi =  (s.beamtime - acclts)/accldt - ia
      vv = take(gapez,ia)*(1. - wi) + take(gapez,ia+1)*wi
      s.ekin = s.ekin + vv*s.gap_len
      s.vzbeam = s.ekintov(s.ekin)
      s.ekinmid = s.ekin[s.nno2]
      s.vzmid = s.vzbeam[s.nno2]
    s.beamtime[:] = s.beamtime + len*s.hlp/s.vzbeam
    s.beam_duration = s.beamtime[-1] - s.beamtime[0]

  #---------------------------------------------------------------------------
  # Adds new drift to lattice. zlast is center of drift just before the new
  # quad.
  def adddrft(s,l):
    s.idrft = s.idrft + 1
    if s.uselattice:
      if s.lattice:
        s.lattice = s.lattice + Drft(l=l,ap=s.aperture(s))
      else:
        s.lattice = Drft(l=l,ap=s.aperture(s))
    else:
      s.drftzs.append(s.latlast)
      s.latlast = s.latlast + l
      s.drftze.append(s.latlast)
      s.drftap.append(s.aperture(s))

  #---------------------------------------------------------------------------
  # Adds new quad to lattice.
  def addquad(s,zlast,hlp,ekin):
    s.iquad = s.iquad + 1
    dedx = 2.*ekin/hlp**2*sqrt(2.*(1.-cos(s.sigma(s)*pi/180.))/
                               (s.occupancy(s)**2*(1.-2./3.*s.occupancy(s))))
    if s.lmagnetic(s):
      quadsign = +(1-2*(s.iquad%2))*s.firstquadsign
      vz = s.ekintov(ekin)
      quad_stren = dedx*quadsign/(vz*s.vtogamma(vz))
    else:
      quadsign = -(1-2*(s.iquad%2))*s.firstquadsign
      quad_stren = dedx*quadsign
    if s.uselattice:
      if s.lmagnetic(s):
        quad = Hele(l=s.occupancy(s)*hlp,nn=[2],vv=[0],am=[quad_stren],
                    ap=s.aperture(s))
      else:
        quad = Hele(l=s.occupancy(s)*hlp,nn=[2],vv=[0],ae=[quad_stren],
                    ap=s.aperture(s))
      s.lattice = s.lattice + quad
    else:
      s.helezs.append(s.latlast)
      s.latlast = s.latlast + s.occupancy(s)*hlp
      s.heleze.append(s.latlast)
      s.heleap.append(s.aperture(s))
      if s.lmagnetic(s):
        s.heleae.append(0.)
        s.heleam.append(quad_stren)
      else:
        s.heleae.append(quad_stren)
        s.heleam.append(0.)

  #---------------------------------------------------------------------------
  # Add new accelerating gap
  def addgap(s,zlast,gapez):
    s.iaccl = s.iaccl + 1
    accl_dt = (s.beamtime[-1] - s.beamtime[0])/s.ngappoints
    accl_ts = s.beamtime[0] - s.nendpoints*accl_dt
    s.geteears(accl_ts,accl_dt)
    if s.uselattice:
      s.accl = Accl(l=s.gap_len,dt=accl_dt,ts=accl_ts,et=gapez+s.eears,
                    ap=s.aperture(s))
      s.lattice = s.lattice + s.accl
    else:
      s.acclzs.append(s.latlast)
      s.latlast = s.latlast + s.gap_len
      s.acclze.append(s.latlast)
      s.acclap.append(s.aperture(s))
      s.acclts.append(accl_ts)
      s.accldt.append(accl_dt)
      s.acclet.append(gapez+s.eears)

  #-----------------------------------------------------------------------
  # --- Add more to the lattice until the endcondition is met
  def generate(s):

    s.moreinit()

    if s.loadandfire:
      s.laflattice()
      return

    while not s.endcondition(s):
      if s.lverbose: s.oneliner()
      s.ihlp = s.ihlp + 1
  
      # --- parameters just after gap
      gapvoltage = s.midpulsegapvoltage(s)
      if gapvoltage is not None:
        s.ekinmid = s.ekin[s.nno2] + gapvoltage*s.hlp
      s.vzmid = s.ekintov(s.ekinmid)

      # --- Calculate hlp length from info at start of hlp
      # --- Uses beam_duration at end of previous hlp and energy after gap
      s.linechgmid = s.chargemid/ \
                    (0.5*(s.beamtime[s.nno2+1]-s.beamtime[s.nno2-1])*s.vzmid)
      s.perveancemid = s.linechgmid/(4.*pi*eps0*s.ekinmid)
      s.hlp = s.amean(s)*sqrt((1.-cos(s.sigma(s)*pi/180.))/(2.*s.perveancemid))

      if gapvoltage is None:
        # --- With no accelerator gap...
        # --- Note that first drift is long enough to cover the missing gap
        s.adddrft((1.-s.occupancy(s))*s.hlp/2.+s.gap_len/2.)
        s.addquad(s.zlast,s.hlp,s.ekinmid)
        s.adddrft((1.-s.occupancy(s))*s.hlp/2.-s.gap_len/2.)

        # --- Advance the beam
        s.advancebeam(gapez=None)

      else:
        # --- beam length at end of this hlp (tilt and vzmid have changed)
        oldbeam_duration = s.beam_duration
        s.beam_duration = s.beam_duration + \
                      s.hlp/s.vzmid*(1./(1.+s.tilt(s)/2)-1./(1.-s.tilt(s)/2.))
        # --- Get time of arrival of beam slices at next gap
        nextbeamtime = s.beamtime[s.nno2]+s.hlp/s.vzmid+(s.itime-0.5)*s.beam_duration
        # --- This expression preserves the time dependence.
        #nextbeamtime = (s.beamtime - s.beamtime[s.nno2])* \
        #               s.beam_duration/oldbeam_duration + \
        #               s.beamtime[s.nno2] + s.hlp/s.vzmid

        # --- Get velocity needed and convert to energy
        nextvz = s.hlp/(nextbeamtime-s.beamtime)
        ekinnext = s.vtoekin(nextvz)
        vgap = ekinnext - s.ekin

        # --- Check for a user supplied vgap
        if s.luservgap: vgap = s.getvgap(s)

        # --- Find params for next gap
        s.gapez = zeros(s.ntaccl+1,'d')
        s.gapez[:s.nendpoints]=vgap[0]
        s.gapez[s.nendpoints:-s.nendpoints] = vgap
        s.gapez[-s.nendpoints:]=vgap[-1]
        s.gapez[:] = s.gapez/s.gap_len

        # --- Scale the accelerating field to both keep it less than the
        # --- maximum gap gradient and to keep it positive. Scaling it (instead
        # --- of simply limiting it) preserves the self-similarity of the beam
        # --- (since otherwise some part of the beam isn't being compressed as
        # --- much as others if the field is simply pinned to the max value
        # --- there).
        # --- This scaling supplies a limit to how fast a tilt can be applied
        # --- as well as how fast the beam can be accelerated.
        if max(s.gapez) > s.maxgapgradient*s.hlp/s.gap_len and \
           min(s.gapez) > 0.:
          s.gapez[:] = s.gapez/max(s.gapez)*s.maxgapgradient*s.hlp/s.gap_len
        elif min(s.gapez) < 0.:
          gapezmax = min(max(s.gapez)-min(s.gapez),
                          s.maxgapgradient*s.hlp/s.gap_len)
          if gapezmax != 0.:
            s.gapez[:] = (s.gapez-min(s.gapez))/(max(s.gapez)-min(s.gapez))* \
                          gapezmax
          else:
            s.gapez[:] = 0.

        # --- Now add in linear ramp before and after pulse. This is done
        # --- checking the limits above so that the linear ramp does not
        # --- effect the scaling.
        s.gapez[:s.nendpoints] = (s.gapez[0] - (s.gapez[1] - s.gapez[0])*
                                               iota(s.nendpoints,1,-1))
        s.gapez[-s.nendpoints:] = (s.gapez[-1] + (s.gapez[-1] - s.gapez[-2])*
                                                 iota(s.nendpoints))

        # --- Generate lattice period
        s.addgap(s.zlast,s.gapez)
        s.adddrft((1.-s.occupancy(s))*s.hlp/2.-s.gap_len/2.)
        s.addquad(s.zlast,s.hlp,s.ekinmid + s.gapez[s.nno2]*s.gap_len/2.)
        s.adddrft((1.-s.occupancy(s))*s.hlp/2.-s.gap_len/2.)

        # --- Advance the beam
        s.advancebeam(s.gapez)


      s.zlast = s.zlast + s.hlp

      # --- Reset beamtime
      # --- This is needed since the beam goes unstable, producing
      # --- growing oscillations in beamtime and then in gapez and vzbeam.
      # --- This is realy only a hack, making the assumption that self-similar
      # --- compression is done.
      #s.beamtime[:] = s.itime*s.beam_duration + s.beamtime[0]

      s.savehist()

  #-----------------------------------------------------------------------
  def savehist(s,linit=0):
    # --- Save histories
    if linit:
      s.hbeam_duration = [s.beam_duration]
      s.hhlp = [s.hlp]
      s.hekinmid = [s.ekinmid]
      s.hvzmid = [s.vzmid]
      s.htiltinput = [s.tilt(s)]
      s.htilt = [(s.vzbeam[-1]-s.vzbeam[0])/s.vzbeam[s.nno2]]
      s.hzlast = [s.zlast]
      s.hlinechgmid = [s.linechgmid]
      s.hbeamtime = [s.beamtime[::s.zhist]+0.]
      s.hvzbeam = [s.vzbeam[::s.zhist]+0.]
    elif s.ihlp%s.nhist == 0:
      s.hbeam_duration.append(s.beam_duration)
      s.hhlp.append(s.hlp)
      s.hekinmid.append(s.ekinmid)
      s.hvzmid.append(s.vzmid)
      s.htiltinput.append(s.tilt(s))
      s.htilt.append((s.vzbeam[-1]-s.vzbeam[0])/s.vzbeam[s.nno2])
      s.hzlast.append(s.zlast)
      s.hlinechgmid.append(s.linechgmid)
      s.hbeamtime.append(s.beamtime[::s.zhist]+0.)
      s.hvzbeam.append(s.vzbeam[::s.zhist]+0.)

  #-----------------------------------------------------------------------
  # --- Do some fine work
  def finalize(s):
    # --- Install the lattice into the WARP database (the Lattice group)
   #if s.uselattice:
   #  madtowarp(s.lattice)
   #else:
    if not s.uselattice:
      top.ndrft = len(s.drftzs) - 1
      top.nhele = len(s.helezs) - 1
      top.nhmlt = 1
      top.naccl = len(s.acclzs) - 1
      top.ntaccl = s.ntaccl
      gchange("Lattice")
      top.drftzs[:] = s.drftzs
      top.drftze[:] = s.drftze
      top.drftap[:] = s.drftap
      top.helezs[:] = s.helezs
      top.heleze[:] = s.heleze
      top.heleap[:] = s.heleap
      top.heleae[0,:] = s.heleae
      top.heleam[0,:] = s.heleam
      top.hele_n[0,:] = 2
      top.hele_v[0,:] = 0
      top.acclzs[:] = s.acclzs
      top.acclze[:] = s.acclze
      top.acclap[:] = s.acclap
      top.acclts[:] = s.acclts
      top.accldt[:] = s.accldt
      top.acclet[:,:] = transpose(array(s.acclet))
    # --- Convert histories into arrays
    s.hbeam_duration = array(s.hbeam_duration)
    s.hhlp = array(s.hhlp)
    s.hekinmid = array(s.hekinmid)
    s.hvzmid = array(s.hvzmid)
    s.htiltinput = array(s.htiltinput)
    s.htilt = array(s.htilt)
    s.hzlast = array(s.hzlast)
    s.hlinechgmid = array(s.hlinechgmid)
    s.hbeamtime = transpose(array(s.hbeamtime))
    s.hvzbeam = transpose(array(s.hvzbeam))
    s.heears = array(s.heears)
    s.hlinechg = array(s.hlinechg)


  #-----------------------------------------------------------------------
  #-----------------------------------------------------------------------
  #-----------------------------------------------------------------------
  #-----------------------------------------------------------------------
  #---    Load and Fire   ------------------------------------------------
  #-----------------------------------------------------------------------

  #-----------------------------------------------------------------------
  # --- Ion energy increases exponentially with z.
  def laftt(s,z):
    result = s.ekinmid*exp(s.accel_gradient/s.ekinmid*z)
    return result

  #-----------------------------------------------------------------------
  # --- Get vz from the position or from the ion energy
  def lafvz(s,z,t=None):
    m = s.ion_mass*amu*clight**2/jperev
    if t is None:
      return sqrt(s.laftt(z)*(2.*m+s.laftt(z)))/(s.laftt(z)+m)*clight
    else:
      return sqrt(t*(2.*m+t))/(t+m)*clight

  #-----------------------------------------------------------------------
  def lafgetvgrad(s,dt=1.e-9,nsteps=20000):
    # --- Integrate to get position as a function of time.
    # --- Step size is hard wired to 1.e-9, small compared to the beam duration.
    # --- Enough points are taken to get beyond the electrostatic section.
    s.gapdt = dt
    zz = zeros(nsteps,'d')
    for i in xrange(nsteps-1):
      zz[i+1] = zz[i] + s.gapdt*s.lafvz(z=zz[i])
    # --- Now calculate voltage gradient at each of the z points, producing
    # --- the gradient as a funtion of time.
    s.vgrad = s.accel_gradient*exp(s.accel_gradient/s.ekinmid*zz)

  #-----------------------------------------------------------------------
  # --- Interpolates from vgrad to gapez
  def lafgetgapez(s):
    accldt = (s.beamtime[-1] - s.beamtime[0])/s.ngappoints
    acclts = s.beamtime[0] - s.nendpoints*accldt
    actime = acclts + iota(0,s.ntaccl)*accldt
    ii = ((actime - s.firetime)/s.gapdt).astype(long)
    ww = ((actime - s.firetime)/s.gapdt)- ii
    ww = where(less(ii,0),0.,ww)
    ii = where(less(ii,0),0,ii)
    gapez = take(s.vgrad,ii)*(1.-ww) + take(s.vgrad,ii+1)*ww
    # --- Set maximum accelerating gradient
    gapez[:] = where(greater(gapez,s.maxgapgradient),s.maxgapgradient,gapez)
    # --- Convert voltage gradient to gap Ez
    gapez[:] = gapez[:]*s.hlp/s.gap_len
    # --- If risetime > 0., add in linear rise
    if s.risetime > 0. and acclts < s.firetime:
      firstgrad = gapez[-1]
      j = 0
      while acclts + j*accldt < s.risetime+s.firetime:
        j = j + 1
      firstgrad = gapez[j]
      j = 0
      while acclts + j*accldt < s.risetime+s.firetime:
        if (acclts + j*accldt) > s.firetime:
          t = (acclts + j*accldt - s.firetime)/s.risetime
          gapez[j] = s.riseprofile(t)*firstgrad
        j = j + 1
      # --- Zero out gap field for times < firetime
      gapez[:] = where(less(actime,s.firetime),0.,gapez[:])
    # --- If falltime is set, add in linear fall
    if s.falltime and acclts+s.ntaccl*accldt > s.stoptime:
      lastgrad = gapez[0]
      j = s.ntaccl
      while acclts + j*accldt > s.stoptime:
        j = j - 1
      lastgrad = gapez[j]
      j = s.ntaccl
      while acclts + j*accldt > s.stoptime:
        if (acclts + j*accldt) < s.stoptime+s.falltime:
          gapez[j]=(s.stoptime+s.falltime-acclts-j*accldt)/s.falltime*lastgrad
        j = j - 1
      # --- Zero out gap field for times > stoptime+falltime
      gapez[:] = where(greater(actime,s.stoptime+s.falltime),0.,gapez[:])
    return gapez

  def laflattice(s):
    #===========================================================================
    # Coming out of the injector, the beam is initially allowed to drift
    # until the whole beam is in the accelerator. At that point the
    # accelerating modules are ramped up (load-and-fire). The lattice is set
    # up so that the lattice period is initially constant (with parameters
    # giving a matched beam at 1.6 MeV). The length of that region is one half
    # of the initial beam length so that when the middle of the beam is at the
    # end of the constant hlp section, the acceleration starts. The changing
    # hlp parameters then match the energy of the beam center.

    # --- First, setup arrays for the accelerating voltage
    s.lafgetvgrad()

    # --- Generate initial drifting region
    while not s.endcondition(s) or \
          (s.unfireandunload and s.beamtime[0] < s.stoptime+s.falltime):
      if s.lverbose: s.oneliner()
      s.ihlp = s.ihlp + 1
      s.gapez = s.lafgetgapez()
      s.perveancemid = s.linechgmid/(4.*pi*eps0*s.ekinmid)
      s.hlp = s.amean(s)*sqrt((1.-cos(s.sigma(s)*pi/180.))/(2.*s.perveancemid))
      s.addgap(s.zlast,s.gapez)
      s.adddrft((1.-s.occupancy(s))*s.hlp/2.-s.gap_len/2.)
      s.addquad(s.zlast,s.hlp,s.ekinmid + s.gapez[s.nno2]*s.gap_len/2.)
      s.adddrft((1.-s.occupancy(s))*s.hlp/2.-s.gap_len/2.)
      s.advancebeam(s.gapez)
      s.zlast = s.zlast + s.hlp
      s.savehist()

    # --- Add a final drift which fills the space of what otherwise would be
    # --- a gap.
    #s.adddrft(s.gap_len/2.)

  #-----------------------------------------------------------------------
  #-----------------------------------------------------------------------
  #-----------------------------------------------------------------------
  #-----------------------------------------------------------------------
  # --- These are the functions that need to be defined by the user
  def endcondition(s,s1):
      raise "The function endcondition needs to be defined"
  def midpulsegapvoltage(s,s1):
      raise "The function midpulsegapvoltage needs to be defined"
  def tilt(s,s1):
      raise "The function tilt needs to be defined"
  def occupancy(s,s1):
      raise "The function occupancy needs to be defined"
  def amean(s,s1):
      raise "The function amean needs to be defined"
  def sigma(s,s1):
      raise "The function sigma needs to be defined"
  def aperture(s,s1):
      raise "The function aperture needs to be defined"
  def lmagnetic(s,s1):
      raise "The function lmagnetic needs to be defined"
  def riseprofile(s,t):
      return t
  def earlength(s,s1):
      return s.hlp

  def __getstate__(s):
    # --- This is needed since the user supplied subroutines cannot be
    # --- pickled properly.
    odict = s.__dict__.copy()
    try: del odict['endcondition']
    except: pass
    try: del odict['midpulsegapvoltage']
    except: pass
    try: del odict['tilt']
    except: pass
    try: del odict['occupancy']
    except: pass
    try: del odict['amean']
    except: pass
    try: del odict['sigma']
    except: pass
    try: del odict['aperture']
    except: pass
    try: del odict['lmagnetic']
    except: pass
    return odict

  ######################################################################
  # --- Routines to get the ear fields
  def earsinit(s):
    if not s.setears: return
    # --- Allocate temporary arrays used by getezbeam
    cir.nittmp = s.nn+1
    gchange("CIRtmp")
    # --- Create arrays used by geteears
    s.cirvar = fones((16,s.nn+1),'d')
    s.heears = []
    s.hlinechg = []
    s.ekintemp = zeros(s.nn+1,'d')
    s.linechg = zeros(s.nn+1,'d')
    s.ezbeam = zeros(s.nn+1,'d')
    s.eears = zeros(s.ntaccl+1,'d')
    # --- Create a charge profile with parabolic falloff in the ends
    s.endlen = int((1.-s.straight)/2.*(s.nn+1))
    s.charge = ones(s.nn+1,'d')
    s.charge[s.endlen::-1] = (1. - (1.*iota(0,s.endlen)/s.endlen)**2)
    s.charge[-s.endlen-1:] = (1. - (1.*iota(0,s.endlen)/s.endlen)**2)
    s.charge[0] = s.charge[1]/2.
    s.charge[-1] = s.charge[-2]/2.
    s.charge[:] = s.charge*s.chargemid

  def geteears(s,accl_ts,accl_dt):
    if not s.setears: return
    if s.icharge == -1:
      s.eears[:] = 0.
      return
    # --- Recalculate beamtime and other params to smooth out glitches
    # --- from the load and fire scheme.
    # --- beamtime is replaced by an array that varies linearly from
    # --- the min of beamtime to the max. The other quantities are then
    # --- appropriately interpolated. The resulting data is dumped directly
    # --- into the cirvar array.
    ii = 1
    beamdt = (s.beamtime[-1]-s.beamtime[0])/s.nn
    s.cirvar[8,:] = s.beamtime[0] + beamdt*iota(0,s.nn)
    # --- Do interpolation using vector notation
    ii = searchsorted(s.beamtime,s.cirvar[8,:]) - 1
    ii[0] = 0
    ii[-1] = s.ngappoints - 1
    ww0 = (s.cirvar[8,:]-take(s.beamtime,ii))/ \
           (take(s.beamtime,ii+1)-take(s.beamtime,ii))
    ww1=1.0 - ww0
    s.cirvar[9,:] = (ww1*take(s.vzbeam,ii) + ww0*take(s.vzbeam,ii+1))/clight
    s.cirvar[13,:] = ww1*take(s.charge,ii) + ww0*take(s.charge,ii+1)
    s.ekintemp[:] = ww1*take(s.ekin,ii) + ww0*take(s.ekin,ii+1)
    #for j in range(s.nn+1):
    #  while s.cirvar[8,:][j] > s.beamtime[ii] and ii < s.nn-1:
    #    ii = min(s.nn,ii + 1)
    #  ww0=(s.cirvar[8,j]-s.beamtime[ii-1])/(s.beamtime[ii]-s.beamtime[ii-1])
    #  ww1=1.0 - ww0
    #  s.cirvar[9,j] = (ww0*s.vzbeam[ii] + ww1*s.vzbeam[ii-1])/clight
    #  s.cirvar[13,j] = ww0*s.charge[ii] + ww1*s.charge[ii-1]
    #  s.ekintemp[j] = ww0*s.ekin[ii] + ww1*s.ekin[ii-1]
    # --- Assume line-charge profile doesn't change
    #s.linechg[:] = s.linechgmid*s.charge/s.chargemid
    currmid = s.linechgmid*s.vzmid
    s.linechg[:] = currmid/s.vzbeam*s.charge/s.chargemid
    # --- Now get perveance
    s.perveance = s.linechg/(4.*pi*eps0*s.ekintemp)
    # --- And from that, the mean beam radius.
    s.cirvar[0,:] = s.hlp*sqrt((2.*s.perveance)/(1.-cos(s.sigma(s)*pi/180.)))
    s.cirvar[2,:] = s.cirvar[0,:]
    # --- Now, get the ez
    lfailed = false
    getezbeam(s.nn+1,s.cirvar,s.ezbeam,s.aperture(s),s.icharge,true,false,
              s.lfixed,lfailed)
    # --- Interpolate ezbeam into eears. This can be done easily since
    # --- the beamtime used here is linearly varying.
    tears = accl_ts + iota(0,s.ntaccl)*accl_dt
    itears = ((tears - s.beamtime[0])/beamdt).astype(long)
    wtears = ((tears - s.beamtime[0])/beamdt) - itears
    wtears[:] = where(less(itears,0),0.,wtears)
    wtears[:] = where(less(s.nn-1,itears),1.,wtears)
    itears[:] = where(less(itears,0),0,itears)
    itears[:] = where(less(s.nn-1,itears),s.nn-1,itears)
    s.eears[:] = take(s.ezbeam,itears  )*(1.-wtears) + \
                 take(s.ezbeam,itears+1)*    wtears
    # --- Find the Ez at or near the end of the beam to replicate through
    # --- the rest of the data.
    #jmin = int((s.beamtime[0] - accl_ts)/accl_dt)
    #jmax = int((s.beamtime[-1] - accl_ts)/accl_dt+1.)
    jmin = argmax(s.eears)
    jmax = argmin(s.eears)
    s.eears[:jmin] = s.eears[jmin]
    s.eears[jmax:] = s.eears[jmax]
    s.eears[:] = -s.eears[:]*s.earlength(s)/s.gap_len
    # --- Do some smoothing and cleaning of the ears. The linearization and
    # --- assumption of none-changing line-charge profile smooth out most of
    # --- glitches from load and fire, but not all. There is still a jagged
    # --- feature in the beams vz profile that leads to some jaggedness in the
    # --- ezbeam.
    # --- Note that the field in the beam center is not zeroed out since
    # --- that component is required because of the long variation in the
    # --- line-charge density resulting from the self-similar current
    # --- scheme of generating the waveforms.
    # --- Now do some smoothing.
    for j in range(s.nendpoints,s.ngappoints+s.nendpoints):
      s.eears[j] = ave(s.eears[j-10:j+10+1])
    s.heears.append(s.eears[::s.zhist]+0.)
    s.hlinechg.append(s.linechg[::s.zhist]+0.)

  ######################################################################
  ######################################################################
  # --- Some extra plotting and data routines
  def plotacclet(s,i=None,t=1,istep=1,acclet=None,acclts=None,accldt=None):
    if not acclet: acclet = top.acclet
    if not acclts: acclts = top.acclts
    if not accldt: accldt = top.accldt
    nta = shape(acclet)[0] - 1
    if not i:
      i1 = 0
      i2 = shape(acclet)[1]
    else:
      i1 = i
      i2 = i + 1
    for i in xrange(i1,i2,istep):
      if t:
        plg(acclet[:,i],acclts[i]+iota(0,nta)*accldt[i])
      else:
        plg(acclet[:,i],iota(-nta/2,nta/2)*accldt[i])

  def getvbeam_z(s,zshift=None):
    if not zshift: zshift = -top.zlatstrt
    # --- gets beam vzmid at each envelope point
    v = zeros(env.nenv+1,'d')
    v[0] = s.hvzmid[0]
    dz = env.dzenv
    ia = 0
    for iz in xrange(1,env.nenv+1):
      zz = env.zl + iz*env.dzenv + zshift
      if zz > s.hzlast[ia+1] and ia+1 < s.ihlp:
        ia = ia + 1
      v[iz] = s.hvzmid[ia]
    return v

  def getibeam_z(s,zshift=None):
    if not zshift: zshift = -top.zlatstrt
    # --- gets beam current at each envelope point
    d = zeros(env.nenv+1,'d')
    d[0] = s.hvzmid[0]*s.hlinechgmid[0]
    dz = env.dzenv
    ia = 0
    for iz in xrange(1,env.nenv+1):
      zz = env.zl + iz*env.dzenv + zshift
      if zz > s.hzlast[ia+1] and ia+1 < s.ihlp:
        ia = ia + 1
      d[iz] = s.hvzmid[ia]*s.hlinechgmid[ia]
    # --- Returns total charge divided by beam duration time the beam length
    # --- scale factor.
    return d


