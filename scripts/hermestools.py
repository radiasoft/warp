from warp import *

def sethermesbeam(lsamecharge = false):
  # Sets initial conditions
  # Before calling this routine, env.zl and env.zu should be set to
  # cover exactly one half lattice period. Also env.dzenv should be preset.
  print "Running sethermesbeam ..."
  sethermesprofile(lsamecharge)
  sethermesenvelope()
  print "done."
  
def sethermesprofile(lsamecharge):
  delzbeam = top.zimax - top.zimin
  def f(z):
    fraction = (1.-top.straight)/2.
    if z <= 0.: return 0.
    if 0. < z <= fraction: return (1.-power(1.-z/fraction,her.iprofile))
    if fraction < z <= 1.-fraction: return 1.
    if 1.-fraction < z <= 1.: return (1.-power((z-1+fraction)/fraction,her.iprofile))
    if z > 1.: return 0.
  def fint(z):
    fraction = (1.-top.straight)/2.
    if z <= 0.:
      return 0.
    if 0. < z <= fraction:
      return (z + (fraction/(her.iprofile+1.))*(power(1-z/fraction,her.iprofile+1.)) - fraction/(her.iprofile+1.))
    if fraction < z <= 1.-fraction:
      return (fraction*her.iprofile/(her.iprofile+1.)+z-fraction)
    if 1. - fraction < z <= 1.:
      return (z - (fraction/(her.iprofile+1.))*(1.+power((z-1.)/fraction+1.,her.iprofile+1.)))
    if z > 1.: return (1.-2*fraction/(her.iprofile+1.))
  if lsamecharge:
    her.var[8,0] = 0.
    her.var[8,-1] = 1.
    for i in range (1,her.niz-1):
      lower = 0.
      upper = 1.
      target = fint(1.)*i/(her.niz-1.)
      while upper - lower > 1.e-12:
        guess = (upper+lower)/2.
        if fint(guess) > target: upper = guess
        else: lower = guess
      her.var[8,i] = (upper+lower)/2.
    her.var[8,:] = delzbeam * her.var[8,:]
    her.var[13,:-1] = top.ibeam*delzbeam*fint(1.)/(top.vbeam*(her.niz-1.))
  else:
    for i in range (her.niz):
      her.var[8,i] = i / (her.niz - 1.)
    if 0 < her.iprofile < 5:
      for i in range (her.niz-1):
        her.var[13,i] = fint(her.var[8,i+1])-fint(her.var[8,i])
    her.var[13,:] = her.var[13,:]*top.ibeam*delzbeam/top.vbeam
    her.var[8,:] = delzbeam * her.var[8,:]
  extebher(0.,her.dther,her.var,her.niz) # Needed to set up rpipe and therefore rmmax
  her.var[9,:] = top.vbeam / top.clight
  her.var[0,:] = wrz.rmmax / 2. # For now
  her.var[1,:] = 0.             # For now
  her.var[2,:] = wrz.rmmax / 2. # For now
  her.var[3,:] = 0.             # For now
  lfail = false
  loadcharge(her.niz,her.var,her.rpipe,her.icharge,lfail)
  if lfail:
    print "Error in loadcharge: Slice positions are not in increasing order"
    return
  getcurrent(her.var,her.niz,her.rpipe,her.icharge,her.lcurgrid,lfail)
  if lfail:
    print "Error in getcurrent: Slice positions are not in increasing order"
    return
  if top.emitn:
    her.var[10:12,:] = top.emitn * her.var[12,:] / top.ibeam
  else:
    gaminv = sqrt(1.-power(her.var[9,:],2))
    her.var[10:12,:] = (top.emit*gaminv/her.var[9,:]) * her.var[12,:] / top.ibeam
  if her.dther < 0.: her.var[8,:] = her.var[8,:] + her.zlher
  else: her.var[8,:] = her.var[8,:] - her.var[8,-1] + her.zlher
                        
def sethermesenvelope (niter = 1000, errorlimit = 1.e-9):
  # --- Finds the matched conditions. Assumes a cigar load.
  her.var[4:8,:] = 0.    # Centroid motion not implemented in HERMES
  her.var[:4,0]  = 0.    # First slice at zero (cigar load)
  her.var[:4,-1] = 0.    # Last slice also at zero (cigar load)
  # --- Save top variables
  ibeam = top.ibeam
  emit = top.emit
  ekin = top.ekin
  a0 = top.a0
  b0 = top.b0
  ap0 = top.ap0
  bp0 = top.bp0
  # --- Load env package
  package ("env"); generate();
  quadindex = searchsorted (top.quadzs,env.zl)
  top.ap0 = sign (top.ap0, -top.quadde[quadindex]+top.quaddb[quadindex]*top.vbeam)
  top.bp0 = - top.ap0
  # --- Calculate lattice half period
  lhp = env.zu - env.zl
  # --- Start loop
  for i in range (1,her.niz-1):
    # --- Load variables for this slice
    top.emit_s = 0.
    top.emitn = sqrt(her.var[10,i]*her.var[11,i])
    top.emit = top.emitn / (her.var[9,i]*top.gammabar)
    top.ibeam = her.var[12,i]
    # --- Recalculate stuff
    derivqty()
    envx()
    for counter in range(niter):
      top.a0  =  0.25*abs(env.aenv[0]+env.aenv[-1]+env.benv[0]+env.benv[-1])
      top.ap0 =  0.25*(env.apenv[0]-env.apenv[-1]-env.bpenv[0]+env.bpenv[-1])
      top.b0  =  top.a0
      top.bp0 = -top.ap0
      envx()
      error = abs(top.a0-env.aenv[-1]) 
      error = error + abs(top.b0-env.benv[-1])
      error = error + abs(top.ap0+env.apenv[-1])
      error = error + abs(top.bp0+env.bpenv[-1])
      if error < errorlimit:
        break
    if (error > errorlimit):
      print "Slice", i, ": After", niter, "steps"
      raise "No convergence in sethermesenvelope"
    # --- Calculate the envelope value at the position of this slice
    zu = env.zu
    dz = env.dzenv
    env.zu = her.var[8,i]
    while env.zu < env.zl: env.zu = env.zu + 2. * lhp
    while env.zu > env.zl + 2. * lhp: env.zu = env.zu - 2. * lhp
    env.dzenv = (env.zu -env.zl) / env.nenv
    envx()
    env.zu = zu
    env.dzenv = dz
    # --- Save the results
    her.var[0,i] = env.aenv[-1]
    her.var[1,i] = env.apenv[-1]
    her.var[2,i] = env.benv[-1]
    her.var[3,i] = env.bpenv[-1]
  top.a0 = a0
  top.b0 = b0
  top.bp0 = bp0
  top.ibeam = ibeam
  top.emit_s = 0.
  top.emit = emit
  derivqty()
  envx()
  package ("her")

