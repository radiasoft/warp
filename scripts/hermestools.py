from warp import *
hermestools_version = "$Id: hermestools.py,v 1.2 2001/03/01 23:48:20 dave Exp $"

def hermestoolsdoc():
  print """
Tools for Hermes
sethermesbeam: 
  """

def sethermesbeam():
  # Sets initial conditions
  # Before calling this routine, env.zl and env.zu should be set to
  # cover exactly one half lattice period. Also env.dzenv should be preset.
  print "Running sethermesbeam ..."
  delzbeam = top.zimax - top.zimin
  her.var[8,:] = arrayrange (her.niz) * delzbeam / (her.niz - 1)
  if her.dther < 0.: her.var[8,:] = her.var[8,:] + her.zlher
  else: her.var[8,:] = her.var[8,:] - her.var[8,-1] + her.zlher
  her.var[9,:] = top.vbeam / top.clight
  sethermescurrent()
  sethermesenvelope()
  print "done."

def sethermescurrent():
  delzbeam = top.zimax - top.zimin
  currise = (1. - top.straight)/2.
  curfall = (1. - top.straight)/2.
  curmax = 0.
  if her.iprofile <= 0:
    part = (her.var[8,:] - her.var[8,0]) /  (her.var[8,-1] - her.var[8,0])
    arg1 = part / currise
    arg2 = (1. - part) / curfall
    her.var[12,:] = tanh(arg1) * tanh(arg2)
  elif her.iprofile < 4:
    def f(x):
      return x - power(x-1,her.iprofile+1.)/(her.iprofile+1.)
    zrise = currise*delzbeam
    zfall = curfall*delzbeam
    z = 0.
    for i in range (her.niz-1):
      dz = her.var[8,i+1] - her.var[8,i]
      if z <= zrise and z + dz <= zrise:
        her.var[13,i] = zrise * (f((z+dz)/zrise) - f(z/zrise))
      if z <= zrise and zrise < z + dz <= delzbeam - zfall:
        her.var[13,i] = zrise * (1. - f(z/zrise)) + (z+dz-zrise)
      if z <= zrise and delzbeam - zfall < z + dz:
        her.var[13,i] = zrise * (1. - f(z/zrise)) + (delzbeam-zrise-zfall) + zfall * (1 - f((delzbeam-z-dz)/zfall))
      if zrise < z <= delzbeam - zfall and z + dz <= delzbeam - zfall:
        her.var[13,i] = dz
      if zrise < z <= delzbeam - zfall and delzbeam - zfall < z + dz:
        her.var[13,i] = (delzbeam-zfall-z) + zfall * (1. - f((delzbeam-z-dz)/zfall))
      if delzbeam - zfall < z and delzbeam - zfall < z + dz:
        her.var[13,i] = zfall * (f((delzbeam-z)/zfall) - f((delzbeam-z-dz)/zfall))
      z = z + dz
  her.var[13,:] = top.ibeam*her.var[13,:]/top.vbeam
  gethertmp(her.var,her.niz)
  if top.emitn:
    her.var[10:12,:] = top.emitn * her.var[12,:] / top.ibeam
  else:
    gaminv = sqrt(1.-power(her.var[9,:],2))
    her.var[10:12,:] = (top.emit*gaminv/her.var[9,:]) * her.var[12,:] / top.ibeam
                        
def sethermesenvelope (niter = 1000, errorlimit = 1.e-9):
  """
This does ...
  - niter=1000: ...
  """
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

