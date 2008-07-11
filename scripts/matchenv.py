from warp import *
import numpy.linalg as linalg
matchenv_version = "$Id: matchenv.py,v 1.2 2008/07/11 21:13:43 dave Exp $"

def matchenvdoc():
  print """
matchenv() Varies 4 quads to match the envelope to the specified final
           values The package env must be generated before calling this routine
"""

def matchenv(iquads,af,bf,apf,bpf,zz=None,maxiter=100,tol=1.e-10,
             usequad=0,usehele=0,useemlt=0,usemmlt=0):
  """
Varies 4 quads to match the envelope to the specified final values
The package env must be generated before calling this routine
  - iquads index of quad elements which are to be varied
  - af final value of a
  - bf final value of b
  - apf final value of a'
  - bpf final value of b'
  - zz=env.zu z location of final values (must be env.zl < zz < env.zu)
  - maxiter=100 maximum number of iterations to perform
  - tol=1.e-10 tolerance to match final values to
  - usequad when true, forces use of quad elements
  - usehele when true, forces use of hele elements
  - useemlt when true, forces use of emlt elements
  - usemmlt when true, forces use of mmlt elements
  """
  # --- First, do some error checking
  if len(iquads) != 4:
    print 'Error: exactly four quads for varying must be specified'
    return
  if zz:
    if zz < env.zl or env.zu < zz:
      print 'Error: the z location speficied must be within zl and zu'
      return
  # Make sure that the indices in iquads are all valid
  try:
    if usequad: dummy = take(top.quadde,iquads)
    if usehele: dummy = take(top.heleae[0,:],iquads)
    if useemlt: dummy = take(top.emltid,iquads)
    if usemmlt: dummy = take(top.mmltid,iquads)
  except IndexError:
    print "Error: iquads index out of bounds"
    return
  # --- Set some constants
  STPMX = 100.0
  TOLX = 1.0e-4   # Minimum relative change in quad strengths for the iteration to continue
  ALF = 1.0e-4    # Ensures sufficient decrease in function value
  EPS = 1.0e-8     # Approximate square root of the machine precision.
  # --- Initialisation
  print '***** Starting the iteration *****'
  le = env.lenvout
  env.lenvout = false
  # --- Save initial quadrupole field data
  elementsfound = 0
  if usequad:
    elementsfound = elementsfound + 1
    savedquadde = take(top.quadde[:],iquads)
    savedquaddb = take(top.quaddb[:],iquads)
  if usehele:
    elementsfound = elementsfound + 1
    savedheleae = zeros((top.nhmlt,4),'d')
    savedheleam = zeros((top.nhmlt,4),'d')
    for i in range(4):
      iq = iquads[i]
      savedheleae[:,i] = top.heleae[:,iq]
      savedheleam[:,i] = top.heleam[:,iq]
  if useemlt:
    elementsfound = elementsfound + 1
    savedesemlt = take(top.esemlt[:,:,:],iquads,2)
  if usemmlt:
    elementsfound = elementsfound + 1
    savedmsmmlt = take(top.msmmlt[:,:,:],iquads,2)
  if elementsfound == 0:
    print "Error: No lattice elements found."
    return
  # --- Make initial envelope calculation
  step()
  if zz and zz < env.zu:
    iz = int((zz - env.zl)/env.dzenv)
    wz =     (zz - env.zl)/env.dzenv - iz
    asave = env.aenv[iz]*(1. - wz) + env.aenv[iz+1]*wz
    bsave = env.benv[iz]*(1. - wz) + env.benv[iz+1]*wz
    apsave = env.apenv[iz]*(1. - wz) + env.apenv[iz+1]*wz
    bpsave = env.bpenv[iz]*(1. - wz) + env.bpenv[iz+1]*wz
  else:
    asave = env.aenv[env.nenv]
    bsave = env.benv[env.nenv]
    apsave = env.apenv[env.nenv]
    bpsave = env.bpenv[env.nenv]
  fvec = array([asave-af,apsave-apf,bsave-bf,bpsave-bpf])
  f = 0.5 * sum (fvec*fvec)
  print 'Initial error is ' + str(max(abs(fvec)))
  scaling = array ([1.0, 1.0, 1.0, 1.0])
  stpmax = 4.0 * STPMX

  # --- Start iterating
  notdone = 1
  iter = 0
  while notdone and iter < maxiter:
    iter = iter + 1
    # --- Vary each quad
    fjac = zeros ((4,4),'d')
    g = zeros ((4),'d') 
    p = zeros ((4),'d')
    xold = zeros ((4),'d')
    for i in range (4):
       temp = scaling[i]
       h = EPS * abs(temp)
       if (h == 0.0): h = EPS
       scaling[i] = temp + h # Trick to reduce finite precision error.
       h = scaling[i] - temp;
       if usequad:
         for j in range (4):
           iq = iquads[j]
           top.quadde[iq] = savedquadde[j] * scaling[j]
           top.quaddb[iq] = savedquaddb[j] * scaling[j]
       if usehele:
         for j in range (4):
           iq = iquads[j]
           top.heleae[:,iq] = savedheleae[:,j] * scaling[j]
           top.heleam[:,iq] = savedheleam[:,j] * scaling[j]
       if useemlt:
         for j in range (4):
           iq = iquads[j]
           top.esemlt[:,:,iq] = savedesemlt[:,:,j] * scaling[j]
       if usemmlt:
         for j in range (4):
           iq = iquads[j]
           top.msmmlt[:,:,iq] = savedmsmmlt[:,:,j] * scaling[j]
       step()
       if zz and zz < env.zu:
         iz = int((zz - env.zl)/env.dzenv)
         wz =     (zz - env.zl)/env.dzenv - iz
         asave = env.aenv[iz]*(1. - wz) + env.aenv[iz+1]*wz
         bsave = env.benv[iz]*(1. - wz) + env.benv[iz+1]*wz
         apsave = env.apenv[iz]*(1. - wz) + env.apenv[iz+1]*wz
         bpsave = env.bpenv[iz]*(1. - wz) + env.bpenv[iz+1]*wz
       else:
         asave = env.aenv[env.nenv]
         bsave = env.benv[env.nenv]
         apsave = env.apenv[env.nenv]
         bpsave = env.bpenv[env.nenv]
       fnew = array([asave-af,apsave-apf,bsave-bf,bpsave-bpf])
       scaling[i] = temp
       fjac[:,i] = (fnew[:]-fvec[:]) / h   # Forward difference formula.
       g[i] = sum (fjac[:,i]*fvec[:])
    xold[:] = scaling[:]
    p[:] = -fvec[:]
    p = linalg.solve(fjac, p);
    # Starting linesearch
    fold = f
    total = sqrt(sum(p[:]*p[:]))
    if (total>stpmax): p[:] = p[:] * stpmax / total # Scale if attempted step is too big.
    slope = sum (g[:]*p[:])
    if (slope >= 0.0): raise 'Roundoff problem in lnsrch'
    junk = zeros ((4),'d')
    for i in range(4): junk[i] = abs (p[i]) / max (abs(xold[i]),1.0)
    test = max (junk[:]) # Compute lambda_min
    alamin = TOLX / test;
    alam = 1.0; # Always try full Newton step first.
    f2 = 0;
    alam2 = 0;
    completed = 0
    while (not completed):         # Start of iteration loop.
      tmplam = 0;
      for i in range(4):
        iq = iquads[i]
        scaling[i] = xold[i] + alam * p[i]
        if usequad:
          top.quadde[iq] = savedquadde[i] * scaling[i]
          top.quaddb[iq] = savedquaddb[i] * scaling[i]
        if usehele:
          top.heleae[:,iq] = savedheleae[:,i] * scaling[i]
          top.heleam[:,iq] = savedheleam[:,i] * scaling[i]
        if useemlt:
          top.esemlt[:,:,iq] = savedesemlt[:,:,i] * scaling[i]
        if usemmlt:
          top.msmmlt[:,:,iq] = savedmsmmlt[:,:,i] * scaling[i]
      step()
      if zz and zz < env.zu:
        iz = int((zz - env.zl)/env.dzenv)
        wz =     (zz - env.zl)/env.dzenv - iz
        asave = env.aenv[iz]*(1. - wz) + env.aenv[iz+1]*wz
        bsave = env.benv[iz]*(1. - wz) + env.benv[iz+1]*wz
        apsave = env.apenv[iz]*(1. - wz) + env.apenv[iz+1]*wz
        bpsave = env.bpenv[iz]*(1. - wz) + env.bpenv[iz+1]*wz
      else:
        asave = env.aenv[env.nenv]
        bsave = env.benv[env.nenv]
        apsave = env.apenv[env.nenv]
        bpsave = env.bpenv[env.nenv]
      fvec = array([asave-af,apsave-apf,bsave-bf,bpsave-bpf])
      total = sum (fvec[:]*fvec[:])
      f = 0.5 * total
      if (alam < alamin): completed = 1 # Convergence on delta x.
      elif (f <= fold + ALF * alam * slope): completed = 1 # Sufficient function decrease
      elif (alam == 1.0): tmplam = - slope / (2.0 * (f-fold-slope)) # Backtrack first time
      else: # Subsequent backtracks.
        rhs1 = f - fold - alam * slope;
        rhs2 = f2 - fold - alam2 * slope;
        a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
        if (a == 0.0): tmplam = - slope / (2.0 * b);
        else:
          disc = b*b-3.0*a*slope;
          if (disc < 0.0): tmplam = 0.5 * alam;
          elif (b <= 0.0): tmplam = (-b + sqrt(disc))/(3.0*a);
          else: tmplam = - slope / (b + sqrt(disc));
        if (tmplam > 0.5 * alam): tmplam = 0.5 * alam;     # lambda <= 0.5 lambda1
      alam2 = alam;
      f2 = f;
      alam = max (tmplam,0.1*alam) # lambda >= 0.1 lambda1
    # End of linesearch
    test = max (abs(fvec));
    print 'Error is ' + str(test)
    if (test < tol):
      print '***** Iteration finished *****'
      notdone = 0
  env.lenvout = le
  for j in range (4):
    iq = iquads[j]
    if useemlt: top.esemltp[:,:,iq] = top.esemltp[:,:,iq] * scaling[j]
    if usemmlt: top.msmmltp[:,:,iq] = top.msmmltp[:,:,iq] * scaling[j]
    if usehele:
      top.heleep[:,iq] = top.heleep[:,iq] * scaling[j]
      top.helemp[:,iq] = top.helemp[:,iq] * scaling[j]




