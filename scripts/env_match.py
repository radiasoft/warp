from warp import *
import LinearAlgebra
env_match_version = "$Id: env_match.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

print 'Envelope matching routines'
print 'match() matches the beam giving the desired value of sigma, varying'
print '        unnormalized emittance'
print 'matchn() matches the beam giving the desired value of sigma, varying'
print '         normalized emittance'
print 'match1() matches the beam by varying a0, b0, ap0, and bp0'
print 'match2() matches the beam by varying a0, b0, ap0, and bp0, forcing'
print '         b0=a0 and bp0=-ap0'
print 'matchenv() modifies the four quads specified to match final values'
print "matchxenv() finds initial x,x',y,y' to match specified final values"
print 'Type "doc(match)" for info on arguments'

# This is a very simple algorithm to match a beam over a lattice period.
# The function sets a0 equal to average of 'a' at the beginning and end of
# the lattice period (assuming that the envelope solve covers only one
# lattice period).  It does the same for b, a', and b' and then
# recalculates the envelope.
def match1(n=1):
  """Matches envelope by setting parameters to average of initial and final values"""
  for i in xrange(n):
    top.a0=0.5*(env.aenv[0]+env.aenv[-1])
    top.b0=0.5*(env.benv[0]+env.benv[-1])
    top.ap0=0.5*(env.apenv[0]+env.apenv[-1]*env.vzenv[-1]/env.vzenv[0])
    top.bp0=0.5*(env.bpenv[0]+env.bpenv[-1]*env.vzenv[-1]/env.vzenv[0])
    step(1)
def match2(n=1):
  """Matches envelope by setting X parameters to average of initial and final values and then setting Y params equal to X params (but b'=-a')"""
  for i in xrange(n):
    top.a0=0.5*(env.aenv[0]+env.aenv[-1])
    top.b0=top.a0
    top.ap0=0.5*(env.apenv[0]+env.apenv[-1]*env.vzenv[-1]/env.vzenv[0])
    top.bp0=-top.ap0
    step(1)

# This routine does a regula-falsi iteration to find the emittance which
# gives the desired sigma.
fff = array([0.,0.])
xxx = array([0.,0.])
def match(n=1,sig_desr=20.,zl=None,zu=None):
  """Matches envelope and varies emittance to obtain specified sigma
     - n=1 number of iterations to perform
     - sig_desr=20. desired value of sigma
     - zl starting z value=env.zl
     - zu ending z value=env.zu"""
  if not zl:
    zl = env.zl
  if not zu:
    zu = env.zu
  zlsave = env.zl
  zusave = env.zu
  tunezssave = env.tunezs
  tunezesave = env.tuneze
  env.zl = zl
  env.zu = zu
  env.tunezs = zl
  env.tuneze = zu
  env.lenvout = false
  xxx[0] = top.emit
  xxx[1] = top.emit*1.05
  top.emit=xxx[0];derivqty();match1(20);fff[0]=top.sigmax - sig_desr
  top.emit=xxx[1];derivqty();match1(20);fff[1]=top.sigmax - sig_desr
  j=1
  while (j <= n and abs(min(fff[0],fff[1])) > 1.e-10):
    j = j + 1
    env.lenvout = false
    nxxx = (xxx[0]*fff[1] - xxx[1]*fff[0])/(fff[1] - fff[0])
    ixxx = 0
    if fff[1] > fff[0]: ixxx = 1
    xxx[ixxx] = nxxx
    top.emit=xxx[ixxx];derivqty();match1(20);fff[ixxx]=top.sigmax - sig_desr
    print "Error = %f" % fff[ixxx]
    env.lenvout = true
    step(1)
  env.zl = zlsave
  env.zu = zusave
  env.tunezs = tunezssave
  env.tuneze = tunezesave

def matchn(n=1,sig_desr=20.,zl=None,zu=None):
  """Matches envelope and varies normalized emittance to obtain specified sigma
     - n=1 number of iterations to perform
     - sig_desr=20. desired value of sigma
     - zl starting z value=env.zl
     - zu ending z value=env.zu"""
  if not zl:
    zl = env.zl
  if not zu:
    zu = env.zu
  zlsave = env.zl
  zusave = env.zu
  tunezssave = env.tunezs
  tunezesave = env.tuneze
  env.zl = zl
  env.zu = zu
  env.tunezs = zl
  env.tuneze = zu
  env.lenvout = false
  xxx[0] = top.emitn
  xxx[1] = top.emitn*1.05
  top.emitn=xxx[0];derivqty();match1(20);fff[0]=top.sigmax - sig_desr
  top.emitn=xxx[1];derivqty();match1(20);fff[1]=top.sigmax - sig_desr
  j=1
  while (j <= n and abs(min(fff[0],fff[1])) > 1.e-10):
    j = j + 1
    env.lenvout = false
    nxxx = (xxx[0]*fff[1] - xxx[1]*fff[0])/(fff[1] - fff[0])
    ixxx = 0
    if fff[1] > fff[0]: ixxx = 1
    xxx[ixxx] = nxxx
    top.emitn=xxx[ixxx];derivqty();match1(20);fff[ixxx]=top.sigmax - sig_desr
    print "Error = %f" % fff[ixxx]
    env.lenvout = true
    step(1)
  env.zl = zlsave
  env.zu = zusave
  env.tunezs = tunezssave
  env.tuneze = tunezesave

# Given a current, this scales the beam parameters appropriately.
def scale(ii=None):
  """Scales beam parameters to match the specified current
     ii = top.ibeam value of current to scale to"""
  if not ii: ii=top.ibeam
  top.a0 = top.a0*sqrt(ii/top.ibeam)
  top.b0 = top.b0*sqrt(ii/top.ibeam)
  top.ap0 = top.ap0*sqrt(ii/top.ibeam)
  top.bp0 = top.bp0*sqrt(ii/top.ibeam)
  top.emit = top.emit*ii/top.ibeam
  top.emitn = top.emitn*ii/top.ibeam
  top.ibeam = ii





########################################################################
def matchenv(quads,af,bf,apf,bpf,zz=None,maxiter=100,tol=1.e-10):
  """Varies 4 quads to match the envelope to the specified final values
     - quads index of quad elements which are to be varied
     - af final value of a
     - bf final value of b
     - apf final value of a'
     - bpf final value of b'
     - zz=env.zu z location of final values (must be env.zl < zz < env.zu)
     - maxiter=100 maximum number of iterations to perform
     - tol=1.e-10 tolerance to match final values to"""
  if len(quads) != 4:
    print 'Error: exactly four quads for varying must be specified'
    return
  if zz:
    if zz < env.zl or env.zu < zz:
      print 'Error: the z location speficied must be within zl and zu'
      return
  mat = zeros((4,4),'d')
  denv = zeros(4,'d')
  vvary = 0.01

  # --- Make initial envelope calculation.
  step()
  if zz:
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
  denv[:] = [af - asave,bf - bsave,apf - apsave, bpf - bpsave]

  notdone = 1
  iter = 0
  while notdone and iter < maxiter:
    iter = iter + 1
    # --- Vary each quad
    le = env.lenvout
    env.lenvout = false
    for iq in xrange(4):

      if top.quadde[quads[iq]] != 0.:
        qorig = top.quadde[quads[iq]]
        top.quadde[quads[iq]] = top.quadde[quads[iq]]*(1. + vvary)
        qorig = qorig/top.vbeam
      elif top.quaddb[quads[iq]] != 0.:
        qorig = top.quaddb[quads[iq]]
        top.quaddb[quads[iq]] = top.quaddb[quads[iq]]*(1. + vvary)

      step()

      if top.quadde[quads[iq]] != 0.:
        top.quadde[quads[iq]] = qorig*top.vbeam
      elif top.quaddb[quads[iq]] != 0.:
        top.quaddb[quads[iq]] = qorig

      if zz:
        ai = env.aenv[iz]*(1. - wz) + env.aenv[iz+1]*wz
        bi = env.benv[iz]*(1. - wz) + env.benv[iz+1]*wz
        api = env.apenv[iz]*(1. - wz) + env.apenv[iz+1]*wz
        bpi = env.bpenv[iz]*(1. - wz) + env.bpenv[iz+1]*wz
      else:
        ai = env.aenv[env.nenv]
        bi = env.benv[env.nenv]
        api = env.apenv[env.nenv]
        bpi = env.bpenv[env.nenv]
      mat[0,iq]=(ai  - asave)/(qorig*vvary)
      mat[1,iq]=(bi  - bsave)/(qorig*vvary)
      mat[2,iq]=(api - apsave)/(qorig*vvary)
      mat[3,iq]=(bpi - bpsave)/(qorig*vvary)

    # --- Invert the matrix
    mati = LinearAlgebra.inverse(mat)

    # --- Get next set of quad voltages.
    for iq in xrange(4):
      delta = sum(mati[iq,:]*denv)
      if top.quadde[quads[iq]] != 0.:
        top.quadde[quads[iq]]=top.quadde[quads[iq]]+delta*top.vbeam
      elif top.quaddb[quads[iq]] != 0.:
        top.quaddb[quads[iq]]=top.quaddb[quads[iq]]+delta

    # --- Calculate envelope with new quad values.
    env.lenvout = le
    step()
    if zz:
      asave = env.aenv[iz]*(1. - wz) + env.aenv[iz+1]*wz
      bsave = env.benv[iz]*(1. - wz) + env.benv[iz+1]*wz
      apsave = env.apenv[iz]*(1. - wz) + env.apenv[iz+1]*wz
      bpsave = env.bpenv[iz]*(1. - wz) + env.bpenv[iz+1]*wz
    else:
      asave = env.aenv[env.nenv]
      bsave = env.benv[env.nenv]
      apsave = env.apenv[env.nenv]
      bpsave = env.bpenv[env.nenv]

    # --- Check for convergence
    denv[:] = [af - asave,bf - bsave,apf - apsave, bpf - bpsave]
    print "error => a = %10.3e b = %10.3e a' = %10.3e b' = %10.3e "%tuple(denv)
    if max(abs(denv)) < tol: notdone = 0

  if iter == maxiter:
    print 'Warning: Maximum number of iterations reached'

########################################################################
def matchxenv(xf=0.,xpf=0.,yf=0.,ypf=0.,zz=None,maxiter=100,tol=1.e-10):
  """Matches single particle orbit to final values
     - xf final value of x
     - xpf final value of x'
     - yf final value of y
     - ypf final value of y'
     - zz=env.zu z location of final values (must be env.zl < zz < env.zu)
     - maxiter=100 maximum number of iterations to perform
     - tol=1.e-10 tolerance to match final values to"""
  if zz:
    if zz < env.zl or env.zu < zz:
      print 'Error: the z location speficied must be within zl and zu'
      return
  matx = zeros((2,2),'d')
  denvx = zeros(2,'d')
  xvary = 0.01
  maty = zeros((2,2),'d')
  denvy = zeros(2,'d')
  yvary = 0.01

  # --- Make initial envelope calculation.
  step()
  if zz:
    iz = int((zz - env.zl)/env.dzenv)
    wz =     (zz - env.zl)/env.dzenv - iz
    xsave = env.xenv[iz]*(1. - wz) + env.xenv[iz+1]*wz
    xpsave = env.xpenv[iz]*(1. - wz) + env.xpenv[iz+1]*wz
    ysave = env.yenv[iz]*(1. - wz) + env.yenv[iz+1]*wz
    ypsave = env.ypenv[iz]*(1. - wz) + env.ypenv[iz+1]*wz
  else:
    xsave = env.xenv[env.nenv]
    xpsave = env.xpenv[env.nenv]
    ysave = env.yenv[env.nenv]
    ypsave = env.ypenv[env.nenv]
  denvx[:] = [xf - xsave,xpf - xpsave]
  denvy[:] = [yf - ysave,ypf - ypsave]

  notdone = 1
  iter = 0
  while notdone and iter < maxiter:
    iter = iter + 1
    # --- Vary each quad
    le = env.lenvout
    env.lenvout = false
    for iq in xrange(2):

      if iq==0:
        xorig = top.x0
        top.x0 = top.x0*(1. + xvary)
        yorig = top.y0
        top.y0 = top.y0*(1. + yvary)
      elif iq==1:
        xorig = top.xp0
        top.xp0 = top.xp0*(1. + xvary)
        yorig = top.yp0
        top.yp0 = top.yp0*(1. + yvary)

      step()

      if iq==0:
        top.x0 = xorig
        top.y0 = yorig
      elif iq==1:
        top.xp0 = xorig
        top.yp0 = yorig

      if zz:
        xi = env.xenv[iz]*(1. - wz) + env.xenv[iz+1]*wz
        xpi = env.xpenv[iz]*(1. - wz) + env.xpenv[iz+1]*wz
        yi = env.yenv[iz]*(1. - wz) + env.yenv[iz+1]*wz
        ypi = env.ypenv[iz]*(1. - wz) + env.ypenv[iz+1]*wz
      else:
        xi = env.xenv[env.nenv]
        xpi = env.xpenv[env.nenv]
        yi = env.yenv[env.nenv]
        ypi = env.ypenv[env.nenv]
      matx[0,iq]=(xi  - xsave)/(xorig*xvary)
      matx[1,iq]=(xpi - xpsave)/(xorig*xvary)
      maty[0,iq]=(yi  - ysave)/(yorig*yvary)
      maty[1,iq]=(ypi - ypsave)/(yorig*yvary)

    # --- Invert the matrix
    matxi = LinearAlgebra.inverse(matx)
    matyi = LinearAlgebra.inverse(maty)

    # --- Get next set of quad voltages.
    top.x0 = top.x0 + sum(matxi[0,:]*denvx)
    top.xp0 = top.xp0 + sum(matxi[1,:]*denvx)
    top.y0 = top.y0 + sum(matyi[0,:]*denvy)
    top.yp0 = top.yp0 + sum(matyi[1,:]*denvy)

    # --- Calculate envelope with new quad values.
    env.lenvout = le
    step()
    if zz:
      xsave = env.xenv[iz]*(1. - wz) + env.xenv[iz+1]*wz
      xpsave = env.xpenv[iz]*(1. - wz) + env.xpenv[iz+1]*wz
      ysave = env.yenv[iz]*(1. - wz) + env.yenv[iz+1]*wz
      ypsave = env.ypenv[iz]*(1. - wz) + env.ypenv[iz+1]*wz
    else:
      xsave = env.xenv[env.nenv]
      xpsave = env.xpenv[env.nenv]
      ysave = env.yenv[env.nenv]
      ypsave = env.ypenv[env.nenv]

    # --- Check for convergence
    denvx[:] = [xf - xsave,xpf - xpsave]
    denvy[:] = [yf - ysave,ypf - ypsave]
    print "error => x = %10.3e x' = %10.3e"%tuple(denvx)
    print "error => y = %10.3e y' = %10.3e"%tuple(denvy)
    if max(abs(denvx)) < tol and max(abs(denvy)) < tol: notdone = 0

  if iter == maxiter:
    print 'Warning: Maximum number of iterations reached'
  else:
    print "top.x0 = %20.15e;top.xp0 = %20.15e"%(top.x0,top.xp0)
    print "top.y0 = %20.15e;top.yp0 = %20.15e"%(top.y0,top.yp0)

