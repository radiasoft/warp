from warp import *
import LinearAlgebra
env_tools_version = "$Id: env_tools.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

print " "
print "This scripts finds matched beam parameters for a periodic lattice."
print "To use, first set the desired value of sigma."
print "  sig_desr = ???"
print "Then type the command ""match(n)"" where n is the number of iterations"
print "which should be of the order 10 to 20."
print "Check to make sure the numbers printed out are what you want."
print "If not, then run more iterations.  If that doesn't work, the more"
print "sophisticated matching algorithm built into the code should be used."
print "Good luck."
print " "

print 'Envelope matching routines'
print 'match(n) matches the beam giving the desired value of sigma'
print 'sig_desr by varying the emittance.'
print 'matchn(n) matches the beam giving the desired value of sigma'
print 'sig_desr by varying the normalized emittance.'
print 'match1(n) matches the beam by varying a0, b0, ap0, and bp0'
print 'matchenv(quads,af,bf,apf,bpf) modifies the four quads specified'
print 'by quads so that the envelope ends with the specified values'

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
    top.ap0=0.5*(env.apenv[0]+env.apenv[-1])
    top.bp0=0.5*(env.bpenv[0]+env.bpenv[-1])
    step(1)
def match2(n=1):
  """Matches envelope by setting X parameters to average of initial and final values and then setting Y params equal to X params (but b'=-a')"""
  for i in xrange(n):
    top.a0=0.5*(env.aenv[0]+env.aenv[-1])
    top.b0=top.a0
    top.ap0=0.5*(env.apenv[0]+env.apenv[-1])
    top.bp0=-top.ap0
    step(1)

# This routine does a regula-falsi iteration to find the emittance which
# gives the desired sigma.
fff = array([0.,0.])
xxx = array([0.,0.])
def match(n=1,sig_desr=20.):
  """Matches envelope and varies emittance to obtain specified sigma
     - n=1 number of iterations to perform
     - sig_desr=20. desired value of sigma"""
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

def matchn(n=1,sig_desr=20.):
  """Matches envelope and varies normalized emittance to obtain specified sigma
     - n=1 number of iterations to perform
     - sig_desr=20. desired value of sigma"""
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
    iz = int((zz - env.zl)/env.nenv)
    wz =     (zz - env.zl)/env.nenv - iz
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
    iz = int((zz - env.zl)/env.nenv)
    wz =     (zz - env.zl)/env.nenv - iz
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


######################################################################
# Find a matched envelope by varying three of the following parameters
#  a0
#  b0
#  ap0
#  bp0
#  ibeam
#  emit or emitn
#  ekin
#  vbeam
#  sigma0
#  sigma0x
#  sigma0y
#  sigma
#  sigmax
#  sigmay
#  dedx or dbdx
#  hlp
#  occupancy
from powell import *
from optimizer import *
#import powell
class Match3:
  def __init__(self,m1,m2,m3,im=None,
               sig0_desr=None,sig_desr=None,
               sig0x_desr=None,sigx_desr=None,sig0y_desr=None,sigy_desr=None,):
    if not im:
      self.im = [m1,m2,m3]
    else:
      self.im = im
    if 'emit' in self.im and 'emitn' in self.im:
      raise "emit and emitn cannot be 2 of the 3 quantities varied"
    if 'dedx' in self.im and 'dbdx' in self.im:
      raise "dedx and dbdx cannot be 2 of the 3 quantities varied"
    if sig0_desr:
      self.sig0x_desr = sig0_desr
      self.sig0y_desr = sig0_desr
    else:
      self.sig0x_desr = sig0x_desr
      self.sig0y_desr = sig0y_desr
    if sig_desr:
      self.sigx_desr = sig_desr
      self.sigy_desr = sig_desr
    else:
      self.sigx_desr = sigx_desr
      self.sigy_desr = sigy_desr
    self.a0scale = max(top.a0,top.b0)
    self.b0scale = self.a0scale
    self.ap0scale = self.a0scale
    self.bp0scale = self.a0scale
    self.ibeamscale = top.ibeam
    self.emitscale = 1.e-6/(top.vbeam/top.clight)
    self.emitnscale = 1.e-6
    self.ekinscale = top.ekin
    self.vbeamscale = top.vbeam
    self.sigma0scale = 90.
    self.sigma0xscale = 90.
    self.sigma0yscale = 90.
    self.sigmascale = 90.
    self.sigmaxscale = 90.
    self.sigmayscale = 90.
    self.dedxscale = max(max(top.quadde),1.e-20)
    self.dbdxscale = max(max(top.quaddb),1.e-20)
    self.hlpscale = (top.quadze[1] - top.quadze[0])
    self.occupancyscale = (top.quadze[0] - top.quadzs[0])/ \
                          (top.quadze[1] - top.quadze[0])
  def loss(self):
    delta = (env.aenv[0] - env.aenv[-1])/self.a0scale
    deltb = (env.benv[0] - env.benv[-1])/self.b0scale
    deltap = (env.apenv[0] - env.apenv[-1])/self.ap0scale
    deltbp = (env.bpenv[0] - env.bpenv[-1])/self.bp0scale
    deltasig0x = 0
    deltasig0y = 0
    #if 'sigma0' not in self.im and 'sigma0x' not in self.im: 
    #  deltasig0x = (self.sig0x_desr - env.sig0x)/self.sigma0xscale
    #else:
    #  deltasig0x = 0
    if 'sigma' not in self.im and 'sigmax' not in self.im: 
      deltasigx = (self.sigx_desr - env.sigx)/self.sigmaxscale/4.
    else:
      deltasigx = 0
    #if 'sigma0' not in self.im and 'sigma0y' not in self.im: 
    #  deltasig0y = (self.sig0y_desr - env.sig0y)/self.sigma0yscale
    #else:
    #  deltasig0y = 0
    if 'sigma' not in self.im and 'sigmay' not in self.im: 
      deltasigy = (self.sigy_desr - env.sigy)/self.sigmayscale/4.
    else:
      deltasigy = 0
    return delta**2 + deltb**2 + deltap**2 + deltbp**2 + \
           deltasig0x**2 + deltasigx**2 + deltasig0y**2 + deltasigy**2
  def printparam(self,p):
    print p+' = ',
    if p == 'a0': print top.a0
    elif p == 'b0': print top.b0
    elif p == 'ap0': print top.ap0
    elif p == 'bp0': print top.bp0
    elif p == 'ibeam': print top.ibeam
    elif p == 'emit': print top.emit
    elif p == 'emitn': print top.emitn
    elif p == 'ekin': print top.ekin
    elif p == 'vbeam': print top.vbeam
    elif p == 'sigma0': print env.sig0x
    elif p == 'sigma0x': print env.sig0x
    elif p == 'sigma0y': print env.sig0y
    elif p == 'sigma': print env.sigx
    elif p == 'sigmax': print env.sigx
    elif p == 'sigmay': print env.sigy
    elif p == 'dedx': print abs(top.quadde[0])
    elif p == 'dbdx': print abs(top.quaddb[0])
    elif p == 'hlp': print (top.quadze[1] - top.quadze[0])
    elif p == 'occupancy':
      print (top.quadze[0]-top.quadzs[0])/(top.quadze[1] - top.quadze[0])
  def getparam(self,p):
    if p == 'a0': return top.a0/self.a0scale
    elif p == 'b0': return top.b0/self.b0scale
    elif p == 'ap0': return top.ap0/self.ap0scale
    elif p == 'bp0': return top.bp0/self.bp0scale
    elif p == 'ibeam': return top.ibeam/self.ibeamscale
    elif p == 'emit': return top.emit/self.emitscale
    elif p == 'emitn': return top.emitn/self.emitnscale
    elif p == 'ekin': return top.ekin/self.ekinscale
    elif p == 'vbeam': return top.vbeam/self.vbeamscale
    elif p == 'sigma0': return env.sig0x/self.sigma0xscale
    elif p == 'sigma0x': return env.sig0x/self.sigma0xscale
    elif p == 'sigma0y': return env.sig0y/self.sigma0yscale
    elif p == 'sigma': return env.sigx/self.sigmaxscale
    elif p == 'sigmax': return env.sigx/self.sigmaxscale
    elif p == 'sigmay': return env.sigy/self.sigmayscale
    elif p == 'dedx': return abs(top.quadde[0])/self.dedxscale
    elif p == 'dbdx': return abs(top.quaddb[0])/self.dbdxscale
    elif p == 'hlp': return (top.quadze[1] - top.quadze[0])/self.hlpscale
    elif p == 'occupancy':
      return (top.quadze[0]-top.quadzs[0])/(top.quadze[1] - top.quadze[0])/ \
             self.occupancyscale
  def setparam(self,p,v):
    if p == 'a0': top.a0 = v*self.a0scale
    elif p == 'b0': top.b0 = v*self.b0scale
    elif p == 'ap0': top.ap0 = v*self.ap0scale
    elif p == 'bp0': top.bp0 = v*self.bp0scale
    elif p == 'ibeam': top.ibeam = v*self.ibeamscale
    elif p == 'emit':
      top.emit = v*self.emitscale
      top.emitx = top.emit
      top.emity = top.emit
    elif p == 'emitn':
      top.emitn = v*self.emitnscale
      top.emitnx = top.emit
      top.emitny = top.emit
    elif p == 'ekin':
      top.ekin = v*self.ekinscale
      top.vbeam = 0
      top.vbeam_s = 0
      derivqty()
    elif p == 'vbeam': top.vbeam = v*self.vbeamscale
    #elif p == 'sigma0':
    #elif p == 'sigma0x':
    #elif p == 'sigma0y':
    #elif p == 'sigma':
    #elif p == 'sigmax':
    #elif p == 'sigmay':
    elif p == 'dedx': top.quadde = v*top.quadde/abs(top.quadde)*self.dedxscale
    elif p == 'dbdx': top.quaddb = v*top.quaddb/abs(top.quaddb)*self.dbdxscale
    elif p == 'hlp':
      hlp = top.quadze[1] - top.quadze[0]
      s = v/hlp*self.hlpscale
      oldc = 0.5*(top.quadzs + top.quadze)
      newc = oldc*s
      top.quadzs = top.quadzs + newc - oldc
      top.quadze = top.quadze + newc - oldc
      env.zu = env.zl + 2.*hlp*s
      env.dzenv = (env.zu - env.zl)/env.nenv
      env.tuneze = env.zu
      top.zlatperi = env.zu
    elif p == 'occupancy':
      qc = 0.5*(top.quadzs + top.quadze)
      hlp = top.quadze[1] - top.quadze[0]
      ql = v*hlp*self.occupancyscale
      top.quadzs = qc - ql/2.
      top.quadze = qc + ql/2.
  def step(self,newv):
    self.setparam(self.im[0],newv[0])
    self.setparam(self.im[1],newv[1])
    self.setparam(self.im[2],newv[2])
    #self.printparam(self.im[0])
    #self.printparam(self.im[1])
    #self.printparam(self.im[2])
    step()
  def func(self,newv):
    self.step(newv)
    return self.loss()
  def go(self,vary=0.01,delta=None,itmax=10,tol=1.e-10):
    env.lenvout = 0
    self.ppp = array([self.getparam(self.im[0]),
                      self.getparam(self.im[1]),
                      self.getparam(self.im[2])])
    if 1:
      self.xxx = identity(3)*1.
      if delta:
        self.delta = ones(3,'d')*delta
      else:
        self.delta = self.ppp*vary
      self.xxx[0,0] = self.delta[0]
      self.xxx[1,1] = self.delta[1]
      self.xxx[2,2] = self.delta[2]
      (self.ppp,self.xxx,self.fret,self.iter) = \
                          powell(self.ppp,self.xxx,self.func,tol,itmax)
    else:
      #bb = Bonehead(self.ppp,self.step,self.loss)
      #bb.minimize(10,2)
      ss = Spsa(3,self.ppp,self.step,self.loss,.01,.01)
      ss.iter(imax=itmax,err=tol)

