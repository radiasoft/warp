"""Functions to find matched beam parameters for a periodic lattice using Circe."""
from warp import *
cir_match_version = "$Id: cir_match.py,v 1.3 2003/08/07 19:00:58 dave Exp $"

print " "
print "This scripts finds matched beam parameters for a periodic lattice."
print "To use, first set the desired value of sigma."
print "  sig_desr = ???"
print "Then type the command 'match(n)' where n is the number of iterations"
print "which should be of the order 10 to 20."
print "Check to make sure the numbers printed out are what you want."
print "If not, then run more iterations.  If that doesn't work, the more"
print "sophisticated matching algorithm built into the code should be used."
print "Good luck."
print " "

print "CIRCE matching routines"
print "match(n) matches the beam giving the desired value of sigma"
print "sig_desr by varying the emittance."
print "matchn(n) matches the beam giving the desired value of sigma"
print "sig_desr by varying the normalized emittance."
print "match1(n) matches the beam by varying a0, b0, ap0, and bp0"
print "matchenv(quads,af,bf,apf,bpf) modifies the four quads specified"
print "by quads so that the envelope ends with the specified values"

# This is a very simple algorithm to match a beam over a lattice period.
# The function sets a0 equal to average of 'a' at the beginning and end of
# the lattice period (assuming that the envelope solve covers only one
# lattice period).  It does the same for b, a', and b' and then
# recalculates the envelope.
def match1(n=1,varsave=None,zl=None,zu=None,savehist=false,errorlimit=1.e-6):
  if not varsave:
    # --- If varsave is not set, then create a beam in the standard circe
    # --- database and make varsave point to it.
    initbeam()
    varsave = cir.var
  if not zl:
    zl = cir.zlcir
  if not zu:
    zu = cir.zucir
  # --- Create a workspace
  varwork = fzeros(shape(varsave),'d')
  # --- Do iterations
  for i in xrange(n):
    # --- Run the beam through to get the final data
    varwork[:,:] = varsave
    cirrun(cir.nit,varwork,zl,zu,cir.dscir,cir.nscir,
           top.aion,top.zion,cir.icharge,false,cir.lperveance,
           cir.lemittance,false,cir.llinear,cir.limage,cir.lendzero,savehist,
	   cir.lfixed)
    # --- Get new values for envelope
    varsave[0,:] = 0.5*(varsave[0,:] + varwork[0,:])
    varsave[1,:] = 0.5*(varsave[1,:] + varwork[1,:])
    varsave[2,:] = 0.5*(varsave[2,:] + varwork[2,:])
    varsave[3,:] = 0.5*(varsave[3,:] + varwork[3,:])
    error = 0.
    for j in range (4): error = error + max(abs(varsave[j,:]-varwork[j,:]))
    if error < errorlimit:
      print i,'steps in match1'
      return
  print 'No convergence in match1 after', n, 'steps.'
  print 'Final error =', error

# --- This is similar to the above routine but forces the beam to beam
# --- to be round (a=b, and a'=-b').
def match2(n=1,varsave=None,zl=None,zu=None,savehist=false,errorlimit=1.e-6):
  if not varsave:
    # --- If varsave is not set, then create a beam in the standard circe
    # --- database and make varsave point to it.
    initbeam()
    varsave = cir.var
  if not zl:
    zl = cir.zlcir
  if not zu:
    zu = cir.zucir
  # --- Create a workspace
  varwork = fzeros(shape(varsave),'d')
  # --- Do iterations
  for i in xrange(n):
    # --- Run the beam through to get the final data
    varwork[:,:] = varsave
    cirrun(cir.nit,varwork,zl,zu,cir.dscir,cir.nscir,
           top.aion,top.zion,cir.icharge,false,cir.lperveance,
           cir.lemittance,false,cir.llinear,cir.limage,cir.lendzero,savehist,
	   cir.lfixed)
    # --- Get new values for envelope
    varsave[0,:] = 0.5*(varsave[0,:] + varwork[0,:])
    varsave[1,:] = 0.5*(varsave[1,:] + varwork[1,:])
    varsave[2,:] =  varsave[0,:]
    varsave[3,:] = -varsave[1,:]
    error = max(abs(varsave[0,:]-varwork[0,:])) + max(abs(varsave[1,:]-varwork[1,:]))
    if error < errorlimit:
      print i,'steps in match2'
      return
  print 'No convergence in match2 after', i, 'steps.'
  print 'Final error =', error

# This routine does a regula-falsi iteration to find the emittance which
# gives the desired sigma.
def match(n=1,_sig_desr=20.,savehist=false):
  fff = zeros((2,cir.nit),'d')
  xxx = zeros((2,cir.nit),'d')
  nxxx = zeros(cir.nit,'d')
  ixxx = zeros(cir.nit)
  sigcir = zeros(cir.nit,'d')
  sig0cir = zeros(cir.nit,'d')
  sig_desr = _sig_desr*ones(cir.nit,'d')
  # --- Get starting values
  xxx[0,:] = cir.var[10,:]
  xxx[1,:] = cir.var[10,:]*1.05
  # --- Calculate f at starting values
  cir.var[10,:] = xxx[0,:];cir.var[11,:] = xxx[0,:]
  getphaset(cir.nit,sig0cir,sigcir,hlp,cir.zlcir,cir.var,top.aion,top.zion,
            cir.lperveance)
  sig0cir = sig0cir*180./top.pi ; sigcir = sigcir*180./top.pi
  fff[0,:] = sigcir - sig_desr
  cir.var[10,:] = xxx[1,:];cir.var[11,:] = xxx[1,:]
  getphaset(cir.nit,sig0cir,sigcir,hlp,cir.zlcir,cir.var,top.aion,top.zion,
            cir.lperveance)
  sig0cir = sig0cir*180./top.pi ; sigcir = sigcir*180./top.pi
  fff[1,:] = sigcir - sig_desr
  j=1
  while (j <= n and abs(min(min(fff))) > 1.e-10):
    j = j + 1
    nxxx = (xxx[0,:]*fff[1,:] - xxx[1,:]*fff[0,:])/(fff[1,:] - fff[0,:])
    for k in range(cir.nit):
      #ixxx[k] = mxx(abs(fff[:,k]))
      ixxx[k] = 1
      if abs(fff[2,k]) > abs(fff[1,k]): ixxx[k] = 2
      xxx[ixxx[k],k] = nxxx[k]
      cir.var[10,:] = xxx[ixxx[k],:];cir.var[11,:] = xxx[ixxx[k],:]
    getphaset(cir.nit,sig0cir,sigcir,hlp,cir.zlcir,cir.var,top.aion,top.zion,
              cir.lperveance)
    sig0cir = sig0cir*180./top.pi ; sigcir = sigcir*180./top.pi
    for k in range(cir.nit):
      fff[ixxx[k],:] = sigcir - sig_desr
    print "Error = ",max(max(fff))


