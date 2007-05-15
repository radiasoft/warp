from warp import *
if lwithnumpy:
  import numpy.oldnumeric.linear_algebra as LinearAlgebra
else:
  import LinearAlgebra
import singleparticle
wxy_match_version = "$Id: wxy_match.py,v 1.8 2007/05/15 18:37:29 dave Exp $"

def wxy_matchdoc():
  print """
Contains tools for matching beams using the WARPxy code.
minit: re-generates the slice beam
match: varies the beam parameters to finds a matched beam (first method)
match1: varies the beam parameters to finds a matched beam (second method)
matchx: varies single particle starting conditions to match final conditions
MatchXY: class to match a beam from a given set of beam parameters to a final
         set of beam parameters by changing the focusing strength on four
         quadrupoles. So far, only linear quadrupole elements can be used.
  """


##############################################################################
# Iterate toward a matched beam using the slice code (instead of the env     #
# package).  It is currently setup to match a beam starting from the center  #
# of a quadrupole.  To match a beam under other conditions, change the lines #
# indicated below.                                                           #
##############################################################################

# --- Set diagnostic parameters to minimize output

# --- This is the critical routine which resets certain variables and
# --- redoes a generate.
zbeam_start = top.zbeam
time_start = top.time
it_start = top.it
np_s_start = None
def minit():
  """Re-inits the slice package"""
  global np_s_start
  if type(np_s_start) == type(None): np_s_start = top.np_s + 0
  top.lprntpara = false
  top.verbosity = 0
  top.zbeam = zbeam_start
  top.zgrid = top.zbeam
  top.zgridprv = top.zbeam
  top.it = it_start
  top.time = time_start
  top.np_s = np_s_start
  top.npmax = sum(np_s_start)
  generate()


# --- This is the function which actually does the iteration
def match(imtch=1,s=128):
  """
Matches the envelope using WARPxy, the beam parameters are set to the
average of the initial and final values
  - imtch=1 the number of iterations to perform
  - s=128 the number of time steps across the region to be matched
  """
  for i in xrange(imtch):
    top.a0=0.5*(top.a0+2.*top.xrms[0,-1])
    top.b0=0.5*(top.b0+2.*top.yrms[0,-1])
    top.ap0=sign(0.5*(abs(top.ap0) +
                 2.*(top.xxpbar[0,-1]-top.xbar[0,-1]*top.xpbar[0,-1])/
                 top.xrms[0,-1]),top.ap0)
    top.bp0=sign(0.5*(abs(top.bp0) +
                 2.*(top.yypbar[0,-1]-top.ybar[0,-1]*top.ypbar[0,-1])/
                 top.yrms[0,-1]),top.bp0)
    top.a0_s = top.a0
    top.b0_s = top.b0
    top.ap0_s = top.ap0
    top.bp0_s = top.bp0
    minit()
    step(s)
    print ("a error = %18.13f a' error = %18.13f"%
         (2.*top.xrms[0,-1]-top.a0,2.*top.vxrms[0,-1]/top.vzbar[0,-1]-top.ap0))
    print ("b error = %18.13f b' error = %18.13f"%
         (2.*top.yrms[0,-1]-top.b0,2.*top.vyrms[0,-1]/top.vzbar[0,-1]+top.bp0))

def match1(imtch=1,s=128):
  """
Matches the envelope using WARPxy, the X beam parameters are set to the
average of the initial and final values, then Y is set equal to X
  - imtch=1 the number of iterations to perform
  - s=128 the number of time steps across the region to be matched
  """
  for i in xrange(imtch):
    top.a0=0.5*(top.a0+2.*top.xrms[0,-1]) #xrms[0,-1]+yrms[0,-1]
    top.b0=top.a0
    top.ap0=sign(0.5*(abs(top.ap0) +
             2.*(top.xxpbar[0,-1]-top.xbar[0,-1]*top.xpbar[0,-1])/top.xrms[0,-1]),top.ap0)
    top.bp0=sign(0.5*(abs(top.bp0) +
             2.*(top.yypbar[0,-1]-top.ybar[0,-1]*top.ypbar[0,-1])/top.yrms[0,-1]),top.bp0)
    top.a0_s = top.a0
    top.b0_s = top.b0
    top.ap0_s = top.ap0
    top.bp0_s = top.bp0
    minit()
    step(s)
    print ("a error = %18.13f a' error = %18.13f"%
         (2.*top.xrms[0,-1]-top.a0,2.*top.vxrms[0,-1]/top.vzbar[0,-1]-top.ap0))
    print ("b error = %18.13f b' error = %18.13f"%
         (2.*top.yrms[0,-1]-top.b0,2.*top.vyrms[0,-1]/top.vzbar[0,-1]+top.bp0))



# --- This is an unrelated function which can be used to calculate sigma0
# --- from the particle in the full applied field.  It uses the calc_sig0
# --- script and resets several variables and turns the moments
# --- calculation off to reduce computation time.
#from calc_sig0 import *
#def go():
#  top.itmomnts = 0
#  top.zbeam = 0.;top.zgrid = 0.;top.it = 0;top.time = 0
#  setlatt()
#  calc_sig0()





########################################################################
def matchx(xf=0.,xpf=0.,yf=0.,ypf=0.,zs=None,ze=None,s=None,
           maxiter=100,tol=1.e-10,vary=0.01,
	   xi=0.,xpi=0.,yi=0.,ypi=0.,xs=1.,xps=1.,ys=1.,yps=1.):
  """Matches a single particle orbit to the specified values
     - xf=0. final value of x
     - xpf=0. final value of x'
     - yf=0. final value of y
     - ypf=0. final value of y'
     - zs=env.zl starting value of z
     - ze=env.zu final value of z
     - s=env.nenv number of time steps across region of interest
     - maxiter=100 maximum number of iterations
     - tol=1.e-10 tolerance position is matched to"""
  mat = zeros((4,4),'d')
  dsp = zeros(4,'d')

  # --- Get number of steps if not specified
  if not s:
    if not zs:
      zs = env.zl
    if not ze:
      ze = env.zu
    s = nint((ze - zs)/wxy.ds)

  notdone = 1
  iter = 0
  while notdone and iter < maxiter:
    iter = iter + 1

    xx = top.x0
    yy = top.y0
    vx = top.xp0*top.vbeam
    vy = top.yp0*top.vbeam
    singleparticle.spinit([xx,xx*(1. + vary),xx,xx,xx],
                          [yy,yy,yy,yy*(1. + vary),yy],
                          [zs,zs,zs,zs,zs],
                          [vx,vx,vx*(1. + vary),vx,vx],
                          [vy,vy,vy,vy,vy*(1. + vary)],
                          5*[top.vbeam],5*[1.],s)
    step(s)

    for ip in xrange(4):
      mat[0,ip]=(top.pgroup.xp[ip+1]  - top.pgroup.xp[0])/((xx-xi)/xs*vary)
      mat[1,ip]=((top.pgroup.uxp[ip+1] - top.pgroup.uxp[0])/
                 ((vx-xpi*top.vbeam)/xps*vary))
      mat[2,ip]=(top.pgroup.yp[ip+1]  - top.pgroup.yp[0])/((yy-yi)/ys*vary)
      mat[3,ip]=((top.pgroup.uyp[ip+1] - top.pgroup.uyp[0])/
                 ((vy-ypi*top.vbeam)/yps*vary))

    # --- Invert the matrix
    mati = LinearAlgebra.inverse(mat)

    # --- Get next starting point
    dsp[:] = [xf - top.pgroup.xp[0],xpf - top.pgroup.uxp[0]/top.pgroup.uzp[0],
              yf - top.pgroup.yp[0],ypf - top.pgroup.uyp[0]/top.pgroup.uzp[0]]
    top.x0 = top.x0 + sum(mati[0,:]*dsp)*xs
    top.xp0 = top.xp0 + sum(mati[1,:]*dsp)*xps
    top.y0 = top.y0 + sum(mati[2,:]*dsp)*ys
    top.yp0 = top.yp0 + sum(mati[3,:]*dsp)*yps

    # --- Check for convergence
    print "error => x = %10.3e x' = %10.3e"%tuple(dsp[:2])
    print "error => y = %10.3e y' = %10.3e"%tuple(dsp[2:])
    if max(abs(dsp)) < tol : notdone = 0

  if iter == maxiter:
    print 'Warning: Maximum number of iterations reached'
  else:
    print "top.x0 = %20.15e;top.xp0 = %20.15e"%(top.x0,top.xp0)
    print "top.y0 = %20.15e;top.yp0 = %20.15e"%(top.y0,top.yp0)






#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

class MatchXY:
  """
Matches a beam to final conditions.
  - afnew: final value of a
  - bfnew: final value of b
  - apfnew: final value of a'
  - bpfnew: final value of b'
  - runsimulation: is a Python function that runs the WarpXY simulation from
                   the starting location to the final location. Before calling
                   iterate, WarpXY should be initialized such that the
                   simulation can run 
  - iquads=[0,1,2,3]: indeces of the quad elements to be varied
  - tolx=1.0e-4: Minimum relative change in quad strengths for the iteration
                 to continue
  - alf=1.0e-4: Ensures sufficient decrease in function value
  - eps=1.0e-8: Approximate square root of the machine precision.
  - stpmax=400.0:
Available functions:
  - iterate(n): perform the iterations
Available data:
  - herror: list of the errors after each iteration
    """
  def __init__(s,afnew,bfnew,apfnew,bpfnew,runsimulation,iquads=[0,1,2,3],
                 tolx=1.0e-4,alf=1.0e-4,eps=1.0e-8,stpmax=400.0):
    # --- Save the beam parameters to be matched to
    s.af = afnew
    s.bf = bfnew
    s.apf = apfnew
    s.bpf = bpfnew
    # --- Find the quads to be varied
    s.iquads = iquads
    s.savedquadde = zeros (4,'d')
    s.savedquaddb = zeros (4,'d')
    try:
      s.savedquadde = take(top.quadde,s.iquads)
      s.savedquaddb = take(top.quaddb,s.iquads)
    except IndexError:
      raise "Error: iquads index out of bounds"
    # --- Save initial aperture and grid size
    s.initialaperture = top.prwall
    s.initialxmmax = w3d.xmmax
    # --- Save stepper function
    s.propagate = runsimulation
    # --- Set some constants
    s.tolx = tolx
    s.alf = alf
    s.eps = eps
    s.stpmax = stpmax
    # --- Initial step has not been done yet
    s.firststepcomplete = 0
    s.niterations = 0
    s.herror = []

  def iterate(s,n=1):
    """Perform the iterations. Optional argument is the number of iterations
to do, defaulting to 1."""
    # --- Do initial iteration if it hasn't been done yet.
    if not s.firststepcomplete:
      error = s.firststep()
      s.herror.append(error)
      s.firststepcomplete = 1
      print 'Initial error is %f'%(error)
    # --- Do more iterations
    for i in xrange(n):
      error = s.nextstep()
      s.herror.append(error)
      s.niterations = s.niterations + 1
      print 'After %d iterations, error = %f'%(s.niterations,error)

  def reset(s,initialaperture,initialxmmax):
    # --- Re-initialize the slice package
    top.lprntpara = false
    top.verbosity = 0
    top.zbeam = 0.
    top.zgrid = 0.
    top.zgridprv = 0.
    top.it = 0.
    top.time = 0.
    top.jhist = -1
    # --- Set grid size in the starting channel
    top.prwallz = initialaperture
    w3d.xmmax =  initialxmmax
    w3d.ymmax =  w3d.xmmax
    w3d.xmmin = -w3d.xmmax
    w3d.ymmin = -w3d.ymmax
    if (w3d.l2symtry):
      w3d.ymmin = 0
    if (w3d.l4symtry):
      w3d.xmmin = 0
      w3d.ymmin = 0
    w3d.dx = (w3d.xmmax - w3d.xmmin)/w3d.nx
    w3d.dy = (w3d.ymmax - w3d.ymmin)/w3d.ny
    # --- Generate the XY code (allocate storage, load ptcls, t=0 plots, etc.).
    w3d.nz = 2
    w3d.zmmin = 0.
    w3d.zmmax = 0. # Needed for parallel running
    package ("wxy"); generate()
    fieldsolxy(0)

  def firststep(s):
    # --- Calculate the initial error
    print '### Calculating the initial error...'
    s.reset(s.initialaperture,s.initialxmmax)
    s.propagate()
    if top.nplive < top.npmax:
      print "Particles lost, remaining number is", top.nplive
    s.fvec = array([2*top.xrms[0,-1]-s.af,
                    2*top.yrms[0,-1]-s.bf,
                    2*top.xxpbar[0,-1]/top.xrms[0,-1]-s.apf,
                    2*top.yypbar[0,-1]/top.yrms[0,-1]-s.bpf])
    s.f = 0.5 * sum (s.fvec*s.fvec)
    error = sqrt(2.*s.f)
    s.scaling = array ([1.0, 1.0, 1.0, 1.0])
    return error

  def nextstep(s):
    """
    Performs one iteration step, and returns the error.
    No arguments are needed.
    """
    # --- Vary each quad
    print '### Calculating the jacobian matrix'
    fjac = zeros ((4,4),'d')
    g = zeros ((4),'d')
    p = zeros ((4),'d')
    xold = zeros ((4),'d')
    for i in range (4):
      temp = s.scaling[i]
      h = s.eps * abs(temp)
      if (h == 0.0): h = s.eps
      # --- Trick to reduce finite precision error.
      s.scaling[i] = temp + h
      h = s.scaling[i] - temp
      for j in range (4):
        iq = s.iquads[j]
        top.quadde[iq] = s.savedquadde[j] * s.scaling[j]
        top.quaddb[iq] = s.savedquaddb[j] * s.scaling[j]
      s.reset(s.initialaperture,s.initialxmmax)
      s.propagate()
      if top.nplive < top.npmax:
        print "Particles lost, remaining number is", top.nplive
      fnew = array([2*top.xrms[0,-1]-s.af,
                    2*top.yrms[0,-1]-s.bf,
                    2*top.xxpbar[0,-1]/top.xrms[0,-1]-s.apf,
                    2*top.yypbar[0,-1]/top.yrms[0,-1]-s.bpf])
      s.scaling[i] = temp
      # --- Forward difference formula.
      fjac[:,i] = (fnew[:]-s.fvec[:]) / h
      g[i] = sum (fjac[:,i]*s.fvec[:])
    xold[:] = s.scaling[:]
    p[:] = -s.fvec[:]
    p = LinearAlgebra.solve_linear_equations (fjac, p)
    # --- Starting linesearch
    fold = s.f
    total = sqrt(sum(p[:]*p[:]))
    # --- Scale if attempted step is too big.
    if (total>s.stpmax): p[:] = p[:] * s.stpmax / total
    slope = sum (g[:]*p[:])
    if (slope >= 0.0): raise 'Roundoff problem in lnsrch'
    junk = zeros ((4),'d')
    for i in range(4): junk[i] = abs (p[i]) / max (abs(xold[i]),1.0)
    # --- Compute lambda_min
    test = max (junk[:])
    alamin = s.tolx / test
    # --- Always try full Newton step first.
    alam = 1.0
    f2 = 0
    alam2 = 0
    completed = 0
    # --- Start of iteration loop.
    while (not completed):
      tmplam = 0
      for i in range(4):
        iq = s.iquads[i]
        s.scaling[i] = xold[i] + alam * p[i]
        top.quadde[iq] = s.savedquadde[i] * s.scaling[i]
        top.quaddb[iq] = s.savedquaddb[i] * s.scaling[i]
      print '### Running with new quad settings'
      s.reset(s.initialaperture,s.initialxmmax)
      s.propagate()
      if top.nplive < top.npmax:
        print "Particles lost, remaining number is", top.nplive
      s.fvec = array([2*top.xrms[0,-1]-s.af,
                      2*top.yrms[0,-1]-s.bf,
                      2*top.xxpbar[0,-1]/top.xrms[0,-1]-s.apf,
                      2*top.yypbar[0,-1]/top.yrms[0,-1]-s.bpf])
      total = sum (s.fvec[:]*s.fvec[:])
      s.f = 0.5 * total
      # --- Convergence on delta x.
      if (alam < alamin): completed = 1
      # --- Sufficient function decrease
      elif (s.f <= fold + s.alf * alam * slope): completed = 1
      # --- Backtrack first time
      elif (alam == 1.0): tmplam = - slope / (2.0 * (s.f-fold-slope))
      # --- Subsequent backtracks.
      else:
        rhs1 = s.f - fold - alam * slope
        rhs2 = f2 - fold - alam2 * slope
        a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2)
        b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2)
        if (a == 0.0): tmplam = - slope / (2.0 * b)
        else:
          disc = b*b-3.0*a*slope
          if (disc < 0.0): tmplam = 0.5 * alam
          elif (b <= 0.0): tmplam = (-b + sqrt(disc))/(3.0*a)
          else: tmplam = - slope / (b + sqrt(disc))
        # --- lambda <= 0.5 lambda1
        if (tmplam > 0.5 * alam): tmplam = 0.5 * alam
      alam2 = alam
      f2 = s.f
      # --- lambda >= 0.1 lambda1
      alam = max (tmplam,0.1*alam)
    # --- End of linesearch
    error = sqrt(2.*s.f)
    return error

