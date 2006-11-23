from warp import *
import time
# FIND_MGPARAM
find_mgparam_version = "$Id: find_mgparam.py,v 1.22 2006/11/23 01:06:35 dave Exp $"
# Author: D. P. Grote, March 1995
# Converted to python: April 1999
# This script optimizes the value of mgparam, the relaxation
# parameter used in the multigrid field solver.  It begins its search with the
# current value of mgparam and moves it initially in increments of .01.  The
# increment is reduced each iteration as the script narrows in on the optimized
# version.
#
# This may not be the most efficient search method, but it was fairly easy
# to write and it works.  It typically takes about ten iterations.
#
# For each iteration, the value of mgparam is printed and then the number
# of field solve iterations.  The last value printed is the optimized value
# and mgparam retains that value upon completion.
#
# An attempt is made to deal with the case where the maximum number of
# field solve iterations is reached.  The search is based of the value of the
# maximum error from the field solve after mgmaxiters iterations instead of the
# number of iterations.
#
# This script was tried for a number of different starting values of mgparam
# and it converged every time.  An effort was made to make the script robust
# and idiot proof.

def find_mgparam(lsavephi=false,resetpasses=0,solver=None,pkg3d=None):
  """
Optimize both mgparam and up and down passes, minimizing the fieldsolve
time.
 - lsavephi=false: When true, the current data in phi is saved and used as the
                   starting point for the field solver, rather than zeroing
                   out phi.
  """
  if solver is None:
    solver = getregisteredsolver()
    if solver is not None:
      solver.find_mgparam()
      return
  if solver is None:
    if(w3d.solvergeom == w3d.XYZgeom and top.fstype not in [7,13]):
      print "The fstype must be set to 7 or 13"
      return
    if(w3d.solvergeom == w3d.Zgeom): return
    if(w3d.solvergeom == w3d.RZgeom or w3d.solvergeom == w3d.XZgeom):
      frz.find_mgparam_rz(lsavephi)
      return
    solver = f3d
    if pkg3d is None: pkg3d = w3d
  else:
    if pkg3d is None: pkg3d = solver
    
  
  # --- Save the cuurrent value of phi to be used as the initial value
  if lsavephi:
    phisave = pkg3d.phi.copy()
  else:
    phisave = None

  if resetpasses:
    solver.downpasses = 1
    solver.uppasses = 1

  # --- Do some error checking first
  if solver.downpasses == 0: solver.downpasses = 1
  if solver.uppasses == 0: solver.uppasses = 1
  # --- Get initial field solve time
  nexttime = field_solve(phisave,solver,pkg3d)
  prevtime = 2*nexttime
  # --- Make sure that the number of iterations is below the maximum
  savemgtol = solver.mgtol
  if solver.mgiters == solver.mgmaxiters:
    try:
      solver.mgtol = 10*solver.mgerror
      nexttime = field_solve(phisave,solver,pkg3d)
      prevtime = 2*nexttime
    except:
      print """Notice: the maximum number of iterations has been reached, so
the values above are unlikely to be optimal. Try increasing the
tolerance, increasing the maximum number of iterations, or making a
better initial guess of mgparam."""
  # --- Loop, increasing the number of passes until the time is minimized.
  while nexttime < prevtime:
    prevparam = solver.mgparam
    prevtime = nexttime
    nexttime = _find_mgparam(phisave,solver,pkg3d)
    print "Field solve time = ",nexttime
    print "f3d.mgparam = ",solver.mgparam
    print "f3d.downpasses = ",solver.downpasses
    print "f3d.uppasses = ",solver.uppasses
    if nexttime < prevtime:
      solver.downpasses = solver.downpasses + 1
      solver.uppasses = solver.uppasses + 1
    else:
      # --- Reset the values to the previous ones (which were the best)
      solver.mgparam = prevparam
      solver.downpasses = solver.downpasses - 1
      solver.uppasses = solver.uppasses - 1
      # --- Do some error checking first
      if solver.downpasses == 0: solver.downpasses = 1
      if solver.uppasses == 0: solver.uppasses = 1
  # --- Print error message if maximum iterations is reached.
  if solver.mgiters == solver.mgmaxiters:
    print """Notice: the maximum number of iterations has been reached, so
the values above are unlikely to be optimal. Try increasing the
tolerance, increasing the maximum number of iterations, or making a
better initial guess of mgparam."""
  else:
    print "-----------------------------------------"
    print "The optimized values:"
    print "Field solve time = ",prevtime
    print "f3d.mgparam = ",solver.mgparam
    print "f3d.downpasses = ",solver.downpasses
    print "f3d.uppasses = ",solver.uppasses
  # --- Restore the value of mgtol
  solver.mgtol = savemgtol
    
def field_solve(phisave,solver,pkg3d):
  if phisave is None:
    ixmin = 1
    ixmax = pkg3d.nx-1
    iymin = 1
    iymax = pkg3d.ny-1
    izmin = 1
    izmax = pkg3d.nz-1
    if (pkg3d.boundxy > 0 or pkg3d.l2symtry or pkg3d.l4symtry): ixmin = 0
    if (pkg3d.boundxy > 0): ixmax = pkg3d.nx
    if (pkg3d.boundxy > 0 or pkg3d.l4symtry): iymin = 0
    if (pkg3d.boundxy > 0): iymax = pkg3d.ny
    if (pkg3d.bound0  > 0): izmin = 0
    if (pkg3d.boundnz > 0): izmax = pkg3d.nz
    pkg3d.phi[ixmin:ixmax+1,iymin:iymax+1,izmin+1:izmax+2] = 0.
  else:
    pkg3d.phi[:,:,:] = phisave
    
  if solver == f3d:
    beforetime = time.time()
    vp3d(-1)
    aftertime = time.time()
  else:
    beforetime = time.time()
    solver.solve(-1)
    aftertime = time.time()

  try:
    nslaves = solver.nslaves
  except AttributeError:
    nslaves = top.nslaves

  if nslaves <= 1:
    return aftertime - beforetime
  else:
    return globalsum(aftertime - beforetime)

def _find_mgparam(phisave,solver,pkg3d):
  icount = 0  # iteration count

# --- Make sure that mgparam is between 0 and 2.
# --- If mgparam is less then zero, put mgparam closer to 2 since the
# --- optimal value is always closer to 2 than to 0.
  if (solver.mgparam <= 0.):
    solver.mgparam = max(1., 2. + solver.mgparam)

# --- If mgparam is greater than two, put it on the other side of two
# --- and reduce the increment.  This keeps mgparam near two.
  if (solver.mgparam > 2.):
    solver.mgparam = max(1., 4. - solver.mgparam)

# --- do initial field solve
  fstime = field_solve(phisave,solver,pkg3d)

# --- set initail values for 'previous' quantities
  mgparam_prev = solver.mgparam
  mgiters_prev = solver.mgiters

# --- set initial increment for mgparam
  sincr = .05

# --- set mgiters to 0 so that while loop executes at least once
  solver.mgiters = 0

# --- increment mgparam so next field solve uses new mgparam
  solver.mgparam = solver.mgparam + sincr

# --- Execute while loop until two iterations give the same number of field
# --- solve iterations or until a maximum number of iterations has been
# --- reached.
  while (mgiters_prev != solver.mgiters and icount < 20):

#   --- print out current value of mgparam
    print "Best parameter so far = %f" % solver.mgparam

#   --- do field solve (which prints out number of field solve iterations)
    up_old = solver.uppasses
    fstime = field_solve(phisave,solver,pkg3d)

#   --- If field solve took more iterations than previous field solve, change
#   --- direction of the increment and reduce its size.  Reducing its size
#   --- removes the possibility of an infinite loop.
    if (solver.mgiters > mgiters_prev):
      sincr = - sincr/2.
      solver.mgparam = mgparam_prev + sincr

#   --- If a smaller number of field solve iterations was returned, then
#   --- reset the previous values and keep changing mgparam in the same
#   --- direction.  The previous number of iterations is saved in mgiters
#   --- temporarily to check if the two iterations had the same number
#   --- of field solver iterations.
    elif (solver.mgiters < mgiters_prev):
      s = mgiters_prev
      mgiters_prev = solver.mgiters
      solver.mgiters = s
      mgparam_prev = solver.mgparam
      solver.mgparam = mgparam_prev + sincr

#   --- Make sure that mgparam stays between 0.01 and 2.  .01 is used instead
#   --- of zero since when mgparam is too close to zero, misleading things
#   --- happen.
#   --- If mgparam is outside the range, start the iterations over at a
#   --- random place near the typical optimum value, 1.9.
    if (solver.mgparam <= 0.01 or 2. < solver.mgparam):
      solver.mgparam = 1.9 + ranf()*.05
      sincr = .01

#   --- increment iteration counter
    icount = icount + 1

# --- print message if an optimal value wasn't found
  if (icount == 20):
    print "Warning: maximum number of iterations reached."
    print "         The value of mgparam may not be optimal."
    print "         Try increasing mgmaxit."


  return fstime

