from warp import *
# FIND_MGPARAM
find_mgparam_version = "$Id: find_mgparam.py,v 1.6 2002/03/14 19:39:46 jlvay Exp $"
# Author: D. P. Grote, March 1995
# Converted to python: April 1999
# Support for RZ multigrid solver added by J.-L. Vay: March 2002
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

def find_mgparam():
  """
Optimize both mgparam and up and down passes, minimizing the fieldsolve
time.
  """
  # --- Do some error checking first
  if f3d.downpasses == 0: f3d.downpasses = 1
  if f3d.uppasses == 0: f3d.uppasses = 1
  # --- Get initial field solve time
  if(w3d.solvergeom == w3d.RZgeom):
    f3d.downpasses = frz.mgridrz_npre 
    f3d.uppasses = frz.mgridrz_npost 
    f3d.mgparam = frz.mgridrz_mgparam 
    f3d.mgmaxiters = frz.mgridrz_ncmax 
  nexttime = field_solve()
  prevtime = 2*nexttime
  # --- Loop, increasing the number of passes until the time is minimized.
  while nexttime < prevtime:
    prevparam = f3d.mgparam
    prevtime = nexttime
    nexttime = _find_mgparam()
    print "Field solve time = ",nexttime
    if(w3d.solvergeom == w3d.XYZgeom):
      print "f3d.mgparam = ",f3d.mgparam
      print "f3d.downpasses = ",f3d.downpasses
      print "f3d.uppasses = ",f3d.uppasses
    else:
      print "frz.mgridrz_mgparam = ",frz.mgridrz_mgparam
      print "frz.mgridrz_npre    = ",frz.mgridrz_npre
      print "frz.mgridrz_npost   = ",frz.mgridrz_npost
      if f3d.mgiters == f3d.mgmaxiters:prevtime=2*nexttime
    if nexttime < prevtime:
      f3d.downpasses = f3d.downpasses + 1
      f3d.uppasses = f3d.uppasses + 1
    else:
      # --- Reset the values to the previous ones (which were the best)
      f3d.mgparam = prevparam
      f3d.downpasses = f3d.downpasses - 1
      f3d.uppasses = f3d.uppasses - 1
      # --- Do some error checking first
      if f3d.downpasses == 0: f3d.downpasses = 1
      if f3d.uppasses == 0: f3d.uppasses = 1
  # --- Print error message if maximum iterations is reached.
  if f3d.mgiters == f3d.mgmaxiters:
    print """Notice: the maximum number of iterations has been reached, so
the values above are unlikely to be optimal. Try increasing the
tolerance, increasing the maximum number of iterations, or making a
better initial guess of mgparam."""
  else:
    if(w3d.solvergeom == w3d.RZgeom):prevtime=findnrecursmin(prevtime)
    print "-----------------------------------------"
    print "The optimized values:"
    print "Field solve time = ",prevtime
    if(w3d.solvergeom == w3d.XYZgeom):
      print "f3d.mgparam = ",f3d.mgparam
      print "f3d.downpasses = ",f3d.downpasses
      print "f3d.uppasses = ",f3d.uppasses
    else:
      print "frz.mgridrz_mgparam     = ",frz.mgridrz_mgparam
      print "frz.mgridrz_npre        = ",frz.mgridrz_npre
      print "frz.mgridrz_npost       = ",frz.mgridrz_npost
      print "frz.mgridrz_nrecurs_min = ",frz.mgridrz_nrecurs_min
    
def field_solve():
  ixmin = 1
  ixmax = w3d.nx-1
  iymin = 1
  iymax = w3d.ny-1
  izmin = 1
  izmax = w3d.nz-1
  if (f3d.boundxy > 0 or w3d.l2symtry or w3d.l4symtry): ixmin = 0
  if (f3d.boundxy > 0): ixmax = w3d.nx
  if (f3d.boundxy > 0 or w3d.l4symtry): iymin = 0
  if (f3d.boundxy > 0): iymax = w3d.ny
  if (f3d.bound0  > 0): izmin = 0
  if (f3d.boundnz > 0): izmax = w3d.ny
  if(w3d.solvergeom == w3d.XYZgeom):
    w3d.phi[ixmin:ixmax+1,iymin:iymax+1,izmin+1:izmax+1] = 0.
  else:
    w3d.phi[ixmin:ixmax+1,:,izmin+1:izmax+1] = 0.
    frz.mgridrz_npre    = f3d.downpasses
    frz.mgridrz_npost   = f3d.uppasses
    frz.mgridrz_mgparam = f3d.mgparam
    frz.mgridrz_ncmax   = f3d.mgmaxiters
    
  beforetime = wtime()
  vp3d(-1)
  aftertime = wtime()
  if(w3d.solvergeom == w3d.RZgeom):
    f3d.mgiters = frz.mgridrz_mgiters
    f3d.downpasses = frz.mgridrz_npre 
    f3d.uppasses = frz.mgridrz_npost
    f3d.mgparam = frz.mgridrz_mgparam
  return aftertime - beforetime

def _find_mgparam():
  icount = 0  # iteration count

# --- Make sure that mgparam is between 0 and 2.
# --- If mgparam is less then zero, put mgparam closer to 2 since the
# --- optimal value is always closer to 2 than to 0.
  if (f3d.mgparam <= 0.):
    f3d.mgparam = max(1., 2. + f3d.mgparam)

  if(w3d.solvergeom == w3d.RZgeom):f3d.mgparam=1.8
 
# --- If mgparam is greater than two, put it on the other side of two
# --- and reduce the increment.  This keeps mgparam near two.
  if (f3d.mgparam > 2.):
    f3d.mgparam = max(1., 4. - f3d.mgparam)

# --- do initial field solve
  fstime = field_solve()

# --- set initail values for 'previous' quantities
  mgparam_prev = f3d.mgparam
  mgiters_prev = f3d.mgiters

# --- set initial increment for mgparam
  sincr = .05

# --- set mgiters to 0 so that while loop executes at least once
  f3d.mgiters = 0

# --- increment mgparam so next field solve uses new mgparam
  f3d.mgparam = f3d.mgparam + sincr

# --- Execute while loop until two iterations give the same number of field
# --- solve iterations or until a maximum number of iterations has been
# --- reached.
  while (mgiters_prev != f3d.mgiters and icount < 200):

#   --- print out current value of mgparam
    print "Best parameter so far = %f" % f3d.mgparam

#   --- do field solve (which prints out number of field solve iterations)
    up_old = f3d.uppasses
    fstime = field_solve()

#   --- If field solve took more iterations than previous field solve, change
#   --- direction of the increment and reduce its size.  Reducing its size
#   --- removes the possibility of an infinite loop.
    if (f3d.mgiters > mgiters_prev):
      sincr = - sincr/2.
      f3d.mgparam = mgparam_prev + sincr

#   --- If a smaller number of field solve iterations was returned, then
#   --- reset the previous values and keep changing mgparam in the same
#   --- direction.  The previous number of iterations is saved in mgiters
#   --- temporarily to check if the two iterations had the same number
#   --- of field solver iterations.
    elif (f3d.mgiters < mgiters_prev):
      s = mgiters_prev
      mgiters_prev = f3d.mgiters
      f3d.mgiters = s
      mgparam_prev = f3d.mgparam
      f3d.mgparam = mgparam_prev + sincr

#   --- Make sure that mgparam stays between 0.01 and 2.  .01 is used instead
#   --- of zero since when mgparam is too close to zero, misleading things
#   --- happen.
#   --- If mgparam is outside the range, start the iterations over at a
#   --- random place near the typical optimum value, 1.9. (1.8 for RZ solver)
    if (f3d.mgparam <= 0.01 or 2. < f3d.mgparam):
      if(w3d.solvergeom == w3d.XYZgeom):
        f3d.mgparam = 1.9 + ranf()*.05
      else:
        f3d.mgparam = 1.8 + ranf()*.05
      sincr = .01

#   --- increment iteration counter
    icount = icount + 1

    if(w3d.solvergeom == w3d.RZgeom and f3d.uppasses <> up_old):
      print "resetting _find_mgparam"
      icount=0
      fstime = field_solve()
      mgparam_prev = 1.8
      mgiters_prev = f3d.mgiters
      sincr = .05
      f3d.mgiters = 0
      f3d.mgparam = f3d.mgparam + sincr


# --- print message if an optimal value wasn't found
  if (icount == 200):
    print "Warning: maximum number of iterations reached."
    print "         The value of mgparam may not be optimal."
    print "         Try increasing mgmaxit."


  return fstime

def findnrecursmin(prevtime):
  """
Optimize nrecurs_min, minimizing the fieldsolve time.
  """
  # --- Get initial field solve time
  nexttime = prevtime
  prevtime = 2*nexttime
  # --- Loop, increasing the number of passes until the time is minimized.
  while nexttime < prevtime and frz.mgridrz_nrecurs_min < frz.mgridrz_nlevels_max:
    prevtime = nexttime
    frz.mgridrz_nrecurs_min = frz.mgridrz_nrecurs_min + 1
    nexttime = field_solve()
    print "Field solve time = ",nexttime
    print "frz.mgridrz_nrecurs_min = ",frz.mgridrz_nrecurs_min
    if f3d.mgiters == f3d.mgmaxiters:prvtime=2*nexttime
    if nexttime > prevtime:
      # --- Reset the values to the previous ones (which were the best)
      frz.mgridrz_nrecurs_min = frz.mgridrz_nrecurs_min - 1
      # --- Do some error checking first
      if frz.mgridrz_nrecurs_min == 0: frz.mgridrz_nrecurs_min = 1
  return prevtime
