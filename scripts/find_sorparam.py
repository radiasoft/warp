from warp import *
find_sorparam_version = "$Id: find_sorparam.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"
# FIND_SORPARAM
# Author: D. P. Grote, March 1995
# Converted to python: April 1999
# This script optimizes the value of sorparam, the successive-overrelaxation
# parameter used in the PSOR field solver.  It begins its search with the
# current value of sorparam and moves it initially in increments of .01.  The
# increment is reduced each iteration as the script narrows in on the optimized
# version.
#
# This may not be the most efficient search method, but it was fairly easy
# to write and it works.  It typically takes about ten iterations.
#
# For each iteration, the value of sorparam is printed and then the number
# of field solve iterations.  The last value printed is the optimized value
# and sorparam retains that value upon completion.
#
# An attempt is made to deal with the case where the maximum number of
# field solve iterations is reached.  The search is based of the value of the
# maximum error from the field solve after sormaxit iterations instead of the
# number of iterations.
#
# This script was tried for a number of different starting values of sorparam
# and it converged every time.  An effort was made to make the script robust
# and idiot proof.

# This function carries out the field solve. Also zeroing out the
# appropriate portion of phi.
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
  w3d.phi[ixmin:ixmax+1,iymin:iymax+1,izmin+1:izmax+1] = 0.
  vp3d(-1)

def find_sorparam():
  icount = 0  # iteration count

  # --- Make sure that sorparam is between 0 and 2.
  # --- If sorparam is less then zero, put sorparam closer to 2 since the
  # --- optimal value is always closer to 2 than to 0.
  if (f3d.sorparam <= 0.):
    f3d.sorparam = max(1., 2. + f3d.sorparam)

  # --- If sorparam is greater than two, put it on the other side of two
  # --- and reduce the increment.  This keeps sorparam near two.
  if (f3d.sorparam > 2.):
    f3d.sorparam = max(1., 4. - f3d.sorparam)

  # --- do initial field solve
  field_solve()

  # --- set initail values for 'previous' quantities
  sorparam_prev = f3d.sorparam
  soriter_prev = f3d.soriter
  sorerror_prev = f3d.sorerror

  # --- set initial increment for sorparam
  sincr = .01

  # --- set soriter to 0 so that while loop executes at least once
  f3d.soriter = 0

  # --- increment sorparam so next field solve uses new sorparam
  f3d.sorparam = f3d.sorparam + sincr

  # --- Execute while loop until two iterations give the same number of field
  # --- solve iterations or until a maximum number of iterations has been
  # --- reached.
  while (soriter_prev != f3d.soriter and icount < 200):

    # --- print out current value of sorparam
    print "Best parameter so far = %f" % f3d.sorparam

    # --- do field solve (which prints out number of field solve iterations)
    field_solve()

    # --- Check if maximum number of field solver iterations was reached.
    # --- If so, base search in value of sorerror.  In this case, the
    # --- increment is continually increased so as not to get stuck in a local
    # --- minimum that is not optimal.
    if (f3d.soriter == f3d.sormaxit and soriter_prev >= f3d.sormaxit):

      # --- Set an artificial value to soriter_prev so that the loop keeps
      # --- iterating even though iterations give the same number of
      # --- field solve iterations.
      soriter_prev = f3d.sormaxit + 1

      # --- If field solve had higher error than the previous, then change the
      # --- sign of the increment and increase its linearly increase its size.
      if (f3d.sorerror > sorerror_prev):
        sincr = - sincr - sign(.01,sincr)
        f3d.sorparam = sorparam_prev + sincr

      # --- If a smaller error was returned, then reset the previous values.
      # --- The increment is increased linearly to try to get to the optimal
      # --- value faster.
      elif (f3d.soriter < soriter_prev):
        s = soriter_prev
        soriter_prev = f3d.soriter
        f3d.soriter = s
        sorerror_prev = f3d.sorerror
        sorparam_prev = f3d.sorparam
        sincr = sincr + sign(.01,sincr)
        f3d.sorparam = sorparam_prev + sincr

    # --- If field solve took more iterations than previous field solve, change
    # --- direction of the increment and reduce its size.  Reducing its size
    # --- removes the possibility of an infinite loop.
    elif (f3d.soriter > soriter_prev):
      sincr = - sincr/2.
      f3d.sorparam = sorparam_prev + sincr

    # --- If a smaller number of field solve iterations was returned, then
    # --- reset the previous values and keep changing sorparam in the same
    # --- direction.  The previous number of iterations is saved in soriter
    # --- temporarily to check if the two iterations had the same number
    # --- of field solver iterations.
    elif (f3d.soriter < soriter_prev):
      s = soriter_prev
      soriter_prev = f3d.soriter
      f3d.soriter = s
      sorerror_prev = f3d.sorerror
      sorparam_prev = f3d.sorparam
      f3d.sorparam = sorparam_prev + sincr

    # --- Make sure that sorparam stays between 0.01 and 2.  .01 is used instead
    # --- of zero since when sorparam is too close to zero, misleading things
    # --- happen.
    # --- If sorparam is outside the range, start the iterations over at a
    # --- random place near the typical optimum value, 1.9.
    if (f3d.sorparam <= 0.01 or 2. < f3d.sorparam):
      f3d.sorparam = 1.9 + ranf()*.05
      sincr = .01

    # --- increment iteration counter
    icount = icount + 1

  # --- print message if an optimal value wasn't found
  if (icount == 200):
    print "Warning: maximum number of iterations reached."
    print "         The value of sorparam may not be optimal."
    print "         Try increasing sormaxit."


