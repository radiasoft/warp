from warp import *
setupvalidation_version = "$Id: setupvalidation.py,v 1.4 2006/11/14 18:21:39 dave Exp $"

def setupvalidationdoc():
  print """
This module provides a convenient way of setting up a validation deck which
compares data from a previous run to a current run. It checks the sums of each
of the particles coordinates and the sums of the charge density and potential
arrays. To use, run the command "setupvalidator()" and follow its
instructions.
  """

##########################################################################
def getvalidationdata():
  """Returns the current values of the quantities to be compared."""
  if lparallel:
    # --- Sum only the non overlapping part of the 3-d arrays. Otherwise
    # --- the answers will be different with differing numbers of processors.
    iz1 = 0
    if me < npes-1:
      iz2 = top.izfsslave[me+1] - top.izfsslave[me]
    else:
      iz2 = iz1 + top.nzfsslave[me] + 1
    return (globalsum(sum(getx(gather=0))),
            globalsum(sum(gety(gather=0))),
            globalsum(sum(getz(gather=0))),
            globalsum(sum(getvx(gather=0))),
            globalsum(sum(getvy(gather=0))),
            globalsum(sum(getvz(gather=0)-top.vbeam)),
            globalsum(sum(1.-getgaminv(gather=0))),
            globalsum(sum(abs(getx(gather=0)))),
            globalsum(sum(abs(gety(gather=0)))),
            globalsum(sum(abs(getz(gather=0)))),
            globalsum(sum(abs(getvx(gather=0)))),
            globalsum(sum(abs(getvy(gather=0)))),
            globalsum(sum(abs(getvz(gather=0)-top.vbeam))),
            globalsum(sum(abs(1.-getgaminv(gather=0)))),
            globalsum(sum(sum(sum(w3d.rho[:,:,iz1:iz2])))),
            globalsum(sum(sum(sum(w3d.phi[:,:,iz1+1:iz2+1])))))
  else:
    return (sum(getx()),
            sum(gety()),
            sum(getz()),
            sum(getvx()),
            sum(getvy()),
            sum(getvz()-top.vbeam),
            sum(1.-getgaminv()),
            sum(abs(getx())),
            sum(abs(gety())),
            sum(abs(getz())),
            sum(abs(getvx())),
            sum(abs(getvy())),
            sum(abs(getvz()-top.vbeam)),
            sum(abs(1.-getgaminv())),
            sum(sum(sum(w3d.rho))),
            sum(sum(sum(w3d.phi))))

##########################################################################
def setupvalidator():
  """Prints out text which is to be inserted into the validation deck. The
text contains the values with which later runs will be compared."""
  output = """
Insert the text below into the place in the validation deck where the check
is to be done.

# ---
from setupvalidation import comparetooriginal
comparetooriginal([
%22.15e,%22.15e,%22.15e,
%22.15e,%22.15e,%22.15e,
%22.15e,
%22.15e,%22.15e,%22.15e,
%22.15e,%22.15e,%22.15e,
%22.15e,
%22.15e,%22.15e]
)
"""
  originaldata = getvalidationdata()
  print output%originaldata

##########################################################################
def comparetooriginal(odata):
  """Prints out the comparison between the original data and the current
  data."""
  fmt1 = "          %22.15e %22.15e"
  fmt2 = "original  %22.15e %22.15e"
  otext = ["x","y","z","vx","vy","vz-vbeam","1-gaminv",
           "x","y","z","vx","vy","vz-vbeam","1-gaminv",
           "rho","phi"]
  vdata = array(getvalidationdata())
  odata = array(odata)
  print "         %10s %-22s %-22s"%(" ","sum()","sum(abs())")
  for i in range(7):
    print "         %10s %22.15e %22.15e"%(otext[i],vdata[i],vdata[i+7])
    print "original %10s %22.15e %22.15e"%(otext[i],odata[i],odata[i+7])
    print

  print "         %10s %-22s"%(" ","sum()")
  print "         %10s %22.15e"%(otext[14],vdata[14])
  print "original %10s %22.15e"%(otext[14],odata[14])
  print
  print "         %10s %22.15e"%(otext[15],vdata[15])
  print "original %10s %22.15e"%(otext[15],odata[15])

# diffs = abs(odata - vdata)
# if max(diffs) > 0.:
#   print "Warning: The following values differ"
#   for i in range(len(vdata)):
#     if abs(odata[i] - vdata[i]) > 0.:
#       print "  ",otext[i]," by ",odata[i] - vdata[i]
# else:
#   print "No differences found"

