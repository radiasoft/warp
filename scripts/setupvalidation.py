from warp import *
setupvalidation_version = "$Id: setupvalidation.py,v 1.2 2001/03/22 23:24:04 dave Exp $"

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
    sumxp = globalsum(sum(abs(getx(gather=0))))
    sumyp = globalsum(sum(abs(gety(gather=0))))
    sumzp = globalsum(sum(abs(getz(gather=0))))
    sumvx = globalsum(sum(abs(getvx(gather=0))))
    sumvy = globalsum(sum(abs(getvy(gather=0))))
    sumvz = globalsum(sum(abs(getvz(gather=0)-top.vbeam)))
    sumgi = globalsum(sum(abs(1.-getgaminv(gather=0))))
    # --- Sum only the non overlapping part of the 3-d arrays. Otherwise
    # --- the answers will be different with differing numbers of processors.
    iz1 = top.izfsslave[me] - top.izslave[me]
    if me < npes-1:
      iz2 = top.izfsslave[me+1] - top.izslave[me]
    else:
      iz2 = iz1 + top.nzfsslave[me] + 1
    sumrho = globalsum(sum(sum(sum(w3d.rho[:,:,iz1:iz2]))))
    sumphi = globalsum(sum(sum(sum(w3d.phi[:,:,iz1+1:iz2+1]))))
  else:
    sumxp = sum(abs(getx()))
    sumyp = sum(abs(gety()))
    sumzp = sum(abs(getz()))
    sumvx = sum(abs(getvx()))
    sumvy = sum(abs(getvy()))
    sumvz = sum(abs(getvz()-top.vbeam))
    sumgi = sum(abs(1.-getgaminv()))
    sumrho = sum(sum(sum(w3d.rho)))
    sumphi = sum(sum(sum(w3d.phi)))
  return (sumxp,sumyp,sumzp,sumvx,sumvy,sumvz,sumgi,sumrho,sumphi)

##########################################################################
def setupvalidator():
  """Prints out text which is to be inserted into the validation deck. The
text contains the values with which later runs will be compared."""
  output = """
Insert the text below into the place in the validation deck where the check
is to be done.

# ---
from setupvalidation import comparetooriginal
comparetooriginal(
%22.15e,%22.15e,%22.15e,
%22.15e,%22.15e,%22.15e,
%22.15e,%22.15e,%22.15e
)
"""
  originaldata = getvalidationdata()
  print output%originaldata

##########################################################################
def comparetooriginal(sumxp,sumyp,sumzp,sumvx,sumvy,sumvz,sumgi,sumrho,sumphi):
  """Prints out the comparison between the original data and the current
  data."""
  fmt0 = "          %-22s %-22s %-22s"
  fmt1 = "          %22.15e %22.15e %22.15e"
  fmt2 = "original  %22.15e %22.15e %22.15e"
  otext = ["sum(x)","sum(y)","sum(z)",
           "sum(vx)","sum(vy)","sum(vz-vbeam)",
           "sum(1-gaminv)","sum(rho)","sum(phi)"]
  vdata = array(getvalidationdata())
  odata = array([sumxp,sumyp,sumzp,sumvx,sumvy,sumvz,sumgi,sumrho,sumphi])
  print fmt0%(otext[0],otext[1],otext[2])
  print fmt1%(vdata[0],vdata[1],vdata[2])
  print fmt2%(odata[0],odata[1],odata[2])
  print
  print fmt0%(otext[3],otext[4],otext[5])
  print fmt1%(vdata[3],vdata[4],vdata[5])
  print fmt2%(odata[3],odata[4],odata[5])
  print
  print fmt0%(otext[6],otext[7],otext[8])
  print fmt1%(vdata[6],vdata[7],vdata[8])
  print fmt2%(odata[6],odata[7],odata[8])
  print

# diffs = abs(odata - vdata)
# if max(diffs) > 0.:
#   print "Warning: The following values differ"
#   for i in range(len(vdata)):
#     if abs(odata[i] - vdata[i]) > 0.:
#       print "  ",otext[i]," by ",odata[i] - vdata[i]
# else:
#   print "No differences found"

