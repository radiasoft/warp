from warp import *
import __main__

def printparameters3d():

  # --- Exit now if output parameters are not to be printed
  if not top.lprntpara: return

  # --- Formats
  f20 = " %s%11.4e%s\n"
  f30 = " %s%8d%s\n"


  textblock = \
      f30%("Number of grid points in x = ",w3d.nx," ") + \
      f30%("Number of grid points in y = ",w3d.ny," ") + \
      f30%("Number of grid points in z = ",w3d.nzfull," ") + \
      f20%("Grid spacing in x = ",w3d.dx," m") + \
      f20%("Grid spacing in y = ",w3d.dy," m") + \
      f20%("Grid spacing in z = ",w3d.dz," m")
  if top.nbend >= 1:
    xrbend = top.bendrc[0]
    xbendlen = top.bendze[0] - top.bendzs[0]
    xstralen = top.bendzs[1] - top.bendze[0]
    xz0bend = top.bendzs[1]
    textblock = textblock + \
      f20%("Bend radius = ",xrbend," m") + \
      f20%("Bend length = ",xbendlen," m") + \
      f20%("Straight section length = ",xstralen," m") + \
      f20%("Z at start of first bend = ",xz0bend," m")
    if top.ndipo == 0:
      xbbend = top.dipoby[0]
      textblock = textblock + f20%("Bending field = ",xbbend," T") + \
  plt(textblock,0.12,0.88,justify="LT")
  fma()

# --- Explicitly add the function to main since it will be called by fortran
__main__.__dict__['printparameters3d'] = printparameters3d

