from warp import *
import __main__

def printparametersrz():

  # --- Exit now if output parameters are not to be printed
  if not top.lprntpara: return

  # --- Formats
  f20 = " %s%11.4e%s\n"
  f40 = " %s%11.4e%s%11.4e%s\n"
  f30 = " %s%8d%s\n"

  textblock = \
      f30%("Number of grid points in r = ",nr," ") + \
      f30%("Number of grid points in z = ",nz," ")
  plt(textblock,0.12,0.88,justify="LT")
  fma()

# --- Explicitly add the function to main since it will be called by fortran
__main__.__dict__['printparametersrz'] = printparametersrz

