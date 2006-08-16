from warp import *
import __main__

def printparameters3d():

  # --- Exit now if output parameters are not to be printed
  if not top.lprntpara: return

  # --- Formats
  f20 = ' %s%11.4e%s\n' # float format
  f30 = ' %s%8d%s\n'    # integer format
  f40 = ' %s%s%s\n'     # string format
  f50 = ' %s%11.4e%s%11.4e%s\n' # format for mesh extent

  textblock = (
    f40%('Particle distribution = ',arraytostr(w3d.distrbtn),'') +
    f30%('Number of grid points in x = ',w3d.nx,'') +
    f30%('Number of grid points in y = ',w3d.ny,'') +
    f30%('Number of grid points in z = ',w3d.nzfull,'') +
    f20%('Grid spacing in x = ',w3d.dx,' m') +
    f20%('Grid spacing in y = ',w3d.dy,' m') +
    f20%('Grid spacing in z = ',w3d.dz,' m') +
    f50%('Grid extends in x from ',w3d.xmmin,' to ',w3d.xmmax,' m') +
    f50%('Grid extends in y from ',w3d.ymmin,' to ',w3d.ymmax,' m') +
    f50%('Grid extends in z from ',w3d.zmminglobal,' to ',w3d.zmmaxglobal,' m')
    )

  if w3d.l2symtry: textblock += ' Two fold symmetry\n'
  if w3d.l4symtry: textblock += ' Four fold symmetry\n'

  if w3d.solvergeom == w3d.XYZgeom: textblock += ' Geometry is 3-D\n'
  if w3d.solvergeom == w3d.RZgeom: textblock += ' Geometry is RZ\n'
  if w3d.solvergeom == w3d.AMRgeom: textblock += ' Geometry is RZ with AMR\n'
  if w3d.solvergeom == w3d.XZgeom: textblock += ' Geometry is XZ\n'
  if w3d.solvergeom == w3d.XYgeom: textblock += ' Geometry is XY\n'
  if w3d.solvergeom == w3d.Zgeom: textblock += ' Geometry is Z\n'
  if w3d.solvergeom == w3d.Rgeom: textblock += ' Geometry is R\n'

  if top.nbend >= 1:
    xrbend = top.bendrc[0]
    xbendlen = top.bendze[0] - top.bendzs[0]
    xstralen = top.bendzs[1] - top.bendze[0]
    xz0bend = top.bendzs[1]
    textblock = textblock + (
      f20%('First Bend radius = ',xrbend,' m') +
      f20%('First Bend length = ',xbendlen,' m') +
      f20%('Next Straight section length = ',xstralen,' m') +
      f20%('Z at start of first bend = ',xz0bend,' m')
      )
    if top.ndipo == 0:
      xbbend = top.dipoby[0]
      textblock = textblock + f20%('First dipole field = ',xbbend,' T')

  plt(textblock,0.12,0.88,justify='LT')
  fma()

# --- Explicitly add the function to main since it will be called by fortran
__main__.__dict__['printparameters3d'] = printparameters3d

