from warp import *

def vtkparticlegrid(xx,yy,zz,gnx=32,gny=32,gnz=32,
                    fname="phasespacegrid.vtk",returngrid=0,
                    comment="phase space grid",cmax=None):
  gg = zeros((gnx+1,gny+1,gnz+1),'d')
  xmin = min(xx)
  xmax = max(xx)
  ymin = min(yy)
  ymax = max(yy)
  zmin = min(zz)
  zmax = max(zz)
  gdx = (xmax - xmin)/(gnx-1)
  gdy = (ymax - ymin)/(gny-1)
  gdz = (zmax - zmin)/(gnz-1)
  xmin = xmin - gdx
  xmax = xmax + gdx
  ymin = ymin - gdy
  ymax = ymax + gdy
  zmin = zmin - gdz
  zmax = zmax + gdz
  setgrid3d(len(xx),xx,yy,zz,gnx,gny,gnz,gg,xmin,xmax,ymin,ymax,zmin,zmax)
  if cmax is not None: gg = where(gg>cmax,cmax,gg)

  ff = open(fname,"w")
  stdoutsave = sys.stdout
  sys.stdout = ff
  print "# vtk DataFile Version 2.0"
  print comment
  print "ASCII"
  print "DATASET RECTILINEAR_GRID"
  print "DIMENSIONS %d %d %d"%(gnx+1,gny+1,gnz+1)
  print "X_COORDINATES %d float"%(gnx+1)
  ss = ""
  for i in range(gnx+1): ss = ss + "%f "%(xmin + i*gdx)
  print ss
  print "Y_COORDINATES %d float"%(gny+1)
  ss = ""
  for i in range(gny+1): ss = ss + "%f "%(ymin + i*gdy)
  print ss
  print "Z_COORDINATES %d float"%(gnz+1)
  ss = ""
  for i in range(gnz+1): ss = ss + "%f "%(zmin + i*gdz)
  print ss
  print "POINT_DATA %d"%((gnx+1)*(gny+1)*(gnz+1))
  print "SCALARS scalars float 1"
  print "LOOKUP_TABLE default"
  for iz in range(gnz+1):
    for iy in range(gny+1):
      for ix in range(gnx+1):
        print "%f"%gg[ix,iy,iz]
  ff.close()
  sys.stdout = stdoutsave

  if returngrid: return gg

