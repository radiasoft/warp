"""Utility and convenience routines used in Warp

iota(): returns a array of sequential integers
remark(): same as print
sign(): emulation of sign function
arrayappend(): appends to multi-dimensional array
exists(): checks if a variable exists
ave(): averages an array of numbers
maxnd(): finds max of multi-dimensional array
minnd(): finds min of multi-dimensional array
getnextfilename(): finds next available file name in a numeric sequence
getmesh2d(): returns tuple of 2 2-d arrays holding the coordinates of the
             two dimensions
getmesh3d(): returns tuple of 3 3-d arrays holding the coordinates of the
             three dimensions
averagezdata(): Does local averaging over the first dimension of the input
                array. Can also do down select on the data.

"""
from warp import *

warputils_version = "$Id: warputils.py,v 1.3 2004/04/16 17:00:43 dave Exp $"

def warputilsdoc():
  import warputils
  print warputils.__doc__

# --- Convenience function modeled after the iota of basis
def iota(low,high=None,step=1):
  if high is None:
    if step > 0:
      return arange(1,low+1,step)
    else:
      return arange(low,0,step)
  else:
    if step > 0:
      return arange(low,high+1,step)
    else:
      return arange(low,high-1,step)

# --- Convenience function to do printing
def remark(s):
  print s

# --- Replicate the sign function with two arguments. If only one is
# --- given, return the value from the Numeric sign function.
# --- This should realy be removed.
numericsign = sign
def sign(x,y=None):
  if y is None: return numericsign(x)
  if type(x) == ArrayType:
    result = where(greater(y,0.),abs(x),-abs(x))
    result = where(equal(y,0.),0.,result)
    return result
  else:
    if y > 0:
      return abs(x)
    elif y < 0:
      return -abs(x)
    else:
      return 0

# --- Convenience function which returns meshes filled with the coordinates.
def getmeshcoordinates(mins,dds,nns):
  """
getmeshcoordinates(mins,dds,nns)
Returns arrays holding the coordinates of the mesh points.
Lenght of list of inputs determines number of dimensions.
  """
  nns = tuple(array(nns) + 1)
  cc = indices(nns,'d')
  for i in xrange(len(mins)): cc[i] = mins[i] + cc[i]*dds[i]
  clist = []
  for i in xrange(len(mins)): clist.append(cc[i])
  return tuple(clist)

def getmesh2d(xmin,dx,nx,ymin,dy,ny):
  """
getmesh2d(xmin,dx,nx,ymin,dy,ny)
Returns 2 2-d arrays holding the coordinates of the mesh points
  """
  return getmeshcoordinates([xmin,ymin],[dx,dy],[nx,ny])

def getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dz,nz):
  """
getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dx,nz)
Returns 3 3-d arrays holding the coordinates of the mesh points
  """
  return getmeshcoordinates([xmin,ymin,zmin],[dx,dy,dz],[nx,ny,nz])

# --- This function appends a new element to the end of an array.
# --- It is not very efficient since it creates a whole new array each time.
def arrayappend(x,a):
  xshape = list(shape(x))
  if type(a) == ArrayType:
    pass
  elif type(a) == ListType:
    a = array(a)
  else:
    a = array([a])
  ashape = list(shape(a))
  if len(xshape)==1 and len(ashape)==1:
    xshape[0] = xshape[0] + ashape[0]
    y = fzeros(xshape,x.typecode())
    y[0:xshape[0]-ashape[0]] = x
    y[xshape[0]-ashape[0]:] = a
  elif len(xshape) == len(ashape)+1 and xshape[:-1] == ashape:
    xshape[-1] = xshape[-1] + 1
    y = fzeros(xshape,x.typecode())
    y[...,0:-1] = x
    y[...,-1] = a
  return y

# Convenience function which returns true if variable exists
def exists(x):
  """
Checks whether or not the variable whose name is specified exists in the
main dictionary.
 - x: Name of variable - must be a string.
  """
  if x in __main__.__dict__.keys(): return true
  if x in locals().keys(): return true
  if x in globals().keys(): return true
  return false

# Returns the average of the input array
def ave(x,index=0):
  if shape(x)[index] > 0:
    return sum(x,index)/shape(x)[index]
  else:
    return 0.

def averagezdata(qty,navg=0,nlines=100,n1=None,n2=None,istep=None,
                 includezeros=false):
  """
Averages data over local region. It also can down select data in the other
dimension.
  - qty: Data to be smoothed. Can be either a 1-D or 2-D array.
  - navg=0: number of data points to average over
  - nlines=100: number of lines from second dimension to choose.
  - n1=shape(qty)[0]-1:
  - n2=shape(qty)[1]-1:
  - istep=max(1,n2/nlines):
  - includezeros=false: by default, only non-zero data is averaged over.
  """
  if navg == 0 or nlines == 0: return qty
  if len(shape(qty)) == 1:
    fixqty = 1
    qty.shape = (len(qty),1)
  else:
    fixqty = 0
  if not n1: n1 = shape(qty)[0] - 1
  if not n2: n2 = shape(qty)[1] - 1
  if istep is None: istep = max(1,n2/nlines)
  hl = qty[:,::istep] + 0.
  hl[navg,:] = sum(qty[navg-navg:navg+navg+1,::istep])
  nn = 2*navg+1 + zeros(shape(hl))
  if not includezeros:
    nn[navg,:] = sum(where(qty[navg-navg:navg+navg+1,::istep]==0.,0,1),0)
  for j in range(navg+1,n1-navg-1):
    hl[j,:] = hl[j-1,:] + (qty[j+navg,::istep] - qty[j-navg-1,::istep])
    nn[j,:] = nn[j-1,:] + (+ where(qty[j+navg,::istep]==0,0,1)
                           - where(qty[j-navg-1,::istep]==0,0,1))
  nn = where(nn==0,1,nn)
  hl = where(qty[:,::istep]==0.,0.,hl)
  hl[navg+1:n1-navg-1,:] = hl[navg+1:n1-navg-1,:]/nn[navg+1:n1-navg-1,:]
  if fixqty: qty.shape = (shape(qty)[0],)
  if n2 > 1: return hl
  else: return hl[:,0]

# --- Returns the max of the multiarray
def maxnd(x):
  """Return the max element of an array of any dimension"""
  xtemp = reshape(x,tuple([product(array(x.shape))]))
  return max(xtemp)
# --- Returns the min of the multiarray
def minnd(x):
  """Return the min element of an array of any dimension"""
  xtemp = reshape(x,tuple([product(array(x.shape))]))
  return min(xtemp)

# Gets next available filename with the format 'root.nnn.suffix'.
def getnextfilename(root,suffix):
  dir = string.join(os.listdir('.'))
  i = 0
  name = root+('.%03d.'%i)+suffix
  while re.search(name,dir):
    i = i + 1
    name = root+('.%03d.'%i)+suffix
  return name

