from warp import *
grid_1d_version = "$Id: grid_1d.py,v 1.4 2001/06/18 20:45:54 dave Exp $"
############################################################################
# This script contains two routines, one to gather particle data onto
# a 1D grid, and the other to scatter data from a 1D grid to the particles.
# DPG April 28, 1997
############################################################################

############################################################################
# This routine gathers a 1-D set of data onto a grid.
# Input:
#   location            1-D array holding the coordinates to be binned
#   data                1-D array holding the actual data (optional)
#   g1dmin              Minimum value of data range (optional)
#   g1dmax              Maximum value of data range (optional)
#   ng1d                Number of grid points minus one
# Output:
#   ng1d                The number of grid points, default is 40.
#   grid_1d(0:ng1d)     The 1-D array the holds the result.
#   grid_1dmesh(0:ng1d) A 1-D array of the mesh values.
#
# Notes on routines:
#   If the argument data is omitted, the grid data is divided by the total
#   number of data elements (to normalize so sum(grid_1d)=1.).  If data is
#   included, the each grid point is divided by the number of particles
#   contributing, so an average of the data is taken.
#
#   The [0, is used to prevent the array passed to sum to have zero length.
#   The way the routine is written, any particles out of the range
#   [g1dmin,g1dmax] are ignored.

def gather_1d(location,data=None,g1dmin=None,g1dmax=None,ng1d=40):
  """
Bins data onto a 1-D mesh. When data is not present, returns normalized density.
  - location: the coordinates of the points to be binned
  - data: the optional data to load onto the mesh
  - g1dmin: the minimum extent of the mesh, defaults to min(location)
  - g1dmax: the minimum extent of the mesh, defaults to max(location)
  - ng1d=40: number of grid points minus one (grid min is at index 0, and grid
             max is at index ng1d)
ouput
  returns a tuple (grid_1d,grid_1dmesh)
  - grid_1d holds the binned data
  - grid_1dmesh holds the bin locations
  """
  grid_1d = zeros(ng1d+1,'d')
  grid_1dmesh = zeros(ng1d+1,'d')

  if (g1dmin is None):
    g1dmin = min(location)
    g1dmin = g1dmin - (max(location)-g1dmin)/ng1d
  if (g1dmax is None):
    g1dmax = max(location)
    g1dmax = g1dmax + (g1dmax - g1dmin)/ng1d
  if (type(data)!=type(location)):
    if (data is None):
      data = 1.
      do_divide = false
  else:
    do_divide = true
  grid_1dmesh = g1dmin + (g1dmax - g1dmin)*iota(0,ng1d)/ng1d
  ig1d = ((location - g1dmin)/((g1dmax - g1dmin)/ng1d)).astype('i')
  wg1d =  (location - g1dmin)/((g1dmax - g1dmin)/ng1d) - ig1d
  if (do_divide):
    for i in xrange(0,ng1d+1):
      grid_1d[i] = (sum(compress(equal(ig1d,i  ),(1. - wg1d))) +
                    sum(compress(equal(ig1d,i-1),(     wg1d))))
    grid_1d[:] = where(greater(grid_1d,0.),grid_1d,ones(len(grid_1d)))
  else:
    grid_1d[:] = len(location)
  for i in xrange(0,ng1d+1):
    grid_1d[i] = (sum(compress(equal(ig1d,i  ),(1. - wg1d)*data)) +
                  sum(compress(equal(ig1d,i-1),(     wg1d)*data)))/grid_1d[i]
  return (grid_1d,grid_1dmesh)


#############################################################################
# This routine scatters 1D data onto particles.
# Input:
#   data                1-D array holding the actual data
#   zmsh                1-D array holding z locations of the data, assumed
#                       to be uniformly spaced
#   location            array holding the particle locations
# Output:
#   (returned)          array holding data at the particle locations
#                       same size as location array
#
# Note that for particles which lie outside of the range of data, a value
# of zero is returned.

def scatter_1d(data,zmsh,location):
  """
Returns the data on the zmsh interpolated to location
  - data input data
  - zmsh mesh locations where data is specified
  - location points to where data is interpolated
ouput
  - returns array same size as location with the interpolated data
  """
  # Perform error checking to make sure the input data makes sense
  assert (len(shape(data)) > 1), \
    ("gather_1d: error - data array is not a 1D array")
  assert (len(data) < 2), \
    ("gather_1d: error - data needs to have more than 1 element")
  assert (len(shape(zmsh)) > 1), \
    ("gather_1d: error - z mesh array is not a 1D array")
  assert (len(zmsh) < 2), \
    ("gather_1d: error - z mesh needs to have more than 1 element")
  assert (len(data) != len(zmsh)), \
    ("gather_1d: error - data array and z mesh array need to be the same size")

  # Get sizes of input arrays
  nz = len(zmsh)
  nd = len(data)

  # Get grid cell size
  dz = zmsh[1] - zmsh[0]

  # Get index and weight of particle data
  iz = ((location-zmsh[0])/dz).astype('i')
  wz =  (location-zmsh[0])/dz  - iz

  # Check for particles outside the range of data.  Set scale 'ww'
  # so that zero is returned for particles outside the range.
  ww = where(logical_and(less_equal(0,iz),less(iz,nd-1)), 1., 0.)
  iz = where(logical_and(less_equal(0,iz),less(iz,nd-1)), iz, 0)

  # Gather data and return it.
  return (ww*(take(data,iz)*(1.-wz) + take(data,iz+1)*wz))

