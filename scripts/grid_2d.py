from warp import *
grid_2d_version = "$Id: grid_2d.py,v 1.2 2001/01/25 18:05:37 dave Exp $"
# This routine lays a 2-D set of data onto a grid.
# Input:
#   gxmin       Minimum x value of data range.
#   gxmax       Maximum x value of data range.
#   gymin       Minimum y value of data range.
#   gymax       Maximum y value of data range.
#   locx        1-D array holding the x location
#   locy        1-D array holding the y location
#   data        1-D array holding the actual data (optional)
# Output:
#   grid_2d(0:ngx,0:ngy)  The 2-D array the holds the result.
#
# The grid dimension must be specified before this script is read in other-
# wise they defaul to 40.
#   ngx                   The number of grid points in x.
#   ngy                   The number of grid points in y.

# Notes on routines:
#   If the argument data is omitted, the grid data is divided by the total
#   number of data elements (to normalize so sum(grid_2d)=1.).  If data is
#   included, the each grid point is divided by the number of particles
#   contributing, so an average of the data is taken.
#
#   The [0 in the sums is used to prevent the array passed to sum to have
#   zero length.
#
#   The way the routine is written, any particles out of range are ignored.
#   The grid does not need to cover all of the data passed.

def setgrid_2d(locx,locy,data=None,ngx=40,ngy=40,
               gxmin=None,gxmax=None,gymin=None,gymax=None):
  """
Bins data onto a 2-D grid. When data is not present, returns normalized density.
  - locx x location of data points
  - locy y location of data points
  - data optional data at locations
  - ngx=40 number of x grid points
  - ngy=40 number of y grid points
  - gxmin minimum x extent of the grid, defaults to min(locx)
  - gxmax maximum x extent of the grid, defaults to max(locx)
  - gxmin minimum y extent of the grid, defaults to min(locy)
  - gxmax maximum y extent of the grid, defaults to max(locy)
ouput
  returns a tuple (grid_2d,grid_2dx,grid_2dy)
  - grid_2d 2-D mesh
  - grid_2dx x values of mesh
  - grid_2dy y values of mesh
  """
  if not gxmin: gxmin = min(locx)
  if not gxmax: gxmax = max(locx)
  if not gymin: gymin = min(locy)
  if not gymax: gymax = max(locy)
  grid_2dx = (gxmin + arange(0,ngx+1)*(gxmax-gxmin)/ngx)*ones(ngy+1)[:,NewAxis]
  grid_2dy = (gymin + arange(0,ngy+1)*(gymax-gymin)/ngy)[:,NewAxis]*ones(ngx+1)
  grid_2d = zeros((ngx+1,ngy+1),'d')
  if not data:
    data = 1.
    do_divide = 0
  else:
    do_divide = 1
  igx = ((locx - gxmin)/((gxmax - gxmin)/ngx)).astype('i')
  wgx =  (locx - gxmin)/((gxmax - gxmin)/ngx) - igx
  igy = ((locy - gymin)/((gymax - gymin)/ngy)).astype('i')
  wgy =  (locy - gymin)/((gymax - gymin)/ngy) - igy
  if do_divide:
    for i in xrange(ngx+1):
      for j in xrange(ngy+1):
        grid_2d[i,j]=(sum(compress(logical_and(equal(igx,i  ),equal(igy,j  )),
                          (1.-wgx)*(1.-wgy)))+
                      sum(compress(logical_and(equal(igx=i-1),equal(igy=j  )),
                          (   wgx)*(1.-wgy)))+
                      sum(compress(logical_and(equal(igx=i  ),equal(igy=j-1)),
                          (1.-wgx)*(   wgy)))+
                      sum(compress(logical_and(equal(igx=i-1),equal(igy=j-1)),
                          (   wgx)*(   wgy))))
    grid_2d = where(greater(grid_2d,0.),grid_2d,1.)
  else:
#   grid_2d = shape(locx)
    grid_2d[:,:] = 1.
  for i in xrange(ngx+1):
    for j in xrange(ngy+1):
       grid_2d[i,j] = (
           (sum(compress(logical_and(equal(igx,i  ),equal(igy,j  )),
                (1.-wgx)*(1.-wgy)*data))+
            sum(compress(logical_and(equal(igx,i-1),equal(igy,j  )),
                (   wgx)*(1.-wgy)*data))+
            sum(compress(logical_and(equal(igx,i  ),equal(igy,j-1)),
                (1.-wgx)*(   wgy)*data))+
            sum(compress(logical_and(equal(igx,i-1),equal(igy,j-1)),
                (   wgx)*(   wgy)*data)))/grid_2d[i,j])
  return (grid_2d,grid_2dx,grid_2dy)







#############################################################################
# This routine scatters 2D data onto particles.
# Input:
#   data                2-D array holding the actual data
#   xmin, ymin          lower range of mesh
#   xmax, ymax          upper range of mesh
#   location            array holding the particle locations
# Output:
#   (returned)          array holding data at the particle locations
#                       same size as location array
#
# Note that for particles which lie outside of the range of data, a value
# of zero is returned.

def scatter_2d(data,xmin,ymin,xmax,ymax,xx,yy):
  """
Returns the data on the zmsh interpolated to location
  - data input data
  - xmin, ymin lower range of mesh
  - xmax, ymax upper range of mesh
  - xx, yy coordinates points to where data is interpolated
ouput
  - returns array same size as coordinates with the interpolated data
  """
  # Perform error checking to make sure the input data makes sense
  if (len(shape(data)) != 2):
    remark("scatter_2d: error - data array is not a 2D array")
    return 0
  if (shape(data)[0] < 2 or shape(data)[1] < 2) :
    remark("scatter_2d: error - data needs to have more than 1 element")
    return 0

  # Get sizes of input arrays
  nx = shape(data)[0]
  ny = shape(data)[1]

  # Get grid cell size
  dx = (xmax - xmin)/(nx-1.)
  dy = (ymax - ymin)/(ny-1.)

  # Get index and weight of particle data
  ix = ((xx-xmin)/dx).astype('i')
  wx =  (xx-xmin)/dx  - ix
  iy = ((yy-ymin)/dy).astype('i')
  wy =  (yy-ymin)/dy  - iy

  # Check for particles outside the range of data.  Set weight to zero
  # so that zero is returned for particles outside the range.
  ww = where(logical_and(logical_and(less_equal(0,ix),less(ix,nx-1)),
                         logical_and(less_equal(0,iy),less(iy,ny-1))), 1., 0.)
  ix = where(logical_and(less_equal(0,ix),less(ix,nx-1)), ix, 0)
  iy = where(logical_and(less_equal(0,iy),less(iy,ny-1)), iy, 0)

  # Convert to 1-D array to allow use of 'take'
  ii = ix*ny + iy
  data1 = reshape(data,tuple([shape(data)[0]*shape(data)[1]]))

  # Gather data and return it.
  return (ww*(take(data1,ii     )*(1.-wx)*(1.-wy) +
              take(data1,ii+1   )*(   wx)*(1.-wy) +
              take(data1,ii+ny  )*(1.-wx)*(   wy) +
              take(data1,ii+ny+1)*(   wx)*(   wy)))


