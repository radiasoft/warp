#
# Python file with some parallel operations
#

from Numeric import *
import mpi
parallel_version = "$Id:"
me = mpi.rank
npes = mpi.procs

# --- By default, set so that only PE0 sends output to the terminal.
mpi.synchronizeQueuedOutput('/dev/null')

def number_of_PE():
  return mpi.procs
def get_rank():
  return mpi.rank

if me == 0:
  from gist import *
else:
  from gistdummy import *

# Enable output from all processors
def EnableAll():
  mpi.synchronizeQueuedOutput(None)

# Disable output from all but processor 0
def DisableAll():
  mpi.synchronizeQueuedOutput('/dev/null')

# Print object on all processors
def pprint(obj):
  mpi.synchronizedWrite(str(obj)+'\n')
  if mpi.rank == 0: print

# Print array (or list) from all processors
def aprint(obj):
  mpi.synchronizedWrite(str(obj))

# Get address of processor
def self_address():
  return mpi.rank

# Copy an array from processor i to processor 0
def getarray(src,v,dest=0):
  if mpi.rank == src:
    mpi.send(v,dest)
  elif mpi.rank == dest:
    return mpi.recv(src)
  return v

# Plot an array that is distributed on all processors
def plotarray(f,z,color="fg",type="solid",marks=0,marker="none",msize=1.0,
              width=1.0):
  if mpi.rank == 0:
    # --- Start list with PE0 data
    x = [f]
    y = [z]
    # --- Append rest of PE's data to list
    for i in range(1,mpi.procs):
      x.append(mpi.recv(i,2))
      y.append(mpi.recv(i,2))
    # --- Get total number of elements in list
    l = 0
    for i in range(0,mpi.procs):
      l = l + len(x[i])
    # --- Create space for total list (flattened out)
    xx = zeros(l,"d")
    yy = zeros(l,"d")
    # --- Copy data to new space
    l = 0
    for i in range(0,mpi.procs):
      xx[l:l+len(x[i])] = x[i]
      yy[l:l+len(y[i])] = y[i]
      l = l + len(y[i])
    # --- Now plot the result
    plg(xx,yy,color=color,marker=marker,msize=msize,width=width)
  else:
    mpi.send(f,0,2)
    mpi.send(z,0,2)
    
# Plot particles that are distributed on all processors
def plotpart(f,z,color="fg",type="none",marker="\1",msize=1.0):
  ff = gatherarray(f)
  zz = gatherarray(z)
  plg(ff,zz,color=color,type="none",marker=marker,msize=msize)
# if mpi.rank == 0:
#   if len(f) > 0 and len(f)==len(z):
#     plg(f,z,color=color,type="none",marker=marker,msize=msize)
#   for i in range(1,mpi.procs):
#     x = mpi.recv(i,3)
#     y = mpi.recv(i,3)
#     if len(x) > 0 and len(y)==len(x):
#       plg(x,y,color=color,type=type,marker=marker,msize=msize)
# else:
#   mpi.send(f,0,3)
#   mpi.send(z,0,3)

# Gather an object from all processors into a list
def gather(obj,dest=0):
  if mpi.rank == dest:
    result = []
    for i in range(mpi.procs):
      if i == dest:
        result.append(obj)
      else:
        result.append(mpi.recv(i))
    return result
  else:
    mpi.send(obj,dest)
    return obj
    
# General gatherarray which returns an array object combining the
# first dimension.
def gatherarray(a,root=0):
  # --- First check if input can be converted to an array
  isinputok = 1
  try:
    if type(a) in [type(0.),type(0)]:
      a = array([a])
    else:
      a = array(a)
  except:
    isinputok = 0
  # --- Make sure the input is ok on all of the processors
  isinputok = globalmin(isinputok)
  # --- If any returned an error, then all exit (to avoid a deadlock)
  if not isinputok:
    print "Object could not be converted to an array"
    return None
  # --- Now, actually gather the array.
  result = gather(a,root)
  # --- All processors but root simply return the input argument
  if me != root: return a
  # --- Processor 1 reshapes the data, removing the first dimension
  # --- Do it bit by bit since the data passed by the other processors may
  # --- not be all the same size.
  newlen = 0
  for i in range(npes):
    newlen = newlen + shape(result[i])[0]
  newshape = list(shape(result[0]))
  newshape[0] = newlen
  newresult = zeros(newshape,a.typecode())
  i1 = 0
  for i in range(npes):
    i2 = i1 + shape(result[i])[0]
    newresult[i1:i2,...] = result[i]
    i1 = i2
  return newresult
  ## --- Old way
  ## --- Its easy if all of the arrays passed from the other processors
  ## --- are the same size.
  #result = array(result)
  #ss = list(shape(result))
  #snew = ss[1:]
  #snew[0] = ss[0]*ss[1]
  #result.shape = snew
  ## --- Return the result
  #return result

# Generic global operation on a distributed array.
def globalop(a,localop,mpiop,defaultval):
  if npes == 0: return localop(a)
  if type(a) in [type(0.),type(0)]:
    local = a
  elif len(a) > 0:
    local = localop(a)
  else:
    local = defaultval
  return mpi.allreduce(local,mpiop)

# Specific operations on a distributed array.
def globalmax(a):
  return globalop(a,max,mpi.MAX,-1.e36)
def globalmin(a):
  return globalop(a,min,mpi.MIN,+1.e36)
def globalsum(a):
  return globalop(a,sum,mpi.SUM,0.)
def globalave(a):
  s = globalop(a,sum,mpi.SUM,0.)
  if type(a) in [type(0.),type(0)]: a = [a]
  n = globalsum(len(a))
  return s/n

# Generic parallel element-by-element operation on a distributed array.
# Note that this would be rather slow on large array. The mpi.allreduce
# routine needs to be fixed to accept arrays.
def parallelop(a,mpiop):
  if npes == 0: return a
  if type(a) == type(array([])):
    a1d = ravel(a) + 0
    for i in range(len(a1d)):
      a1d[i] = mpi.allreduce(a1d[i],mpiop)
    a1d.shape = shape(a)
    return a1d
  else:
    return mpi.allreduce(a,mpiop)

# Specific parallel element-by-element operations on a distributed array.
def parallelmax(a):
  return parallelop(a,mpi.MAX)
def parallelmin(a):
  return parallelop(a,mpi.MIN)
def parallelsum(a):
  return parallelop(a,mpi.SUM)


