#
# Python file with some parallel operations
#
parallel_version = "$Id: parallel.py,v 1.13 2002/07/25 21:59:22 dave Exp $"

from Numeric import *
# --- Try import mpi - if not found, then run in serial mode
try:
  import mpi
  me = mpi.rank
  npes = mpi.procs
except ImportError:
  me = 0
  npes = 0

if npes > 0: lparallel = 1
else:        lparallel = 0

# --- The interface has changed some in the newest version of pyMPI.
try:
  mpiallreduce = mpi.allreduce
  def mpirecv(pe=0,ms=0):
    return mpi.recv(pe,ms)
except (AttributeError,NameError):
  pass

try:
  mpiallreduce = mpi.nativeAllreduce
  def mpirecv(pe=0,ms=0):
    result,stat = mpi.recv(pe,ms)
    return result
except (AttributeError,NameError):
  pass

# ---------------------------------------------------------------------------
# --- By default, set so that only PE0 sends output to the terminal.
if lparallel:
  mpi.synchronizeQueuedOutput('/dev/null')

def number_of_PE():
  if not lparallel: return 0
  return mpi.procs
def get_rank():
  if not lparallel: return 0
  return mpi.rank

if me == 0:
  try:
    from gist import *
  except ImportError:
    #from gistdummy import *
    pass
else:
  try:
    from gistdummy import *
  except ImportError:
    pass

# ---------------------------------------------------------------------------
# Enable output from all processors
def EnableAll():
  if not lparallel: return
  mpi.synchronizeQueuedOutput(None)

# ---------------------------------------------------------------------------
# Disable output from all but processor 0
def DisableAll():
  if not lparallel: return
  mpi.synchronizeQueuedOutput('/dev/null')

# ---------------------------------------------------------------------------
# Print object on all processors
def pprint(obj):
  if not lparallel:
    print str(obj)
    return
  mpi.synchronizedWrite(str(obj)+'\n')
  if mpi.rank == 0: print

# ---------------------------------------------------------------------------
# Print array (or list) from all processors
def aprint(obj):
  if not lparallel:
    print str(obj)
    return
  mpi.synchronizedWrite(str(obj))

# ---------------------------------------------------------------------------
# Get address of processor
def self_address():
  if not lparallel: return 0
  return mpi.rank

# ---------------------------------------------------------------------------
# Copy an array from processor i to processor 0
def getarray(src,v,dest=0):
  if not lparallel: return v
  if mpi.rank == src:
    mpi.send(v,dest)
  elif mpi.rank == dest:
    return mpirecv(src)
  return v

# ---------------------------------------------------------------------------
# Plot an array that is distributed on all processors
def plotarray(f,z,color="fg",linetype="solid",marks=0,marker="none",msize=1.0,
              width=1.0):
  ff = gatherarray(f)
  zz = gatherarray(z)
  if len(ff) > 0 and len(zz) == len(ff):
    plg(ff,zz,color=color,marker=marker,msize=msize,width=width,type=linetype)
# if mpi.rank == 0:
#   # --- Start list with PE0 data
#   x = [f]
#   y = [z]
#   # --- Append rest of PE's data to list
#   for i in range(1,mpi.procs):
#     x.append(mpirecv(i,2))
#     y.append(mpirecv(i,2))
#   # --- Get total number of elements in list
#   l = 0
#   for i in range(0,mpi.procs):
#     l = l + len(x[i])
#   # --- Create space for total list (flattened out)
#   xx = zeros(l,"d")
#   yy = zeros(l,"d")
#   # --- Copy data to new space
#   l = 0
#   for i in range(0,mpi.procs):
#     xx[l:l+len(x[i])] = x[i]
#     yy[l:l+len(y[i])] = y[i]
#     l = l + len(y[i])
#   # --- Now plot the result
#   plg(xx,yy,color=color,marker=marker,msize=msize,width=width,type=linetype)
# else:
#   mpi.send(f,0,2)
#   mpi.send(z,0,2)
    
# ---------------------------------------------------------------------------
# Plot particles that are distributed on all processors
def plotpart(f,z,color="fg",linetype="none",marker="\1",msize=1.0,**kw):
# ff = gatherarray(f)
# zz = gatherarray(z)
# if len(ff) > 0 and len(zz) == len(ff):
#   plg(ff,zz,color=color,type=linetype,marker=marker,msize=msize)
  kw['color'] = color
  kw['type'] = linetype
  kw['marker'] = marker
  kw['msize'] = msize
  if not lparallel:
    apply(plg,(f,z),kw)
    return
  # --- This way is preferred since for large data sets, it reduces the
  # --- risk of running out of memory since only part of the data is stored
  # --- on PE0 at a time.
  if mpi.rank == 0:
    if len(f) > 0 and len(f)==len(z):
      apply(plg,(f,z),kw)
    for i in range(1,mpi.procs):
      x = mpirecv(i,3)
      y = mpirecv(i,3)
      if len(x) > 0 and len(y)==len(x):
        apply(plg,(x,y),kw)
  else:
    mpi.send(f,0,3)
    mpi.send(z,0,3)

# ---------------------------------------------------------------------------
# Gather an object from all processors into a list
def gather(obj,dest=0):
  if not lparallel: return obj
  if mpi.rank == dest:
    result = []
    for i in range(mpi.procs):
      if i == dest:
        result.append(obj)
      else:
        result.append(mpirecv(i))
    return result
  else:
    mpi.send(obj,dest)
    return obj

# ---------------------------------------------------------------------------
# Define a barrier
def barrier():
  if not lparallel: return
  mpi.barrier()

# ---------------------------------------------------------------------------
# Broadcast an object to all procerssors
def broadcast(obj,root=0):
  if not lparallel: return obj
  return mpi.bcast(obj,root)

# ---------------------------------------------------------------------------
# Gather an object from all processors into a list and scatter back to all
def gatherall(obj):
  if not lparallel: return obj
  obj = gather(obj)
  return mpi.bcast(obj)
    
# ---------------------------------------------------------------------------
# General gatherarray which returns an array object combining the
# first dimension.
def gatherarray(a,root=0,othersempty=0):
  if not lparallel: return a
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
  # --- All processors but root simply return either the input argument
  # --- or an empty array
  if me != root:
    if othersempty: return zeros(len(shape(a))*[0],a.typecode())
    else: return a
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

# ---------------------------------------------------------------------------
# Generic global operation on a distributed array.
def globalop(a,localop,mpiop,defaultval):
  if type(a) in [FloatType,IntType]:
    local = a
  elif len(a) > 0:
    local = localop(a)
  else:
    local = defaultval
  if not lparallel: return local
  return mpiallreduce(local,eval("mpi."+mpiop))

# ---------------------------------------------------------------------------
# Specific operations on a distributed array.
def globalmax(a):
  return globalop(a,max,"MAX",-1.e36)
def globalmin(a):
  return globalop(a,min,"MIN",+1.e36)
def globalsum(a):
  return globalop(a,sum,"SUM",0.)
def globalave(a):
  s = globalop(a,sum,"SUM",0.)
  if type(a) in [FloatType,IntType]: a = [a]
  n = globalsum(len(a))
  if n > 0: return s/n
  else:     return 0.

# ---------------------------------------------------------------------------
# Generic parallel element-by-element operation on a distributed array.
# Note that this would be rather slow on large array. The mpi.allreduce
# routine needs to be fixed to accept arrays.
def parallelop(a,mpiop):
  if not lparallel: return a
  if type(a) == type(array([])):
    a1d = ravel(a) + 0
    for i in range(len(a1d)):
      a1d[i] = mpiallreduce(a1d[i],eval("mpi."+mpiop))
    a1d.shape = shape(a)
    return a1d
  else:
    return mpiallreduce(a,eval("mpi."+mpiop))

# ---------------------------------------------------------------------------
# Specific parallel element-by-element operations on a distributed array.
def parallelmax(a):
  return parallelop(a,"MAX")
def parallelmin(a):
  return parallelop(a,"MIN")
def parallelsum(a):
  return parallelop(a,"SUM")


