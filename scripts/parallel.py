#
# Python file with some parallel operations
#
parallel_version = "$Id: parallel.py,v 1.41 2010/07/01 21:56:43 dave Exp $"

from numpy import *
# --- Try import mpi - if not found, then run in serial mode
# --- Note that:
# --- mpi.COMM_WORLD is same as MPI_COMM_WORLD
# --- mpi.WORLD is a duplicate of MPI_COMM_WORLD
# --- comm_world is used for most communications (and defaults to mpi.WORLD)
try:
  import mpi
  mpi.synchronizeQueuedOutput(None)
  me = mpi.rank
  npes = mpi.procs
  comm_world = mpi.WORLD
except ImportError:
  me = 0
  npes = 1
  comm_world = None

lparallel = (npes > 1)

def setdefaultcomm_world(comm):
  global comm_world,me,npes
  if not lparallel: return
  comm_world = comm
  me = comm_world.rank
  npes = comm_world.procs

"""
# --- This check is old enough that it is no longer needed.

# --- The interface has changed some in the newest version of pyMPI.
# --- Check the interface to the mpi.recv command. The newer versions
# --- return a tuple instead of just the data itself.
# --- Is there a better way of doing this?
if lparallel:
  mpi.send(me,me)
  _i = mpi.recv(me)
  if type(_i) == TupleType: _newpympi = 1
  else:                     _newpympi = 0
else:
  _newpympi = 1
"""
_newpympi = 1

if _newpympi:
  def mpirecv(pe=0,ms=0,comm=None):
    if comm is None: comm = comm_world
    result,stat = comm.recv(pe,ms)
    return result
else:
  def mpirecv(pe=0,ms=0,comm=None):
    if comm is None: comm = comm_world
    return comm.recv(pe,ms)

# ---------------------------------------------------------------------------
# --- By default, set so that only PE0 sends output to the terminal.
if lparallel:
  mpi.synchronizeQueuedOutput('/dev/null')

def number_of_PE(comm=None):
  if not lparallel: return 1
  if comm is None: comm = comm_world
  return comm.procs
def get_rank(comm=None):
  if not lparallel: return 0
  if comm is None: comm = comm_world
  return comm.rank

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
  # Ignore all exceptions to make sure that there is not a lock up.
  try:
    ss = str(obj)
  except:
    ss = ''
  mpi.synchronizedWrite(ss+'\n')
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
def self_address(comm=None):
  if not lparallel: return 0
  if comm is None: comm = comm_world
  return comm.rank

# ---------------------------------------------------------------------------
# Copy an array from processor i to processor 0
def getarray(src,v,dest=0,comm=None):
  if not lparallel: return v
  if comm is None: comm = comm_world
  if comm.rank == src:
    comm.send(v,dest)
  elif comm.rank == dest:
    return mpirecv(src,comm=comm)
  return v

# ---------------------------------------------------------------------------
# Gather an object from all processors in a communicator (default comm_world) into a list
def gather(obj,dest=0,comm=None):
  if not lparallel: return [obj]
  if comm is None: comm = comm_world
  if comm.rank == dest:
    result = []
    for i in range(comm.procs):
      if i == dest:
        result.append(obj)
      else:
        result.append(mpirecv(i,comm=comm))
    return result
  else:
    comm.send(obj,dest)
    return [obj]

# ---------------------------------------------------------------------------
# Gather an object from a list of processors into a list on destination processor,
# eventually broadcasting result. 
def gatherlist(obj,dest=0,procs=None,bcast=0,comm=None):
  if not lparallel: return [obj]
  if comm is None: comm = comm_world
  if procs is None:
    procs=range(comm.procs)
  else:
    procs=list(procs)
  result = []
  if comm.rank == dest:
    for i in procs:
      if i == dest:
        result.append(obj)
      else:
        result.append(mpirecv(i,comm=comm))
  else:
    if comm.rank in procs:
      comm.send(obj,dest)

  if bcast:
    result = comm.bcast(result,dest)
      
  return result

# ---------------------------------------------------------------------------
# Gather an object from a list of processors into a list on destination processor,
# eventually broadcasting result. The send is decomposed into log2(n) sends 
# between processors.
def gatherlog(obj,dest=0,procs=None,bcast=0,comm=None):
  if not lparallel: return [obj]
  if comm is None: comm = comm_world
  if procs is None:
    procs=range(comm.procs)
  else:
    procs=list(procs)
  n = int(ceil(log2(len(procs))))
  obj = [obj]
  for i in range(n):
    st = 2**i
    stp = 2**(i+1)
    listsend = procs[st::stp]
    listrecv = procs[::stp][:len(listsend)]
    if comm.rank in listsend:
      ip = listsend.index(comm.rank)
      comm.send(obj,listrecv[ip])
    if comm.rank in listrecv:
      ip = listrecv.index(comm.rank)
      obj+=mpirecv(listsend[ip],comm=comm)

  if bcast:
    obj = comm.bcast(obj,procs[0])
  else:
    if dest<>procs[0]:
      if comm.rank==procs[0]:
        comm.send(obj,dest)
      if comm.rank==dest:
        obj=mpirecv(procs[0],comm=comm)
      
  return obj

# ---------------------------------------------------------------------------
# Define a barrier
def barrier(comm=None):
  if not lparallel: return
  if comm is None: comm = comm_world
  comm.barrier()

# ---------------------------------------------------------------------------
# Broadcast an object to all processors
def broadcast(obj,root=0,comm=None):
  if not lparallel: return obj
  if comm is None: comm = comm_world
  return comm.bcast(obj,root)

# ---------------------------------------------------------------------------
# Gather an object from all processors into a list and scatter back to all
def gatherall(obj,comm=None):
  if not lparallel: return obj
  if comm is None: comm = comm_world
  obj = gather(obj,comm=comm)
  return comm.bcast(obj)
    
# ---------------------------------------------------------------------------
# General gatherarray which returns an array object combining the
# first dimension.
def gatherarray(a,root=0,othersempty=0,bcast=0,comm=None):
  if not lparallel: return a
  if comm is None: comm = comm_world
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
  isinputok = globalmin(isinputok,comm=comm)
  # --- If any returned an error, then all exit (to avoid a deadlock)
  if not isinputok:
    print "Object could not be converted to an array"
    return None
  # --- Now, actually gather the array.
  # --- The check of whether the result is ok may not be needed.
  try:
    result = gather(a,root,comm=comm)
    isinputok = 1
  except:
    isinputok = 0
  # --- Make sure again that the input is ok on all of the processors
  isinputok = globalmin(isinputok,comm=comm)
  if not isinputok:
    print "Error in gather object"
    try:
      "Object has shape ",shape(a)
    except NameError:
      pass
    return None
  # --- All processors but root simply return either the input argument
  # --- or an empty array unless the result is to be broadcast
  if comm.rank != root and not bcast:
    if othersempty: return zeros(len(shape(a))*[0],a.dtype.char)
    else: return a
  # --- Root processor reshapes the data, removing the first dimension
  # --- Do it bit by bit since the data passed by the other processors may
  # --- not be all the same size.
  if comm.rank == root:
    newlen = 0
    for i in range(comm.size):
      newlen = newlen + shape(result[i])[0]
    newshape = list(shape(result[0]))
    newshape[0] = newlen
    newresult = zeros(newshape,a.dtype.char)
    i1 = 0
    for i in range(comm.size):
      i2 = i1 + shape(result[i])[0]
      newresult[i1:i2,...] = result[i]
      i1 = i2
  else:
    newresult = 0
  if bcast: newresult = comm.bcast(newresult,root)
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
# Find the nonzero value of array over all processors. This assumes that the
# non-zero values for each index are the same for all processors.
# Resulting data is broadcast to all processors.
def parallelnonzeroarray(a,comm=None):
  if comm is None: comm = comm_world
  dmax = parallelmax(a,comm=comm)
  dmin = parallelmin(a,comm=comm)
  result = where(not_equal(dmax,0),dmax,dmin)
  return result

# ---------------------------------------------------------------------------
# Generic global operation on a distributed array.
def globalop(a,localop,mpiop,defaultval,comm=None):
  if comm is None: comm = comm_world
  if len(shape(a)) == 0:
    local = a
  elif len(a) > 0:
    try:
      # --- Protect application of localop since the data may not be
      # --- appropriate on all processors, for example it may be an empty array.
      local = localop(a)
    except:
      local = defaultval
  else:
    local = defaultval
  if not lparallel: return local
  return comm.allreduce(local,getattr(mpi,mpiop))

# ---------------------------------------------------------------------------
# Specific operations on a distributed array.
def globalmax(a,comm=None):
  if comm is None: comm = comm_world
  return globalop(a,max,"MAX",-1.e36,comm=comm)
def globalmin(a,comm=None):
  if comm is None: comm = comm_world
  return globalop(a,min,"MIN",+1.e36,comm=comm)
def globalsum(a,comm=None):
  if comm is None: comm = comm_world
  return globalop(a,sum,"SUM",0.,comm=comm)
def globalave(a,comm=None):
  if comm is None: comm = comm_world
  s = globalop(a,sum,"SUM",0.,comm=comm)
  if len(shape(a)) == 0: a = [a]
  n = globalsum(len(a),comm=comm)
  if n > 0: return s/n
  else:     return 0.

# ---------------------------------------------------------------------------
# Generic parallel element-by-element operation on a distributed array.
# --- Note that this is no long needed
def parallelop(a,mpiop,comm=None):
  if not lparallel: return a
  if comm is None: comm = comm_world
  if type(a) == type(array([])):
    a1d = ravel(a) + 0
    for i in range(len(a1d)):
      a1d[i] = comm.allreduce(a1d[i],getattr(mpi,mpiop))
    a1d.shape = shape(a)
    return a1d
  else:
    return comm.allreduce(a,getattr(mpi,mpiop))

# ---------------------------------------------------------------------------
# Specific parallel element-by-element operations on a distributed array.
def parallelmax(a,comm=None):
  #return parallelop(a,"MAX")
  if not lparallel: return a
  if comm is None: comm = comm_world
  return comm.allreduce(a,maximum)
def parallelmin(a,comm=None):
  #return parallelop(a,"MIN")
  if not lparallel: return a
  if comm is None: comm = comm_world
  return comm.allreduce(a,minimum)
def parallelsum(a,comm=None):
  if not lparallel: return a
  if comm is None: comm = comm_world
  return comm.allreduce(a,mpi.SUM)


