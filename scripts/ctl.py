# Control module
ctl_version = "$Id: ctl.py,v 1.3 2002/11/27 01:17:54 dave Exp $"
from warp import *

def generate():
  "Generates the current package"
  for p in package():
    try:
      exec 'command = '+p+'gen'
      break
    except:
      pass
  #try:
  command(1,1)
  #except:
    #pass
  # --- Get generate time
  top.gentime = wtime() - top.starttime

beforestepfuncs = []
afterstepfuncs = []
def step(n=1,maxcalls=None):
  for p in package():
    try:
      exec 'command = '+p+'exe'
      break
    except:
      pass
  #try:
  if maxcalls is None: maxcalls = n
  top.maxcalls = maxcalls
  ncalls = n
  top.ncall = 0
  while top.ncall < ncalls:
    top.ncall = top.ncall + 1
    for f in beforestepfuncs: f()
    command(1,1)
    ruthere()
    for f in afterstepfuncs: f()
    # --- Get step time
    top.steptime = wtime() - top.starttime - top.gentime
    # --- Flush the stdout buffer
    sys.stdout.flush()
  #except:
    #pass

def finish():
  for p in package():
    try:
      exec 'command = '+p+'fin'
      break
    except:
      pass
  try:
    command(1,1)
  except:
    pass

########################################################################
def stepz(zstep=0.):
  """
Runs for the specified z distance.
  """
  zfinal = top.zbeam + zstep
  # --- Step until the beam is just before zfinal, calling step in a way
  # --- so that split leap-frog advances can be done.
  while top.zbeam + top.vbeamfrm*top.dt < zfinal:
    step(1,10)
  # --- Step until the final value is reached, synchronizing the advance
  while top.zbeam < zfinal:
    step(1)

########################################################################
def stept(tstep=0.):
  """
Runs for the specified time
  """
  tfinal = top.time + tstep
  # --- Step until the beam is just before tfinal, calling step in a way
  # --- so that split leap-frog advances can be done.
  while top.time + top.dt < tfinal:
    step(1,10)
  # --- Step until the final value is reached, synchronizing the advance
  while top.time < tfinal:
    step(1)

