# Control module
ctl_version = "$Id: ctl.py,v 1.2 2001/10/29 17:27:19 dave Exp $"
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
def step(n=1):
  for p in package():
    try:
      exec 'command = '+p+'exe'
      break
    except:
      pass
  #try:
  top.maxcalls = n
  top.ncall = 0
  while top.ncall < top.maxcalls:
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
