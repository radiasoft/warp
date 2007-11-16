# Control module
ctl_version = "$Id: ctl.py,v 1.15 2007/11/16 01:09:56 dave Exp $"
from warp import *
import controllers
import signal

#############################################################################
# --- Setup signal handler to capture Control-C
# --- When a interrupt request is received, all this handler does is set a
# --- flag to that effect. Then, a subsequent call to ruthere will check
# --- that flag, and if set, raise an exception. This allows a graceful
# --- stop with the current time step completed.
# --- Set the following two in case ruthere is called before setinterrupt.
_defaultcontrolC = signal.getsignal(signal.SIGINT)
_controlCrecieved = 0
def _handlecontrolC(signum,frame):
  global _controlCrecieved
  _controlCrecieved = 1
def ruthere():
  """
Checks if an interrupt was requested (usually control-C). If so, then raise
an exception. This always restores the original interrupt handler so that the
calling code does not have to, and so that, if there is an exception, it gets
restored (since the calling code is not returned to).
  """
  global _controlCrecieved
  oldsignal = signal.signal(signal.SIGINT,_defaultcontrolC)
  if _controlCrecieved:
    _controlCrecieved = 0
    raise "Interrupt requested"
def setinterrupt():
  global _controlCrecieved,_defaultcontrolC
  _controlCrecieved = 0
  _defaultcontrolC = signal.signal(signal.SIGINT,_handlecontrolC)

def _getcommand(ext1,ext2):
  """
This retreives the appropriate command, searching through the packages
until it finds one which has the command. Note that the package can either
be a fortran package (looking for ext1) or a instance of a class which
inheritted from PackageClass (looking for ext2).
  """
  for p in package():
    pkg = packageobject(p)
    try:
      command = getattr(pkg,p+ext1)
      break
    except:
      pass
    try:
      command = getattr(pkg,ext2)
      break
    except:
      pass
  return command

#############################################################################
def generate(command=None):
  "Generates the current package"
  #setinterrupt()
  a = wtime()
  if command is None:
    command = _getcommand('gen','generate')
  command()
  controllers.aftergenerate()
  # --- Get generate time
  top.gentime = wtime() - a
  #ruthere()

def step(n=1,maxcalls=None,command=None):
  if command is None:
    command = _getcommand('exe','step')
  if maxcalls is None: maxcalls = n
  top.maxcalls = int(maxcalls)
  ncalls = n
  top.ncall = 0
  while top.ncall < ncalls:
    starttime = wtime()
    #setinterrupt()
    top.ncall = top.ncall + 1

    bb = wtime()
    controllers.callbeforestepfuncs()
    aa = wtime()
    try: step.beforetime = step.beforetime + (aa - bb)
    except: step.beforetime = 0.

    command()

    bb = wtime()
    controllers.callafterstepfuncs()
    aa = wtime()
    try: step.aftertime = step.aftertime + (aa - bb)
    except: step.aftertime = 0.

    # --- Get step time
    #top.steptime = wtime() - top.starttime - top.gentime
    # --- Flush the stdout buffer
    sys.stdout.flush()

    #ruthere()

    # --- Accumulate step time
    endtime = wtime()
    top.steptime += (endtime - starttime)

def finish(command=None):
  #setinterrupt()
  if command is None:
    command = _getcommand('fin','finish')
  try:
    command()
  except:
    pass
  #ruthere()

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

