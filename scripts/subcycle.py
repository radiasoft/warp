from warp import *
subcycle_version = "$Id: subcycle.py,v 1.3 2001/10/25 22:32:24 dave Exp $"

def subcycledoc():
  print """
Setups up subcycling
  """

class Subcycle:
  """
Setups up the code to run with subcycling, the field solve is only done
periodically. Init function takes one optional argument, the number of time
steps between field solves, defaulting to 25.
The functions disable and enable are available to turn the subcycling off and
on. The enable function is automatically called initially.
  """
  # --------------------------------------------------------------------
  def __init__(s,nsubcycle=25):
    s.nsubcycle = nsubcycle
    s.enabled = 0
    s.enable()
  # --------------------------------------------------------------------
  def enable(s):
    if not s.enabled:
      s.enabled = 1
      installbeforestep(s.dobeforestep)
      installafterstep(s.doafterstep)
  # --------------------------------------------------------------------
  def disable(s):
    if s.enabled:
      s.enabled = 0
      uninstallbeforestep(s.dobeforestep)
      uninstallafterstep(s.doafterstep)
  # --------------------------------------------------------------------
  def dobeforestep(s):
    "Turns off charge deposition for optimization"
    if (top.it+1)%s.nsubcycle == 0:
      top.depos = "vector"
    else:
      top.depos = "none"
  # --------------------------------------------------------------------
  def doafterstep(s):
    """Turns off the field solver for subcycled steps and do a field
solve when needed"""
    if top.it%s.nsubcycle != 0:
      # --- Done this was so that the correct field solve is done, even in
      # --- cases where the user changes fstype midrun.
      if top.fstype > -1: s.fstypesave = top.fstype
      top.fstype = -1
    if top.it%s.nsubcycle == 0:
      top.fstype = s.fstypesave
      fieldsol(-1)
      top.fstype = -1

