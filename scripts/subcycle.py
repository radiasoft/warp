"""
Defines class Subcycle, which does periodic fieldsolves.
"""
from warp import *
subcycle_version = "$Id: subcycle.py,v 1.8 2002/07/25 23:02:14 dave Exp $"

def subcycledoc():
  import subcycle
  print subcycle.__doc__

class Subcycle:
  """
Setups up the code to run with subcycling, the field solve is only done
periodically. Init function takes one optional argument, the number of time
steps between field solves, defaulting to 25.
The functions disable and enable are available to turn the subcycling off and
on. The enable function is automatically called initially.
The function changecycle can be called to change the number of steps between
fieldsolves.
  """
  # --------------------------------------------------------------------
  def __init__(s,nsubcycle=25):
    if nsubcycle < 1:
      print "Warning, number of steps must be greater than zero."
      print "         Default value of 25 will be used."
      nsubcycle = 25
    s.nsubcycle = nsubcycle
    s.enabled = 0
    s.enable()
    s.fstypesave = top.fstype
    if top.fstype < 0:
      print "Warning, no fieldsolve type defined. Subcycling will continue"
      print "as normal, except that no fieldsolves will be done."
      print "Set top.fstype as appropriate."
    s.depossave = top.depos
    s.laccumulate_rhosave = top.laccumulate_rho
  # --------------------------------------------------------------------
  def enable(s):
    "Enable the subcycling"
    if not s.enabled:
      s.enabled = 1
      installbeforestep(s.dobeforestep)
      installafterstep(s.doafterstep)
  # --------------------------------------------------------------------
  def disable(s):
    """
Disable the subcycling. Returns code to normal operation where a
fieldsolve is done every timestep.
    """
    if s.enabled:
      s.enabled = 0
      uninstallbeforestep(s.dobeforestep)
      uninstallafterstep(s.doafterstep)
  # --------------------------------------------------------------------
  def changecycle(s,newcycle):
    assert newcycle > 0,"Number of steps must be greater than zero"
    s.nsubcycle = newcycle

  # --------------------------------------------------------------------
  def turnoffdeposition(s):
    """
Turns off charge deposition (and zeroing of rho). Saves user's values.
    """
    try:
      s.depossave = top.depos.copy()
    except:
      s.depossave = top.depos + ''
    s.laccumulate_rhosave = top.laccumulate_rho
    top.depos = "none"
    top.laccumulate_rho = true
  # --------------------------------------------------------------------
  def turnondeposition(s):
    """
Turns on charge deposition by restoring the values.
    """
    top.depos = s.depossave
    top.laccumulate_rho = s.laccumulate_rhosave
  # --------------------------------------------------------------------
  def turnofffieldsolve(s):
    """
Turns off the field solve (save the field solve type).
    """
    s.fstypesave = top.fstype
    top.fstype = -1
  # --------------------------------------------------------------------
  def turnonfieldsolve(s):
    """
Turns on the field solver by restoring the field solve type.
    """
    top.fstype = s.fstypesave

  # --------------------------------------------------------------------
  def dobeforestep(s):
    """
Turns off charge deposition when not needed for optimization. Turns off field
solve so that it is only explicitly done after the step.
    """
    if (top.it+1)%s.nsubcycle == 0: s.turnondeposition()
    else:                           s.turnoffdeposition()
    s.turnofffieldsolve()
  # --------------------------------------------------------------------
  def doafterstep(s):
    """
Does a fieldsolve as needed. Restores quantities to the values expected
by the user.
    """
    s.turnonfieldsolve()
    s.turnondeposition()
    if top.it%s.nsubcycle == 0: fieldsol(-1)

