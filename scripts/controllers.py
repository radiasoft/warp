"""Controller operations

For each time, the following three functions are defined.
install___: Installs a function to be called at that specified time
uninstall___: Uninstalls the function (so it won't be called anymore)
isinstalled___: Checks if the function is installed

The functions all take a function or instance method as an argument. Note that
if an instance method is used, a reference to the instace must be kept,
otherwise the method will not be called. If the instance is deleted, the
method will no longer be called.

Functions can be called at the following times:
beforefs: before the field solve
afterfs: after the field solve
beforestep: before the time step
afterstep: after the time step
particlescraper: at the time that the particle boundary conditions are applied
addconductor: at the start of the multigrid solver (to install conductors)
beforeplot: before a plot (actually after a frame advance)
afterplot: after a plot (acutally just before a frame advance)
plseldom: during a special time step, when position and velocity are
          synchronized, specified by itplseldom or zzplseldom
plalways: during a special time step, when position and velocity are
          synchronized, specified by itplalways or zzplalways

Here is the complete list of functions:
installbeforefs, uninstallbeforefs, isinstalledbeforefs

installafterfs, uninstallafterfs, isinstalledafterfs

installbeforestep, uninstallbeforestep, isinstalledbeforestep

installafterstep, uninstallafterstep, isinstalledafterstep

installparticlescraper, uninstallparticlescraper, isinstalledparticlescraper

installaddconductor, uninstalladdconductor, isinstalledaddconductor

installbeforeplot, uninstallbeforeplot, isinstalledbeforeplot

installafterplot, uninstallafterplot, isinstalledafterplot

installplseldom, uninstallplseldom, isinstalledplseldom

installplalways, uninstallplalways, isinstalledplalways

"""
from __future__ import generators
controllers_version = "$Id: controllers.py,v 1.7 2004/09/09 19:43:12 dave Exp $"
def controllersdoc():
  import controllers
  print controllers.__doc__

from warp import *
from types import *
import weakref
import copy
import time

class ControllerFunction:
  """
# --- Class to handle the function lists.
# --- Note that for functions passed in that are methods of a class instance,
# --- a weak reference is saved. If this is not done, the the reference to the
# --- instance's method in the function list will preserve a reference to the
# --- instance. That means that the instance would not be deleted when all
# --- other references are deleted. If the user deletes an instance that
# --- has a method referred to in a function list, then that method will also
# --- be removed from the list.
  """

  def __init__(self,name=None):
    self.funcs = []
    self.time = 0.
    self.name = name

  def __call__(self):
    "Call all of the functions in the list"
    tt = self.callfuncsinlist()
    self.time = self.time + tt

  def __getstate__(self):
    """
    The instance is picklable. Only functions in the list are preserved, and
    only by their names.
    """
    dict = self.__dict__.copy()
    del dict['funcs']
    funcnamelist = []
    for f in self.controllerfunclist():
      if type(f) != ListType:
        funcnamelist.append(f.__name__)
    dict['funcnamelist'] = funcnamelist
    return dict

  def __setstate__(self,dict):
    """
    The instance is picklable. Only functions in the list are preserved, and
    only by their names.
    """
    self.__dict__.update(dict)
    import __main__
    self.funcs = []
    for fname in dict['funcnamelist']:
      if fname in __main__.__dict__:
        func = __main__.__dict__[fname]
        self.installfuncinlist(func)
      else:
        # --- If the function is not saved in main, then keep the name in case
        # --- it will be added later.
        self.installfuncinlist(fname)
    # --- When in instance is unpickled, make sure it replaces whatever
    # --- copy was in the controllers dict. This must be done since the
    # --- places in the code which refer to these instances, refer to the
    # --- ones in the controllers dictionary.
    # --- Also put is into main, since some instances are called from
    # --- fortran which can only access main.
    if self.name is not None:
      import controllers
      controllers.__dict__[self.name] = self
      import __main__
      __main__.__dict__[self.name] = self

  def hasfuncsinstalled(self):
    "Checks if there are any functions installed"
    return len(self.funcs) > 0

  def controllerfunclist(self):
    funclistcopy = copy.copy(self.funcs)
    for f in funclistcopy:
      if type(f) == ListType:
        object = f[0]()
        if object is None:
          self.funcs.remove(f)
          continue
        result = [object,f[1]]
      elif type(f) == StringType:
        import __main__
        if f in __main__.__dict__:
          result = __main__.__dict__[f]
          # --- If the function with the name is found, then replace the
          # --- name in the list with the function.
          self.funcs[self.funcs.index(f)] = result
        else:
          continue
      else:
        result = f
      yield result

  def installfuncinlist(self,f):
    if type(f) == MethodType:
      # --- If the function is a method of a class instance, then save a weak
      # --- reference to that instance and the method name.
      finstance = weakref.ref(f.im_self)
      fname = f.__name__
      self.funcs.append([finstance,fname])
    else:
      self.funcs.append(f)

  def uninstallfuncinlist(self,f):
    # --- An element by element search is needed
    funclistcopy = copy.copy(self.funcs)
    for func in funclistcopy:
      if f == func:
        self.funcs.remove(f)
        return
      elif type(func) == ListType and type(f) == MethodType:
        object = func[0]()
        if f.im_self is object and f.__name__ == func[1]:
          self.funcs.remove(func)
          return
      elif type(func) == StringType:
        if f.__name__ == func:
          self.funcs.remove(func)
          return
    raise 'Warning: no such function had been installed'

  def isinstalledfuncinlist(self,f):
    # --- An element by element search is needed
    funclistcopy = copy.copy(self.funcs)
    for func in funclistcopy:
      if f == func:
        return 1
      elif type(func) == ListType and type(f) == MethodType:
        object = func[0]()
        if f.im_self is object and f.__name__ == func[1]:
          return 1
      elif type(func) == StringType:
        if f.__name__ == func:
          return 1
    return 0

  def callfuncsinlist(self):
    bb = time.time()
    for f in self.controllerfunclist():
      if type(f) == ListType:
        f = getattr(f[0],f[1])
      f()
    aa = time.time()
    return aa - bb

#=============================================================================

# --- This is primarily needed by warp.py so that these objects can be removed
# --- from the list of python objects which are not written out.
controllerfunctionlist = ['beforefs','afterfs',
                          'callscraper','calladdconductor',
                          'callbeforestepfuncs','callafterstepfuncs',
                          'callbeforeplotfuncs','callafterplotfuncs',
                          'callplseldomfuncs','callplalwaysfuncs']

# --- Now create the actual instances.
beforefs = ControllerFunction('beforefs')
afterfs = ControllerFunction('afterfs')
callscraper = ControllerFunction('callscraper')
calladdconductor = ControllerFunction('calladdconductor')
callbeforestepfuncs = ControllerFunction('callbeforestepfuncs')
callafterstepfuncs = ControllerFunction('callafterstepfuncs')
callbeforeplotfuncs = ControllerFunction('callbeforeplotfuncs')
callafterplotfuncs = ControllerFunction('callafterplotfuncs')
callplseldomfuncs = ControllerFunction('callplseldomfuncs')
callplalwaysfuncs = ControllerFunction('callplalwaysfuncs')

# ----------------------------------------------------------------------------
def installbeforefs(f):
  "Adds a function to the list of functions called before a field-solve"
  beforefs.installfuncinlist(f)
  w3d.lbeforefs = true
def uninstallbeforefs(f):
  "Removes the function from the list of functions called before a field-solve"
  beforefs.uninstallfuncinlist(f)
  if not beforefs.hasfuncsinstalled(): w3d.lbeforefs = false
def isinstalledbeforefs(f):
  return beforefs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterfs(f):
  "Adds a function to the list of functions called after a field-solve"
  afterfs.installfuncinlist(f)
  w3d.lafterfs = true
def uninstallafterfs(f):
  "Removes the function from the list of functions called after a field-solve"
  afterfs.uninstallfuncinlist(f)
  if not afterfs.hasfuncsinstalled(): w3d.lafterfs = false
def isinstalledafterfs(f):
  return afterfs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installparticlescraper(f):
  "Adds a function to the list of functions called to scrape particles"
  callscraper.installfuncinlist(f)
  w3d.lcallscraper = true
def uninstallparticlescraper(f):
  "Removes the function from the list of functions called to scrape particles"
  callscraper.uninstallfuncinlist(f)
  if not callscraper.hasfuncsinstalled(): w3d.lcallscraper = false
def isinstalledparticlescraper(f):
  return callscraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installaddconductor(f):
  "Adds a function to the list of functions called to add conductors"
  calladdconductor.installfuncinlist(f)
  f3d.laddconductor = true
def uninstalladdconductor(f):
  "Removes the function from the list of functions called to add conductors"
  calladdconductor.uninstallfuncinlist(f)
  if not calladdconductor.hasfuncsinstalled(): f3d.laddconductor = false
def isinstalledaddconductor(f):
  return calladdconductor.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installbeforestep(f):
  "Adds a function to the list of functions called before a step"
  callbeforestepfuncs.installfuncinlist(f)
def uninstallbeforestep(f):
  "Removes the function from the list of functions called before a step"
  callbeforestepfuncs.uninstallfuncinlist(f)
def isinstalledbeforestep(f):
  return callbeforestepfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterstep(f):
  "Adds a function to the list of functions called after a step"
  callafterstepfuncs.installfuncinlist(f)
def uninstallafterstep(f):
  "Removes the function from the list of functions called after a step"
  callafterstepfuncs.uninstallfuncinlist(f)
def isinstalledafterstep(f):
  return callafterstepfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installbeforeplot(f):
  "Adds a function to the list of functions called before a plot"
  beforeplotfuncs.installfuncinlist(f)
def uninstallbeforeplot(f):
  "Removes the function from the list of functions called before a plot"
  beforeplotfuncs.uninstallfuncinlist(f)
def isinstalledbeforeplot(f):
  return beforeplotfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterplot(f):
  "Adds a function to the list of functions called after a plot"
  callafterplotfuncs.installfuncinlist(f)
def uninstallafterplot(f):
  "Removes the function from the list of functions called after a plot"
  callafterplotfuncs.uninstallfuncinlist(f)
def isinstalledafterplot(f):
  return callafterplotfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installplseldom(f):
  "Adds a function to the list of functions controlled by itplseldom and zzplseldom"
  callplseldomfuncs.installfuncinlist(f)
def uninstallplseldom(f):
  "Removes the function from the list of functions controlled by itplseldom and zzplseldom"
  callplseldomfuncs.uninstallfuncinlist(f)
def isinstalledplseldom(f):
  return callplseldomfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installplalways(f):
  "Adds a function to the list of functions controlled by itplalways and zzplalways"
  callplalwaysfuncs.installfuncinlist(f)
def uninstallplalways(f):
  "Removes the function from the list of functions controlled by itplalways and zzplalways"
  callplalwaysfuncs.uninstallfuncinlist(f)
def isinstalledplalways(f):
  return callplalwaysfuncs.isinstalledfuncinlist(f)

