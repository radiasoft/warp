"""Controller operations

For each time, the following three functions are defined.
install___: Installs a function to be called at that specified time
uninstall___: Uninstalls the function (so it won't be called anymore)
isinstalled___: Checks if the function is installed

The functions all take a function or instance method as an argument. Note that
if an instance method is used, the user must not delete the instance,
otherwise the method will not be called and it will be removed from the list.

Functions can be called at the following times:
aftergenerate: immediately after the generate is complete
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

installaftergenerate, uninstallaftergenerate, isinstalledaftergenerate

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
controllers_version = "$Id: controllers.py,v 1.12 2006/02/16 22:17:46 dave Exp $"
def controllersdoc():
  import controllers
  print controllers.__doc__

#from warp import *
import warp
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
#
# --- This class also provides what is effectively a picklable function
# --- reference. Though it is not complete, since in some cases, functions
# --- won't be restorable.
  """

  def __init__(self,name=None,lcallonce=0):
    self.funcs = []
    self.time = 0.
    self.name = name
    self.lcallonce = lcallonce

  def __call__(self):
    "Call all of the functions in the list"
    tt = self.callfuncsinlist()
    self.time = self.time + tt
    if self.lcallonce: self.funcs = []

  def __getstate__(self):
    """
    The instance is picklable. Only functions in the list are preserved, and
    only by their names. The list of functions names replaces the funcs
    attribute in the dictionary returned. Note that nothing special is
    needed on a restore since the function names will automatically be
    converted back into functions the first time they are called
    (so there is no __setstate__).
    The ControllerFunctionContainer class below ensures that top level
    controllers are restored properly.
    """
    dict = self.__dict__.copy()
    del dict['funcs']
    funcnamelist = []
    for f in self.controllerfuncnames():
      funcnamelist.append(f)
    dict['funcs'] = funcnamelist
    return dict

  def hasfuncsinstalled(self):
    "Checks if there are any functions installed"
    return len(self.funcs) > 0

  def controllerfuncnames(self):
    "Returns the names of the functions in the list, skipping methods"
    for f in self.funcs:
      if type(f) == ListType:
        continue
      elif type(f) == StringType:
        import __main__
        if f in __main__.__dict__:
          result = f
        else:
          continue
      else:
        result = f.__name__
      yield result

  def controllerfunclist(self):
    funclistcopy = copy.copy(self.funcs)
    for f in funclistcopy:
      if type(f) == ListType:
        object = f[0]()
        if object is None:
          self.funcs.remove(f)
          continue
        result = getattr(object,f[1])
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
      if not callable(result):
        print "\n\nWarning: a controller was found that is not callable."
        print "Only callable objects can be installed."
        print "It is possible that the callable's name has been overwritten"
        print "by something not callable. This can happen during restart"
        print "if a function name had later been used as a variable name."
        if type(f) == StringType:
          print "The name of the controller is ",f
        print "\n\n"
        continue
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

  def callfuncsinlist(self,*args,**kw):
    bb = time.time()
    for f in self.controllerfunclist():
      f(*args,**kw)
    aa = time.time()
    return aa - bb

#=============================================================================

# --- Now create the actual instances.
aftergenerate = ControllerFunction('aftergenerate')
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
callafterrestartfuncs = ControllerFunction('callafterrestartfuncs',lcallonce=1)

#=============================================================================
class ControllerFunctionContainer:
  """
This is a somewhat kludgy fix to how to get any saved functions restored.
A single instance of this class is created and this instance is what is save
in a dump. This instance will have a list of the controllers, so the
controllers will be saved, but not as top level python variables.
Upon restoration, this container will go through each of the saved controllers
and reinstall the functions saved therein. This installs the functions in the
original set of controllers created when this module was first imported.
Anything that may have already been installed will therefore be unaffected.
  """
  def __init__(self,clist):
    self.clist = clist
  def __setstate__(self,dict):
    import controllers
    import __main__
    self.__dict__.update(dict)
    for c in self.clist:
      for f in c.funcs:
        # --- Check if f is already in the original list of functions,
        # --- and skip it if it is. Both the function name (f) and the
        # --- actual function in main are checked.
        # --- This will be the case if, for example, the user execs the
        # --- original input file, which sets up some functions, before
        # --- doing the restart.
        origfuncs = controllers.__dict__[c.name].funcs
        try:
          ffunc = __main__.__dict__[f]
        except KeyError:
          ffunc = None
        if (f not in origfuncs and ffunc not in origfuncs):
          controllers.__dict__[c.name].installfuncinlist(f)
    # --- The clist is obtained from the original instance in the controllers
    # --- module so that the list contains references to the original
    # --- controller instances. This is needed, since in the next dump,
    # --- this instance will be written out and must contain an updated
    # --- list of controllers.
    self.clist = controllers.__dict__['controllerfunctioncontainer'].clist

# --- This is primarily needed by warp.py so that these objects can be removed
# --- from the list of python objects which are not written out.
controllerfunctioncontainer = ControllerFunctionContainer(
                               [aftergenerate,beforefs,afterfs,
                                callscraper,calladdconductor,
                                callbeforestepfuncs,callafterstepfuncs,
                                callbeforeplotfuncs,callafterplotfuncs,
                                callplseldomfuncs,callplalwaysfuncs,
                                callafterrestartfuncs])


#=============================================================================
# ----------------------------------------------------------------------------
def installaftergenerate(f):
  "Adds a function to the list of functions called after a generate"
  aftergenerate.installfuncinlist(f)
def uninstallaftergenerate(f):
  "Removes the function from the list of functions called after a generate"
  aftergenerate.uninstallfuncinlist(f)
def isinstalledaftergenerate(f):
  return aftergenerate.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installbeforefs(f):
  "Adds a function to the list of functions called before a field-solve"
  beforefs.installfuncinlist(f)
  warp.w3d.lbeforefs = warp.true
def uninstallbeforefs(f):
  "Removes the function from the list of functions called before a field-solve"
  beforefs.uninstallfuncinlist(f)
  if not beforefs.hasfuncsinstalled(): warp.w3d.lbeforefs = warp.false
def isinstalledbeforefs(f):
  return beforefs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installafterfs(f):
  "Adds a function to the list of functions called after a field-solve"
  afterfs.installfuncinlist(f)
  warp.w3d.lafterfs = warp.true
def uninstallafterfs(f):
  "Removes the function from the list of functions called after a field-solve"
  afterfs.uninstallfuncinlist(f)
  if not afterfs.hasfuncsinstalled(): warp.w3d.lafterfs = warp.false
def isinstalledafterfs(f):
  return afterfs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installparticlescraper(f):
  "Adds a function to the list of functions called to scrape particles"
  callscraper.installfuncinlist(f)
  warp.w3d.lcallscraper = warp.true
def uninstallparticlescraper(f):
  "Removes the function from the list of functions called to scrape particles"
  callscraper.uninstallfuncinlist(f)
  if not callscraper.hasfuncsinstalled(): warp.w3d.lcallscraper = warp.false
def isinstalledparticlescraper(f):
  return callscraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def installaddconductor(f):
  "Adds a function to the list of functions called to add conductors"
  calladdconductor.installfuncinlist(f)
  warp.f3d.laddconductor = warp.true
def uninstalladdconductor(f):
  "Removes the function from the list of functions called to add conductors"
  calladdconductor.uninstallfuncinlist(f)
  if not calladdconductor.hasfuncsinstalled(): warp.f3d.laddconductor = warp.false
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

# ----------------------------------------------------------------------------
def installafterrestart(f):
  "Adds a function to the list of functions called immediately after a restart"
  callafterrestartfuncs.installfuncinlist(f)
def uninstallafterrestart(f):
  "Removes the function from the list of functions called immediately after a restart"
  callafterrestartfuncs.uninstallfuncinlist(f)
def isinstalledafterrestart(f):
  return callafterrestartfuncs.isinstalledfuncinlist(f)


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def fixcontrollersfromolddump():
  import __main__
  import controllers
  controllernames = ['aftergenerate','beforefs','afterfs','callscraper',
                     'calladdconductor','callbeforestepfuncs',
                     'callafterstepfuncs','callbeforeplotfuncs',
                     'callafterplotfuncs','callplseldomfuncs',
                     'callplalwaysfuncs']
  for cname in controllernames:
    if cname in __main__.__dict__:
      controller = __main__.__dict__[cname]
      if 'funcnamelist' in controller.__dict__:
        controllers.__dict__[controller.name] = controller
        controller.funcs = controller.funcnamelist
        del controller.funcnamelist
    else:
      print "Controller ",cname," not found"

