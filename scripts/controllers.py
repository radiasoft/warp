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
controllers_version = "$Id: controllers.py,v 1.3 2004/07/29 17:31:46 dave Exp $"
def controllersdoc():
  import controllers
  print controllers.__doc__

from warp import *
from types import *
import weakref
import copy

# --- Functions to handle the function lists.
# --- Note that for functions passed in that are methods of a class instance,
# --- a weak reference is saved. If this is not done, the the reference to the
# --- instance's method in the function list will preserve a reference to the
# --- instance. That means that the instance would not be deleted when all
# --- other references are deleted. If the user deletes an instance that
# --- has a method referred to in a function list, then that method will also
# --- be removed from the list.

def _controllerfunclist(flist):
  i = 0
  while i < len(flist):
    f = flist[i]
    if type(f) == ListType:
      object = f[0]()
      if object is None:
        del flist[i]
        continue
      result = [object,f[1]]
    else:
      result = flist[i]
    i = i + 1
    yield i-1,result

def _installfuncinlist(flist,f):
  if type(f) == MethodType:
    # --- If the function is a method of a class instance, then save a weak
    # --- reference to that instance and the method name.
    finstance = weakref.ref(f.im_self)
    fname = f.__name__
    flist.append([finstance,fname])
  else:
    flist.append(f)
def _uninstallfuncinlist(flist,f):
  if type(f) == MethodType:
    # --- If the function is a method of a class instance, then an element
    # --- by element search is needed to find the weak reference to the
    # --- method's instance.
    flistcopy = copy.copy(flist)
    for fins in flistcopy:
      if type(fins) == ListType:
        object = fins[0]()
        if f.im_self is object and f.__name__ == fins[1]:
          flist.remove(fins)
          return
  else:
    if f in flist:
      flist.remove(f)
      return
  raise 'Warning: no such function had been installed'
def _isinstalledfuncinlist(flist,f):
  if type(f) == MethodType:
    # --- If the function is a method of a class instance, then an element
    # --- by element search is needed to find the weak reference to the
    # --- method's instance.
    for fins in flist:
      if type(fins) == ListType:
        object = fins[0]()
        if f.im_self is object and f.__name__ == fins[1]:
          return 1
  else:
    if f in flist:
      return 1
  return 0
def _callfuncsinlist(flist):
  bb = wtime()
  flistcopy = copy.copy(flist)
  for f in flistcopy:
    # --- If the function is a method of a class instance, then an element
    # --- by element search is needed to find the weak reference to the
    # --- method's instance.
    if type(f) == ListType:
      object = f[0]()
      if object is not None:
        getattr(object,f[1])()
      else:
        # --- If the instance has been deleted (the weak reference returned
        # --- None, then remove this function from the list.
        flist.remove(f)
    else:
      f()
  aa = wtime()
  return aa - bb

#=============================================================================
# --- Setup mechanism for "before" and "after" python scripts
beforefsfuncs = []
afterfsfuncs = []
callscraperfuncs = []
addconductorfuncs = []
beforestepfuncs = []
afterstepfuncs = []
beforeplotfuncs = []
afterplotfuncs = []
plseldomfuncs = []
plalwaysfuncs = []
_controllerfuncs = {'beforefs':beforefsfuncs,
                   'afterfs':afterfsfuncs,
                   'callscraper':callscraperfuncs,
                   'addconductor':addconductorfuncs,
                   'beforestep':beforestepfuncs,
                   'afterstep':afterstepfuncs,
                   'beforeplot':beforeplotfuncs,
                   'afterplot':afterplotfuncs,
                   'plseldom':plseldomfuncs,
                   'plalways':plalwaysfuncs}

def beforefs():
  tt = _callfuncsinlist(beforefsfuncs)
  try: beforefs.time = beforefs.time + tt
  except: beforefs.time = tt
def afterfs():
  tt = _callfuncsinlist(afterfsfuncs)
  try: afterfs.time = afterfs.time + tt
  except: afterfs.time = tt
def callscraper():
  tt = _callfuncsinlist(callscraperfuncs)
  try: callscraper.time = callscraper.time + tt
  except: callscraper.time = tt
def calladdconductor():
  tt = _callfuncsinlist(addconductorfuncs)
  try: calladdconductor.time = calladdconductor.time + tt
  except: calladdconductor.time = tt
def callbeforestepfuncs():
  tt = _callfuncsinlist(beforestepfuncs)
  try: callbeforestepfuncs.time = callbeforestepfuncs.time + tt
  except: callbeforestepfuncs.time = tt
def callafterstepfuncs():
  tt = _callfuncsinlist(afterstepfuncs)
  try: callafterstepfuncs.time = callafterstepfuncs.time + tt
  except: callafterstepfuncs.time = tt
def callbeforeplotfuncs():
  tt = _callfuncsinlist(beforeplotfuncs)
  try: callbeforeplotfuncs.time = callbeforeplotfuncs.time + tt
  except: callbeforeplotfuncs.time = tt
def callafterplotfuncs():
  tt = _callfuncsinlist(afterplotfuncs)
  try: callafterplotfuncs.time = callafterplotfuncs.time + tt
  except: callafterplotfuncs.time = tt
def callplseldomfuncs():
  tt = _callfuncsinlist(plseldomfuncs)
  try: callplseldomfuncs.time = callplseldomfuncs.time + tt
  except: callplseldomfuncs.time = tt
def callplalwaysfuncs():
  tt = _callfuncsinlist(plalwaysfuncs)
  try: callplalwaysfuncs.time = callplalwaysfuncs.time + tt
  except: callplalwaysfuncs.time = tt

# ----------------------------------------------------------------------------
def installbeforefs(f):
  "Adds a function to the list of functions called before a field-solve"
  _installfuncinlist(beforefsfuncs,f)
  w3d.lbeforefs = true
def uninstallbeforefs(f):
  "Removes the function from the list of functions called before a field-solve"
  _uninstallfuncinlist(beforefsfuncs,f)
  if len(beforefsfuncs) == 0: w3d.lbeforefs = false
def isinstalledbeforefs(f):
  return _isinstalledfuncinlist(beforefsfuncs,f)

# ----------------------------------------------------------------------------
def installafterfs(f):
  "Adds a function to the list of functions called after a field-solve"
  _installfuncinlist(afterfsfuncs,f)
  w3d.lafterfs = true
def uninstallafterfs(f):
  "Removes the function from the list of functions called after a field-solve"
  _uninstallfuncinlist(afterfsfuncs,f)
  if len(afterfsfuncs) == 0: w3d.lafterfs = false
def isinstalledafterfs(f):
  return _isinstalledfuncinlist(afterfsfuncs,f)

# ----------------------------------------------------------------------------
def installparticlescraper(f):
  "Adds a function to the list of functions called to scrape particles"
  _installfuncinlist(callscraperfuncs,f)
  w3d.lcallscraper = true
def uninstallparticlescraper(f):
  "Removes the function from the list of functions called to scrape particles"
  _uninstallfuncinlist(callscraperfuncs,f)
  if len(callscraperfuncs) == 0: w3d.lcallscraper = false
def isinstalledparticlescraper(f):
  return _isinstalledfuncinlist(callscraperfuncs,f)

# ----------------------------------------------------------------------------
def installaddconductor(f):
  "Adds a function to the list of functions called to add conductors"
  _installfuncinlist(addconductorfuncs,f)
  f3d.laddconductor = true
def uninstalladdconductor(f):
  "Removes the function from the list of functions called to add conductors"
  _uninstallfuncinlist(addconductorfuncs,f)
  if len(addconductorfuncs) == 0: f3d.laddconductor = false
def isinstalledaddconductor(f):
  return _isinstalledfuncinlist(addconductorfuncs,f)

# ----------------------------------------------------------------------------
def installbeforestep(f):
  "Adds a function to the list of functions called before a step"
  _installfuncinlist(beforestepfuncs,f)
def uninstallbeforestep(f):
  "Removes the function from the list of functions called before a step"
  _uninstallfuncinlist(beforestepfuncs,f)
def isinstalledbeforestep(f):
  return _isinstalledfuncinlist(beforestepfuncs,f)

# ----------------------------------------------------------------------------
def installafterstep(f):
  "Adds a function to the list of functions called after a step"
  _installfuncinlist(afterstepfuncs,f)
def uninstallafterstep(f):
  "Removes the function from the list of functions called after a step"
  _uninstallfuncinlist(afterstepfuncs,f)
def isinstalledafterstep(f):
  return _isinstalledfuncinlist(afterstepfuncs,f)

# ----------------------------------------------------------------------------
def installbeforeplot(f):
  "Adds a function to the list of functions called before a plot"
  _installfuncinlist(beforeplotfuncs,f)
def uninstallbeforeplot(f):
  "Removes the function from the list of functions called before a plot"
  _uninstallfuncinlist(beforeplotfuncs,f)
def isinstalledbeforeplot(f):
  return _isinstalledfuncinlist(beforeplotfuncs,f)

# ----------------------------------------------------------------------------
def installafterplot(f):
  "Adds a function to the list of functions called after a plot"
  _installfuncinlist(afterplotfuncs,f)
def uninstallafterplot(f):
  "Removes the function from the list of functions called after a plot"
  _uninstallfuncinlist(afterplotfuncs,f)
def isinstalledafterplot(f):
  return _isinstalledfuncinlist(afterplotfuncs,f)

# ----------------------------------------------------------------------------
def installplseldom(f):
  "Adds a function to the list of functions controlled by itplseldom and zzplseldom"
  _installfuncinlist(plseldomfuncs,f)
def uninstallplseldom(f):
  "Removes the function from the list of functions controlled by itplseldom and zzplseldom"
  _uninstallfuncinlist(plseldomfuncs,f)
def isinstalledplseldom(f):
  return _isinstalledfuncinlist(plseldomfuncs,f)

# ----------------------------------------------------------------------------
def installplalways(f):
  "Adds a function to the list of functions controlled by itplalways and zzplalways"
  _installfuncinlist(plalwaysfuncs,f)
def uninstallplalways(f):
  "Removes the function from the list of functions controlled by itplalways and zzplalways"
  _uninstallfuncinlist(plalwaysfuncs,f)
def isinstalledplalways(f):
  return _isinstalledfuncinlist(plalwaysfuncs,f)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# --- Functions handling dumping and restoring of controllers
# ----------------------------------------------------------------------------
def controllerspreparefordump():
  # --- Convert control functions to their names so they can be written.
  # --- For methods, store the object and the method name.
  # --- Note that functions defined interactively will need to be redefined
  # --- in the restarted run since the source is not available.
  import __main__
  for n,flist in _controllerfuncs.iteritems():
    __main__.__dict__['controllercount%s'%n] = 0
    for i,f in _controllerfunclist(flist):
      __main__.__dict__['controllercount%s'%n] += 1
      if type(f) is ListType:
        # --- This won't work since the dump routine will write out a separate
        # --- copy of the object, independent of the original
        #__main__.__dict__['controller%s_%d_ref'%(n,i)] = f[0]
        #__main__.__dict__['controller%s_%d_name'%(n,i)] = f[1]
        pass
      else:
        __main__.__dict__['controller%s_%d'%(n,i)] = f.__name__
def controllerscleanafterdump():
  import __main__
  for n,flist in _controllerfuncs.iteritems():
    count = __main__.__dict__['controllercount%s'%n]
    del __main__.__dict__['controllercount%s'%n]
    for i in range(count):
      try:
        # --- This won't work since the dump routine will write out a separate
        # --- copy of the object, independent of the original
        #del __main__.__dict__['controller%s_%d_ref'%(n,i)]
        #del __main__.__dict__['controller%s_%d_name'%(n,i)]
        pass
      except KeyError:
        del __main__.__dict__['controller%s_%d'%(n,i)]
def controllersrecreatelists():
  import __main__
  for n,flist in controllerfuncs.iteritems():
    count = __main__.__dict__['controllercount%s'%n]
    for i in range(count):
      if 'controller%s_%d_ref'%(n,i) in __main__.__dict__:
        # --- This won't work since the dump routine will write out a separate
        # --- copy of the object, independent of the original
        #obj = __main__.__dict__['controller%s_%d_ref'%(n,i)]
        #name = __main__.__dict__['controller%s_%d_name'%(n,i)]
        #meth = getattr(obj,name)
        #if not _isinstalledfuncinlist(flist,meth):
        #  _installfuncinlist(flist,meth)
        pass
      else:
        fname = __main__.__dict__['controller%s_%d'%(n,i)]
        func = __main__.__dict__[fname]
        if not _isinstalledfuncinlist(flist,func):
          _installfuncinlist(flist,func)
  controllerscleanafterdump()

