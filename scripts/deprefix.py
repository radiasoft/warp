from warp import *
import __main__
deprefix_version = "$Id: deprefix.py,v 1.1 2000/10/16 18:34:19 dave Exp $"

# --- Creates python object for each variable in each of the specified
# --- packages. Does all packages by default.
def deprefix(pkgs=None):
  if not pkgs:
     pkgs = package()
  if type(pkgs) != type([]):
     pkgs = [pkgs]

  # --- For each package, get a list of all of the variables
  for p in pkgs:
    vlist = eval(p+'.varlist("")',globals())

    # --- For each variable, create python object in the global dictionary
    for v in vlist:
      if eval(p+'.getpyobject("'+v+'")',__main__.__dict__)!=[]:
        exec('__main__.__dict__["'+v+'"] = '+p+'.'+v,globals())

def reprefix(pkgs=None):
  if not pkgs:
     pkgs = package()
  if type(pkgs) != type([]):
     pkgs = [pkgs]

  # --- Get list of global variables
  vlist = __main__.__dict__.keys()

  # --- For each package, set the variables in the list of globals which
  # --- match a variable in the package. The ordering of the loops is done
  # --- in such a way so that package precedence controls which package
  # --- variable is set from the global variable.
  for p in pkgs:
    plist = eval(p+'.varlist("")',globals())
    for v in vlist:
      if v in plist:
	try:
	  vlist.remove(v)
          exec(p+'.'+v+' = __main__.__dict__["'+v+'"]',globals())
	except:
	  pass

