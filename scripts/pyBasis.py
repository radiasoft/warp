# import all of the neccesary packages
from Numeric import *
from types import *
import RNG
import RandomArray
from pybasisC import *
import pybasisC
import os
import string
import re
try:
  import PW
  import PR
except ImportError:
  pass
try:
  import PWpyt
  import PRpyt
except ImportError:
  pass
import __main__
import sys
import cPickle
try:
  import inspect
except ImportError:
  pass
# --- Add line completion capability
try:
  import readline
except ImportError:
  pass
else:
  import rlcompleter
  readline.parse_and_bind("tab: complete")

Basis_version = "$Id: pyBasis.py,v 1.48 2004/01/27 17:25:14 dave Exp $"

if sys.platform in ['sn960510','linux-i386']:
  true = -1
  false = 0
else:
  true = 1
  false = 0

# --- Convenience function modeled after the iota of basis
def iota(low,high=None,step=1):
  if high is None:
    if step > 0:
      return arange(1,low+1,step)
    else:
      return arange(low,0,step)
  else:
    if step > 0:
      return arange(low,high+1,step)
    else:
      return arange(low,high-1,step)

# --- Converts an array of characters into a string.
def arraytostr(a,strip=true):
  a = array(a)
  if len(shape(a)) == 1:
    result = ''
    for c in a:
      result = result + c
    if strip: result = string.strip(result)
  elif len(shape(a)) == 2:
    result = []
    for i in xrange(shape(a)[1]):
      result.append(arraytostr(a[:,i]))
  return result

# --- Convenience function to do printing
def remark(s):
  print s

# --- Allows int operation on arrrays
builtinint = int
def int(x):
  if type(x) == ArrayType:
    return x.astype(Int)
  else:
    return builtinint(x)

# --- Return the nearest integer
def nint(x):
  if type(x) == ArrayType:
    return where(greater(x,0),int(x+0.5),-int(abs(x)+0.5))
  else:
    if x >= 0: return int(x+0.5)
    else: return -int(abs(x)+0.5)

# --- Replicate the sign function with two arguments. If only one is
# --- given, return the value from the Numeric sign function.
# --- This should realy be removed.
numericsign = sign
def sign(x,y=None):
  if y is None: return numericsign(x)
  if type(x) == ArrayType:
    result = where(greater(y,0.),abs(x),-abs(x))
    result = where(equal(y,0.),0.,result)
    return result
  else:
    if y > 0:
      return abs(x)
    elif y < 0:
      return -abs(x)
    else:
      return 0

# --- These are replacements for array creation routines which create
# --- arrays which have the proper ordering for fortran. When arrays created
# --- with these commands are passed to a fortran subroutine, no copies are
# --- needed to get the data into the proper order for fortran.
def fones(shape,typecode=Int):
  try:
    s = list(shape)
  except TypeError:
    s = list([shape])
  s.reverse()
  return transpose(ones(s,typecode))
def fzeros(shape,typecode=Int):
  try:
    s = list(shape)
  except TypeError:
    s = list([shape])
  s.reverse()
  return transpose(zeros(s,typecode))

# --- Convenience function which returns meshes filled with the coordinates.
def getmeshcoordinates(mins,dds,nns):
  """
getmeshcoordinates(mins,dds,nns)
Returns arrays holding the coordinates of the mesh points.
Lenght of list of inputs determines number of dimensions.
  """
  nns = tuple(array(nns) + 1)
  cc = indices(nns,'d')
  for i in xrange(len(mins)): cc[i] = mins[i] + cc[i]*dds[i]
  clist = []
  for i in xrange(len(mins)): clist.append(cc[i])
  return tuple(clist)

def getmesh2d(xmin,dx,nx,ymin,dy,ny):
  """
getmesh2d(xmin,dx,nx,ymin,dy,ny)
Returns 2 2-d arrays holding the coordinates of the mesh points
  """
  return getmeshcoordinates([xmin,ymin],[dx,dy],[nx,ny])

def getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dz,nz):
  """
getmesh3d(xmin,dx,nx,ymin,dy,ny,zmin,dx,nz)
Returns 3 3-d arrays holding the coordinates of the mesh points
  """
  return getmeshcoordinates([xmin,ymin,zmin],[dx,dy,dz],[nx,ny,nz])

# --- This function appends a new element to the end of an array.
# --- It is not very efficient since it creates a whole new array each time.
def arrayappend(x,a):
  xshape = list(shape(x))
  if type(a) == ArrayType:
    pass
  elif type(a) == ListType:
    a = array(a)
  else:
    a = array([a])
  ashape = list(shape(a))
  if len(xshape)==1 and len(ashape)==1:
    xshape[0] = xshape[0] + ashape[0]
    y = zeros(xshape,x.typecode())
    y[0:xshape[0]-ashape[0]] = x
    y[xshape[0]-ashape[0]:] = a
  elif len(xshape)>1 and len(ashape)==1 and xshape[-2]==ashape[0]:
    xshape[-1] = xshape[-1] + 1
    y = zeros(xshape,x.typecode())
    y[...,0:-1] = x
    y[...,-1] = a
  return y

# Convenience function which returns true if variable exists
def exists(x):
  """
Checks whether or not the variable whose name is specified exists in the
main dictionary.
 - x: Name of variable - must be a string.
  """
  if x in __main__.__dict__.keys(): return true
  if x in locals().keys(): return true
  if x in globals().keys(): return true
  return false

# Returns the average of the input array
def ave(x,index=0):
  if shape(x)[index] > 0:
    return sum(x,index)/shape(x)[index]
  else:
    return 0.

def averagezdata(qty,navg=0,nlines=100,n1=None,n2=None,istep=None,
                 includezeros=false):
  """
Averages data over local region. It also can down select data in the other
dimension.
  - qty: Data to be smoothed. Can be either a 1-D or 2-D array.
  - navg=0: number of data points to average over
  - nlines=100: number of lines from second dimension to choose.
  - n1=shape(qty)[0]-1:
  - n2=shape(qty)[1]-1:
  - istep=max(1,n2/nlines):
  - includezeros=false: by default, only non-zero data is averaged over.
  """
  if navg == 0 or nlines == 0: return qty
  if len(shape(qty)) == 1:
    fixqty = 1
    qty.shape = (len(qty),1)
  else:
    fixqty = 0
  if not n1: n1 = shape(qty)[0] - 1
  if not n2: n2 = shape(qty)[1] - 1
  if istep is None: istep = max(1,n2/nlines)
  hl = qty[:,::istep] + 0.
  hl[navg,:] = sum(qty[navg-navg:navg+navg+1,::istep])
  nn = 2*navg+1 + zeros(shape(hl))
  if not includezeros:
    nn[navg,:] = sum(where(qty[navg-navg:navg+navg+1,::istep]==0.,0,1),0)
  for j in range(navg+1,n1-navg-1):
    hl[j,:] = hl[j-1,:] + (qty[j+navg,::istep] - qty[j-navg-1,::istep])
    nn[j,:] = nn[j-1,:] + (+ where(qty[j+navg,::istep]==0,0,1)
                           - where(qty[j-navg-1,::istep]==0,0,1))
  nn = where(nn==0,1,nn)
  hl = where(qty[:,::istep]==0.,0.,hl)
  hl[navg+1:n1-navg-1,:] = hl[navg+1:n1-navg-1,:]/nn[navg+1:n1-navg-1,:]
  if fixqty: qty.shape = shape(qty)[0]
  if shape(qty)[1] > 1: return hl
  else: return hl[:,0]

# --- Returns the max of the multiarray
def maxnd(x):
  """Return the max element of an array of any dimension"""
  xtemp = reshape(x,tuple([product(array(x.shape))]))
  return max(xtemp)
# --- Returns the min of the multiarray
def minnd(x):
  """Return the min element of an array of any dimension"""
  xtemp = reshape(x,tuple([product(array(x.shape))]))
  return min(xtemp)

# Gets next available filename with the format 'root.nnn.suffix'.
def getnextfilename(root,suffix):
  dir = string.join(os.listdir('.'))
  i = 0
  name = root+('.%03d.'%i)+suffix
  while re.search(name,dir):
    i = i + 1
    name = root+('.%03d.'%i)+suffix
  return name

# --- Prints out the documentation of the subroutine or variable.
def doc(f,printit=1):
  # --- The for loop only gives the code something to break out of. There's
  # --- probably a better way of doing this.
  for i in range(1):
    if type(f) == StringType:
        # --- Check if it is a WARP variable
        try:
          d = listvar(f)
          break
        except NameError:
          pass
        # --- Check if it is a module name
        try:
          m = __import__(f)
          try:
            d = m.__dict__[f+'doc']()
            if d is None: d = ''
          except KeyError:
            d = m.__doc__
          break
        except ImportError:
          pass
        # --- Try to get the actual value of the object
        try:
          v = __main__.__dict__[f]
          d = v.__doc__
          break
        except KeyError:
          d = "Name not found"
        except AttributeError:
          d = "No documentation found"
    else:
      # --- Check if it has a doc string
      try:
        d = f.__doc__
      except AttributeError:
        d = "No documentation found"
  if printit: print d
  else:       return d

# --- Get size of all variables in a group
def getgroupsize(pkg,grp):
  ll = pkg.varlist(grp)
  ss = 0
  for v in ll:
    vv = pkg.getpyobject(v)
    if type(vv) == type(array([1])):
      ss = ss + product(array(shape(vv)))
    else:
      ss = ss + 1
  return ss

# --- Print out all variables in a group
def printgroup(pkg,group='',maxelements=10):
  """
Print out all variables in a group or with an attribute
  - pkg: package name
  - group: group name
  - maxelements=10: only up to this many elements of arrays are printed
  """
  if type(pkg) == StringType: pkg = __main__.__dict__[pkg]
  vlist = pkg.varlist(group)
  if not vlist:
    print "Unknown group name "+group
    return
  for vname in vlist:
    v = pkg.getpyobject(vname)
    if v is None:
      print vname+' is not allocated'
    elif type(v) != ArrayType:
      print vname+' = '+str(v)
    else:
      if v.typecode() == 'c':
        print vname+' = "'+str(arraytostr(v))+'"'
      elif size(v) <= maxelements:
        print vname+' = '+str(v)
      else:
        if rank(v) == 1:
          print vname+' = '+str(v[:maxelements])[:-1]+" ..."
        else:
          if shape(v)[0] <= maxelements:
            if rank(v) == 2:
              print vname+' = ['+str(v[:,0])+"] ..."
            elif rank(v) == 3:
              print vname+' = [['+str(v[:,0,0])+"]] ..."
            elif rank(v) == 4:
              print vname+' = [[['+str(v[:,0,0,0])+"]]] ..."
            elif rank(v) == 5:
              print vname+' = [[[['+str(v[:,0,0,0,0])+"]]]] ..."
            elif rank(v) == 6:
              print vname+' = [[[[['+str(v[:,0,0,0,0,0])+"]]]]] ..."
          else:
            if rank(v) == 2:
              print vname+' = ['+str(v[:maxelements,0])[:-1]+" ..."
            elif rank(v) == 3:
              print vname+' = [['+str(v[:maxelements,0,0])[:-1]+" ..."
            elif rank(v) == 4:
              print vname+' = [[['+str(v[:maxelements,0,0,0])[:-1]+" ..."
            elif rank(v) == 5:
              print vname+' = [[[['+str(v[:maxelements,0,0,0,0])[:-1]+" ..."
            elif rank(v) == 6:
              print vname+' = [[[[['+str(v[:maxelements,0,0,0,0,0])[:-1]+" ..."
  
##############################################################################
##############################################################################
def pydumpbasisobject(ff,attr,objname,obj,varsuffix,writtenvars,fobjlist,
                      serial,verbose):
  # --- General work of this object
  if verbose: print "object "+objname+" being written"
  # --- Write out the value of fobj so that in restore, any links to this
  # --- object can be restored. Only do this if fobj != 0, which means that
  # --- it is not a top level package, but a variable of fortran derived type.
  fobj = obj.getfobject()
  if fobj != 0:
    ff.write('FOBJ'+varsuffix,fobj)
    ff.write('TYPENAME'+varsuffix,obj.gettypename())
    # --- If this object has already be written out, then return.
    if fobj in fobjlist: return
    # --- Add this object to the list of object already written out.
    fobjlist.append(fobj)
  # --- Get variables in this package which have attribute attr.
  vlist = []
  for a in attr:
    if type(a) == StringType: vlist = vlist + obj.varlist(a)
  # --- Loop over list of variables
  for vname in vlist:
    # --- Check if object is available (i.e. check if dynamic array is
    # --- allocated).
    v = obj.getpyobject(vname)
    if v is None: continue
    # --- If serial flag is set, get attributes and if has the parallel
    # --- attribute, don't write it.
    if serial:
      a = obj.getvarattr(vname)
      if re.search('parallel',a):
        if verbose: print "variable "+vname+varsuffix+" skipped since it is a parallel variable"
        continue
    # --- Check if variable is a complex array. Currently, these
    # --- can not be written out.
    if type(v) == ArrayType and v.typecode() == Complex:
      if verbose: print "variable "+vname+varsuffix+" skipped since it is a complex array"
      continue
    # --- Check if variable with same name has already been written out.
    # --- This only matters when the variable is being written out as
    # --- a plane python variable.
    if '@' not in varsuffix:
      if vname in writtenvars:
        if verbose: print "variable "+objname+"."+vname+" skipped since other variable would have same name in the file"
        continue
      writtenvars.append(vname)
    # --- Check if variable is a PyBasisType, if so, recursively call this
    # --- function.
    if type(v) == PyBasisType:
      # --- Note that the attribute passed in is blank, since all components
      # --- are to be written out to the file.
      pydumpbasisobject(ff,[''],vname,v,'@'+vname+varsuffix,writtenvars,
                        fobjlist,serial,verbose)
      continue
    # --- If this point is reached, then variable is written out to file
    if verbose: print "writing "+objname+"."+vname+" as "+vname+varsuffix
    ff.write(vname+varsuffix,v)

##############################################################################
# Python version of the dump routine. This uses the varlist command to
# list of all of the variables in each package which have the
# attribute attr (and actually attr could be a group name too). It then
# checks on the state of the python object, making sure that unallocated
# arrays are not written out.  Finally, the variable is written out to the
# file with the name in the format vame@pkg.  Additionally, python
# variables can be written to the file by passing in a list of the names
# through vars. The '@' sign is used between the package name and the
# variable name so that no python variable names can be clobbered ('@'
# is not a valid character in python names). The 'ff.write' command is
# used, allowing names with an '@' in them. The writing of python variables
# is put into a 'try' command since some variables cannot be written to
# a pdb file.
def pydump(fname=None,attr=["dump"],vars=[],serial=0,ff=None,varsuffix=None,
           verbose=false,hdf=0):
  """
Dump data into a pdb file
  - fname: dump file name
  - attr=["dump"]: attribute or list of attributes of variables to dump
       Any items that are not strings are skipped. To write no variables,
       use attr=None.
  - vars=[]: list of python variables to dump
  - serial=0: switch between parallel and serial versions
  - ff=None: Allows passing in of a file object so that pydump can be called
       multiple times to pass data into the same file. Note that
       the file must be explicitly closed by the user.
  - varsuffix=None: Suffix to add to the variable names. If none is specified,
       the suffix '@pkg' is used, where pkg is the package name that the
       variable is in. Note that if varsuffix is specified, the simulation
       cannot be restarted from the dump file.
  - verbose=false: When true, prints out the names of the variables as they are
       written to the dump file
  - hdf=0: when true, dump into an HDF file rather than a PDB.
  """
  assert fname is not None or ff is not None,\
         "Either a filename must be specified or a pdb file pointer"
  # --- Open the file if the file object was not passed in.
  # --- If the file object was passed in, then don't close it.
  if ff is None:
    if not hdf:
      # --- Try to open file with PDB format as requested.
      try:
        ff = PW.PW(fname)
        # --- With PDB, pickle dumps can only be done in ascii.
        dumpsmode = 0
      except:
        pass
    if hdf or ff is None:
      # --- If HDF requested or PDB not available, try HDF.
      try:
        ff = PWpyt.PW(fname)
        # --- An advantage of HDF is that pickle dumps can be done in binary
        dumpsmode = 1
      except:
        pass
    if hdf and ff is None:
      # --- If HDF was requested and didn't work, try PDB anyway.
      try:
        ff = PW.PW(fname)
        # --- With PDB, pickle dumps can only be done in ascii.
        dumpsmode = 0
      except:
        pass
    assert ff is not None,"Dump file cannot be opened, no data formats available"
    closefile = 1
  else:
    try:
      if ff.file_type == "HDF":
        dumpsmode = 1
      else:
        dumpsmode = 0
    except:
      dumpsmode = 0
    closefile = 0
  # --- Make sure the file has a file_type. Older versions of the pdb
  # --- wrappers did not define a file type.
  try:
    ff.file_type
  except:
    ff.__dict__["file_type"] = "oldPDB"
  # --- Convert attr into a list if needed
  if not (type(attr) == ListType): attr = [attr]
  # --- Loop through all of the packages (getting pkg object).
  # --- When varsuffix is specified, the list of variables already written
  # --- is created. This solves two problems. It gives proper precedence to
  # --- variables of the same name in different packages. It also fixes
  # --- an obscure bug in the pdb package - writing two different arrays with
  # --- the same name causes a problem and the pdb file header is not
  # --- properly written. The pdb code should really be fixed.
  pkgsuffix = varsuffix
  packagelist = package()
  writtenvars = []
  fobjlist = []
  for pname in packagelist:
    pkg = __main__.__dict__[pname]
    if varsuffix is None: pkgsuffix = '@' + pname
    pydumpbasisobject(ff,attr,pname,pkg,pkgsuffix,writtenvars,fobjlist,
                      serial,verbose)

  # --- Now, write out the python variables (that can be written out).
  # --- If supplied, the varsuffix is append to the names here too.
  if varsuffix is None: varsuffix = ''
  for vname in vars:
    # --- Skip python variables that would overwrite fortran variables.
    if len(writtenvars) > 0:
      if vname in writtenvars:
        if verbose: print "variable "+vname+" skipped since other variable would have same name in the file"
        continue
    # --- Get the value of the variable.
    vval = __main__.__dict__[vname]
    # --- Write out the source of functions. Note that the source of functions
    # --- typed in interatively is not retrieveable - inspect.getsource
    # --- returns an IOError.
    if type(vval) in [FunctionType]:
      try:
        source = inspect.getsource(vval)
        #if verbose:
        if verbose: print "writing python function "+vname+" as "+vname+varsuffix+'@function'
        ff.write(vname+varsuffix+'@function',source)
      except (IOError,NameError):
        if verbose: print "could not write python function "+vname
      continue
    # --- Zero length arrays cannot by written out.
    if type(vval) == ArrayType and product(array(shape(vval))) == 0:
      continue
    # --- Check if variable is a PyBasisType
    if type(vval) == PyBasisType:
      pydumpbasisobject(ff,attr,vname,vval,'@'+vname+varsuffix,writtenvars,
                        fobjlist,serial,verbose)
      continue
    # --- Try writing as normal variable.
    # --- The docontinue temporary is needed since python1.5.2 doesn't
    # --- seem to like continue statements inside of try statements.
    docontinue = 0
    try:
      if verbose: print "writing python variable "+vname+" as "+vname+varsuffix
      ff.write(vname+varsuffix,vval)
      docontinue = 1
    except:
      pass
    if docontinue: continue
    # --- If that didn't work, try writing as a pickled object. This is only
    # --- needed for old pdb files.
    if ff.file_type == 'oldPDB':
      try:
        if verbose:
          print "writing python variable "+vname+" as "+vname+varsuffix+'@pickle'
        ff.write(vname+varsuffix+'@pickle',cPickle.dumps(vval,dumpsmode))
        docontinue = 1
      except (cPickle.PicklingError,TypeError):
        pass
      if docontinue: continue
    # --- All attempts failed so write warning message
    if verbose: print "cannot write python variable "+vname
  if closefile: ff.close()


# Python version of the restore routine. It restores all of the variables
# in the pdb file, including any that are not part of a pybasis package.
# An '@' in the name distinguishes between the two. The 'ff.__getattr__' is
# used so that variables with an '@' in the name can be read. The reading
# in of python variables is put in a 'try' command to make it idiot proof.
# More fancy foot work is done to get new variables read in into the
# global dictionary.
def pyrestore(filename=None,fname=None,verbose=0,skip=[],ff=None,varsuffix=None,ls=0):
  """
Restores all of the variables in the specified file.
  - filename: file to read in from (assumes PDB format)
  - verbose=0: When true, prints out the names of variables which are read in
  - skip=[]: list of variables to skip
  - ff=None: Allows passing in of a file object so that pydump can be called
       multiple times to pass data into the same file. Note that
       the file must be explicitly closed by the user.
  - varsuffix: when set, all variables read in will be given the suffix
               Note that fortran variables are then read into python vars
  - ls=0: when true, prints a list of the variables in the file
          when 1 prints as tuple
          when 2 prints in a column
Note that it will automatically detect whether the file is PDB or HDF.
  """
  assert filename is not None or fname is not None or ff is not None,\
         "Either a filename must be specified or a pdb file pointer"
  if ff is None:
    # --- The original had fname, but changed to filename to be consistent
    # --- with restart and dump.
    if filename is None: filename = fname
    # --- Make sure a filename was input.
    assert filename is not None,"A filename must be specified"
    # --- open pdb file
    try:
      ff = PR.PR(filename)
    except:
      ff = PRpyt.PR(filename)
    closefile = 1
  else:
    closefile = 0
  try:
    ff.file_type
  except:
    ff.__dict__["file_type"] = "oldPDB"
  # --- Get a list of all of the variables in the file, loop over that list
  vlist = ff.inquire_names()
  # --- Print list of variables
  if ls:
    if ls == 1:
      print vlist
    else:
      for l in vlist: print l

  # --- First, sort out the list of variables
  groups = _sortvarsbysuffix(vlist,skip)
  fobjdict = {}

  # --- Read in the variables with the standard suffices.

  # --- These would be interpreter variables written to the file
  # --- from python (or other sources). A simple assignment is done and
  # --- the variable in put in the main dictionary.
  if groups.has_key(''):
    plist = groups['']
    del groups['']
    for vname in plist:
      pyname = vname
      if varsuffix is not None: pyname = pyname + str(varsuffix)
      try:
        if verbose: print "reading in python variable "+vname
        __main__.__dict__[pyname] = ff.__getattr__(vname)
      except:
        if verbose: print "error with variable "+vname

  # --- These would be interpreter variables written to the file
  # --- as pickled objects. The data is unpickled and the variable
  # --- in put in the main dictionary. This is only needed for
  # --- the old pdb wrapper.
  if groups.has_key('pickle') and ff.file_type == "oldPDB":
    picklelist = groups['pickle']
    del groups['pickle']
    for vname in picklelist:
      pyname = vname
      if varsuffix is not None: pyname = pyname + str(varsuffix)
      try:
        if verbose: print "reading in pickled variable "+vname
        __main__.__dict__[pyname]=cPickle.loads(ff.__getattr__(vname+'@pickle'))
      except:
        if verbose: print "error with variable "+vname

  # --- These would be interpreter variables written to the file
  # --- from Basis. A simple assignment is done and the variable
  # --- in put in the main dictionary.
  if groups.has_key('global'):
    globallist = groups['global']
    del groups['global']
    for vname in globallist:
      pyname = vname
      if varsuffix is not None: pyname = pyname + str(varsuffix)
      try:
        if verbose: print "reading in Basis variable "+vname
        __main__.__dict__[pyname] = ff.__getattr__(vname+'@global')
      except:
        if verbose: print "error with variable "+vname

  # --- User defined Python functions
  if groups.has_key('function'):
    functionlist = groups['function']
    del groups['function']
    for vname in functionlist:
      # --- Skip functions which have already been defined in case the user
      # --- has made source updates since the dump was made.
      if __main__.__dict__.has_key(vname): 
        if verbose:
          print "skipping python function %s since it already is defined"%vname
      else:
        try:
          if verbose: print "reading in python function"+vname
          source = ff.__getattr__(vname+'@function')
          exec(source,__main__.__dict__)
        except:
          if verbose: print "error with function "+vname

  # --- Ignore variables with suffix @parallel
  if groups.has_key('parallel'):
    del groups['parallel']

  for gname in groups.keys():
    pyrestorepybssisobject(ff,gname,groups[gname],fobjdict,varsuffix,
                           verbose,doarrays=0)
  for gname in groups.keys():
    pyrestorepybssisobject(ff,gname,groups[gname],fobjdict,varsuffix,
                           verbose,doarrays=1)

  if closefile: ff.close()

def _sortvarsbysuffix(vlist,skip):
  # --- Sort the variables, collecting them in groups based on there suffix.
  groups = {}
  for v in vlist:
    if '@' in v:
      i = string.rfind(v,'@')
      vname = v[:i]
      gname = v[i+1:]
    else:
      # --- Otherwise, variable is plain python variable.
      vname = v
      gname = ''

    # --- If variable is in the skip list, then skip
    if (vname in skip or
        (len(v) > 4 and v[-4]=='@' and v[-3:]+'.'+v[:-4] in skip)):
#     if verbose: print "skipping "+v
      continue

    # --- Now add the variable to the appropriate group list.
    groups.setdefault(gname,[]).append(vname)

  return groups

# This is a list of variables that have been renamed at some point. Note that
# this should be cleaned up periodically.
renamed = {'f3d.bound0':'w3d.bound0',
           'f3d.boundnz':'w3d.boundnz',
           'f3d.boundxy':'w3d.boundxy',
           'top.prwelips':'top.prwelipz'}

def pyrestorepybssisobject(ff,gname,vlist,fobjdict,varsuffix,verbose,doarrays):

  # --- Convert gname in pdb-style name
  gsplit = string.split(gname,'.')
  gsplit.reverse()
  gpdbname = string.join(gsplit,'@')

  # --- Check is the variable gname exists or is allocated.
  # --- If not, create a new variable.
  neednew = 0
  try:
    v = eval(gname,__main__.__dict__)
    if v is None: neednew = 1
  except:
    neednew = 1

  if neednew:
    # --- A new variable needs to be created.
    try:
      fobj = ff.read("FOBJ@"+gpdbname)
    except:
      return
    # --- First, check if the object has already be restored.
    if fobj in fobjdict:
      # --- If so, then point new variable to existing object
      exec("%s = %s"%(gname,fobjdict[fobj]),__main__.__dict__)
      # return ???
    else:
      # --- Otherwise, create a new instance of the appropriate type,
      # --- and add it to the list of objects.
      typename = ff.read("TYPENAME@"+gpdbname)
      exec("%s = %s()"%(gname,typename),__main__.__dict__)
      fobjdict[fobj] = gname

  # --- Sort out the list of variables
  groups = _sortvarsbysuffix(vlist,[])

  # --- Get "leaf" variables
  if groups.has_key(''):
    leafvars = groups['']
    del groups['']
  else:
    leafvars = []

  # --- Read in leafs.
  for vname in leafvars:
    if vname == 'FOBJ' or vname == 'TYPENAME': continue
    fullname = gname + '.' + vname
    vpdbname = vname + '@' + gpdbname

    if fullname in renamed.keys(): fullname = renamed[fullname]

    # --- Add suffix to name if given.
    # --- varsuffix is wrapped in str in case a nonstring was passed in.
    if varsuffix is not None: fullname = vname + str(varsuffix)

    try:
      if type(ff.__getattr__(vpdbname)) != ArrayType and not doarrays:
        # --- Simple assignment is done for scalars, using the exec command
        if verbose: print "reading in "+fullname
        exec(fullname+'=ff.__getattr__(vpdbname)',__main__.__dict__,locals())
      elif type(ff.__getattr__(vpdbname)) == ArrayType and doarrays:
        pkg = eval(gname,__main__.__dict__)
        # --- forceassign is used, allowing the array read in to have a
        # --- different size than the current size of the warp array.
        if verbose: print "reading in "+gname+"."+fullname
        pkg.forceassign(vname,ff.__getattr__(vpdbname))
    except:
      # --- The catches errors in cases where the variable is not an
      # --- actual warp variable, for example if it had been deleted
      # --- after the dump was originally made.
      print "Warning: There was problem restoring %s"% (fullname)

  # --- Read in rest of groups.
  for g,v in groups.items():
    pyrestorepybssisobject(ff,gname+'.'+g,v,fobjdict,varsuffix,verbose,doarrays)


# --- create an alias for pyrestore
restore = pyrestore

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

def Basisdoc():
  print """
iota(): returns a array of sequential integers
arraytostr(): converts an array of chars to a string
remark(): same as print
int(): converts data to integer
nint(): converts data to nearest integer
sign(): emulation of sign function
fones(): returns multi-dimensional array with fortran ordering
fzeros(): returns multi-dimensional array with fortran ordering
arrayappend(): appends to multi-dimensional array
exists(): checks if a variable exists
ave(): averages an array of numbers
maxnd(): finds max of multi-dimensional array
minnd(): finds min of multi-dimensional array
getnextfilename(): finds next available file name in a numeric sequence
doc(): prints info about variables and functions
printgroup(): prints all variables in the group or with an attribute
pydump(): dumps data into pdb format file
pyrestore(): reads data from pdb format file
restore(): equivalent to pyrestore
"""
