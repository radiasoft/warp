__all__ = ['pickledump','picklerestore']
import __main__
import os
import types
import inspect
import Forthon
try:
  import numpy
except:
  pass

#try:
#  import tables
#  import hdf5pickle
#except:
if 1:
  tables = None
  import cPickle

def pickledump(filename=None,attr=['dump'],vars=[],serial=0,ff=None,
               varsuffix=None,verbose=0):
  """
Dump data into a pdb file
  - filename: dump file name
  - attr=['dump']: attribute or list of attributes of variables to dump
       Any items that are not strings are skipped. To write no variables,
       use attr=None.
  - vars=[]: list of python variables to dump
  - serial=0: switch between parallel and serial versions
  - ff=None: Allows passing in of a file object so that pydump can be called
       multiple times to pass data into the same file. Note that
       the file must be explicitly closed by the user.
  - varsuffix=None: Suffix to add to the variable names.
  - verbose=0: When true, prints out the names of the variables as they are
       written to the dump file
  - hdf=0: when true, dump into an HDF file rather than a PDB.
  """
  assert filename is not None or ff is not None,\
         'Either a filename must be specified or a tables file pointer'

  # --- Open the file if the file object was not passed in.
  # --- If the file object was passed in, then don't close it.
  if ff is None:
    if tables is not None:
      ff = tables.openFile(filename, 'w')
    else:
      ff = open(filename, 'w')
    closefile = 1
  else:
    closefile = 0

  # --- A single pickler instance is used so that references among objects
  # --- are preserved.
  if tables is not None:
    pickler = hdf5pickle.Pickler(file=ff)
  else:
    pickler = cPickle.Pickler(ff,cPickle.HIGHEST_PROTOCOL)

  # --- It turns out however, the writing multiple pickles to the same file
  # --- is broken when using arrays. So now everything is put into a single
  # --- list and it is pickled all at once.
  picklelist = []

  # --- Put the varsuffix in the file if it was specified
  if varsuffix is not None:
    if tables is not None:
      ff.setNodeAttr('/','varsuffix',varsuffix)
    else:
      #pickler.dump('varsuffix')
      picklelist.append('varsuffix')

  # --- Convert attr into a list if needed
  if not (type(attr) == types.ListType): attr = [attr]

  # --- Loop through all of the packages.
  for pname in Forthon.package():

    # --- Get the package object from the name.
    pkg = Forthon.packageobject(pname)

    # --- PackageBase instances will be written out with the python variables.
    if isinstance(pkg,Forthon.PackageBase): continue

    # --- Get variables in this package which have attribute attr.
    vlist = []
    for a in attr:
      if type(a) == types.StringType: vlist = vlist + pkg.varlist(a)

    # --- Create an empty dictionary to put everything into.
    pkgdict = {}

    # --- Loop over list of variables
    for vname in vlist:

      # --- Check if object is available (i.e. check if dynamic array is
      # --- allocated, and if not, skip it).
      v = pkg.getpyobject(vname)
      if v is None: continue

      # --- If serial flag is set, get attributes, and if it has the parallel
      # --- attribute then don't write it.
      if serial:
        a = pkg.getvarattr(vname)
        if re.search('parallel',a):
          if verbose:
            print 'variable '+vname+' skipped since it is a parallel variable'
          continue

      # --- Add the variable to the dictionary.
      pkgdict[vname] = v

    # --- Now the package dictionary can be written out
    if tables is not None:
      pickler.dump('/'+pname,pkgdict)
    else:
      #pickler.dump((pname,pkgdict))
      picklelist.append((pname,pkgdict))

  # --- Create an empty dictionary to put pytrhon variables into.
  pkgdict = {}

  # --- Put an empty dict of functions into the dict. They are handled
  # --- differently upon restore.
  pkgdict['_functions'] = {}

  # --- Now, write out the python variables (that can be written out).
  for vname in vars:

    # --- Get the value of the variable.
    vval = __main__.__dict__[vname]

    # --- Write out the source of functions. Note that the source of functions
    # --- typed in interactively is not retrieveable - inspect.getsource
    # --- returns an IOError.
    if isinstance(vval,types.FunctionType):

      source = None
      try:
        # --- Check if the source had been saved as an attribute of itself.
        # --- This allows functions to be saved that would otherwise
        # --- not be because inspect.getsource can't find them.
        source = vval._source
      except AttributeError:
        pass

      if source is None:
        # --- If the source wasn't found, try using inspect to get it.
        try:
          source = inspect.getsource(vval)
        except (IOError,NameError):
          pass

      if source is not None:
        if verbose: print 'writing python function '+vname
        # --- Clean up any indentation in case the function was defined in
        # --- an indented block of code
        while source[0] == ' ': source = source[1:]
        # --- Now put it into the dict of functions in pkgdict
        pkgdict['_functions'][vname] = source
        # --- Save the source of a function as an attribute of itself to make
        # --- retreival easier the next time.
        setattr(vval,'_source',source)
      else:
        if verbose: print 'could not write python function '+vname
      continue

    # --- Zero length arrays cannot by written out so they are skipped.
    if isinstance(vval,numpy.ndarray) and numpy.product(vval.shape) == 0:
      continue

    # --- Make sure that the object is picklable. There is no easy way of
    # --- doing this except by trying to pickle the object. This could be
    # --- very inefficient.
    try:
      testpickle = cPickle.dumps(vval)
    except:
      if verbose: print 'python variable '+vname+' could no be pickled'
      continue
    try:
      testunpickle = cPickle.loads(testpickle)
    except:
      print 'python variable '+vname+' could no be unpickled'
      continue

    if verbose: print 'writing python variable '+vname
    pkgdict[vname] = vval

  # --- Now the dictionary can be written out
  if tables is not None:
    pickler.dump('/__main__',pkgdict)
  else:
    #pickler.dump(('__main__',pkgdict))
    picklelist.append(('__main__',pkgdict))
    pickler.dump(picklelist)

  if closefile: ff.close()

#############################################################################
def picklerestore(filename,verbose=0,skip=[],
                  varsuffix=None):
  """
Restores all of the variables in the specified pickle file.
  - filename: file to read in from
  - verbose=0: When true, prints out the names of variables which are read in
  - skip=[]: list of variables to skip (now ignored)
  - varsuffix: when set, all variables read in will be given the suffix
               Note that fortran variables are then read into python vars
  """
  ff = open(filename, 'rb')

  picklelist = cPickle.load(ff)

# while 1:
#   try:
#     object = cPickle.load(ff)
#   except EOFError:
#     break

  for object in picklelist:

    # --- Check if it is varsuffix. If varsuffix was input, the varsuffix
    # --- from the file is ignored.
    if isinstance(object,types.StringType):
      if varsuffix is None: varsuffix = object
      continue

    # --- The object is otherwise a tuple (pname,dict).
    pname,dict = object

    # --- Add the varsuffix to all names, if given.
    if isinstance(varsuffix,types.StringType):
      newdict = {}
      for key,value in dict.iteritems():
        newdict[varsuffix+key] = value
      dict = newdict

    # --- Check for the __main__ dictionary. Also, if a varsuffix is
    # --- given, then put everything into main.
    if pname == '__main__' or isinstance(varsuffix,types.StringType):
      __main__.__dict__.update(dict)
      continue

    # --- The rest of the objects for be Forthon package objects

    # --- Get the package object from the name.
    try:
      pkg = __main__.__dict__[pname]
    except KeyError:
      # --- Import the package if it isn't already present.
      pkg = __import__(pname,globals(),locals())
      try:
        # --- Check if it is the right kind of object
        pkg.getfobject
      except AttributeError:
        # --- There is a module with the same name as the package. See if it
        # --- contains the package. Note that if it doesn't, skip this package
        # --- since it can't be found. Or should an exception be raised?
        try:
          pkg = getattr(pkg,pname)
        except AttributeError:
          continue

    # --- Now the dictionary can be updated.
    pkg.setdict(dict)

  # --- User defined Python functions need to be handle specially.
  try:
    functiondict = __main__._functions
    del __main__._functions
  except:
    functiondict = {}
  for vname,source in functiondict.iteritems():
    # --- Skip functions which have already been defined in case the user
    # --- has made source updates since the dump was made.
    if vname in __main__.__dict__:
      if verbose:
        print "skipping python function %s since it already is defined"%vname
    else:
      try:
        exec(source,__main__.__dict__)
        # --- Save the source of the function as an attribute of itself
        # --- so that it can be latter saved in a dump file again.
        # --- This is needed since for any functions defined here,
        # --- inspect.getsource cannot get the source.
        setattr(__main__.__dict__[vname],'_source',source)
      except:
        if verbose: print "error with function "+vname

  ff.close()

