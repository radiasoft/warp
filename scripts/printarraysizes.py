from warp import *
printarraysizes_version = "$Id: printarraysizes.py,v 1.2 2001/06/18 20:45:55 dave Exp $"

def printarraysizesdoc():
  print "printarraysizes: prints sizes of all allocated arrays"

def printarraysizes(filename=None,threshold=0):
  """Print the sizes of all allocated arrays
  - filename=None: when specified, the data will be written to that file
  - threshold=0: only arrays with sizes greater than this will be listed
  """
  if filename is not None:
    ff = open(filename,'w')
    printtofile = 1
  else:
    printtofile = 0
  for pkg in package():
    p = eval(pkg)
    vlist = p.varlist("")
    for vname in vlist:
      try:
        v = eval(pkg+'.'+vname)
        s = product(array(shape(v)))
        if s > threshold:
          if printtofile:
            ff.write("%s.%s %d\n"%(pkg,vname,s))
          else:
            print "%s.%s %d"%(pkg,vname,s)
      except:
        pass
  if printtofile: ff.close()

