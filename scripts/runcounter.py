# This sets up a runcounter that can keep track of a runnumber from run to
# run. The run number is save in a file in between runs.
from warp import *
import time
import string
runcounter_version = "$Id: runcounter.py,v 1.8 2004/01/08 21:09:26 dave Exp $"

def runcounter(init=0,delta=1,suffix=None,sleep=0,ensambles=None):
  if not suffix: suffix = arraytostr(top.runid)
  
  try:
    # --- Try to open the runnumber file
    runnumberfile = open(suffix+"_runnumber","r")
    runnumber = string.atoi(runnumberfile.readline()) + delta
    runnumberfile.close()
  except IOError:
    # --- If no such file, then set to initial runnumber
    runnumber = init

  # --- Make sure that every processor has read in runnumber already
  if npes>0: mpi.barrier()

  # --- PE0 (or serial job) can now write out the next value
  if me == 0:
    runnumberfile = open(suffix+"_runnumber","w")
    runnumberfile.write("%d\n"%(runnumber))
    runnumberfile.close()

  # --- Make this job wait a few seconds to make sure it is not running
  # --- the same time as another job in this series.
  # --- This is an attempt to prevent simultaneaous writes to a data file.
  if sleep > 0: time.sleep((runnumber-init)*sleep)

  if ensambles is None:
    return runnumber
  else:
    r = []
    for e in ensambles:
      r.append(runnumber % e)
      runnumber = int(runnumber/e)
    r.append(runnumber)
    return tuple(r)

def accumulatedata(filename,datadict,scalardict={},pe0only=1):
  """
Accumulates data from multiple runs in a single file.
 - filename: name of the data file
 - datadict: dictionary of data to be written to the file.
             "keys" are names to be used for the data.
             "values" are the data.
             Note that the shape of the data need not be the same for each run.
             If it is not, the data will be stored as a list, otherwise as an
             array with one extra dimension.
 - scalardict={}: dictionary of scalars that will be written to the file.
                  These are not accumulated but are written directly to the
                  file. The values in the file will be from the most recent
                  call to this function.
 - pe0only=1: When true, only processor zero writes any data. Only effects
              parallel version.
  """
  if pe0only and me > 0: return
  # --- Dictionary of accumulated data to write to file.
  accumulateddata = {}

  # --- If file already exist, then open it.
  if os.access(filename,os.F_OK):
    f_in = PR.PR(filename)
    varlist = list(f_in.inquire_names())
  else:
    f_in = None
    varlist = []

  # --- Remove names in scalardict from varlist.
  for name,data in map(None,scalardict.keys(),scalardict.values()):
    if name in varlist: varlist.remove(name)

  # --- Loop over the items in the dictionary.
  for name,data in map(None,datadict.keys(),datadict.values()):
    if name in varlist: varlist.remove(name)

    # --- Read in data if file exists.
    if f_in is not None:
      try:
        originaldata = f_in.__getattr__(name)
      except:
        # --- File exists but data could not be found in it.
        originaldata = None
    else:
      originaldata = None

    if originaldata is None:
      # --- If there was no original data, then create new array.
      # --- Assume that the data will always have the same shape so
      # --- setup the array appropriately.
      dshape = shape(array(data))
      newshape = tuple(list(dshape) + [0])
      newarray = zeros(newshape,'d')
      accumulateddata[name] = arrayappend(newarray,data)
    else:
      # --- Append new data to existing data.
      if type(originaldata) is not ListType:
        if shape(originaldata)[:-1] == shape(data):
          # --- If the shape is the same as the previous data, then
          # --- continue appending the data to the array.
          accumulateddata[name] = arrayappend(originaldata,data)
        else:
          # --- The shape of the data has changed, so need to convert data into
          # --- a list. New data is append in next if block.
          dlist = []
          for i in range(shape(originaldata)[-1]):
            dlist.append(originaldata[...,i])
          accumulateddata[name] = dlist
          accumulateddata[name].append(data)
      if type(originaldata) is ListType:
        accumulateddata[name] = originaldata
        accumulateddata[name].append(data)

  # --- Read in any left over data from f_in.
  f_indict = {}
  for name in varlist:
    f_indict[name] = f_in.__getattr__(name)

  if f_in is not None: f_in.close()

  f_out = PW.PW(filename)
  for name,data in map(None,accumulateddata.keys(),accumulateddata.values()):
    f_out.write(name,data)
  for name,data in map(None,scalardict.keys(),scalardict.values()):
    f_out.write(name,data)
    if name in varlist: varlist.remove(name)
  for name in varlist:
    # --- Write out any data in the file that is not in datadict or scalardict.
    f_out.write(name,f_indict[name])
  f_out.close()



