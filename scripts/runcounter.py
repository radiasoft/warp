# This sets up a runcounter that can keep track of a counter from run to
# run. The run number is save in a file in between runs.
from warp import *
import time
import string
import sys, copy
runcounter_version = "$Id: runcounter.py,v 1.12 2004/05/24 21:38:26 dave Exp $"

def runcounter(init=0,delta=1,ensambles=[],prefix=None,suffix="_runcounter",
               sleep=0):
  """
Each time runcounter is called, it will return a number incremented from the
previous result by a given delta, defaulting to 1. The state is stored
in filed so that incrementing can be done across multiple runs. This
allows for example parameter scans.
 - init=0: initial value which is returned the first time runcounter is
           called.
 - delta=1: value number if incremented by
 - ensambles=None: size of multiple ensambles
 - prefix=runid: prefix for file name where state is stored.
 - suffix="_runcounter": suffix for file name

Using ensambles, multiple counters can be returned. The number of counters
returned is one greater than the number of ensambles. When a counter reaches
its ensamble value, it is reset and the next counter is incremented. init and
delta can be lists specifying the initial value and delta for each ensamble.
They default to 0 and 1 respectively.
For example, if ensambles = [10,15], the returned value will be a tuple of
three numbers. The first between 0 and 9, the second between 0 and 14, and the
third ever increasing. Each run, the first number will be incremented, every
10 runs the second will be incremented and the first reset, etc.
  """

  # --- Set default prefix
  if not prefix: prefix = arraytostr(top.runid)

  # --- Handle ensambles
  assert type(ensambles) in [IntType,ListType,TupleType,ArrayType],\
         'ensambles must either be an integer or one of a list, tuple, or array'
  # --- Make sure it is a list and make sure that it is a copy since it is
  # --- appended to.
  if type(ensambles) == IntType: ensambles = [ensambles]
  ensambles = copy.copy(list(ensambles))
  # --- Add the maximum value of integers to the last value. This simplifies
  # --- the code below.
  ensambles.append(sys.maxint)
  
  # --- Handle init values
  assert type(init) in [IntType,ListType,TupleType,ArrayType],\
         'init must either be an integer or one of a list, tuple, or array'
  # --- Make sure it is a list
  if type(init) == IntType: init = [init]
  init = list(init)
  # --- Make sure that init is the same length as ensambles by appending zeroes
  while len(init) < len(ensambles): init.append(0)

  # --- Handle delta values
  assert type(delta) in [IntType,ListType,TupleType,ArrayType],\
         'delta must either be an integer or one of a list, tuple, or array'
  # --- Make sure it is a list
  if type(delta) == IntType: delta = [delta]
  delta = list(delta)
  # --- Make sure that delta is the same length as ensambles by appending ones
  while len(delta) < len(ensambles): delta.append(1)

  try:
    # --- Try to open the runcounter file
    runcounterfile = open(prefix+suffix,"r")
    # --- Read in the state, converting each number into an integer
    counter = map(string.atoi,string.split(runcounterfile.readline()))
    runcounterfile.close()
  except IOError:
    # --- If no such file, then flag counter
    counter = None

  # --- Increment the counter (if not the initial call)
  if counter is not None:
    # --- Add the delta to the first element
    counter[0] = counter[0] + delta[0]
    # --- Loop over elements, checking if each is greater than the
    # --- size of the ensambles
    # --- Note that because sys.maxint is in the last element, i can
    # --- never go beyond the end of the lists.
    i = 0
    while counter[i] >= ensambles[i]:
      counter[i] = init[i]
      i = i + 1
      counter[i] = counter[i] + delta[i]
  else:
    counter = init

  # --- Make sure that every processor has read in counter already
  if npes>0: mpi.barrier()

  # --- PE0 (or serial job) can now write out the next value
  if me == 0:
    runcounterfile = open(prefix+suffix,"w")
    runcounterfile.write(string.join(map(repr,counter)) + '\n')
    runcounterfile.close()

  # --- Make this job wait a few seconds to make sure it is not running
  # --- the same time as another job in this series.
  # --- This is an attempt to prevent simultaneaous writes to a data file.
  if sleep > 0: time.sleep((counter[0]-init[0])*sleep)

  # --- Return the counter as a tuple so it can be returned to multiple
  # --- variables.
  if len(ensambles) > 1:
    return tuple(counter)
  else:
    return counter[0]

def accumulatedata(filename,datadict,scalardict={},pe0only=1,
                   sortname=None):
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
 - sortname=None: If specified, the data will be sorted based on this
                  quantity.
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

  # --- Sort data if requested.
  if sortname is not None:
    assert accumulateddata.has_key(sortname),"The quantity %s must be included for sorting."%sortname
    s = accumulateddata[sortname]
    i = argsort(s)
    sorteddata = {}
    for name,data in map(None,accumulateddata.keys(),accumulateddata.values()):
      if len(transpose(data)) == len(s):
        sorteddata[name] = transpose(take(transpose(data),i))
      else:
        print "Could not sort %s since its length is different than %s"%(name,sortname)
        sorteddata[name] = data
    accumulateddata = sorteddata

  # --- Now write out the data.
  # --- File is first written into a temporary file and then copied into
  # --- the original. This is done is case there is a write error. If there
  # --- is an error, the file being written to will be corrupted and all
  # --- of the data lost.
  f_out = PW.PW(filename+'temp')
  for name,data in map(None,accumulateddata.keys(),accumulateddata.values()):
    f_out.write(name,data)
  for name,data in map(None,scalardict.keys(),scalardict.values()):
    f_out.write(name,data)
    if name in varlist: varlist.remove(name)
  for name in varlist:
    # --- Write out any data in the file that is not in datadict or scalardict.
    f_out.write(name,f_indict[name])
  f_out.close()
  os.system('mv '+filename+'temp '+filename)



