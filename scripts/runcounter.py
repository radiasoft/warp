# This sets up a runcounter that can keep track of a runnumber from run to
# run. The run number is save in a file in between runs.
from warp import *
import time
import string
runcounter_version = "$Id: runcounter.py,v 1.2 2002/07/17 16:54:51 dave Exp $"

def runcounter(init=0,delta=1,suffix=None,sleep=0):
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

  return runnumber

def accumulatedata(filename,datadict):
  actualdata = {}
  try:
    ff = PR.PR(filename)
    for k,v in map(None,datadict.keys(),datadict.values())
      d = ff.read(k)
      actualdata[k] = arrayappend(d,eval(v,globals()))
  except IOError:
    for k,v in map(None,datadict.keys(),datadict.values())
      actualdata[k] = eval(v,globals())
  ff = PW.PW(filename)
  for k,v in map(None,actualdata.keys(),actualdata.values())
    ff.write(k,v)
  ff.close()
