# This sets up a runcounter that can keep track of a runnumber from run to
# run. The run number is save in a file in between runs.
import warp
import time
import string
runcounter_version = "$Id: runcounter.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

def runcounter(init=0,delta=1,suffix=None,sleep=0):
  if not suffix: suffix = warp.arraytostr(warp.top.runid)
  
  try:
    # --- Try to open the runnumber file
    runnumberfile = open(suffix+"_runnumber","r")
    runnumber = string.atoi(runnumberfile.readline()) + delta
    runnumberfile.close()
  except IOError:
    # --- If no such file, then set to initial runnumber
    runnumber = init

  # --- Make sure that every processor has read in runnumber already
  if warp.npes>0: warp.mpi.barrier()

  # --- PE0 (or serial job) can now write out the next value
  if warp.me == 0:
    runnumberfile = open(suffix+"_runnumber","w")
    runnumberfile.write("%d\n"%(runnumber))
    runnumberfile.close()

  # --- Make this job wait a few seconds to make sure it is not running
  # --- the same time as another job in this series.
  # --- This is an attempt to prevent simultaneaous writes to a data file.
  if sleep > 0: time.sleep((runnumber-init)*sleep)

  return runnumber

