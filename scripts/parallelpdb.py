from warp import *
import PW
parallelpdb_version = "$Id: parallelpdb.py,v 1.6 2007/12/20 00:36:41 dave Exp $"

# Routines to allow all processors to write data to one file.
# The PE's sequentially append data to the file. When a file is open,
# PE0 keeps it open except during a write. Then, PE0 first writes its data
# to the file and the temporarily closes the file. The other processors
# open it, append to it, and then close it. When the last processor is
# done, PE0 reopens it.  This allows the user to access data in the file
# at any time (via PE0).

class ParallelPW(PW.PW):

    def __init__ (self, filename='', mode="w", verbose = 1):
        "PW(filename='', verbose=1) creates filename if given"
        self.__dict__['_handle'] = None
        self.set_verbosity (verbose)
        if filename:
          if me == 0:
            self.open (filename, mode)
        self.__dict__['filename'] = filename
        self.__dict__['mode'] = mode

    def write (self, name, quantity, record = 0, pe=None):
        """Write quantity to file as 'name'"""
	# --- Only one processor is to write to the file.
        if pe is not None:
	  # --- If it is not pe0, then pe0 should close the file and the other
	  # --- should open it.
	  if pe != 0:
	    if me == 0:
	      self.close()
	    mpi.barrier()
	    if me == pe:
              self.open(self.filename,'a')
          if me == pe:
            self.inquire_handle().write(name,quantity,record,None)
	  # --- Now, that pe should close it and pe0 should open it again.
	  if pe != 0:
	    if me == pe:
	      self.close()
	    mpi.barrier()
	    if me == 0:
              self.open(self.filename,'a')
          return
	# -------------------------------------------------
	# --- Coding for all pe's writing to the file.
	nn = gatherall(list(shape(quantity)))
	ndims = len(nn[0])
	# --- Check for consistency of the data shape
	for i in range(npes):
	  if len(nn[i]) != ndims:
	    raise "number of dimensions must be the same on all processors"
	  for id in range(ndims-1):
	    if nn[i][id] != nn[0][id]:
	      raise "all dimensions must be of the same length except the last"
	# --- Get partial sum of last dimension
	nlast = zeros(npes+1,'l')
	for i in range(npes):
	  nlast[i+1] = nlast[i] + nn[i][-1]
	# --- Now, PE0 creates space in the file and closes it.
        if me == 0:
          if not self.is_open(): self.check_open()
          nn[0][-1] = nlast[-1]
	  self.defent(name,quantity,tuple(nn[0]))
	  self.close()
	# --- Now, all open file at once and write out data simultaneously.
	mpi.barrier()
        self.open(self.filename,'a')
	index = ndims*[0]
        index[-1] = nlast[me]
        self.inquire_handle().write(name,quantity,0,tuple(index))
	# --- All but PE0 close the file.
	if me != 0:
	  self.close()

  
# ----------------------------------------------------------------------
# --- The old version
class ParallelPWold(PW.PW):

    def __init__ (self, filename='', mode="w", verbose = 1):
        "PW(filename='', verbose=1) creates filename if given"
        self.__dict__['_handle'] = None
        self.set_verbosity (verbose)
        if filename:
          if me == 0:
            self.open (filename, mode)
        self.__dict__['filename'] = filename
        self.__dict__['mode'] = mode

    def write (self, name, quantity, record = 0, pe=None):
        """Write quantity to file as 'name'"""
        if pe is not None:
          if me == pe:
            PW.PW.write(self,name,quantity,record)
          return
        if me == 0:
          if not self.is_open():
            mpi.send(0,1,1)
            self.check_open()
        else:
          i = mpi.recv(me-1,1)
          if not i:
            if me < npes-1:
              mpi.send(1,(me+1) % npes,1)
            return
          record = 2
          self.open(self.filename,'a')
        if self.inquire_verbosity () > 1:
            if record == 0:
                print "PW::write writing", name
            else:
                print "PW::write writing record", record, \
                      "of", name
        self.inquire_handle ().write (name, quantity, record)
        self.close()
        mpi.send(1,(me+1) % npes,1)
        if me == 0:
          i = mpi.recv(npes-1,1)
          self.open(self.filename,'a')

  
