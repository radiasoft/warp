from warp import *
import PW
parallelpdb_version = "$Id:"

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
        if pe != None:
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

  
