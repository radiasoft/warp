# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.

"""
HDF basic writer class PW by David Grote, LLNL
Modified from PW.py originally written by Paul Dubois, LLNL, to use
PDB files.
$Id: PWpyt.py,v 1.3 2003/09/16 00:09:20 dave Exp $
"""
import tables
import cPickle
import re
from Numeric import *
from types import *

class IntScalar(tables.IsDescription):
  name = tables.StringCol(50)
  value = tables.IntCol()
class FloatScalar(tables.IsDescription):
  name = tables.StringCol(50)
  value = tables.FloatCol()

class PW:
    "HDF file writer class."

    no_file_message = '(PW object not open on any file)'

    Error = 'PW error'

    def __del__(self):
        "Close any file open when this object disappears."
        self.close()

    def __init__(self, filename=None, mode="w", verbose = 1, compress=0,complib='zlib'):
        "PW(filename='', verbose=1, compress=0,complib='zlib') creates filename if given" 
        self.__dict__['file'] = None
        self.set_verbosity(verbose)
        self.set_compress(compress)
        self.set_group('/')
        if filename is not None:
            self.open(filename, mode)
        self.__dict__['fixamp'] = re.compile('@')
        self.__dict__['amprepl'] = 'XXX'

    def __setattr__(self, name, value):
        self.write(name, value)

    def __repr__(self):
        if self.is_open():
            current_mode = 'opened for writing'
            return 'HDF file %s %s.' % \
               (self.inquire_filename(), current_mode)
        else:
            return PW.no_file_message

    __str__ = __repr__

    def check_open(self):
        "check_open(): raise exception if not open for write."
        if not self.is_open():
            raise PW.Error, 'PW object not open for write.'

    def close(self):
        "close(): close the file."
        h = self.inquire_file()
        if h is not None:
            if self.inquire_verbosity():
                print "Closing HDF file being written:",self.inquire_filename()
            h.close()
        self.__dict__['file'] = None
        self.__dict__['ints'] = None
        self.__dict__['floats'] = None

    def inquire_filename(self):
        "inquire_filename() = name of this file."
        if self.is_open():
            return self.inquire_file().filename
        else:
            return ''

    def inquire_file(self):
        "inquire_file() = object open on this file."
        return self.file

    def inquire_mode(self):
        "inquire_mode() = mode('w', or 'a') of this file."
        self.check_open()
        return self.mode

    def inquire_group(self):
        "inquire_group() = present HDF group"
        return self.__dict__['group']

    def inquire_verbosity(self):
        "inquire_verbosity() = current value of verbose flag."
        return self._verbose_flag

    def inquire_compress(self):
        "inquire_compress() = current value of compress flag."
        return self._compress_flag

    def is_open(self):
        "is_open() = true if file is open"
        if self.inquire_file() is None: return 0
        return self.inquire_file().isopen

    def open(self, filename, mode = "w"):
        "open(filename, 'w')"
        self.close()
        assert mode in ['w','a'],"Improper mode: " + mode
        self.__dict__['mode'] = mode
        self.__dict__['file'] = tables.openFile(filename,mode=mode,
                                                rootUEP=self.inquire_group())
        if mode == 'a':
          # --- If mode is append, check if the ints and floats tables
          # --- have been created. If not, catch the error.
          try:
            self.__dict__['ints'] = self.inquire_file().root.ints
            self.__dict__['floats'] = self.inquire_file().root.floats
          except LookupError:
            pass
        if self.__dict__['ints'] is None:
          # --- If mode is 'w' or if the mode is 'a' but the tables have
          # --- not yet been writte, then create the tables.
          self.__dict__['ints'] = self.inquire_file().createTable(
                                        self.inquire_group(),
                                        'ints',IntScalar,"Scalar Ints")
          self.__dict__['floats'] = self.inquire_file().createTable(
                                        self.inquire_group(),
                                        'floats',FloatScalar,"Scalar Floats")

    def make_group(self, name):
        """make_group(name) 
        -- create a new HDF group, return status"""
        self.check_open()
        self.inquire_file().createGroup(self.inquire_group(),name)
        return 1

    def make_link(self, var, link):
        """make_link(var, link) 
        -- make a link, return status"""
        #self.check_open()
        #return self.inquire_file().ln(name)
        raise "Links unsupported"

    def set_group(self, name):
        """set_group(name) 
        -- change HDF group to name, return status"""
        if name[0] == '/':
          group = name
        elif len(self.inquire_group()) == 1:
          group = '/' + name
        else:
          group = self.inquire_group() + '/' + name
        self.__dict__['group'] = group

    def set_verbosity(self, flag):
        """set_verbosity(flag) sets verbosity level to flag.
        0 for quiet operation, 
        1 to report closing only, 
        2 to report access to data."""
        if 0 <= flag <= 2:
            self.__dict__['_verbose_flag'] = flag
        else:
            self.__dict__['_verbose_flag'] = 2

    def set_compress(self,flag):
        """Set level of compress when writing the HDF file"""
        self.__dict__['_compress_flag'] = flag

    def write(self, name, quantity, record = 0, indx = None):
        """Write quantity to file as 'name'"""
        self.check_open()
        if self.inquire_verbosity() > 1: 
            if record == 0:
                print "PW::write writing", name
            else:
                print "PW::write writing record", record,"of", name
        h = self.inquire_file()
        name = self.fixamp.sub(self.amprepl,name)
        if type(quantity) == IntType:
          self.ints.row['name'] = name
          self.ints.row['value'] = quantity
          self.ints.row.append()
        elif type(quantity) == FloatType:
          self.floats.row['name'] = name
          self.floats.row['value'] = quantity
          self.floats.row.append()
        elif type(quantity) in [ListType, TupleType, StringType, ArrayType]:
          if type(quantity) == ArrayType and min(array(shape(quantity))) == 0:
            raise "Array must not have a dimension of length zero"
          h.createArray(self.inquire_group(),name,quantity)
        else:
          h.createArray(self.inquire_group(),name,cPickel.dumps(quantity,bin=1),"Pickled")

    def defent(self, name, quantity, indx):
        """Define entry for quantity in file as 'name'"""
        self.check_open()
        if self.inquire_verbosity() > 1: 
            print "PW::defining entry for", name
        raise "defent not supported"


if __name__ == "__main__":
    f=PW("foo.pdb")
    a = 1
    b = 2.0
    c = "Hello world"
    from multiarray import *
    d = array([1.,2., 3.])
    e = array([1,2,3])
    g = array(["hello", "world", "array"])
    h = array([[1.,2.,3.], [4.,5.,6]])
    k = 3
    f.a = a
    f.b = b
    f.c = c
    f.d = d
    f.e = e
    f.g = g
    f.h = h
    f.close()
    f.open("foo.pdb", "a")
    f.k = k
    f.close()
# read-back test
    from PR import PR
    f = PR('foo.pdb')
    for x in f.inquire_names():
        print x, "is", eval(x), ", in file it is", eval('f.'+x)
    f.close()
# record-writing
    g = PW('goo.pdb')
    g.set_verbosity(1)
    xh = array([0.]*4)
    for i in range(len(xh)):
        x = i / 10.
        g.write('xh', x, i + 1)
        xh [i] = x
    g.close()
    g = PR('goo.pdb')
    print "xh is", xh, ", file it is ", g.xh
    g.close()










