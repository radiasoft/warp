"""
HDF self-describing file format writer class PW
by David Grote, LLNL
Modified from PW.py originally written by Paul Dubois, LLNL, to use
PDB files.
$Id: PWhdf.py,v 1.3 2003/08/07 19:00:57 dave Exp $
"""
import _pyhl
from Numeric import *
from types import *

class PW:
    "HDF file writer class."

    no_file_message = '(PW object not open on any file)'

    Error = 'PW error'

    type_dict = {'l':'int','d':'double','c':'string',
                 IntType:'int',FloatType:'double',StringType:'string'}

    def __del__ (self):
        "Close any file open when this object disappears."
        self.close()

    def __init__ (self, filename='', mode="w", verbose = 1, compression = 6):
        "PW(filename='', verbose=1, compression=6) creates filename if given" 
        self.__dict__['_nodelist'] = None
        self.set_verbosity (verbose)
        self.set_compression (compression)
        if filename:
            self.open (filename, mode)

    def __setattr__ (self, name, value):
        self.write (name, value)

    def __repr__ (self):
        if self.is_open ():
            current_mode = 'opened for writing'
            return 'HDF file %s %s.' % \
               (self.inquire_filename(), current_mode)
        else:
            return PW.no_file_message

    __str__ = __repr__

    def check_open (self):
        "check_open (): raise exception if not open for write."
        if not self.is_open ():
            raise PW.Error, 'PW object not open for write.'

    def close (self):
        "close(): close the file."
        h = self.inquire_nodelist ()
        if not (h is None):
            if self.inquire_verbosity():
                print "Closing HDF file being written:", \
                      self.inquire_filename()
            h.write(self.inquire_filename(),self.inquire_compression())
        self.__dict__['_nodelist'] = None

    def inquire_filename (self):
        "inquire_filename() = name of this file."
        if self.is_open ():
            return self._filename
        else:
            return ''

    def inquire_nodelist (self):
        "inquire_nodelist () = nodelist object open on this file."
        return self._nodelist

    def inquire_mode (self):
        "inquire_mode () = mode ('w', or 'a') of this file. Only 'w' supported"
        self.check_open()
        return 'w'

    def inquire_pwd(self):
        "inquire_pwd () = present HDF directory"
        self.check_open ()
        return self._pwd

    def inquire_verbosity (self):
        "inquire_verbosity () = current value of verbose flag."
        return self._verbose_flag

    def inquire_compression (self):
        "inquire_compression () = current value of compression flag."
        return self._compression_flag

    def is_open (self):
        "is_open () = true if file is open"
        if self.inquire_nodelist () == None:
            return 0
        return 1

    def open (self, filename, mode = "w"):
        "open (filename, 'w')"
        self.close ()
        self.__dict__['_filename'] = filename
        if mode == "w":
            self.__dict__['_nodelist'] = _pyhl.nodelist()
            self.set_directory("/data")
            self.make_directory(self._pwd)
        #elif mode == "a":
            #self.__dict__['_nodelist'] = pypdb.append (filename)
        else:
            raise ValueError, "Improper mode: " + mode


    def new_directory_name(self,name):
        if name[0] == '/': return name
        else:              return self._pwd + '/' + name

    def make_directory(self, name):
        """make_directory(name) 
        -- create a new HDF directory, return status"""
        self.check_open()
        self.inquire_nodelist().addNode(_pyhl.node(_pyhl.GROUP_ID,
                                                self.new_directory_name(name)))
        return 1

    def make_link(self, var, link):
        """make_directory(var, link) 
        -- make a link, return status"""
        #self.check_open()
        #return self.inquire_nodelist().ln(name)
        raise "Links unsupported"

    def set_directory (self, name):
        """set_directory(name) 
        -- change HDF directory to name, return status"""
        self.check_open()
        self.__dict__['_pwd'] = self.new_directory_name(name)
        return 1

    def set_verbosity (self, flag):
        """set_verbosity (flag) sets verbosity level to flag.
        0 for quiet operation, 
        1 to report closing only, 
        2 to report access to data."""
        if 0 <= flag <= 2:
            self.__dict__['_verbose_flag'] = flag
        else:
            self.__dict__['_verbose_flag'] = 2

    def set_compression(self,flag):
        """Set level of compression when writing the HDF file"""
        self.__dict__['_compression_flag'] = flag

    def write (self, name, quantity, record = 0, indx = None):
        """Write quantity to file as 'name'"""
        self.check_open()
        if self.inquire_verbosity () > 1: 
            if record == 0:
                print "PW::write writing", name
            else:
                print "PW::write writing record", record, \
                      "of", name
        if type(quantity) is ArrayType:
          anode = _pyhl.node(_pyhl.DATASET_ID,self._pwd+'/'+name)
          anode.setArrayValue(-1,shape(quantity),quantity,
                              self.type_dict[quantity.typecode()],-1)
          self.inquire_nodelist().addNode(anode)
        elif type(quantity) in [FloatType,IntType,StringType]:
          anode = _pyhl.node(_pyhl.ATTRIBUTE_ID,self._pwd+'/'+name)
          anode.setScalarValue(-1,quantity,
                               self.type_dict[type(quantity)],-1)
          self.inquire_nodelist().addNode(anode)
        else:
          # TODO: write as pickled object as a string
          raise "type unsupported"

    def defent (self, name, quantity, indx):
        """Define entry for quantity in file as 'name'"""
        self.check_open()
        if self.inquire_verbosity () > 1: 
            print "PW::defining entry for", name
        raise "defent not supported"


if __name__ == "__main__":
    f=PW("foo.pdb")
    a = 1
    b = 2.0
    c = "Hello world"
    from multiarray import *
    d = array ([1.,2., 3.])
    e = array ([1,2,3])
    g = array (["hello", "world", "array"])
    h = array ([[1.,2.,3.], [4.,5.,6]])
    k = 3
    f.a = a
    f.b = b
    f.c = c
    f.d = d
    f.e = e
    f.g = g
    f.h = h
    f.close ()
    f.open ("foo.pdb", "a")
    f.k = k
    f.close ()
# read-back test
    from PR import PR
    f = PR('foo.pdb')
    for x in f.inquire_names ():
        print x, "is", eval(x), ", in file it is", eval('f.'+x)
    f.close ()
# record-writing
    g = PW ('goo.pdb')
    g.set_verbosity (1)
    xh = array([0.]*4)
    for i in range (len(xh)):
        x = i / 10.
        g.write ('xh', x, i + 1)
        xh [i] = x
    g.close ()
    g = PR ('goo.pdb')
    print "xh is", xh, ", file it is ", g.xh
    g.close ()










