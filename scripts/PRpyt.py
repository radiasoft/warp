# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

"""
HDF basic read-access class PR.py by Dave Grote, LLNL
"""
import pprint
import tables
import sys
import re

class PR:
    "HDF file read-access class."

    no_file_message = '(PR object not open on any file)'

    Error = 'PR error'

    def __del__(self):
        "Close any file open when this object disappears."
        self.close()

    def __init__(self, filename='', group='/', verbose = 1):
        """PR(filename='', group='', verbose=1)
        opens filename if given"""
        self.__dict__['file'] = None
        self.__dict__['tabledict'] = {}
        self.set_verbosity(verbose)
        self.set_group(group)
        self.__dict__['fixamp'] = re.compile('@')
        self.__dict__['amprepl'] = 'XXX'
        self.__dict__['fixXXX'] = re.compile('XXX')
        self.__dict__['XXXrepl'] = '@'
        if filename:
            self.open(filename)


    def __getattr__(self, name):
        return self.read(name)

    def __repr__(self):
        if self.is_open():
            current_mode = 'opened for reading'
        else:
            current_mode = PR.no_file_message
        return 'HDF file %s %s.' % (self.inquire_filename(), current_mode)

    __str__ = __repr__

    def check_open(self):
        "check_open(): raise exception if not open for read."
        if not self.is_open():
            raise PR.Error, 'PR object not open for read.'

    def close(self):
        "close(): close the file."
        if self.is_open():
            if self.inquire_verbosity():
                print "Closing HDF file",self.inquire_filename()
            self.inquire_file().close()
            d = self.__dict__
            v = self.inquire_verbosity()

            for n in d.keys():
                del d[n]

            d['file'] = None
            self.__dict__['tabledict'] = {}
            self.set_verbosity(v)

    def inquire_filename(self):
        "inquire_filename() = name of this file."
        if self.is_open ():
            return self.inquire_file().filename
        else:
            return ''

    def inquire_file(self):
        "inquire_file() = object open on this file."
        return self.file

    def inquire_high(self, name):
        "inquire_high(name) = high indices of name."
        v = self.__getattr__(name)
        return shape(v)

    def inquire_low(self, name):
        "inquire_low(name) = low indices of the name."
        v = self.__getattr__(name)
        return tuple(zeros(len(shape(v))))

    def inquire_ls(self,group=None):
        """inquire_ls(group=None) = list of those hdf names in current group or group"""
        self.check_open()
        if group is None:
          ll = self.inquire_file().listNodes(self.inquire_group(),classname='Leaf')
        else:
          ll = self.inquire_file().listNodes(group,classname='Leaf')
        ll = map(lambda l:l.name,ll)
        ll = ll + self.tabledict.keys()
        ll = map(lambda l:self.fixXXX.sub(self.XXXrepl,l),ll)
        return ll

    def inquire_mode(self):
        "inquire_mode() = mode ('r', 'w', or 'a') of this file."
        self.check_open()
        return self.inquire_file().mode

    def inquire_names(self):
        "inquire_names() = sorted list of names"
        r = self.inquire_ls()
        r.sort()
        return r

    def inquire_group(self):
        "inquire_group() = present HDF group"
        return self.__dict__['group']

    def inquire_shape(self,name):
        "inquire_shape(name) = shape of name."
        return self.inquire_high(name)

    def inquire_type(self, name):
        "inquire_type(name) = type of name."
        v = self.read(name)
        return type(v)

    def inquire_verbosity(self):
        "inquire_verbosity() = current value of verbose flag."
        return self._verbose

    def is_column_major(self):
        "is_column_major() = true if file was written by Fortran."
        self.check_open()
        return 0

    def is_open(self):
        "is_open() = true if file is open for read"
        if self.inquire_file() is None: return 0
        return self.inquire_file().isopen

    def print_names(self, file=sys.stdout):
        """print_list(file=sys.stdout): pretty print the list of names to file."""
        pprint.pprint(self.inquire_names(), file)

    def open(self, filename, group=None, mode='r'):
        """open(name, group=None): open file, optionally starting at group"""
        self.close()
        if group is None: group = self.inquire_group()
        self.__dict__['file'] = tables.openFile(filename,mode=mode,rootUEP=group)
        self.set_group('/')
        self.__dict__['ints'] = self.file.getNode(self.inquire_group(),'ints')
        self.readtable(self.ints)
        self.__dict__['floats'] = self.file.getNode(self.inquire_group(),'floats')
        self.readtable(self.floats)

    def readtable(self,table):
        for i in table.iterrows():
          self.tabledict[i['name']] = i['value']

    def read(self,name):
        "read(name) = the value of name as a Python object."
        self.check_open()
        name = self.fixamp.sub(self.amprepl,name)
        if self.tabledict.has_key(name): return self.tabledict[name]
        return self.inquire_file().getNode(self.inquire_group()+name).read()

    def read_part(self, name, triples):
        """read_part(name, triples) =  part of the value of name 
        as a Python object. Triples must be an array of 3 integers
        per subscript giving the low/high/stride values desired. 
        Note that these subscripts are relative to the values stored
        in the HDF file, and do not follow Python conventions."""
        assert ((len(triples)%3) == 0),"length of triples must be a multiple of three"
        v = self.__getattr__[name]
        ndims = len(triples)/3
        ii = '['
        for d in range(ndims):
          ii = ii + '%d:%d+1:%d'%tuple(triples[d*3:d*3+3])
          if d < ndims-1: ii = ii + ','
        ii = ii + ']'
        return eval('v'+ii,locals())

    def set_group(self, name):
        """set_group(name) 
        -- change HDF group to name, return status"""
        if name[0] == '/':
          group = name
        elif len(self.inquire_group()) == 1:
          group = '/' + name
        else:
          group = self.inquire_group() + '/' + name
        if group[-1] != '/': group = group + '/'
        self.__dict__['group'] = group

    def set_verbosity(self, flag):
        """verbose(flag) sets verbosity level to flag.
        0 for quiet operation, 
        1 to report closing only, 
        2 to report access to data."""
        if 0 <= flag <= 2:
            self.__dict__['_verbose'] = flag
        else:
            raise PR.Error, 'Illegal value for verbosity: '+`flag`

