"""Array type which can be appended to in an efficient way.
"""
from warp import *
# Class which allows an appendable array.
# DPG 8/19/99
appendablearray_version = "$Id: appendablearray.py,v 1.8 2005/02/15 19:55:35 dave Exp $"

class AppendableArray:
  """
Creates an array which can be appended to in an efficient manner. The object
keeps an internal array which is bigger than the actual data. When new data is
appended, it is only copied to fill in the extra space. More space is only
allocated when that space fills up. This saves both the allocation time, and
the time to copy the existing data to the new space.
 - initlen=1: The initial size of the array. The most efficiency is gained
              when initlen is close to the total expected length.
 - shape=None: The appendable unit can be an array. This gives the shape of
               the unit. The full shape of the array then would be
               [n]+shape, where n is the number of units appended.
 - type='i': Typecode of the array. Uses the same default as the standard
             array creation routines.
 - autobump=100: The size of the increment used when additional extra space
                 is needed.

Create an instance like so
>>> a = AppendableArray(initlen=100,type='d')

Append a single unit like this
>>> a.append(7.)

or multiple units like this
>>> a.append(ones(5,'d'))

The data can be obtained by directly indexing the instance, like this
>>> print a[:4]
[ 7., 1., 1., 1.,]
will give the first four number appended

Other methods include len, data, setautobump, cleardata, reshape
  """
  def __init__(self,initlen=1,shape=None,type='i',autobump=100):
    self._maxlen = initlen
    self._type = type
    self._shape = shape
    self._datalen = 0
    self._autobump = autobump
    self._allocatearray()
  def _extend(self,deltalen):
    # --- Only increase of the size of the array if the extra space fills up
    if len(self) + deltalen > self._maxlen:
      self._maxlen = self._maxlen + max(deltalen,self._autobump)
      a = self._array[:len(self),...] + 0
      self._allocatearray()
      self._array[:len(self),...] = a
  def _allocatearray(self):
    if self._shape is None:
      self._array = zeros(self._maxlen,self._type)
    else:
      self._array = zeros([self._maxlen]+list(self._shape),self._type)
  def append(self,data):
    if self._shape is None:
      # --- If data is just a scalar, then set length to one. Otherwise
      # --- get length of data to add.
      try:
        lendata = shape(data)[0]
      except (TypeError,IndexError):
        lendata = 1
    else:
      # --- Data must be an array in this case.
      # --- If the shape of data is the same as the original shape,
      # --- then only one unit is added. Otherwise, get the number
      # --- of units to add. The length is always added to the first
      # --- dimension.
      if len(shape(data)) == len(self._shape): lendata = 1
      else:                                    lendata = shape(data)[0]
    self._extend(lendata)
    newlen = self._datalen + lendata
    self._array[self._datalen:newlen,...] = data
    self._datalen = newlen
  def data(self):
    """
Return the data.
    """
    return self._array[:len(self),...]
  def setautobump(self,a):
    """
Set the autobump attribute to the value specified.
    """
    self._autobump = a
  def cleardata(self):
    """
Reset the array so it has a length of zero.
    """
    self._datalen = 0
  def reshape(self,newshape):
    """
Change the shape of the appendable unit. Can only be used if a shape was
specified on creation.
 - newshape: must have the same number of dimensions as the original shape
    """
    assert self._shape is not None,\
           'Only an array with a specified shape can be reshaped'
    assert len(newshape) == len(self._shape),\
           'New shape must have the same number of dimensions as original shape'
    # --- Save old data
    oldshape = self._shape
    oldarray = self._array
    # --- Create new array
    self._shape = newshape
    self._allocatearray()
    # --- Copy data from old to new
    ii = [None] + list(minimum(oldshape,newshape))
    ss = map(slice,ii)
    self._array[ss] = oldarray[ss]

  def __len__(self):
    return self._datalen
  def __getitem__(self,key):
    return self.data()[key]
  def __setitem__(self,key,value):
    self.data()[key] = value

if sys.version < "2.0":
  def _appendablearray__getslice__(self,i,j):
    return self.data()[i:j,...]
  AppendableArray.__getslice__ = _appendablearray__getslice__

