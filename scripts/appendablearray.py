from warp import *
# Class which allows an appendable array.
# DPG 8/19/99
appendablearray_version = "$Id: appendablearray.py,v 1.4 2002/01/29 22:29:38 dave Exp $"

class AppendableArray:
  def __init__(self,initlen,shape=None,type='i',autobump=100):
    self._maxlen = initlen
    self._type = type
    self._shape = shape
    self._datalen = 0
    self._autobump = autobump
    self.allocatearray()
  def append(self,data):
    if self._shape is None:
      try:
        lendata = shape(data)[0]
      except TypeError:
        lendata = 1
    else:
      if len(shape(data)) == len(self._shape): lendata = 1
      else:                                    lendata = shape(data)[-1]
    self.extend(lendata)
    newlen = self._datalen + lendata
    self._array[self._datalen:newlen,...] = data
    self._datalen = newlen
  def data(self):
    return self._array[:len(self),...]
  def setautobump(self,a):
    self._autobump = a
  def cleardata(self):
    self._datalen = 0
  def extend(self,deltalen):
    if len(self) + deltalen > self._maxlen:
      self._maxlen = self._maxlen + max(deltalen,self._autobump)
      a = self._array[:len(self),...] + 0
      self.allocatearray()
      self._array[:len(self),...] = a
  def allocatearray(self):
    if self._shape is None:
      self._array = zeros(self._maxlen,self._type)
    else:
      self._array = zeros([self._maxlen]+list(self._shape),self._type)
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

