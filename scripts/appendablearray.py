from warp import *
# Class which allows an appendable array.
# DPG 8/19/99
appendablearray_version = "$Id: appendablearray.py,v 1.3 2001/07/19 21:32:22 dave Exp $"

class AppendableArray:
  def __init__(self,initlen,type='i',autobump=100):
    self._maxlen = initlen
    self._type = type
    self._array = zeros(self._maxlen,self._type)
    self._datalen = 0
    self._autobump = autobump
  def append(self,data):
    try:
      lendata = len(data)
    except TypeError:
      lendata = 1
    self.extend(lendata)
    newlen = self._datalen + lendata
    self._array[self._datalen:newlen] = data
    self._datalen = newlen
  def data(self):
    return self._array[:len(self)]
  def setautobump(self,a):
    self._autobump = a
  def cleardata(self):
    self._datalen = 0
  def extend(self,deltalen):
    if len(self) + deltalen > self._maxlen:
      self._maxlen = self._maxlen + max(deltalen,self._autobump)
      a = zeros(self._maxlen,self._type)
      a[:len(self)] = self.data()
      self._array = a
  def __len__(self):
    return self._datalen
  def __getitem__(self,key):
    return self.data()[key]
  def __setitem__(self,key,value):
    self.data()[key] = value

if sys.version < "2.0":
  def _appendablearray__getslice__(self,i,j):
    return self.data()[i:j]
  AppendableArray.__getslice__ = _appendablearray__getslice__

