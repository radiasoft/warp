from warp import *
# Class which allows an appendable array.
# DPG 8/19/99
appendablearray_version = "$Id: appendablearray.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

class AppendableArray:
  def __init__(self,initlen,type='i',autobump=100):
    self._maxlen = initlen
    self._type = type
    self._array = zeros(maxlen,type)
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
    return self._array[:self._datalen]
  def setautobump(self,a):
    self._autobump = a
  def len(self):
    return self._datalen
  def cleardata(self):
    self._datalen = 0
  def extend(self,deltalen):
    if self._datalen + deltalen > self._maxlen:
      self._maxlen = self._maxlen + max(deltalen,self._autobump)
      a = zeros(self._maxlen,self._type)
      a[:self.len()] = self.data()
      self._array = a
