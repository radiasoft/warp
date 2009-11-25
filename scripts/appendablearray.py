"""Array type which can be appended to in an efficient way.
"""
__all__ = ['AppendableArray','DynamicHistogram','DynamicHistogramIntersect']
import sys
import numpy
from toppy import deposeintersect,setgrid1d,deposgrid1d
# Class which allows an appendable array.
# DPG 8/19/99
appendablearray_version = "$Id: appendablearray.py,v 1.17 2009/11/25 22:28:25 jlvay Exp $"

class AppendableArray:
  """
Creates an array which can be appended to in an efficient manner. The object
keeps an internal array which is bigger than the actual data. When new data is
appended, it is only copied to fill in the extra space. More space is only
allocated when that space fills up. This saves both the allocation time, and
the time to copy the existing data to the new space.
 - initlen=1: The initial size of the array. The most efficiency is gained
              when initlen is close to the total expected length.
 - unitshape=None: The appendable unit can be an array. This gives the shape of
               the unit. The full shape of the array then would be
               [n]+unitshape, where n is the number of units appended.
 - typecode='i': Typecode of the array. Uses the same default as the standard
                 array creation routines.
 - autobump=100: The size of the increment used when additional extra space
                 is needed.
 - initunit=None: When given, the unitshape and the typecode are taken from
                  it. Also, this unit is make the first unit in the array.
 - aggressivebumping=1.: Whenever new space is added, autobump will be
                         increased by this factor in an attempt to
                         minimize the number of reallocations. Its new
                         value will be the max of the size of the appended data
                         and this times the old autobump size. A good value is
                         1.5 - this can greatly reduce the amount of
                         rallocation withoug a significant amount of wasted
                         space.

Create an instance like so
>>> a = AppendableArray(initlen=100,typecode='d')

Append a single unit like this
>>> a.append(7.)

or multiple units like this
>>> a.append(numpy.ones(5,'d'))

The data can be obtained by directly indexing the instance, like this
>>> print a[:4]
[ 7., 1., 1., 1.,]
will give the first four number appended

Other methods include len, data, setautobump, cleardata, reshape
  """
  def __init__(self,initlen=1,unitshape=None,typecode=None,autobump=100,
               initunit=None,aggressivebumping=0):
    if typecode is None: typecode = numpy.zeros(1).dtype.char
    self._maxlen = initlen
    if initunit is None:
      self._typecode = typecode
      self._unitshape = unitshape
    else:
      # --- Get typecode and unitshape from initunit
      if isinstance(initunit,numpy.ndarray):
        self._typecode = initunit.dtype.char
        self._unitshape = initunit.shape
      else:
        if type(initunit) == IntType: self._typecode = 'i'
        else:                         self._typecode = 'd'
        self._unitshape = None
    self._datalen = 0
    self._autobump = autobump
    self.aggressivebumping = aggressivebumping
    self._allocatearray()
    if initunit is not None: self.append(initunit)
  def checkautobumpsize(self,deltalen):
    # --- The factor of 1.5 gives nearly the same amount of savings as 2,
    # --- but doesn't waste quite as much space.
    newautobump = max(deltalen,int(1.5*self.getautobump()))
    self.setautobump(newautobump)
  def _extend(self,deltalen):
    # --- Only increase of the size of the array if the extra space fills up
    if len(self) + deltalen > self._maxlen:
      self.checkautobumpsize(deltalen)
      self._maxlen = self._maxlen + max(deltalen,self._autobump)
      a = self._array[:len(self),...] + 0
      self._allocatearray()
      self._array[:len(self),...] = a
  def _allocatearray(self):
    if self._unitshape is None:
      self._array = numpy.zeros(self._maxlen,self._typecode)
    else:
      self._array = numpy.zeros([self._maxlen]+list(self._unitshape),self._typecode)
  def append(self,data):
    if self._unitshape is None:
      # --- If data is just a scalar, then set length to one. Otherwise
      # --- get length of data to add.
      try:
        lendata = len(data)
      except (TypeError,IndexError):
        lendata = 1
    else:
      # --- Data must be an array in this case.
      # --- If the shape of data is the same as the original shape,
      # --- then only one unit is added. Otherwise, get the number
      # --- of units to add. The length is always added to the first
      # --- dimension.
      if len(data.shape) == len(self._unitshape): lendata = 1
      else:                                       lendata = data.shape[0]
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
  def getautobump(self):
    """
Get the autobump attribute
    """
    return self._autobump
  def cleardata(self):
    """
Reset the array so it has a length of zero.
    """
    self._datalen = 0
  def resetdata(self,data):
    """
Resets the data to be the input values - all of the original data is thrown
away. The unit shape of data must be the same.
    """
    self.cleardata()
    self.append(data)
  def unitshape(self):
    if self._unitshape is None:
      return (1,)
    else:
      return self._unitshape
  def reshape(self,newunitshape):
    """
Change the shape of the appendable unit. Can only be used if a unitshape was
specified on creation.
 - newunitshape: must have the same number of dimensions as the original
                 unitshape
    """
    assert self._unitshape is not None,\
           'Only an array with a specified unitshape can be reshaped'
    assert len(newunitshape) == len(self._unitshape),\
           ('New unitshape must have the same number of dimensions as original'+
            'unitshape')
    # --- Save old data
    oldunitshape = self._unitshape
    oldarray = self._array
    # --- Create new array
    self._unitshape = newunitshape
    self._allocatearray()
    # --- Copy data from old to new
    ii = [None] + list(numpy.minimum(oldunitshape,newunitshape))
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

class DynamicHistogram:
  """
  """
  def __init__(self,n=100,min=None,max=None,overfrac=0.2):
    self.n=n
    self.min = min
    self.max = max
    self.overfrac=overfrac
    self.data=numpy.zeros(n,'d')
    self.bins=None
    self.binsdx=None
    self.binsmin=None
    self.binsmax=None
    
  def checkbounds(self,d):
    dmin = min(d)
    dmax = max(d)
    if self.min is None: self.min=dmin
    if self.max is None: self.max=dmax
    if self.bins is None:
      self.binsdx = (self.max-self.min)/(self.n-1)
      self.bins = self.min+numpy.arange(self.n)*self.binsdx
      self.binsminmax = numpy.zeros(self.n*2)
      self.binsminmax[0::2] = self.min+numpy.arange(self.n)*self.binsdx
      self.binsminmax[1::2] = self.min+numpy.arange(self.n)*self.binsdx+self.binsdx
    l_rescale_array=0
    if dmin<self.min:
      minold = self.min
      self.min=dmin-self.overfrac*(self.max-self.min)
      l_rescale_array=1
    if dmax>self.max:
      maxold = self.max
      self.max=dmax+self.overfrac*(self.max-self.min)
      l_rescale_array=1
    if l_rescale_array:
      newdata  = numpy.zeros(self.n,'d')
      tmpcount = numpy.zeros(self.n,'d')
      deposgrid1d(1,self.n,self.bins,self.data,self.n-1,newdata,tmpcount,self.min,self.max)
      self.binsdx = (self.max-self.min)/(self.n-1)
      self.bins = self.min+numpy.arange(self.n)*self.binsdx
      self.binsminmax = numpy.zeros(self.n*2)
      self.binsminmax[0::2] = self.min+numpy.arange(self.n)*self.binsdx
      self.binsminmax[1::2] = self.min+numpy.arange(self.n)*self.binsdx+self.binsdx
      self.data = newdata

  def accumulate(self,d,weights=None):
    if type(d) is type(0.):d=array([d])
    self.checkbounds(d)
    if weights is None:     
      setgrid1d(numpy.shape(d)[0],d,self.n-1,self.data,self.min,self.max)
    else:
      setgrid1dw(numpy.shape(d)[0],d,weights,self.n-1,self.data,self.min,self.max)

class DynamicHistogramIntersect(DynamicHistogram):
  """
  """
  def __init__(self,n=100,min=None,max=None,overfrac=0.2):
     DynamicHistogram.__init__(self,n,min,max,overfrac)
     
  def accumulate(self,z1,d1,ssn1,w1,z2,d2,ssn2,w2,z0):
     if self.min is None:self.min=1.e36
     if self.max is None:self.max=1.e36
     minold = self.min
     maxold = self.max
     i1 = numpy.argsort(ssn1)
     i2 = numpy.argsort(ssn2)
     d1 = numpy.take(d1,i1)
     d2 = numpy.take(d2,i2)
     z1 = numpy.take(z1,i1)
     z2 = numpy.take(z2,i2)
     ssn1 = numpy.take(ssn1,i1)
     ssn2 = numpy.take(ssn2,i2)
     n1 = len(i1)
     n2 = len(i2)
     fminmax = numpy.array([self.min,self.max])
     deposeintersect(z1,d1,ssn1,w1,n1,z2,d2,ssn2,w2,n2,z0,self.data,fminmax,self.n-1,True,self.overfrac)
     self.min = fminmax[0]
     self.max = fminmax[1]
     if self.min<>minold or self.max<>maxold:
       self.bins = self.min+numpy.arange(self.n)*(self.max-self.min)/(self.n-1)

  def getave(self):
    return numpy.sum(self.data[:]*self.bins[:])/numpy.sum(self.data[:])
