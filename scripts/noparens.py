# Creates a class representation of functions that allow
# the function to be called without the added '()' when there are no
# arguments. Note that in that case, an extra blank line is output.
noparens_version = "$Id: noparens.py,v 1.3 2007/06/21 16:47:36 dave Exp $"


class Noparens:
  def __init__(self,func):
    self.func = func
    self.__doc__ = func.__doc__
  def __repr__(self):
    self.func()
    return ''
  __str__ = __repr__
  def __call__(self,*k,**kw):
    return apply(self.func,k,kw)
