# Creates a class representation of functions that allow
# the function to be called without the added '()' when there are no
# arguments. Note that in that case, an extra blank line is output.
noparens_version = "$Id: noparens.py,v 1.2 2001/01/11 00:21:48 dave Exp $"


class Noparens:
  def __init__(self,func):
    self.func = func
  def __repr__(self):
    self.func()
    return ''
  __str__ = __repr__
  def __call__(self,*k,**kw):
    return apply(self.func,k,kw)
