# Creates a class representation of functions that allow
# the function to be called without the added '()' when there are no
# arguments. Note that in that case, an extra blank line is output.
noparens_version = "$Id: noparens.py,v 1.1 2000/10/16 18:34:19 dave Exp $"


class Noparens:
  def __init__(self,func):
    self.func = func
  def __repr__(self):
    self.func()
    return ''
  __str__ = __repr__
  def __call__(self,*k,**kw):
    apply(self.func,k,kw)
