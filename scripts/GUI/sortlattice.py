from warp import *
import string

class ElemObj:
  def __init__(self,type,id):
    # --- Must refer directly to the dict since setattr does not allow new
    # --- attributes.
    self.__dict__['type'] = type
    self.__dict__['id'] = id
  def getpyobj(self,suffix):
    return top.getpyobject(string.lower(self.type)+suffix)
  def _getparam(self,suffix):
    topvar = self.getpyobj(suffix)
    return topvar[self.id]
  def _setparam(self,suffix,v):
    topvar = self.getpyobj(suffix)
    topvar[self.id] = v
  def __getattr__(self,name):
    if name == 'enabled': return self.getpyobj('s')
    if name in self.__dict__.keys(): return self.__dict__[name]
    try:
      v = self.getpyobj(name)
      return v[self.id]
    except:
      raise AttributeError
  def getattr(self,name):
    return self.__getattr__(name)
  def __setattr__(self,name,value):
    #if name in self.__dict__.keys(): self.__dict__[name] = value
    try:
      v = self.getpyobj(name)
      v[self.id] = value
    except:
      raise AttributeError
  def setattr(self,name,val):
    self.__setattr__(name,val)

def sortlattice():

  # --- First, get list of all of the existing elements in one nice place.
  elems = []
  if top.bends:
    for ii in xrange(top.nbend+1): elems.append(ElemObj("Bend",ii))
  if top.dipos:
    for ii in xrange(top.ndipo+1): elems.append(ElemObj("Dipo",ii))
  if top.quads:
    for ii in xrange(top.nquad+1): elems.append(ElemObj("Quad",ii))
  if top.sexts:
    for ii in xrange(top.nsext+1): elems.append(ElemObj("Sext",ii))
  if top.heles:
    for ii in xrange(top.nhele+1): elems.append(ElemObj("Hele",ii))
  if top.accls:
    for ii in xrange(top.naccl+1): elems.append(ElemObj("Accl",ii))
  if top.emlts:
    for ii in xrange(top.nemlt+1): elems.append(ElemObj("Emlt",ii))
  if top.mmlts:
    for ii in xrange(top.nmmlt+1): elems.append(ElemObj("Mmlt",ii))
  if top.bgrds:
    for ii in xrange(top.nbgrd+1): elems.append(ElemObj("Bgrd",ii))
  if top.pgrds:
    for ii in xrange(top.npgrd+1): elems.append(ElemObj("Pgrd",ii))
  if top.drfts:
    for ii in xrange(top.ndrft+1): elems.append(ElemObj("Drft",ii))

  elemszs = []
  for e in elems: elemszs.append(e.zs)
  elemszs = array(elemszs)
  sortindex = argsort(elemszs)
  sortedelems = []
  for i in sortindex: sortedelems.append(elems[i])

  return sortedelems

