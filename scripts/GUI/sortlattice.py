from warp import *

class ElemObj:
  def __init__(self,type,id):
    # --- Must refer directly to the dict since setattr does not allow new
    # --- attributes.
    self.__dict__['type'] = type
    self.__dict__['id'] = id
  def _getparam(self,suffix):
    topvar = top.getpyobject(self.type+suffix)
    return topvar[self.id]
  def _setparam(self,suffix,v):
    topvar = top.getpyobject(self.type+suffix)
    topvar[self.id] = v
  def __getattr__(self,name):
    if name == 'enabled': return top.getpyobject(self.type+'s')
    if name in self.__dict__.keys(): return self.__dict__[name]
    try:
      v = top.getpyobject(self.type+name)
      return v[self.id]
    except:
      raise AttributeError
  def __setattr__(self,name,value):
    #if name in self.__dict__.keys(): self.__dict__[name] = value
    try:
      v = top.getpyobject(self.type+name)
      v[self.id] = value
    except:
      raise AttributeError

def sortlattice():

  # --- First, get list of all of the existing elements in one nice place.
  elems = []
  if top.bends:
    for ii in xrange(top.nbend+1): elems.append(ElemObj("bend",ii))
  if top.dipos:
    for ii in xrange(top.ndipo+1): elems.append(ElemObj("dipo",ii))
  if top.quads:
    for ii in xrange(top.nquad+1): elems.append(ElemObj("quad",ii))
  if top.sexts:
    for ii in xrange(top.nsext+1): elems.append(ElemObj("sext",ii))
  if top.heles:
    for ii in xrange(top.nhele+1): elems.append(ElemObj("hele",ii))
  if top.accls:
    for ii in xrange(top.naccl+1): elems.append(ElemObj("accl",ii))
  if top.emlts:
    for ii in xrange(top.nemlt+1): elems.append(ElemObj("emlt",ii))
  if top.mmlts:
    for ii in xrange(top.nmmlt+1): elems.append(ElemObj("mmlt",ii))
  if top.bgrds:
    for ii in xrange(top.nbgrd+1): elems.append(ElemObj("bgrd",ii))
  if top.pgrds:
    for ii in xrange(top.npgrd+1): elems.append(ElemObj("pgrd",ii))
  if top.drfts:
    for ii in xrange(top.ndrft+1): elems.append(ElemObj("drft",ii))

  elemszs = []
  for e in elems: elemszs.append(e.zs)
  elemszs = array(elemszs)
  sortindex = argsort(elemszs)
  sortedelems = []
  for i in sortindex: sortedelems.append(elems[i])

  return sortedelems

