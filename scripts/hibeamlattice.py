# Setup classes to mirror lattice element type in hibeam.
# DPG 7/22/99
import re
import string
hibeamlattice_version = "$Id: hibeamlattice.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

class LINE:
  def __init__(self,*elems):
    self.type = 'LINE'
    # --- Unravel any imbedded lists.
    i = 0
    elems = list(elems)
    while i < len(elems):
      if type(elems[i]) == type([]):
        elems[i:i+1] = elems[i]
      else:
        i = i + 1
    # --- save list of elements
    self.elems = elems
    # --- initialize expanded list as a blank
    self.elemslist = []
  def expand(self):
    if self.elemslist: return self.elemslist
    for e in self.elems:
      self.elemslist.append(e.expand())
    # --- Unravel any imbedded lists.
    i = 0
    self.elemslist = list(self.elemslist)
    while i < len(self.elemslist):
      if type(self.elemslist[i]) == type([]):
        self.elemslist[i:i+1] = self.elemslist[i]
      else:
        i = i + 1
    return self.elemslist
  def __mul__(self,other):
    # --- This allows multiplication of elements
    return LINE(other*[self])
  def __rmul__(self,other):
    # --- This allows multiplication of elements
    return LINE(other*[self])

class Elem:
  def __init__(self,l=0,length=0):
    self.type = ''
    self.length = max(l,length)
  #def isin(self,z=top.zbeam):
    #if self._zs < z and z < self._ze: return 1
    #return 0
  def length(self):
    return self._length
  def expand(self):
    return self
  def __mul__(self,other):
    return LINE(other*[self])
  def __rmul__(self,other):
    return LINE(other*[self])

class drift(Elem):
  def __init__(self,l=0,length=0,aperture=0,error_type='',
                    offset_x=0,offset_y=0,i_cap_pointer=0,n_cap_nodes=0):
    self.type = 'drift'
    Elem.__init__(self,l=l,length=length)
    self.aperture = aperture
    self.error_type = error_type
    self.offset_x = offset_x
    self.offset_y = offset_y
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes

class box(Elem):
  def __init__(self,l=0,length=0,width_x=0,width_y=0,error_type='',
                    offset_x=0,offset_y=0,i_cap_pointer=0,n_cap_nodes=0):
    self.type = 'box'
    Elem.__init__(self,l=l,length=length)
    self.width_x = width_x
    self.width_y = width_y
    self.error_type = error_type
    self.offset_x = offset_x
    self.offset_y = offset_y
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes

class quad(Elem):
  def __init__(self,l=0,length=0,aperture=0,voltage=0,gradient=0,r_elem=0,
                    width_x=0,width_y=0,error_type='',
                    offset_x=0,offset_y=0,i_cap_pointer=0,n_cap_nodes=0):
    self.type = 'quad'
    Elem.__init__(self,l=l,length=length)
    self.aperture = aperture
    self.voltage = voltage
    self.gradient = gradient
    self.r_elem = r_elem
    self.error_type = error_type
    self.offset_x = offset_x
    self.offset_y = offset_y
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes
    if gradient != 0:
      self.voltage = gradient*aperture**2
    else:
      self.gradient = voltage/aperture**2

class hyperb(Elem):
  def __init__(self,l=0,length=0,aperture=0,voltage=0,r_elem=0,error_type='',
                    i_cap_pointer=0,n_cap_nodes=0):
    self.type = 'hyperb'
    Elem.__init__(self,l=l,length=length)
    self.aperture = aperture
    self.voltage = voltage
    self.r_elem = r_elem
    self.error_type = error_type
    self.i_cap_pointer = i_cap_pointer
    self.n_cap_nodes = n_cap_nodes

class wire(Elem):
  def __init__(self,l=0,length=0,aperture=0,):
    self.type = 'wire'
    Elem.__init__(self,l=l,length=length)
    self.aperture = aperture

class child:
  def __init__(self,parent,**changes):
    self.parent = eval(parent,globals())
    self.__dict__.update(changes)
    if self.__dict__.has_key('l'):
      self.__dict__['length'] = self.__dict__['l']
  def __getattr__(self,name):
    return eval('self.parent.'+name,locals())


#########################################################################
# --- Now, using above classes, read in and parse a hibeam lattice file.
def hibeamlattice(file):
  ff = open(file,'r')
  data = ff.readlines()
  ff.close()
  # --- Massage the data removing carriage returns, comments and blank lines
  i = 0
  while i < len(data):
    # --- Remove carriage return at end of line
    data[i] = data[i][:-1]
    # --- Remove all white space
    data[i] = re.sub('\s','',data[i])
    # --- Remove comments
    data[i] = re.sub('!.*','',data[i])
    # --- Delete empty lines
    if data[i]:
      i = i + 1
    else:
      del data[i]
  # --- Massage the data removing line continuations.
  i = 0
  while i < len(data):
    mextend = re.search('&',data[i])
    while mextend:
      data[i] = re.sub('&',data[i+1],data[i])
      del data[i+1]
      mextend = re.search('&',data[i])
    i = i + 1
  # --- Massage the data into python syntax
  for i in xrange(len(data)):
    # --- Replace '=' with '(' in LINE command
    # --- Replace first ',' with '(' in all except LINE commands
    mline = re.search('LINE',data[i])
    if mline:
      data[i] = re.sub('=','(',data[i])
    else:
      data[i] = re.sub(',','(',data[i],1)
    # --- Replace ':' with '=' in all commands
    data[i] = re.sub(':','=',data[i])
    # --- Append ')' to the end of the line
    data[i] = data[i] + ')'
  # --- Now the data is ready to be processed
  for d in data:
    mfinish = re.search(r'(?P<u>USE\((?P<use>\w*)\))|(?P<end>END)',d)
    if mfinish:
      if mfinish.group('u'):
        return eval(mfinish.group('use'),globals())
    else:
      exec(d,globals())

