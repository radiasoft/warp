"""
This module contains convenience functions making use of the python wrapper of
OpenDX.
ppxxpy: draws a plot of x-xp-y with coloring based on yp
ppxyz: draws a plot of x-y-z with coloring based on xp
viewisosurface1: example of a isosurface plot
viewisosurface: example of a isosurface plot with a plane of data
viewparticles: views particles
DXImage: views a OpenDX object using a nice interactor handler
"""
from warp import *
from wxPython.wx import *
import __main__
selectbox = 1

try:
  # --- Try importing. If not found, still define the functions and classes.
  # --- But any calls require the pyOpenDX attributes will raise an exception.
  from pyDXObject import *
except:
  pass

pyOpenDX_version = "$Id: pyOpenDX.py,v 1.16 2004/11/24 22:52:31 jlvay Exp $"
def pyOpenDXdoc():
  import pyOpenDX
  print pyOpenDX.__doc__

###########################################################################
def ppxxpy(iw = 0,labels=1,display=1,**kw):
  """Plots X-Y-Z"""
  checkparticleplotarguments(kw)
  ii = selectparticles(iw=iw,kwdict=kw)
  if labels == 1: labels = ['X',"X'",'Y']
  return viewparticles(take(top.xp,ii),(take(top.uxp,ii)/take(top.uzp,ii)),
                       take(top.yp,ii),(take(top.uyp,ii)/take(top.uzp,ii)),
                       labels,name='WARP viz',display=display)

def ppxyz(iw = 0,cc=None,labels=1,display=1,rscale=None,zscale=None,scale=1.,ratio=1.,**kw):
  """Plots X-Y-Z"""
  checkparticleplotarguments(kw)
  ii = selectparticles(iw=iw,kwdict=kw)
  if labels == 1: labels = ['X','Y','Z']
  xx = take(top.xp,ii)
  yy = take(top.yp,ii)
  zz = take(top.zp,ii)
  if cc is None: cc = top.uxp
  cc = take(cc,ii)
  if rscale is not None:
    xx = xx*rscale
    yy = yy*rscale
  if zscale is not None:
    zz = zz*zscale
  return viewparticles(xx,yy,zz,cc,
                       labels,name='WARP viz',display=display,scale=scale,ratio=ratio)

###########################################################################
def viewisosurface1(data,isovalue,origins=None,deltas=None,name='WARP viz'):
  if origins is None: origins = [0.,0.,0.]
  if deltas is None: deltas = [1.,1.,1.]
  f = DXObject_fromarray(data,origins,deltas)

  minput = {'data':f,'value':isovalue}
  moutput = ['surface']
  (isosurface,) = DXCallModule('Isosurface',minput,moutput)

  DXReference(isosurface)
  minput = {'object':isosurface}
  moutput = ['camera']
  (camera,) = DXCallModule('AutoCamera',minput,moutput)

  DXImage(isosurface,camera,name)

def viewisosurface(data,isovalue,name='WARP viz'):
  f = DXObject_fromarray(data)
  DXReference(f)

  minput = {'data':f,'value':isovalue}
  moutput = ['surface']
  (isosurface,) = DXCallModule('Isosurface',minput,moutput)

  DXReference(isosurface)
  minput = {'object':isosurface}
  moutput = ['camera']
  (camera,) = DXCallModule('AutoCamera',minput,moutput)

  minput = {'input':f,'dimension':2,'position':128}
  moutput = ['output']
  (slab,) = DXCallModule('Slab',minput,moutput)

  minput = {'data':slab}
  moutput = ['mapped']
  (slabcolor,) = DXCallModule('AutoColor',minput,moutput)

  minput = {'object':isosurface,'object1':slabcolor}
  moutput = ['group']
  (group,) = DXCallModule('Collect',minput,moutput)

  DXImage(group,camera,name)

###########################################################################
def viewparticles(x,y,z,v,labels=None,name='WARP viz',
                  display=1,scale=1.,ratio=1.):
  x = gatherarray(x)
  y = gatherarray(y)
  z = gatherarray(z)
  v = gatherarray(v)
  if me > 0: return
  # --- First combine particle data and create a DX array
  n = len(x)
  p = zeros((n,3),'d')
  p[:,0] = x
  p[:,1] = y
  p[:,2] = z
  dxp = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
  DXAddArrayData(dxp,0,n,p.astype(Float32))
  # --- Create a DX array for the data
  dxd = DXNewArray(TYPE_DOUBLE,CATEGORY_REAL,0)
  DXAddArrayData(dxd,0,n,v)
  DXSetStringAttribute(dxd,'dep','positions')
  # --- Create the field
  dxf = DXNewField()
  DXSetComponentValue(dxf,'positions',dxp)
  DXSetComponentValue(dxf,'data',dxd)
  DXEndField(dxf)

  # --- Create glyphs to render data
  minput = {'data':dxf,'type':'speedy','ratio':ratio,'scale':scale}
  moutput = ['glyphs']
  (glyphs,) = DXCallModule('AutoGlyph',minput,moutput)

  minput = {'data':glyphs,'opacity':.5}
  moutput = ['mapped']
  (dxobject,) = DXCallModule('AutoColor',minput,moutput)

  if display:
    DXImage(dxobject,name=name,labels=labels)
  else:
    return dxobject

###########################################################################
def viewboundingbox(xmin,xmax,ymin,ymax,zmin,zmax,color='yellow'):
  """
Create a box
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of the box
  """
  # --- Create a field containing two points on opposite side of the
  # --- mesh and use it to find the boinding box.
  n = 2
  p = zeros((n,3),'d')
  p[:,0] = [xmin,xmax]
  p[:,1] = [ymin,ymax]
  p[:,2] = [zmin,zmax]
  dxp = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
  DXAddArrayData(dxp,0,n,p.astype(Float32))
  # --- Create a DX array for the data
  dxd = DXNewArray(TYPE_DOUBLE,CATEGORY_REAL,0)
  DXAddArrayData(dxd,0,n,array([1.,1.]))
  DXSetStringAttribute(dxd,'dep','positions')
  # --- Create the field
  dxf = DXNewField()
  DXSetComponentValue(dxf,'positions',dxp)
  DXSetComponentValue(dxf,'data',dxd)
  DXEndField(dxf)

  minput = {'input':dxf}
  moutput = ['box']
  (box,) = DXCallModule('ShowBox',minput,moutput)
  minput = {'input':box,'color':color}
  moutput = ['colored']
  (box,) = DXCallModule('Color',minput,moutput)
  return box

###########################################################################
class Visualizable:
  """
Virtual class which defines the API for objects that can be visualized.
Note that any class inheriting this class must define the createdxobject
method.
  """
  def createdxobject(self,kwdict={},**kw):
    """
Creates the dxobject. Note that this method must define the dxobject
attribute.
    """
    raise 'createdxobject is a virtual function'
  def getdxobject(self,kwdict={},**kw):
    """
Returns the dxobject, creating it is necessary. Note that the local
reference to the dxobject is deleted.
    """
    try:
      self.dxobject
    except AttributeError:
      kw.update(kwdict)
      self.createdxobject(kwdict=kw)
    dxobject = self.dxobject
    del self.dxobject
    return dxobject
  def visualize(self,kwdict={},**kw):
    """
Displays the dxobject, creating it if necessary.
    """
    kw.update(kwdict)
    dxobject = self.getdxobject(kwdict=kw)
    DXImage(dxobject)

###########################################################################
class DXCollection(Visualizable):
  global selectbox
  def __init__(self,*dxobjects,**kw):
    self.labels = kw.get('labels',None)
    self.preserve = kw.get('preserve',0)
    self.collection = None
    for o in dxobjects:
      if __main__.selectbox:
        try:
          o.SelectBox()
        except:
          pass
      self.addobject(o)
  def createdxobject(self,kwdict={},**kw):
    assert self.collection is not None,"Empty collection can not be visualized"
    self.dxobject = self.collection
    if not self.preserve: self.collection = None
  def extractdxobject(self,dxobject):
    # --- This is somewhat kludgy, but works. This keeps looking until
    # --- it finds a actual dxobject.
    # --- Assumes that that will be an actual dxobject. This will fail
    # --- miserably if an objects dxobject reference is recursive.
    # --- Please, let that never happen!
    while not isDXObject(dxobject) and dxobject is not None:
      dxobject = dxobject.getdxobject()
    return dxobject
  def initializedxobject(self,object):
    dxobject = self.extractdxobject(object)
    if dxobject is None: return
    minput = {'object':dxobject}
    moutput = ['group']
    (self.collection,) = DXCallModule('Collect',minput,moutput)
    if self.preserve: DXReference(self.collection)
    self.len = 1
  def checkmaxlegth(self):
    if self.len == 21:
      if self.preserve: oldcollection = self.collection
      self.initializedxobject(self.collection)
      if self.preserve: DXDelete(oldcollection)
  def addobject(self,object):
    if self.collection is None:
      self.initializedxobject(object)
    else:
      dxobject = self.extractdxobject(object)
      if dxobject is None: return
      self.checkmaxlegth()
      if self.preserve: oldcollection = self.collection
      minput = {'input':self.collection,'object':dxobject}
      moutput = ['group']
      (self.collection,) = DXCallModule('Append',minput,moutput)
      if self.preserve:
        DXDelete(oldcollection)
        DXReference(self.collection)
  def reset(self):
    if self.collection is None: return
    if self.preserve: DXDelete(self.collection)
    self.collection = None
  def __add__(self,right):
    return DXCollection(self,right)

###########################################################################
#==========================================================================
interactor = 0
def interactor_handler():
  global interactor
  getchar = sys.stdin.readline()[:-1]
  if len(getchar) > 0: interactor = eval(getchar)
  else:                interactor = -1

def DXImage(object,camera=None,name='WARP viz',labels=None):
  """
Displays an image of the input object, allowing mouse controls for moving the
image. Default mode is rotation. Press 1 for panning, 2 for zooming.
 - object: object to be imaged. Note that it will be deleted afterward
           (to be consistent with the OpenDX standard).
 - camera: optional camera used to view the image. Default uses AutoCamera.
 - name='WARP viz': Title given to X window. When given a new name, a new X
                    windows will be opened.
 - labels=None: when specified, axis with labels are included. labels must be
                of list containing three strings.
  """
  global interactor

  group = DXCollection(object)
  dxobject = group.getdxobject()
  assert isDXObject(dxobject),'Object to display is not a DXObject'

  # --- Turn off the caching of images - caching can use a lot of memory
  # --- which is not freed.
  minput = {'input':dxobject,'attribute':'cache','value':0}
  moutput = ['output']
  (dxobject,) = DXCallModule('Options',minput,moutput)

  if camera is None:
    DXReference(dxobject)
    minput = {'object':dxobject}
    moutput = ['camera']
    (camera,) = DXCallModule('AutoCamera',minput,moutput)

  if labels is not None:
    assert (len(labels) == 3),"Length of lables list must be three"
    DXReference(camera)
    labels = DXMakeStringList(labels)
    minput = {'input':dxobject,'camera':camera,'labels':labels,
              'ticks':3,'colors':'grey'}#'yellow'}
    moutput = ['axes']
    (dxobject,) = DXCallModule('AutoAxes',minput,moutput)

  DXRegisterInputHandler(interactor_handler)

  def dxinter(l_init=0,name=None,dxobject=None):
    if dxobject is None: dxobject = __main__.wgui.dxobject
    if name is None: name = __main__.wgui.dxname

    DXCheckRIH(0)

    minput = {'name':name}
    moutput = ['where','size','events']
    (wwhere,wsize,wevents,) = DXCallModule('SuperviseWindow',minput,moutput)

    if wevents.isnull() and not l_init:
      DXDelete(wwhere)
      DXDelete(wsize)
    else:
      # --- This is needed and sometimes not.
      DXReference(dxobject)
      minput = {'where':wwhere,'size':wsize,'events':wevents,
                'object':dxobject,'mode':interactor,'resetObject':1}
      if l_init:
        minput['defaultCamera'] = camera
        minput['resetCamera'] = 1
      moutput = ['object','camera','where']
      (dobject,dcamera,dwhere,) = DXCallModule('SuperviseState',minput,moutput)

      minput = {'object':dobject,'camera':dcamera,'where':dwhere}
      moutput = []
      DXCallModule('Display',minput,moutput)

  try:
    __main__.wgui.dxinter = dxinter
    __main__.wgui.dxname  = name
    __main__.wgui.dxobject = dxobject
    dxinter(1)
    __main__.wgui.dx_timer = wxPyTimer(__main__.wgui.dxinter)
    __main__.wgui.dx_timer.Start(100)
  except:
    dxinter(1,name,dxobject)
    while interactor >=0:
      dxinter(0,name,dxobject)
    DXDelete(dxobject)
    interactor = 0

def DXNewImage():
  _group.reset()

# This may not be the best thing to do, but it works - it gives access to this
# function without having to explicitly import this module. Note that this
# only works for interactive sessions.
__main__.__dict__['DXNewImage'] = DXNewImage


#==========================================================================
def DXWriteImage(filename='image',object=None,camera=None,labels=None,
                 format=None):
  """
Writes an image of the object to a file.
 - filename='image': file name. Its suffix can specify the format.
 - object: the object to image - it must be specified
 - camera: camera viewing the object, the default to view down the z-axis
 - labels: list of plot labels for the three axis. If given, axis are plotted.
 - format: file format, defaults to 'rgb'. Other options include
           'tiff', 'ps', 'eps'. See OpenDX docs for more format options.
  """

  assert object is not None,"object must be specified"

  group = DXCollection(object)
  dxobject = group.getdxobject()
  assert isDXObject(dxobject),'Object to display is not a DXObject'
  if labels is None: labels = group.labels

  if camera is None:
    DXReference(dxobject)
    minput = {'object':dxobject}
    moutput = ['camera']
    (camera,) = DXCallModule('AutoCamera',minput,moutput)

  if labels is not None:
    assert (len(labels) == 3),"Length of lables list must be three"
    labels = DXMakeStringList(labels)
    DXReference(camera)
    minput = {'input':dxobject,'camera':camera,'labels':labels,
              'ticks':3,'colors':'yellow'}
    moutput = ['axes']
    (dxobject,) = DXCallModule('AutoAxes',minput,moutput)

  minput = {'object':dxobject,'camera':camera}
  moutput = ['image']
  (image,) = DXCallModule('Render',minput,moutput)

  minput = {'image':image,'name':filename}
  if format is not None: minput['format'] = format
  moutput = []
  DXCallModule('WriteImage',minput,moutput)

