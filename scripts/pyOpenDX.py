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
from pyDXObject import *
import __main__

pyOpenDX_version = "$Id: pyOpenDX.py,v 1.8 2004/05/22 01:48:35 dave Exp $"
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

def ppxyz(iw = 0,labels=1,display=1,**kw):
  """Plots X-Y-Z"""
  checkparticleplotarguments(kw)
  ii = selectparticles(iw=iw,kwdict=kw)
  if labels == 1: labels = ['X','Y','Z']
  return viewparticles(take(top.xp,ii),take(top.yp,ii),take(top.zp,ii),
                       take(top.uxp,ii),
                       labels,name='WARP viz',display=display)

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

  DXReference(camera)
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

  DXReference(camera)
  DXReference(group)
  DXImage(group,camera,name)

###########################################################################
def viewparticles(x,y,z,v,labels=None,name='WARP particles',
                  display=1):
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
  minput = {'data':dxf,'type':'speedy','ratio':1.0}
  moutput = ['glyphs']
  (glyphs,) = DXCallModule('AutoGlyph',minput,moutput)

  minput = {'data':glyphs,'opacity':.5}
  moutput = ['mapped']
  (dxobject,) = DXCallModule('AutoColor',minput,moutput)
  DXReference(dxobject)

  if display:
    DXImage(dxobject,name=name,labels=labels)
    DXDelete(dxobject)
  else:
    return dxobject

###########################################################################
class DXCollection:
  def __init__(self,*dxobjects,**kw):
    labels = kw.get('labels',None)
    self.labels = labels
    self.dxobject = None
    for o in dxobjects: self.addobject(o)
  def getdxobject(self):
    return self.dxobject
  def extractdxobject(self,object):
    dxobject = object
    try:
      # --- This is horribly kludgy, but works. This keeps looking until
      # --- it finds an object which doesn't have dxobject as an attribute.
      # --- Assumes that that will be an actual dxobject. This will fail
      # --- miserably if an objects dxobject reference is recursive.
      # --- Please, let that never happen!
      while 1:
        dxobject = dxobject.getdxobject()
    except AttributeError:
      pass
    return dxobject
  def createdxobject(self,object):
    dxobject = self.extractdxobject(object)
    minput = {'object':dxobject}
    moutput = ['group']
    (self.dxobject,) = DXCallModule('Collect',minput,moutput)
    DXReference(self.dxobject)
  def addobject(self,object):
    if self.dxobject is None:
      self.createdxobject(object)
    else:
      dxobject = self.extractdxobject(object)
      minput = {'input':self.dxobject,'object':dxobject}
      moutput = ['group']
      (g,) = DXCallModule('Append',minput,moutput)
      DXDelete(self.dxobject)
      self.dxobject = g
      DXReference(self.dxobject)
  def reset(self):
    DXDelete(self.dxobject)
    self.dxobject = None
  def reference(self):
    DXReference(self.dxobject)

###########################################################################
#==========================================================================
interactor = 0
def interactor_handler():
  global interactor
  getchar = sys.stdin.readline()[:-1]
  if len(getchar) > 0: interactor = eval(getchar)
  else:                interactor = -1

_group = DXCollection()
def DXImage(object,camera=None,name='WARP viz',labels=None):
  global interactor

  _group.addobject(object)
  dxobject = _group.dxobject
  if labels is None: labels = _group.labels

  if camera is None:
    cameraadded = 1
    DXReference(dxobject)
    minput = {'object':dxobject}
    moutput = ['camera']
    (camera,) = DXCallModule('AutoCamera',minput,moutput)
    DXReference(camera)
  else:
    cameraadded = 0

  if labels is not None:
    assert (len(labels) == 3),"Length of lables list must be three"
    labelsadded = 1
    labels = DXMakeStringList(labels)
    minput = {'input':dxobject,'camera':camera,'labels':labels,
              'ticks':3,'colors':'yellow'}
    moutput = ['axes']
    olddxobject = dxobject
    (dxobject,) = DXCallModule('AutoAxes',minput,moutput)
    DXReference(dxobject)
  else:
    labelsadded = 0

  i = 0
  DXRegisterInputHandler(interactor_handler)
  while interactor >= 0:
    i = i + 1
    DXCheckRIH(0)

    minput = {'name':name}
    moutput = ['where','size','events']
    (wwhere,wsize,wevents,) = DXCallModule('SuperviseWindow',minput,moutput)

    if wevents.isnull() and i != 1:
      DXDelete(wwhere)
      DXDelete(wsize)
    else:
      minput = {'where':wwhere,'size':wsize,'events':wevents,
                'object':dxobject,'mode':interactor,'resetObject':1}
      if i == 1:
        minput['defaultCamera'] = camera
        minput['resetCamera'] = 1
      moutput = ['object','camera','where']
      (dobject,dcamera,dwhere,) = DXCallModule('SuperviseState',minput,moutput)

      minput = {'object':dobject,'camera':dcamera,'where':dwhere}
      moutput = ['where']
      (dwhere,) = DXCallModule('Display',minput,moutput)

  DXDelete(dwhere)
  if cameraadded: DXDelete(camera)
  if labelsadded: DXDelete(olddxobject)
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
  if labels is None: labels = group.labels

  if camera is None:
    cameraadded = 1
    DXReference(dxobject)
    minput = {'object':dxobject}
    moutput = ['camera']
    (camera,) = DXCallModule('AutoCamera',minput,moutput)
    DXReference(camera)
  else:
    cameraadded = 0

  if labels is not None:
    assert (len(labels) == 3),"Length of lables list must be three"
    labelsadded = 1
    labels = DXMakeStringList(labels)
    minput = {'input':dxobject,'camera':camera,'labels':labels,
              'ticks':3,'colors':'yellow'}
    moutput = ['axes']
    olddxobject = dxobject
    (dxobject,) = DXCallModule('AutoAxes',minput,moutput)
    DXReference(dxobject)
  else:
    labelsadded = 0

  minput = {'object':dxobject,'camera':camera}
  moutput = ['image']
  (image,) = DXCallModule('Render',minput,moutput)

  minput = {'image':image,'name':filename}
  if format is not None: minput['format'] = format
  moutput = []
  DXCallModule('WriteImage',minput,moutput)

  if cameraadded: DXDelete(camera)
  if labelsadded: DXDelete(olddxobject)
