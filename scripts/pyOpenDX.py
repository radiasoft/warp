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
try:
  from wxPython.wx import *
except:
  pass
import __main__
__main__.selectbox = 1
__main__.l_hardware_acceleration = 1
__main__.l_dxforceupdate = 0
__main__.l_dxupdate_all_windows = 1
__main__.l_dxnewwindow = 0
__main__.DXWindows={}

try:
  # --- Try importing. If not found, still define the functions and classes.
  # --- But any calls require the pyOpenDX attributes will raise an exception.
  from pyDXObject import *
except:
  pass

pyOpenDX_version = "$Id: pyOpenDX.py,v 1.24 2004/12/22 21:53:20 jlvay Exp $"
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

def ppxyz(iw = 0,cc=None,labels=1,display=1,rscale=None,zscale=None,size=1.,ratio=1.,stride=1,type='speedy',scale=None,**kw):
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
                       labels,name=None,display=display,size=size,ratio=ratio,stride=stride,type=type,scale=scale)

def ppxyzvxvyvz(iw = 0,labels=1,display=1,rscale=None,zscale=None,size=3.,ratio=0.,stride=1,type='standard',scale=None,**kw):
  """Vector Plot X-Y-Z-Vx-Vy-Vz"""
  checkparticleplotarguments(kw)
  ii = selectparticles(iw=iw,kwdict=kw)
  if labels == 1: labels = ['X','Y','Z']
  xx = take(top.xp,ii)
  yy = take(top.yp,ii)
  zz = take(top.zp,ii)
  vx = take(top.uxp,ii)
  vy = take(top.uyp,ii)
  vz = take(top.uzp,ii)
  gg = take(top.gaminv,ii)
  vx=vx*gg
  vy=vy*gg
  vz=vz*gg
  if rscale is not None:
    xx = xx*rscale
    yy = yy*rscale
  if zscale is not None:
    zz = zz*zscale
  return viewvparticles(xx,yy,zz,vx,vy,vz,
                       labels,name=None,display=display,size=size,ratio=ratio,stride=stride,type=type,scale=scale,normalize=1)

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
def viewparticles(x,y,z,v,labels=None,name=None,
                  display=1,size=1.,ratio=1.,stride=1,type='standard',scale=None):
  if stride==1:
    x = gatherarray(x)
    y = gatherarray(y)
    z = gatherarray(z)
    v = gatherarray(v)
  else:
    x = gatherarray(x[::stride])
    y = gatherarray(y[::stride])
    z = gatherarray(z[::stride])
    v = gatherarray(v[::stride])
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
  if scale is not None:
    origin = array([0.,0.,0.]).astype(Float32)
    dxorigin = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
    DXAddArrayData(dxorigin,0,1,origin)

    dxdata = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,0)
    DXAddArrayData(dxdata,0,1,array(1.,Float32))

    minput = {'origin':dxorigin,'data':dxdata}
    moutput = ['output']
    (glyph_data,) = DXCallModule('Construct',minput,moutput)

    minput = {'data':glyph_data,'type':type,'scale':1.}
    moutput = ['glyphs']
    (glyph,) = DXCallModule('Glyph',minput,moutput)
    
    dxscale = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
    DXAddArrayData(dxscale,0,1,array([1./scale[0],1./scale[1],1./scale[2]],Float32))
    minput  = {'input':glyph,'scaling':dxscale}
    moutput = ['output']
    (type,) = DXCallModule('Scale',minput,moutput)

  minput = {'data':dxf,'type':type,'ratio':ratio,'scale':size}
  moutput = ['glyphs']
  (glyphs,) = DXCallModule('AutoGlyph',minput,moutput)

  minput = {'data':glyphs}#,'opacity':.5} # significant slow down with opacity<1.
  moutput = ['mapped']
  (dxobject,) = DXCallModule('AutoColor',minput,moutput)

  if display:
    DXImage(dxobject,name=name,labels=labels,scale=scale)
  else:
    return dxobject

###########################################################################
def viewvparticles(x,y,z,vx,vy,vz,labels=None,name=None,
                  display=1,size=1.,ratio=0.,stride=1,type='standard',scale=None,normalize=0):
  if stride==1:
    x = gatherarray(x)
    y = gatherarray(y)
    z = gatherarray(z)
    vx = gatherarray(vx)
    vy = gatherarray(vy)
    vz = gatherarray(vz)
  else:
    x = gatherarray(x[::stride])
    y = gatherarray(y[::stride])
    z = gatherarray(z[::stride])
    vx = gatherarray(vx[::stride])
    vy = gatherarray(vy[::stride])
    vz = gatherarray(vz[::stride])
  if me > 0: return
  if normalize:
    vx = vx/ave(abs(vx))
    vy = vy/ave(abs(vy))
    vz = vz/ave(abs(vz))

  # --- First combine particle data and create a DX array
  n = len(x)
  p = zeros((n,3),'d')
  p[:,0] = x
  p[:,1] = y
  p[:,2] = z
  v = zeros((n,3),'d')
  v[:,0] = vx
  v[:,1] = vy
  v[:,2] = vz
  dxp = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
  DXAddArrayData(dxp,0,n,p.astype(Float32))
  dxv = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
  DXAddArrayData(dxv,0,n,v.astype(Float32))

  # --- Create the field
  dxf = DXNewField()
  DXSetComponentValue(dxf,'positions',dxp)
  DXSetComponentValue(dxf,'data',dxv)
  DXEndField(dxf)

  # --- Create glyphs to render data
  minput = {'data':dxf,'type':type,'ratio':ratio,'scale':size}
  moutput = ['glyphs']
  (glyphs,) = DXCallModule('AutoGlyph',minput,moutput)

  minput = {'data':glyphs}#,'opacity':.5}
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

def DXMountainPlot(f,xmin=0.,ymin=0.,dx=1.,dy=1.,scale=1.,display=1,labels=['x','y','z'],name = None,perspective=1):
  fshape = shape(f)
  nx = fshape[0]
  ny = fshape[1]
  dxdata = DXNewArray(TYPE_DOUBLE,CATEGORY_REAL,0)
  DXAddArrayData(dxdata,0,nx*ny,f)

  origin = array([xmin,ymin]).astype(Float32)
  dxorigin = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,2)
  DXAddArrayData(dxorigin,0,1,origin)

  deltas = array([dx,dy]).astype(Float32)
  dxdeltas = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,2)
  DXAddArrayData(dxdeltas,0,1,deltas)

  counts = array([nx,ny])
  dxcounts = DXNewArray(TYPE_INT,CATEGORY_REAL,1,2)
  DXAddArrayData(dxcounts,0,1,counts)

  minput = {'origin':dxorigin,'deltas':dxdeltas,'counts':dxcounts,'data':dxdata}
  moutput = ['output']
  (dxf2,) = DXCallModule('Construct',minput,moutput)

  minput = {'data':dxf2}
  moutput = ['mapped']
  (colored,) = DXCallModule('AutoColor',minput,moutput)
    
  frange = maxnd(f)-minnd(f)
  if frange == 0.:
    scale = 1.
  else:
    xrange = sqrt((nx*dx)**2+(ny*dy)**2)
    scale = 2.*scale*(0.1*xrange/frange)
  minput = {'data':colored,'scale':scale}
  moutput = ['graph']
  (dxobject,) = DXCallModule('RubberSheet',minput,moutput)

  up = array([0,0,1])
  dxup = DXNewArray(TYPE_INT,CATEGORY_REAL,1,3)
  DXAddArrayData(dxup,0,1,up)
  DXReference(dxobject)
  minput = {'object':dxobject,'direction':'diagonal','up':dxup,'perspective':perspective}
  moutput = ['camera']
  (camera,) = DXCallModule('AutoCamera',minput,moutput)

  if display:
    DXImage(dxobject,name=name,labels=labels,camera=camera)
  else:
    return dxobject
 

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

class DXWindow:
  def __init__(self,name):
    self.name = name

###########################################################################
#==========================================================================
def interactor_handler():
  getchar = sys.stdin.readline()[:-1]
  if len(getchar) > 0: __main__.dxwindow.interactor = eval(getchar)
  else:                __main__.dxwindow.interactor = -1

def DXSuperviseWindow(name,visibility=1):
    minput = {'name':name,'visibility':visibility}
    moutput = ['where','size','events']
    return DXCallModule('SuperviseWindow',minput,moutput)

def DXScale(dxobject,scale):
    dxscale = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
    DXAddArrayData(dxscale,0,1,array(scale,Float32))
    minput  = {'input':dxobject,'scaling':dxscale}
    moutput = ['output']
    (dxobject_out,) = DXCallModule('Scale',minput,moutput)
    return dxobject_out

def DXAutoAxes(dxobject,camera,labels):
    assert (len(labels) == 3),"Length of labels list must be three"
    DXReference(camera)
    labels = DXMakeStringList(labels)
    minput = {'input':dxobject,'camera':camera,'labels':labels,
              'ticks':3,'adjust':1,'colors':'grey'}#'yellow'}
    moutput = ['axes']
    (dxobject_out,) = DXCallModule('AutoAxes',minput,moutput)
    return dxobject_out

def DXImage(object,camera=None,name=None,labels=None,hardware=1,scale=None):
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

  if name is None:
    if __main__.l_dxnewwindow or (__main__.DXWindows=={}): 
      name='WARP viz %g'%(len(__main__.DXWindows.keys())+1)
    else:
      name=__main__.dxwindow.name

  group = DXCollection(object)
  dxobject = group.getdxobject()
  assert isDXObject(dxobject),'Object to display is not a DXObject'

  # --- Turn off the caching of images - caching can use a lot of memory
  # --- which is not freed.
  minput = {'input':dxobject,'attribute':'cache','value':0}
  moutput = ['output']
  (dxobject,) = DXCallModule('Options',minput,moutput)
  if not __main__.DXWindows.has_key('name'):
    __main__.DXWindows[name]=DXWindow(name)
  __main__.dxwindow=__main__.DXWindows[name]
  
  __main__.dxwindow.dxobject_init = dxobject
  __main__.dxwindow.dxobject      = dxobject
  __main__.dxwindow.interactor    = 0
  __main__.dxwindow.previousmode  = 0

  if hardware:
    __main__.dxwindow.l_hardware_acceleration=1
  else:
    __main__.dxwindow.l_hardware_acceleration=0

  if scale is None:scale=[1.,1.,1.]
  __main__.dxwindow.dxscale=scale
  __main__.dxwindow.dxobject = DXScale(__main__.dxwindow.dxobject_init,__main__.dxwindow.dxscale)

  if camera is None:
    DXReference(__main__.dxwindow.dxobject)
    minput = {'object':__main__.dxwindow.dxobject}
    moutput = ['camera']
    (camera,) = DXCallModule('AutoCamera',minput,moutput)

  try:
    __main__.gui_is_on = __main__.wgui.initialized
  except:
    __main__.gui_is_on = 0
  if  not __main__.gui_is_on:
    DXRegisterInputHandler(interactor_handler)

  def dxinter_loop(dxwindow,l_init=0,forceupdate=0):
    if dxwindow.wevents.isnull() and not l_init and not forceupdate:
      DXDelete(dxwindow.wwhere)
      DXDelete(dxwindow.wsize)
    else:
      if dxwindow.l_dxhaslabels:
        if dxwindow.l_dxrescaled:
          dxwindow.dxobject = DXAutoAxes(dxwindow.dxobject,camera,labels) 
          if dxwindow.l_hardware_acceleration:
            DXRendering(dxwindow,rendering='hardware')
          else:
            DXRendering(dxwindow,rendering='software')
          dxwindow.l_dxrescaled=0
      dxobject = dxwindow.dxobject

      # --- This is needed and sometimes not.
      DXReference(dxobject)
      minput = {'where':dxwindow.wwhere,'size':dxwindow.wsize,'events':dxwindow.wevents,
                'object':dxobject,'mode':dxwindow.interactor,'resetObject':1}
      if l_init:
        minput['defaultCamera'] = camera
        minput['resetCamera'] = 1
      moutput = ['object','camera','where']
      (dxwindow.dobject,dxwindow.dcamera,dxwindow.dwhere,) = DXCallModule('SuperviseState',minput,moutput)

      minput = {'object':dxwindow.dobject,'camera':dxwindow.dcamera,'where':dxwindow.dwhere}
      moutput = []
      DXCallModule('Display',minput,moutput)

  def dxinter(l_init=0):
    DXCheckRIH(0)
    winit = __main__.dxwindow

    (winit.wwhere,winit.wsize,winit.wevents,) = DXSuperviseWindow(winit.name)

    # --- Special coding is needed when hardware acceleration is used. The
    # --- hardware seems to be capturing all events so that the events
    # --- returned by SuperviseWindow is always null. Because of this, the
    # --- code below is never called. This causes requests to change the
    # --- interactor mode to be ignored. Now, if the interactor mode is
    # --- changed, the code is forced to reexecute the display.
    # --- This qualifies as a kludge since it is not understood exactly the
    # --- events is alway null when using hardware acceleration.
    try: dxinter.previousmode
    except: dxinter.previousmode = 0
    forceupdate = (((winit.interactor != winit.previousmode) and
                   (hardware or dxwindow.l_hardware_acceleration)) or 
                   __main__.l_dxforceupdate)
    if __main__.l_dxupdate_all_windows:
      for w in __main__.DXWindows.itervalues():
        if w is not winit:
          (w.wwhere,w.wsize,w.wevents,) = DXSuperviseWindow(w.name)
        w.interactor=winit.interactor
        w.l_hardware_acceleration=winit.l_hardware_acceleration
        dxinter_loop(w,l_init,forceupdate)
    else:
      dxinter_loop(winit,l_init,forceupdate)

    __main__.l_dxforceupdate = 0
    winit.previousmode = winit.interactor

  if labels is not None:
    __main__.dxwindow.l_dxrescaled=1
    __main__.dxwindow.l_dxhaslabels=1
    __main__.dxwindow.dxlabels = labels
  else:
    __main__.dxwindow.l_dxhaslabels=0
  if __main__.gui_is_on:
    __main__.wgui.dxinter = dxinter
    __main__.wgui.dxname   = name
    if __main__.dxwindow.l_hardware_acceleration:
      DXRendering(__main__.dxwindow,rendering='hardware')
    tmpvar=__main__.l_dxupdate_all_windows
    __main__.l_dxupdate_all_windows=0
    dxinter(1)
    __main__.l_dxupdate_all_windows=tmpvar
    try:
      timer = __main__.wgui.dx_timer
    except:
      __main__.wgui.dx_timer = wxPyTimer(__main__.wgui.dxinter)
      __main__.wgui.dx_timer.Start(100)
  else:
    dxinter(1)
    print "Press 0 for rotation, 1 for panning, 2 for zooming, return to exit."
    print "Default mode is rotation."
    while __main__.dxwindow.interactor >=0:
      dxinter()
    DXDelete(dxobject)
    __main__.interactor = 0

def DXRendering(dxwindow,rendering='software'):
  minput = {'input':dxwindow.dxobject,'attribute':'rendering mode','value':rendering}
  moutput = ['output']
  (dxwindow.dxobject,) = DXCallModule('Options',minput,moutput)
  __main__.l_dxforceupdate=1
  if rendering=='hardware':
    dxwindow.l_hardware_acceleration=1
  else:
    dxwindow.l_hardware_acceleration=0

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

  if format=="dx" or format=="vrml":
    minput = {'object':dxobject,'name':filename,'format':format}
    moutput = []
    DXCallModule('Export',minput,moutput)
  else:
    minput = {'object':dxobject,'camera':camera}
    moutput = ['image']
    (image,) = DXCallModule('Render',minput,moutput)

    minput = {'image':image,'name':filename}
    if format is not None: minput['format'] = format
    moutput = []
    DXCallModule('WriteImage',minput,moutput)

