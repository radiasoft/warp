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
  import wx
except:
  pass
import __main__
__main__.selectbox        = 0
__main__.selectbox_remove = 1
__main__.selectbox_xmin   = None
__main__.selectbox_xmax   = None
__main__.selectbox_ymin   = None
__main__.selectbox_ymax   = None
__main__.selectbox_zmin   = None
__main__.selectbox_zmax   = None
__main__.l_hardware_acceleration = 1
__main__.l_dxforceupdate = 0
__main__.l_dxupdate_all_windows = 1
__main__.l_dxnewwindow = 0
__main__.l_dxresetcamera=0 
__main__.dxperspective = 0
__main__.dxangle = 30
__main__.dxtimerconstant = 10#0 # time between each interruption in ms
__main__.DXWindows={}

try:
  # --- Try importing. If not found, still define the functions and classes.
  # --- But any calls require the pyOpenDX attributes will raise an exception.
  from pyDXObject import *
except:
  pass

pyOpenDX_version = "$Id: pyOpenDX.py,v 1.34 2007/12/20 01:34:43 dave Exp $"
def pyOpenDXdoc():
  import pyOpenDX
  print pyOpenDX.__doc__

###########################################################################
def ppxxpy(iw = 0,labels=1,display=1,**kw):
  """Plots X-Y-Z"""
  checkparticleplotarguments(kw)
  pgroup = kw.get('pgroup',top.pgroup)
  if labels == 1: labels = ['X',"X'",'Y']
  ii = selectparticles(iw=iw,kwdict=kw)
  x  = getx(ii=ii,**kw)
  xp = getxp(ii=ii,**kw)
  y  = gety(ii=ii,**kw)
  yp = getyp(ii=ii,**kw)
  return viewparticles(x,xp,y,yp,
                       labels,name='WARP viz',display=display)

def ppxyz(iw = 0,cc=None,labels=1,display=1,rscale=None,zscale=None,size=1.,ratio=1.,stride=1,
          type='speedy',color=None,scale=None,**kw):
  """Plots X-Y-Z"""
  checkparticleplotarguments(kw)
  pgroup = kw.get('pgroup',top.pgroup)
  ii = selectparticles(iw=iw,kwdict=kw)
  if labels == 1: labels = ['X','Y','Z']
  xx = getx(ii=ii,**kw)
  yy = gety(ii=ii,**kw)
  zz = getz(ii=ii,**kw)
  if cc is None:
    cc = getux(ii=ii,**kw)
  else:
    if isinstance(ii,slice): cc = cc[ii]
    else:                    cc = take(cc,ii)
  if rscale is not None:
    xx = xx*rscale
    yy = yy*rscale
  if zscale is not None:
    zz = zz*zscale
  return viewparticles(xx,yy,zz,cc,
                       labels,name=None,display=display,size=size,ratio=ratio,
                       stride=stride,type=type,scale=scale,color=color)

def ppxyzvxvyvz(iw = 0,labels=1,display=1,rscale=None,zscale=None,size=3.,ratio=0.,stride=1,type='standard',scale=None,**kw):
  """Vector Plot X-Y-Z-Vx-Vy-Vz"""
  checkparticleplotarguments(kw)
  pgroup = kw.get('pgroup',top.pgroup)
  if labels == 1: labels = ['X','Y','Z']
  ii = selectparticles(iw=iw,kwdict=kw)
  xx = getx(ii=ii,**kw)
  yy = gety(ii=ii,**kw)
  zz = getz(ii=ii,**kw)
  vx = getvx(ii=ii,**kw)
  vy = getvy(ii=ii,**kw)
  vz = getvz(ii=ii,**kw)
  if rscale is not None:
    xx = xx*rscale
    yy = yy*rscale
  if zscale is not None:
    zz = zz*zscale
  return viewvparticles(xx,yy,zz,vx,vy,vz,
                       labels,name=None,display=display,size=size,ratio=ratio,stride=stride,type=type,scale=scale,normalize=1)

def pprzvrvtvz(iw = 0,labels=1,display=1,rscale=None,zscale=None,size=3.,ratio=0.,stride=1,type='standard',scale=None,**kw):
  """Vector Plot X-Y-Z-Vx-Vy-Vz"""
  checkparticleplotarguments(kw)
  if labels == 1: labels = ['X','Y','Z']
  ii = selectparticles(iw=iw,kwdict=kw)
  xx = getx(ii=ii,**kw)
  yy = gety(ii=ii,**kw)
  zz = getz(ii=ii,**kw)
  r  = getr(ii=ii,**kw)
  vx = getvx(ii=ii,**kw)
  vy = getvy(ii=ii,**kw)
  vz = getvz(ii=ii,**kw)
  vr =  vx*xx/r+vy*yy/r
  vt = -vx*yy/r+vy*xx/r
  if rscale is not None:
    xx = xx*rscale
    yy = yy*rscale
  if zscale is not None:
    zz = zz*zscale
  return viewvparticles(r,zeros(shape(r)[0],'d'),zz,vr,vt,vz,
                       labels,name=None,display=display,size=size,ratio=ratio,stride=stride,type=type,scale=scale,normalize=1)

###########################################################################
def viewisosurface1(data,isovalue,origins=None,deltas=None,name='WARP viz',display=1,color=None,intensity=1.,
                    opacity=0.5,colorbar=1,colormap=None,include_mins=None,include_maxs=None,exclude_mins=None,exclude_maxs=None):
  if origins is None: origins = [0.,0.,0.]
  if deltas is None: deltas = [1.,1.,1.]
  f = DXObject_fromarray(data,origins,deltas)

  if exclude_mins is not None and exclude_maxs is None:
    exclude_maxs=[w3d.xmmax,w3d.ymmax,w3d.zmmax]
  if exclude_maxs is not None and exclude_mins is None:
    exclude_mins=[w3d.xmmin,w3d.ymmin,w3d.zmmin]
  if exclude_mins is not None:   
    f=DXExcavate(f,exclude_mins,exclude_maxs)

  if include_mins is not None and include_maxs is None:
    include_maxs=[w3d.xmmax,w3d.ymmax,w3d.zmmax]
  if include_maxs is not None and include_mins is None:
    include_mins=[w3d.xmmin,w3d.ymmin,w3d.zmmin]
  if include_mins is not None:   
    f=DXUnmark( DXInclude( DXMark(f,'positions') ,include_mins,include_maxs) ,'positions')

  if type(isovalue) is type([]) or type(isovalue) is type(array(0)):
    isovalues=array(isovalue).astype(Float32)
    isovalue=DXNewArray(TYPE_FLOAT,CATEGORY_REAL,0)
    DXAddArrayData(isovalue,0,shape(isovalues)[0],isovalues)

  minput = {'data':f,'value':isovalue}
  moutput = ['surface']
  (isosurface,) = DXCallModule('Isosurface',minput,moutput)

  if color is not None:
    if color=='auto' or colormap is not None:
      if colormap is None:
        isosurface,colormap = DXAutoColor(isosurface,opacity=opacity,intensity=intensity)
      else:
        isosurface = DXColor(isosurface,color=colormap,opacity=opacity)        
    elif color=='greyscale':
      isosurface = DXAutoGrayScale(isosurface,opacity=opacity,saturation=intensity)
    else:
      isosurface = DXColor(isosurface,color=color,opacity=opacity)

  if color is not None and colorbar:
    dxcolorbar = DXColorBar(colormap)

  if display:
    DXReference(isosurface)
    minput = {'object':isosurface}
    moutput = ['camera']
    (camera,) = DXCallModule('AutoCamera',minput,moutput)
    if color is not None and colorbar:
      DXImage(DXCollect([isosurface,dxcolorbar]),camera,name)
    else:
      DXImage(isosurface,camera,name)
  else:
    if color is not None and colorbar:
      return isosurface,dxcolorbar
    else:
      return isosurface

def viewisosurface(data,isovalue,name='WARP viz',display=1):
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

  if display:
    DXImage(group,camera,name)
  else:
    return group

def viewcoloredvolume(data,display=1,origins=None,deltas=None,name='WARP viz',opacity=0.5,intensity=2.):
  if origins is None: origins = [0.,0.,0.]
  if deltas is None: deltas = [1.,1.,1.]
  f = DXObject_fromarray(data,origins,deltas)

  dxobject,colormap = DXAutoColor(f,opacity=opacity,intensity=intensity)
  if display:
    DXImage(dxobject,camera=None,name=name)
  else:
    return dxobject
  
def viewgreyscalevolume(data,display=1,origins=None,deltas=None,name='WARP viz',opacity=0.5,hue=0.,saturation=0.):
  if origins is None: origins = [0.,0.,0.]
  if deltas is None: deltas = [1.,1.,1.]
  f = DXObject_fromarray(data,origins,deltas)

#  if color is not None:
#    f = DXColor(f,color)
  dxobject = DXAutoGrayScale(f,opacity=opacity,hue=hue,saturation=saturation)
  if display:
    DXImage(dxobject,camera=None,name=name)
  else:
    return fxobject
  
###########################################################################
def viewparticles(x,y,z,v,labels=None,name=None,
#                  display=1,size=1.,ratio=1.,stride=1,type='standard',scale=None,auto=0,shape=None,
                  display=1,size=1.,ratio=None,stride=1,type='standard',scale=None,auto=0,shape=None,min=None,max=None,
                  color=None,intensity=1.,opacity=0.5,colorbar=1,colormap=None,cmin=None,cmax=None,
                  include_mins=None,include_maxs=None,exclude_mins=None,exclude_maxs=None):
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

  if 1: # alternative to method below
  # include/exclude particles
   if exclude_mins is not None and exclude_maxs is None:
    exclude_maxs=[w3d.xmmax,w3d.ymmax,w3d.zmmax]
   if exclude_maxs is not None and exclude_mins is None:
    exclude_mins=[w3d.xmmin,w3d.ymmin,w3d.zmmin]
   if exclude_mins is not None:   
     n = len(x);
     ii = compress( 1- array((x>=exclude_mins[0]) & (x<=exclude_maxs[0]) & \
                             (y>=exclude_mins[1]) & (y<=exclude_maxs[1]) & \
                             (z>=exclude_mins[2]) & (z<=exclude_maxs[2])),arange(n))
     x = take(x,ii)
     y = take(y,ii)
     z = take(z,ii)
     v = take(v,ii)
     
   if include_mins is not None and include_maxs is None:
    include_maxs=[w3d.xmmax,w3d.ymmax,w3d.zmmax]
   if include_maxs is not None and include_mins is None:
    include_mins=[w3d.xmmin,w3d.ymmin,w3d.zmmin]
   if include_mins is not None:   
     n = len(x)
     ii = compress((x>=include_mins[0]) & (x<=include_maxs[0]) & \
                   (y>=include_mins[1]) & (y<=include_maxs[1]) & \
                   (z>=include_mins[2]) & (z<=include_maxs[2]),arange(n))
     x = take(x,ii)
     y = take(y,ii)
     z = take(z,ii)
     v = take(v,ii)

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
  
  if 0: #for some reason, it does not work
  # include/exclude particles
   if exclude_mins is not None and exclude_maxs is None:
    exclude_maxs=[w3d.xmmax,w3d.ymmax,w3d.zmmax]
   if exclude_maxs is not None and exclude_mins is None:
    exclude_mins=[w3d.xmmin,w3d.ymmin,w3d.zmmin]
   if exclude_mins is not None:   
    dxf=DXExcavate(dxf,exclude_mins,exclude_maxs)

   if include_mins is not None and include_maxs is None:
    include_maxs=[w3d.xmmax,w3d.ymmax,w3d.zmmax]
   if include_maxs is not None and include_mins is None:
    include_mins=[w3d.xmmin,w3d.ymmin,w3d.zmmin]
   if include_mins is not None:   
    dxf=DXUnmark( DXInclude( DXMark(dxf,'positions') ,include_mins,include_maxs) ,'positions')

  # --- Create glyphs to render data
  if scale is not None:
    dxorigin = DXVector([0.,0.,0.])

    dxdata = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,0)
    DXAddArrayData(dxdata,0,1,array(1.,Float32))

    minput = {'origin':dxorigin,'data':dxdata}
    moutput = ['output']
    (glyph_data,) = DXCallModule('Construct',minput,moutput)

    minput = {'data':glyph_data,'type':type,'scale':1.}
    moutput = ['glyphs']
    (glyph,) = DXCallModule('Glyph',minput,moutput)
    
    dxscale = DXVector([1./scale[0],1./scale[1],1./scale[2]])
    minput  = {'input':glyph,'scaling':dxscale}
    moutput = ['output']
    (type,) = DXCallModule('Scale',minput,moutput)

#  minput = {'data':dxf,'type':type,'ratio':ratio,'scale':size}
#  moutput = ['glyphs']
#  (glyphs,) = DXCallModule('Glyph',minput,moutput)
  glyphs = DXGlyph(data=dxf,type=type,scale=size,ratio=ratio,auto=auto,min=min,max=max,shape=shape)

  if color is not None:
    if color=='auto' or colormap is not None:
      if colormap is None:
        glyphs,colormap = DXAutoColor(glyphs,opacity=opacity,intensity=intensity,min=cmin,max=cmax)
      else:
        glyphs = DXColor(glyphs,color=colormap,opacity=opacity)        
        colorbar=0
    elif color=='greyscale':
      glyphs = DXAutoGrayScale(glyphs,opacity=opacity,saturation=intensity)
    else:
      glyphs = DXColor(glyphs,color=color,opacity=opacity)
      colorbar=0

  if color is not None and colorbar and display:
    dxcolorbar = DXColorBar(colormap)

  if display:
    if color is not None and colorbar:
      DXImage(DXCollect([glyphs,dxcolorbar]),name=name,labels=labels,scale=scale)
    else:
      DXImage(glyphs,name=name,labels=labels,scale=scale)
  else:
    if color is not None and colorbar:
      return glyphs,colormap
    else:
      return glyphs

###########################################################################
def viewvparticles(x,y,z,vx,vy,vz,labels=None,name=None,
                  color=None,intensity=1.,opacity=1.,colorbar=1,colormap=None,
                  display=1,size=1.,ratio=0.,stride=1,type='standard',scale=None,shape=1.,normalize=0):
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
    vmax=max(sqrt(vx*vx+vy*vy+vz*vz))
    vx = vx/vmax
    vy = vy/vmax
    vz = vz/vmax

  dxf = DXParticles(x,y,z,vx,vy,vz)

  # --- Create glyphs to render data
  minput = {'data':dxf,'type':type,'ratio':ratio,'scale':size,'shape':shape}
  moutput = ['glyphs']
  (glyphs,) = DXCallModule('Glyph',minput,moutput)

  if color is not None:
    if color=='auto' or colormap is not None:
      if colormap is None:
        glyphs,colormap = DXAutoColor(glyphs,opacity=opacity,intensity=intensity)
      else:
        glyphs = DXColor(glyphs,color=colormap,opacity=opacity)        
    elif color=='greyscale':
      glyphs = DXAutoGrayScale(glyphs,opacity=opacity,saturation=intensity)
    else:
      glyphs = DXColor(glyphs,color=color,opacity=opacity)

  if color is not None and colorbar:
    dxcolorbar = DXColorBar(colormap)

  if display:
    if color is not None and colorbar:
      DXImage(DXCollect([glyphs,dxcolorbar]),name=name,labels=labels,scale=scale)
    else:
      DXImage(glyphs,name=name,labels=labels,scale=scale)
  else:
    if color is not None and colorbar:
      return glyphs,dxcolorbar
    else:
      return glyphs

###########################################################################
def viewboundingbox(xmin,xmax,ymin,ymax,zmin,zmax,color='yellow'):
  """
Create a box
  - xmin,xmax,ymin,ymax,zmin,zmax: extent of the box
  """
  # --- Create a field containing two points on opposite side of the
  # --- mesh and use it to find the bounding box.
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

  counts = array([nx,ny]).astype(Int32)
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

  up = array([0,0,1]).astype(Int32)
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
    if __main__.selectbox:
      try:
        dxobject.SelectBox(__main__.selectbox_xmin,
                           __main__.selectbox_xmax,
                           __main__.selectbox_ymin,
                           __main__.selectbox_ymax,
                           __main__.selectbox_zmin,
                           __main__.selectbox_zmax,
                           l_remove=__main__.selectbox_remove)
      except:
        pass
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
  def __init__(self,*dxobjects,**kw):
    self.labels = kw.get('labels',None)
    self.preserve = kw.get('preserve',0)
    self.collection = None
    for o in dxobjects:
      if __main__.selectbox:
        try:
          o.SelectBox(__main__.selectbox_xmin,
                           __main__.selectbox_xmax,
                           __main__.selectbox_ymin,
                           __main__.selectbox_ymax,
                           __main__.selectbox_zmin,
                           __main__.selectbox_zmax,
                           l_remove=__main__.selectbox_remove)
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
    minput  = {'input':dxobject,'scaling':DXVector(scale)}
    moutput = ['output']
    (dxobject_out,) = DXCallModule('Scale',minput,moutput)
    return dxobject_out

def DXAppend(group,dxobject):
    minput = {'input':group,'object':dxobject}
    moutput = ['group']
    (group,) = DXCallModule('Append',minput,moutput)
    return group
    
def DXAutoAxes(dxobject,camera,labels,ticks=3,frame=0,adjust=1,grid=0,colors='grey',
               annotation="all",labelscale=1.,font="variable",corners=None,
               xticklocations=None,yticklocations=None,zticklocations=None,
               xticklabels=None,yticklabels=None,zticklabels=None,):
    assert (len(labels) == 3),"Length of labels list must be three"
    DXReference(camera)
    labels = DXMakeStringList(labels)
    if type(ticks)==type([]):ticks=DXIntegerList(ticks)
    if type(annotation)==type([]):annotation=DXMakeStringList(annotation)
    if type(colors)==type([]):colors=DXMakeStringList(colors)
    minput = {'input':dxobject,
              'camera':camera,
              'labels':labels,
              'ticks':ticks,
              'frame':frame,
              'adjust':adjust,
              'grid':grid,
              'colors':colors,
              'annotation':annotation,
              'labelscale':labelscale,
              'font':font}
    if corners is not None:minput['corners']=DXVector(corners)
    if xticklocations is not None:minput['xticklocations']=DXIntegerList(xticklocations)
    if yticklocations is not None:minput['yticklocations']=DXIntegerList(yticklocations)
    if zticklocations is not None:minput['zticklocations']=DXIntegerList(zticklocations)
    if xticklabels is not None:minput['xticklabels']=DXMakeStringList(xticklabels)
    if yticklabels is not None:minput['yticklabels']=DXMakeStringList(yticklabels)
    if zticklabels is not None:minput['zticklabels']=DXMakeStringList(zticklabels)
    moutput = ['axes']
    (dxobject_out,) = DXCallModule('AutoAxes',minput,moutput)
    return dxobject_out

def DXAutoColor(dxobject,opacity=None,intensity=None,start=None,range=None,saturation=None,
                min=None,max=None,delayed=None,outofrange=None):
    minput = {'data':dxobject}
    moutput = ['mapped','colormap']
    if opacity is not None:minput['opacity']=opacity
    if intensity is not None:minput['intensity']=intensity
    if start is not None:minput['start']=start
    if range is not None:minput['range']=range
    if saturation is not None:minput['saturation']=saturation
    if min is not None:minput['min']=min
    if max is not None:minput['max']=max
    if delayed is not None:minput['delayed']=delayed
    if outofrange is not None:minput['outofrange']=outofrange
    (dxobject_out,colormap) = DXCallModule('AutoColor',minput,moutput)
    return dxobject_out,colormap

def DXAutoGrayScale(dxobject,opacity=0.5,hue=0.,saturation=0.):
    minput = {'data':dxobject,'opacity':opacity,'hue':hue,'saturation':saturation}
    moutput = ['mapped']
    (dxobject_out,) = DXCallModule('AutoGrayScale',minput,moutput)
    return dxobject_out

def DXCaption(string,position=None,flag=None,reference=None,alignment=None,
              height=None,font=None,direction=None,up=None):
    minput = {'string':string}
    if position is not None:minput['position']=DXVector(position)
    if flag is not None:minput['flag']=flag
    if reference is not None:minput['reference']=DXVector(reference)
    if alignment is not None:minput['alignment']=alignment
    if height is not None:minput['height']=height
    if font is not None:minput['font']=font
    if direction is not None:minput['direction']=DXVector(direction)
    if up is not None:minput['up']=DXVector(up)
    moutput = ['caption']
    (dxobject,) = DXCallModule('Caption',minput,moutput)
    return dxobject

def DXColor(dxobject,color=None,opacity=0.5,component='colors',delayed=0):
    if color is None:
      return dxobject
    if type(color) in [type([]),type(()),type(ones(0))]:color=DXVector(color)
    minput = {'input':dxobject,'color':color,'opacity':opacity,'delayed':delayed}
    moutput = ['colored']
    (dxobject_out,) = DXCallModule('Color',minput,moutput)
    return dxobject_out

def DXColorBar(colormap,position=None,shape=None,horizontal=None,ticks=None,min=None,max=None,label=None,colors=None,
               annotation=None,labelscale=None,font=None,ticklocations=None,ticklabels=None,background=None):
    minput = {'colormap':colormap}
    moutput = ['colorbar']
    if position is not None:minput['position']=DXVector(position)
    if shape is not None:minput['shape']=DXVector(shape)
    if horizontal is not None:minput['horizontal']=horizontal
    if ticks is not None:minput['ticks']=ticks
    if min is not None:minput['min']=min
    if max is not None:minput['max']=max
    if label is not None:minput['label']=label
    if type(colors)==type([]):colors=DXMakeStringList(colors)
    if colors is not None:minput['colors']=colors
    if type(annotation)==type([]):annotation=DXMakeStringList(annotation)
    if annotation is not None:minput['annotation']=annotation
    if labelscale is not None:minput['labelscale']=labelscale
    if font is not None:minput['font']=font
    if ticklocations is not None:minput['ticklocations']=ticklocations
    if ticklabels is not None:minput['ticklabels']=ticklabels
    (colorbar1,) = DXCallModule('ColorBar',minput,moutput)
    if background is not None:
      if position is not None:minput['position']=DXVector(position)
      if shape is None:
        minput['shape']=DXVector([320,110])
      else:
        minput['shape']=DXVector([shape[0]+20,shape[1]+85])
      minput['colormap']=DXColormap(min=0.,max=1.,ncolors=2,colors=[[0.,0.,0.],[0.,0.,0.]])[0]
      minput['labelscale']=0.
      CBbackground = DXColor(DXCallModule('ColorBar',minput,moutput)[0],background,opacity=1.)
      colorbar=DXCollect([colorbar1,CBbackground])
    else:
      colorbar=colorbar1
      
    return colorbar

def DXColormap(data=None,min=None,max=None,ncolors=None,colorstart=None,colorend=None,colors=None,opacitystart=None,opacityend=None,opacities=None):
    if data is not None:
      if min is None:min=minnd(data)
      if max is None:max=maxnd(data)
    if ncolors is None:ncolors=256
    if colors is None:
      if colorstart is None:colorstart=[0.,0.,1.]
      if colorend is None:colorend=[1.,0.,0.]
      colorstart=array(colorstart)
      colorend=array(colorend)
    if opacities is None:
      if opacitystart is None: opacitystart=1.
      if opacityend is None: opacityend=1.
    dxpositions=[]
    dxcolors = []
    dxopacities = []
    dc = 1./float(ncolors-1)
    dp = (max-min)/(ncolors)
    do = (opacityend-opacitystart)/(ncolors)
    for i in range(ncolors+1):
      if i<ncolors:
        if colors is None:
          alpha = i*dc
          dxcolors.append((1.-alpha)*colorstart+alpha*colorend)
        else:
          dxcolors.append(colors[i])
      dxpositions.append([min+i*dp])
      if opacities is None:
        dxopacities.append([opacitystart+i*do])
      else:
        dxopacities.append([opacities[i]])
    dxpositions=DXVector(dxpositions)
    dxcolors=DXVector(dxcolors)
    dxopacities=DXVector(dxopacities)
    DXReference(dxpositions)
 
    colormap  = DXConstruct(origin=dxpositions,data=dxcolors)
    opacities = DXConstruct(origin=dxpositions,data=dxopacities)

    DXDelete(dxpositions)
    return colormap,opacities
    
def DXConstruct(origin=None,deltas=None,counts=None,data=None):
    minput = {}
    moutput = ['output']
    if origin is not None:minput['origin']=origin
    if deltas is not None:minput['deltas']=deltas
    if counts is not None:minput['counts']=counts
    if data is not None:minput['data']=data
    (dxobject,) = DXCallModule('Construct',minput,moutput)
    return dxobject


def DXGlyph(data,type=None,shape=None,scale=None,ratio=None,min=None,max=None,auto=0):
    minput = {'data':data}
    moutput = ['glyphs']
    if type is not None:minput['type']=type
    if shape is not None:minput['shape']=shape
    if scale is not None:minput['scale']=scale
    if ratio is not None:minput['ratio']=ratio
    if min is not None:minput['min']=min
    if max is not None:minput['max']=max
    if auto:
      (dxobject_out,) = DXCallModule('AutoGlyph',minput,moutput)
    else:
      (dxobject_out,) = DXCallModule('Glyph',minput,moutput)
    return dxobject_out

def DXAutoGlyph(data,type=None,shape=None,scale=None,ratio=None,min=None,max=None,auto=0):
    return DXGlyph(data,type,shape,scale,ratio,min,max,auto=1)

def DXLight(where=[0.,0.,1.],color=[1.,1.,1.],camera=0):
    if not camera:
      where=DXVector(where)
    minput = {'where':where,'color':DXVector(color),'camera':camera}
    moutput = ['light']
    (dxobject_out,) = DXCallModule('Light',minput,moutput)
    return dxobject_out

def DXAmbientLight(color):
    minput = {'color':DXVector(color)}
    moutput = ['light']
    (dxobject_out,) = DXCallModule('AmbientLight',minput,moutput)
    return dxobject_out

def DXMark(dxobject,name):
    minput = {'input':dxobject,'name':name}
    moutput = ['output']
    (dxobject_out,) = DXCallModule('Mark',minput,moutput)
    return dxobject_out

def DXUnmark(dxobject,name):
    minput = {'input':dxobject,'name':name}
    moutput = ['output']
    (dxobject_out,) = DXCallModule('Unmark',minput,moutput)
    return dxobject_out

def DXInclude(dxobject,min,max,exclude=0,cull=1,pointwise=0):
    dxmin = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
    DXAddArrayData(dxmin,0,1,array(min,Float32))
    dxmax = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
    DXAddArrayData(dxmax,0,1,array(max,Float32))
    minput = {'data':dxobject,'min':dxmin,'max':dxmax,'exclude':exclude,'cull':cull,'pointwise':pointwise}
    moutput = ['output']
    (dxobject_out,) = DXCallModule('Include',minput,moutput)
    return dxobject_out

def DXExcavate(dxobject,min,max):
    return DXUnmark( DXInclude( DXMark(dxobject,'positions') ,min,max,exclude=1) ,'positions')

def DXTranslate(dxobject,translation):
    minput = {'input':dxobject,'translation':DXVector(translation)}
    moutput = ['output']
    (dxobject_out,) = DXCallModule('Translate',minput,moutput)
    return dxobject_out    

def DXShade(dxobject,shade=None,how=None,specular=None,shininess=None,diffuse=None,ambient=None):
    minput = {'input':dxobject}
    moutput = ['output']
    if shade is not None:minput['shade'] = shade
    if how is not None:minput['how'] = how
    if specular is not None:minput['specular'] = specular
    if shininess is not None:minput['shininess'] = shininess
    if diffuse is not None:minput['diffuse'] = diffuse
    if ambient is not None:minput['ambient'] = ambient
    (dxobject_out,) = DXCallModule('Shade',minput,moutput)
    return dxobject_out    
    
def DXSlab(dxobject,dimension,position,thickness=0):
    minput = {'input':dxobject,'dimension':dimension,'position':position,'thickness':thickness}
    moutput = ['output']
    (dxobject_out,) = DXCallModule('Slab',minput,moutput)
    return dxobject_out    

def DXSlice(dxobject,dimension,position):
    minput = {'input':dxobject,'dimension':dimension,'position':position}
    moutput = ['output']
    (dxobject_out,) = DXCallModule('Slice',minput,moutput)
    return dxobject_out    

def DXDescribe(dxobject,options="all"):
    minput={'object':dxobject,'options':options}
    moutput=[]
    DXCallModule('Describe',minput,moutput)

def DXPrint(dxobject):
    DXReference(dxobject)
    minput={'object':dxobject}
    moutput=[]
    DXCallModule('Print',minput,moutput)

def DXVector(v,vtype=Float32):
    """
Returns DX vector from Numeric array v. The shape of v is [nx,ny,...,r], 
where r is the rank of the vector and nx,ny,... are the dimensions of the 
array of vectors. For a scalar vector, just pass a list [vx,vy,vz] or a 
1-d array. 
    """
    if type(v) is type([]):v=array(v)
    dxv = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,shape(v)[-1])
    DXAddArrayData(dxv,0,product(shape(v)[:-1]),v.astype(vtype))
    return dxv


def DXIntegerList(i):
  if type(i) is type([]):i=array(i)
  dxi = DXNewArray(TYPE_INT,CATEGORY_REAL,0,1)
  DXAddArrayData(dxi,0,shape(i)[0],i.astype(Int32))
  return dxi

def DXParticles(x,y,z,vx,vy,vz):
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

  return dxf

def DXTube(x,y,z,c,diameter=None,ngon=None):
  # --- First combine particle data and create a DX array
  n = len(x)
  p = zeros((n,3),'d')
  p[:,0] = array(x)
  p[:,1] = array(y)
  p[:,2] = array(z)
  dxp = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3)
  DXAddArrayData(dxp,0,n,p.astype(Float32))
  dxc = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,1)
  DXAddArrayData(dxc,0,n,array(c).astype(Float32))

  # --- Create the field
  dxf = DXNewField()
  DXSetComponentValue(dxf,'positions',dxp)
  DXSetComponentValue(dxf,'data',dxc)
  DXEndField(dxf)

  dxf2=DXConstruct(DXVector(p),counts=n,data=dxc)
  minput = {'line':dxf2}
  if diameter is not None:minput['diameter']=diameter
  if ngon is not None:minput['ngon']=ngon
  moutput = ['tube']
  (dxobject,) = DXCallModule('Tube',minput,moutput)
  return dxobject  

def DXCamera(lookto=[0.,0.,0.],lookfrom=[0.,0.,1.],width=100.,resolution=640,aspect=0.75,up=[0,1,0],perspective=0,angle=30.,background='black'):
    minput = {'to':DXVector(lookto),
             'from':DXVector(lookfrom),
              'width':width,       
              'resolution':resolution,
              'aspect':aspect,
              'up':DXVector(up),
              'perspective':perspective,
              'angle':angle,
              'background':background}
    moutput = ['camera']
    (camera,) = DXCallModule('Camera',minput,moutput)
    return camera

def DXAutocamera(dxobject,direction='front',width=100.,resolution=640,aspect=0.75,up=[0,1,0],perspective=0,angle=30.,background='black'):
    if type(dxobject) is type([]):
       dxobject=DXVector(direction)
    else:
       DXReference(dxobject)
       width=dxobject
    if type(direction) is type([]):
       direction=DXVector(direction)
    minput = {'object':dxobject,
              'direction':direction,
              'width':width,       
              'resolution':resolution,
              'aspect':aspect,
              'up':DXVector(up),
              'perspective':perspective,
              'angle':angle,
              'background':background}
    moutput = ['camera']
    (camera,) = DXCallModule('AutoCamera',minput,moutput)
    return camera

def DXCollect(dxobjects):
  minput = {'object':dxobjects[0]}
  n=len(dxobjects)
  moutput=['group']
  (group,) = DXCallModule('Collect',minput,moutput)
  if n>1:
    for i in range(1,n):
      group = DXAppend(group,dxobjects[i])
  return group

def DXCollectold(dxobjects):
  minput = {'object':dxobjects[0]}
  n=len(dxobjects)
  if n>1:
    for i in range(1,n):
      minput['object%g'%i]=dxobjects[i]
  moutput=['group']
  (group,) = DXCallModule('Collect',minput,moutput)
  return group

def DXImage(object,camera=None,name=None,labels=None,hardware=1,scale=None,ticks=3,frame=0,
            adjust=1,grid=0,colors='grey',annotation="all",labelscale=1.,direction='front'):
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
#    minput = {'object':__main__.dxwindow.dxobject,'perspective':__main__.dxperspective,
#              'angle':__main__.dxangle,'direction':direction}
#    moutput = ['camera']
#    (camera,) = DXCallModule('AutoCamera',minput,moutput)
    camera=DXAutocamera(dxobject=__main__.dxwindow.dxobject,
                        perspective=__main__.dxperspective,
                        angle=__main__.dxangle,
                        direction=direction)
  __main__.dxcamera=camera

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
          dxwindow.dxobject = DXAutoAxes(dxwindow.dxobject,__main__.dxcamera,labels,
                                         ticks,frame,adjust,grid,colors,annotation,labelscale) 
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
      if l_init or __main__.l_dxresetcamera:
        if not l_init and __main__.l_dxresetcamera:camera=__main__.dxcamera
        minput['defaultCamera'] = __main__.dxcamera
        minput['resetCamera'] = 1
        __main__.l_dxresetcamera=0 
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
                   (hardware or __main__.dxwindow.l_hardware_acceleration)) or 
                   __main__.l_dxforceupdate or __main__.l_dxresetcamera)
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
      __main__.wgui.dx_timer = wx.wxPyTimer(__main__.wgui.dxinter)
      __main__.wgui.dx_timer.Start(__main__.dxtimerconstant)
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
    camera

  if labels is not None:
    assert (len(labels) == 3),"Length of lables list must be three"
    labels = DXMakeStringList(labels)
    l_delete_camera=1
    DXReference(camera)
    minput = {'input':dxobject,'camera':camera,'labels':labels,
              'ticks':3,'colors':'yellow'}
    moutput = ['axes']
    (dxobject,) = DXCallModule('AutoAxes',minput,moutput)
  else:
    l_delete_camera=0

  if format=="dx" or format=="vrml":
    minput = {'object':dxobject,'name':filename,'format':format}
    moutput = []
    DXCallModule('Export',minput,moutput)
  else:
    minput = {'object':dxobject,'camera':camera}
    moutput = ['image']
    (image,) = DXCallModule('Render',minput,moutput)

    minput = {'image':image,'name':filename}
    if format is not None: minput['format'] = format+' gamma=1.'
    moutput = []
    DXCallModule('WriteImage',minput,moutput)

  if l_delete_camera:
    DXDelete(camera)

def DXPrintCamera():
  DXPrint(__main__.dxwindow.dcamera)
