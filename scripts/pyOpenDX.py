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

pyOpenDX_version = "$Id: pyOpenDX.py,v 1.3 2003/02/24 16:13:41 dave Exp $"
def pyOpenDXdoc():
  import pyOpenDX
  print pyOpenDX.__doc__

###########################################################################
def ppxxpy(iw = 0,**kw):
  """Plots X-Y-Z"""
  checkparticleplotarguments(kw)
  ii = selectparticles(iw=iw,kwdict=kw)
  viewparticles(take(top.xp,ii),(take(top.uxp,ii)/take(top.uzp,ii)),
                take(top.yp,ii),(take(top.uyp,ii)/take(top.uzp,ii)),
                l1='X',l2="X'",l3='Y',name='WARP viz')
  return

def ppxyz(iw = 0,**kw):
  """Plots X-Y-Z"""
  checkparticleplotarguments(kw)
  ii = selectparticles(iw=iw,kwdict=kw)
  viewparticles(take(top.xp,ii),take(top.yp,ii),take(top.zp,ii),
                take(top.uxp,ii),
                l1='X',l2='Y',l3='Z',name='WARP viz')
  return

###########################################################################
def viewisosurface1(data,isovalue,name='WARP viz'):
  f = DXObject_fromarray(data)

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
def viewparticles(x,y,z,v,l1=None,l2=None,l3=None,name='WARP particles'):
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
  (colorglyphs,) = DXCallModule('AutoColor',minput,moutput)
  DXReference(colorglyphs)

  minput = {'object':colorglyphs}
  moutput = ['camera']
  (camera,) = DXCallModule('AutoCamera',minput,moutput)
  DXReference(camera)

  if l1 is not None and l2 is not None and l3 is not None:
    labelsadded = 1
    labels = DXMakeStringList([l1,l2,l3])
    minput = {'input':colorglyphs,'camera':camera,'labels':labels,
              'ticks':3,'colors':'yellow'}
    moutput = ['axes']
    (ob,) = DXCallModule('AutoAxes',minput,moutput)
    DXReference(ob)
  else:
    labelsadded = 0
    ob = colorglyphs

  DXImage(ob,camera,name)
  DXDelete(ob)
  DXDelete(camera)
  if labelsadded:
    DXDelete(colorglyphs)

###########################################################################
#==========================================================================
interactor = 0
def interactor_handler():
  global interactor
  getchar = sys.stdin.readline()[:-1]
  if len(getchar) > 0: interactor = eval(getchar)
  else:                interactor = -1

_group = [None]
_groupn = [0]
def DXImage(object,camera,name='WARP viz'):
  global interactor

  if _group[0] is None:
    minput = {'object':object}
    moutput = ['group']
    (_group[0],) = DXCallModule('Collect',minput,moutput)
    DXReference(_group[0])
  else:
    minput = {'input':_group[0],'object':object}
    moutput = ['group']
    (g,) = DXCallModule('Append',minput,moutput)
    DXDelete(_group[0])
    _group[0] = g
    DXReference(_group[0])

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
                'object':_group[0],'mode':interactor,'resetObject':1}
      if i == 1:
        minput['defaultCamera'] = camera
        minput['resetCamera'] = 1
      moutput = ['object','camera','where']
      (dobject,dcamera,dwhere,) = DXCallModule('SuperviseState',minput,moutput)

      minput = {'object':dobject,'camera':dcamera,'where':dwhere}
      moutput = ['where']
      (dwhere,) = DXCallModule('Display',minput,moutput)

  DXDelete(dwhere)
  interactor = 0

def DXNewImage():
  DXDelete(_group[0])
  _group[0] = None
  _groupn[0] = 0

# This may not be the best thing to do, but it works - it gives access to this
# function without having to explicitly import this module. Note that this
# only works for interactive sessions.
__main__.__dict__['DXNewImage'] = DXNewImage

