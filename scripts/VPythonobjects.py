"""Utility routines for doing 3-d visualization using the VPython module.
Partly taken from VPython demos.
Modified by DPG

VisualMesh: can plot 3-D surfaces corresponding to meshed data.
"""
from warp import *
VPythonobjects_version = "$Id: VPythonobjects.py,v 1.6 2004/04/16 22:59:29 dave Exp $"

def VPythonobjectsdoc():
  import VPythonobjects
  print VPythonobjects.__doc__

##########################################################################

class VisualModel:
  def __init__(self,twoSided=1,scene=None,title='Visualization',vrange=None,
                    viewer='OpenDX'):
    self.viewer = viewer
    self.triangles = []
    self.colors = []
    self.normals = []
    self.connections = None
    self.twoSided = twoSided  # add every face twice with opposite normals
    if vrange is not None:
      self.vrange = vrange
      self.autoscale = 0
      self.uniform = 1
    else:
      self.vrange = [10.,10.,10.]
      self.autoscale = 1
      self.uniform = 1
    self.title = title
    self.scene = scene

  def Display(self,showgrid=0,showaxes=0):

    if self.viewer == 'VPython':
      import visual

      if self.scene is None:
        self.scene = visual.display(exit=0,width=500,height=500,
                                    uniform=self.uniform,title=self.title,
                                    autoscale=self.autoscale,range=self.vrange)

      self.frame = visual.frame(display=self.scene)
      if not self.colors: self.colors=None
      self.model = visual.faces(frame=self.frame,pos=self.triangles,
                                normal=self.normals,color=self.colors,
                                display=self.scene)
    else:
      import pyOpenDX

      n = len(self.triangles)
      p = array(self.triangles).astype(Float32)
      ss = pyOpenDX.DXNewArray(pyOpenDX.TYPE_FLOAT,pyOpenDX.CATEGORY_REAL,1,3)
      pyOpenDX.DXAddArrayData(ss,0,n,p)

      normals = array(self.normals)
      nn = pyOpenDX.DXNewArray(pyOpenDX.TYPE_FLOAT,pyOpenDX.CATEGORY_REAL,1,3)
      pyOpenDX.DXAddArrayData(nn,0,n,normals.astype(Float32))
      pyOpenDX.DXSetStringAttribute(nn,'dep','positions')

      if not self.colors: colors = array(n*[[0.5,0.7,1.0]])
      else:               colors = array(self.colors)
      co = pyOpenDX.DXNewArray(pyOpenDX.TYPE_FLOAT,pyOpenDX.CATEGORY_REAL,1,3)
      pyOpenDX.DXAddArrayData(co,0,n,colors.astype(Float32))
      pyOpenDX.DXSetStringAttribute(co,'dep','positions')

      if self.connections is None:
        self.connections = arange(n)
        self.connections.shape = (nint(n/3),3)
      cc = pyOpenDX.DXNewArray(pyOpenDX.TYPE_INT,pyOpenDX.CATEGORY_REAL,1,3)
      pyOpenDX.DXAddArrayData(cc,0,nint(n/3),self.connections.astype(Int))
      pyOpenDX.DXSetStringAttribute(cc,'ref','positions')
      pyOpenDX.DXSetStringAttribute(cc,'element type','triangles')
      pyOpenDX.DXSetStringAttribute(cc,'dep','connections')

      ff = pyOpenDX.DXNewField()
      pyOpenDX.DXSetComponentValue(ff,'positions',ss)
      pyOpenDX.DXSetComponentValue(ff,'colors',co)
      pyOpenDX.DXSetComponentValue(ff,'normals',nn)
      pyOpenDX.DXSetComponentValue(ff,'connections',cc)
      pyOpenDX.DXEndField(ff)

      if showgrid:

        # --- Create a field containing two points on opposite side of the
        # --- mesh and use it to find the boinding box.
        n = 2
        p = zeros((n,3),'d')
        p[:,0] = [w3d.xmmin,w3d.xmmax]
        p[:,1] = [w3d.xmmin,w3d.xmmax]
        p[:,2] = [w3d.xmmin,w3d.xmmax]
        dxp = pyOpenDX.DXNewArray(pyOpenDX.TYPE_FLOAT,pyOpenDX.CATEGORY_REAL,1,3)
        pyOpenDX.DXAddArrayData(dxp,0,n,p.astype(Float32))
        # --- Create a DX array for the data
        dxd = pyOpenDX.DXNewArray(pyOpenDX.TYPE_DOUBLE,pyOpenDX.CATEGORY_REAL,0)
        pyOpenDX.DXAddArrayData(dxd,0,n,array([1.,1.]))
        pyOpenDX.DXSetStringAttribute(dxd,'dep','positions')
        # --- Create the field
        dxf = pyOpenDX.DXNewField()
        pyOpenDX.DXSetComponentValue(dxf,'positions',dxp)
        pyOpenDX.DXSetComponentValue(dxf,'data',dxd)
        pyOpenDX.DXEndField(dxf)

        minput = {'input':dxf}
        moutput = ['box']
        (box,) = pyOpenDX.DXCallModule('ShowBox',minput,moutput)

        minput = {'object':ff,'object1':box}
        moutput = ['group']
        (group,) = pyOpenDX.DXCallModule('Collect',minput,moutput)

      else:
        group = ff

      pyOpenDX.DXReference(group)
      minput = {'object':group}
      moutput = ['camera']
      (camera,) = pyOpenDX.DXCallModule('AutoCamera',minput,moutput)
      pyOpenDX.DXReference(camera)

      if showaxes:
        labels = pyOpenDX.DXMakeStringList(['x','y','z'])
        minput = {'input':group,'camera':camera,'labels':labels,
                  'ticks':3,'colors':'yellow'}
        moutput = ['axes']
        (axes,) = pyOpenDX.DXCallModule('AutoAxes',minput,moutput)
        pyOpenDX.DXReference(axes)
      else:
        axes = group

      pyOpenDX.DXImage(axes,camera,self.title)

      pyOpenDX.DXDelete(camera)
      pyOpenDX.DXDelete(axes)

  def FacetedTriangle(self, v1, v2, v3, color=None):
    """Add a triangle to the model, apply faceted shading automatically"""
    normal = self.Norm( self.Cross(v2-v1, v3-v1) )
    for v in (v1,v2,v3):
      self.triangles.append(v)
      if color is not None: self.colors.append(color)
      self.normals.append(normal)
      #self.model.append( pos=v, color=color, normal=normal )
    if self.twoSided:
      for v in (v1,v3,v2):
        #self.model.append( pos=v, color=color, normal=-normal )
        self.triangles.append(v)
        if color is not None: self.colors.append(color)
        self.normals.append(-normal)

  def FacetedPolygon(self, v, color=None):
    """Appends a planar polygon of any number of vertices to the model,
       applying faceted shading automatically."""
    for t in range(len(v)-2):
      self.FacetedTriangle( v[0], v[t+1], v[t+2] ,color=color)

  def DoSmoothShading(self):
    rsq = sum(self.triangles**2,1)
    ii = argsort(rsq)
    # to be completed later

  def DoSmoothShading1(self):
    """Change a faceted model to smooth shaded, by averaging normals at
    coinciding vertices.
    
    This is a very slow and simple smooth shading
    implementation which has to figure out the connectivity of the
    model and does not attempt to detect sharp edges.

    It attempts to work even in two-sided mode where there are two
    opposite normals at each vertex.  It may fail somehow in pathological
    cases. """

    pos = self.model.pos
    normal = self.model.normal

    vertex_map = {}  # vertex position -> vertex normal
    vertex_map_backface = {}
    for i in range( len(pos) ):
      tp = tuple(pos[i])
      old_normal = vertex_map.get( tp, (0,0,0) )
      if dot(old_normal, normal[i]) >= 0:
        vertex_map[tp] = normal[i] + old_normal
      else:
        vertex_map_backface[tp] = normal[i] + vertex_map_backface.get(tp, (0,0,0))

    for i in range( len(pos) ):
      tp = tuple(pos[i])
      if dot(vertex_map[tp], normal[i]) >= 0:
        normal[i] = vertex_map[tp] and self.Norm( vertex_map[ tp ] )
      else:
        normal[i] = vertex_map_backface[tp] and self.Norm(vertex_map_backface[tp] )

  def DrawNormal(self, scale):
    pos = self.model.pos
    normal = self.model.normal
    for i in range(len(pos)):
      arrow(pos=pos[i], axis=normal[i]*scale)

  def Cross(self,v1,v2):
    return array([v1[1]*v2[2] - v1[2]*v2[1],
                  v1[2]*v2[0] - v1[0]*v2[2],
                  v1[0]*v2[1] - v1[1]*v2[0]])

  def Norm(self,v):
    magv = sqrt(sum(v**2))
    if magv != 0.: return v/magv
    else:          return v

########################################################################
class VisualMesh (VisualModel):
  """
xvalues, yvalues, zvalues: 2-D arrays containing the coordinates and data
twoSided=1: when true, surface is two sided
color=None: can be specified as an [r,g,b] list
scene=None: an already existing display scene. When None, create a new one.
title='Mesh': Display title - only used when new scene created.
  """
  def __init__(self, xvalues=None, yvalues=None, zvalues=None,
               xscaled=0,zscaled=1,
               twoSided=1,color=None,color1=None,color2=None,
               scene=None,title=None,vrange=None,viewer=None):
    if not title: title = 'Mesh'
    VisualModel.__init__(self,twoSided=twoSided,scene=scene,title=title,
                              vrange=vrange,viewer=viewer)

    assert zvalues is not None,"zvalues must be specified"

    s = shape(zvalues)
    if len(s) != 2:
      print 'First argument must be a 2-Dimensional array'
      return
    if xvalues is None:
      xvalues = arange(s[0])[:,NewAxis]*ones(s[1],'d')
    elif len(shape(xvalues))==1:
      xvalues = xvalues[:,NewAxis]*ones(s[1],'d')
    if yvalues is None:
      yvalues = arange(s[1])*ones(s[0],'d')[:,NewAxis]
    elif len(shape(yvalues))==1:
      yvalues = yvalues*ones(s[0],'d')[:,NewAxis]

    if xscaled:
      xrange = maxnd(xvalues) - minnd(xvalues)
      xvalues = xvalues/xrange
      yvalues = yvalues/xrange
    if zscaled:
      xrange = maxnd(xvalues) - minnd(xvalues)
      zrange = maxnd(zvalues) - minnd(zvalues)
      zvalues = zvalues/zrange*xrange/2.

    if color1 is not None:
      getcolor = 1
      zmin = minnd(zvalues)
      zmax = maxnd(zvalues)
      if color2 is None: color2 = zeros(3,'d')
    else:
      getcolor = 0

    points = zeros( xvalues.shape + (3,), Float )
    points[...,0] = xvalues
    points[...,1] = yvalues
    points[...,2] = zvalues

    for i in range(zvalues.shape[0]-1):
      for j in range(zvalues.shape[1]-1):
        if getcolor:
          color = color1 + color2*(points[i,j,2] - zmin)/(zmax - zmin)
        self.FacetedPolygon([points[i,j], points[i,j+1],
                             points[i+1,j+1], points[i+1,j]],
                            color=color)
    self.Display()

########################################################################
class VisualRevolution(VisualModel):
  """
Visualize surface of revolution
  """
  def __init__(self,srfrv,zzmin,zzmax,nz=20,nth=20,xoff=0,yoff=0,zoff=0,
                    twoSided=1,color=None,color1=None,color2=None,
                    scene=None,title=None,vrange=None):
    if not title: title = 'Surface of revolution'
    VisualModel.__init__(self,twoSided=twoSided,scene=scene,title=title,
                              vrange=vrange)

    zz = arange(0,nz+1)*(zzmax - zzmin)/nz + zzmin
    rr = ones(nz+1,'d')
    for i in range(nz+1):
      warp.f3d.srfrv_z = zz[i]
      srfrv()
      rr[i] = warp.f3d.srfrv_r

    xx = cos(2.*pi*arange(0,nth+1)/nth) + xoff
    yy = sin(2.*pi*arange(0,nth+1)/nth) + yoff

    if color1 is not None:
      getcolor = 1
      zmin = minnd(zvalues)
      zmax = maxnd(zvalues)
      if color2 is None: color2 = zeros(3,'d')
    else:
      getcolor = 0

    for i in xrange(nz-1):
      for j in xrange(len(xx)-1):
        p1 = array([rr[i  ]*xx[j  ], rr[i  ]*yy[j  ], zz[i  ]+zoff])
        p2 = array([rr[i  ]*xx[j+1], rr[i  ]*yy[j+1], zz[i  ]+zoff])
        p3 = array([rr[i+1]*xx[j+1], rr[i+1]*yy[j+1], zz[i+1]+zoff])
        p4 = array([rr[i+1]*xx[j  ], rr[i+1]*yy[j  ], zz[i+1]+zoff])
        self.FacetedPolygon([p1,p2,p3,p4],color=color)
    self.Display()


