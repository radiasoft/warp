"""
Partly taken from VPython demos.
Modified by DPG

VisualMesh: can plot 3-D surfaces corresponding to meshed data.
"""
VPythonobjects_version = "$Id: VPythonobjects.py,v 1.2 2003/01/16 20:21:00 dave Exp $"

def VPythonobjectsdoc():
  import VPythonobjects
  print VPythonobjects.__doc__

##########################################################################
from visual import *
import warp

class VisualModel:
    def __init__(self,twoSided=true,scene=None,title='VPython',vrange=None):
        self.triangles = []
        self.colors = []
        self.normals = []
        self.twoSided = twoSided  # add every face twice with opposite normals
        if vrange is not None:
          autoscale = 0
          uniform = 1
        else:
          vrange = [10.,10.,10.]
          autoscale = 1
          uniform = 1
        if scene is None:
          self.scene = display(exit=0,width=500,height=500,
                               uniform=uniform,title=title,
                               autoscale=autoscale,range=vrange)
        else:
          self.scene = scene

    def Display(self):
        self.frame = frame(display=self.scene)
        if not self.colors: self.colors=None
        self.model = faces(frame=self.frame,pos=self.triangles,
                           normal=self.normals,color=self.colors,
                           display=self.scene)

    def FacetedTriangle(self, v1, v2, v3, color=None):
        """Add a triangle to the model, apply faceted shading automatically"""
        try:
            normal = norm( cross(v2-v1, v3-v1) )
        except:
            normal = vector(0,0,0)
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
                normal[i] = vertex_map[tp] and norm( vertex_map[ tp ] )
            else:
                normal[i] = vertex_map_backface[tp] and norm(vertex_map_backface[tp] )

    def DrawNormal(self, scale):
        pos = self.model.pos
        normal = self.model.normal
        for i in range(len(pos)):
            arrow(pos=pos[i], axis=normal[i]*scale)

########################################################################
class VisualMesh (VisualModel):
    """
xvalues, yvalues, zvalues: 2-D arrays containing the coordinates and data
twoSided=true: when true, surface is two sided
color=None: can be specified as an [r,g,b] list
scene=None: an already existing display scene. When None, create a new one.
title='Mesh': Display title - only used when new scene created.
    """
    def __init__(self, xvalues, yvalues, zvalues,
                 twoSided=true,color=None,
                 scene=None,title=None,vrange=None):
        if not title: title = 'Mesh'
        VisualModel.__init__(self,twoSided=twoSided,scene=scene,title=title,
                                  vrange=vrange)

        points = zeros( xvalues.shape + (3,), Float )
        points[...,0] = xvalues
        points[...,1] = yvalues
        points[...,2] = zvalues

        for i in range(zvalues.shape[0]-1):
            for j in range(zvalues.shape[1]-1):
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
                    twoSided=true,color=None,
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

    for i in xrange(nz-1):
      for j in xrange(len(xx)-1):
        p1 = array([rr[i  ]*xx[j  ], rr[i  ]*yy[j  ], zz[i  ]+zoff])
        p2 = array([rr[i  ]*xx[j+1], rr[i  ]*yy[j+1], zz[i  ]+zoff])
        p3 = array([rr[i+1]*xx[j+1], rr[i+1]*yy[j+1], zz[i+1]+zoff])
        p4 = array([rr[i+1]*xx[j  ], rr[i+1]*yy[j  ], zz[i+1]+zoff])
        self.FacetedPolygon([p1,p2,p3,p4],color=color)
    self.Display()


