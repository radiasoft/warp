"""This module contains a class Plarr3d with methods to plot a 3-D array of points connected or not by lines.  Also included is a method to make a color-separated stereoscopic plot. For more info, import Plarr3d and type doc(Plarr3d)."""
plarr3d_version = "$Id: plarr3d.py,v 1.3 2009/01/08 19:27:35 dave Exp $"
def plarr3ddoc():
   import plarr3d
   print plarr3d.__doc__
from numpy import *
from gist import *
true=1; false=0
defaultxoffset=0.
defaultyoffset=0.
defaultstereooffset=.1
# a dummy array for type testing:
dumarr=arange(2.)
lcolor4blk=[200,0,0]
rcolor4blk=[0,0,255]

class Plarr3d:
    """
Class Plarr3d plots 3-D array of points connected or not connected
by lines.
 Usage: create instance of Plarr3d, e.g. p=Plarr3d()
 Optional arguments:
    winz, the distance from the observer to the window onto
      which the 3D object is projected.  This, relative to
      the object's dimensions, will determine the perspective
    objz is z position of origin of local coordinates
    theta is inclination of local z axis relative to plotting z
      (rotation w.r.t. lab x axis)
    phi is rotation of local coord. about local y axis, relative
      to plotting window's x,y.

 Basic plotting function is then p.plot3d, and the (optional)
   auxiliary functions p.makeframe(xmin,xmax,ymin,ymax,zmin,zmax)
   and p.plotstereo(x,y,z,lcolor,rcolor,autoframe,linetype,marker,
     stereooffset):

 p.plot3d(x,y,z,color,autoframe,linetype,marker)

  where x,y,z are arrays of x, y, and z values of a set of points
  (in a coordinate system local to the points).  The default
  is to plot lines connecting these points.  Optional arguments:
    color (default "black") of points and connecting lines
    autoframe (default true): construct a 3-D frame around the
      object which just encloses the object (automatic scaling).
      If set to false, plot3d will use a previously calculated
      frame, if one exists; otherwise autoframe is set to true.
      The scaling factor can be calculated from arbitrary
      user-supplied limits using p.getscalefac.
    linetype (default 1): connect the points if 1; display markers at
      points if 0.
    marker (default "\1"): marker to display at points if they are
      not being connected.  "\1" plots dots at points; see gist
      documentation for other possibilties

 p.makeframe(xmin,xmax,ymin,ymax,zmin,zmax):
    Generates a frame with minimum and maximum values as specified.
    If it has never been called, autoframe in plotframe will be set
    to "true".

 p.plotstereo(x,y,z,lcolor,rcolor,autoframe,linetype,marker,
     stereooffset):
    Makes a color-separated stereo plot.
    lcolor is color for left eye (default "cyan")
    rcolor is color for right eye (default "red")
    lcolor and rcolor can be strings, or 3-tuples of rgb values
      (between 0 and 255).
    autoframe, linetype, marker, as for plot3d.
    stereooffset (default 0.1) is how far right eye is displaced
     horizontally (to right, if positive) as fraction of lab-frame
     x extent of object.
    Note the above color defaults work (more or less) for red-blue
     3-D glasses and plots on a white background.  Colors
     should be reversed, or stereooffset set negative, for
     black background.  lcolor4blk and rcolor4blk are provided
     in class definitions to provide some defaults in this case.

    """

    def __init__(self,winz=5.,objz=15.,theta=10.,phi=0.):
        """
plots to window 1 (scaled) unit by 1 unit.
winz is z position of window, relative to eye
objz is z position of origin of local coordinates
theta is inclination of local z axis relative to plotting z
    (rotation w.r.t. lab x axis)
These become object attributes and can be reset at any time
   (e.g. p.phi = [newvalue])
phi is rotation of local coord. about local y axis, relative
    to plotting window's x,y.
        """
        self.winz=winz
        self.objz=objz
        self.theta=theta
        self.phi=phi
        self.framelab=zeros([12,2,3],"d")
        self.frame=zeros([12,2,3],"d")
        self.frame0=zeros([12,2,3],"d")
        self.scaleframe=zeros([12,2,3],"d")
        self.framecolor="black"
        self.calledmakeframe = false

    def setorientation(self):
        " sets orientation from theta and phi "
        thetarad=pi*self.theta/180.
        phirad=self.phi*pi/180.
        self.costh = cos(thetarad)
        self.sinth = sin(thetarad)
        self.cosphi=cos(phirad)
        self.sinphi = sin(phirad)

    def convertolab(self,x,y,z):
        "converts local coordinates to lab coordinates"
        # first do rotation by phi about y axis:
        x1=x*self.cosphi-z*self.sinphi
        z1=x*self.sinphi+z*self.cosphi
        y1=y
        # now rotate by theta about x axis:
        ylab = y1*self.costh+z1*self.sinth
        zlab = z1*self.costh - y1*self.sinth + self.objz
        xlab=x1
        return (xlab,ylab,zlab)

    def makeframe(self,xmin,xmax,ymin,ymax,zmin,zmax):
        """
make a collection of twelve lines which are the border of the frame in
the local space.  And get the lab frame limits of the frame.
        """
        self.calledmakeframe = "true"
        # set orientation when a frame is made.
        self.setorientation()
        self.frame=zeros([12,2,3],"d")
        self.frame[0,0,:]=array([xmin,ymin,zmin])
        self.frame[0,1,:]=array([xmax,ymin,zmin])
        self.frame[1,0,:]=array([xmin,ymin,zmin])
        self.frame[1,1,:]=array([xmin,ymax,zmin])
        self.frame[2,0,:]=array([xmin,ymin,zmin])
        self.frame[2,1,:]=array([xmin,ymin,zmax])
        self.frame[3,0,:]=array([xmax,ymax,zmax])
        self.frame[3,1,:]=array([xmin,ymax,zmax])
        self.frame[4,0,:]=array([xmax,ymax,zmax])
        self.frame[4,1,:]=array([xmax,ymin,zmax])
        self.frame[5,0,:]=array([xmax,ymax,zmax])
        self.frame[5,1,:]=array([xmax,ymax,zmin])
        self.frame[6,0,:]=array([xmax,ymin,zmin])
        self.frame[6,1,:]=array([xmax,ymax,zmin])
        self.frame[7,0,:]=array([xmin,ymin,zmax])
        self.frame[7,1,:]=array([xmin,ymax,zmax])
        self.frame[8,0,:]=array([xmin,ymin,zmax])
        self.frame[8,1,:]=array([xmax,ymin,zmax])
        self.frame[9,0,:]=array([xmax,ymin,zmin])
        self.frame[9,1,:]=array([xmax,ymin,zmax])
        self.frame[10,0,:]=array([xmin,ymax,zmin])
        self.frame[10,1,:]=array([xmin,ymax,zmax])
        self.frame[11,0,:]=array([xmin,ymax,zmin])
        self.frame[11,1,:]=array([xmax,ymax,zmin])
        self.makelabframe()
        self.getlablims(xmin,xmax,ymin,ymax,zmin,zmax)

    def makelabframe(self):
        "convert frame to lab coordinates"
        for i in range(12):
          (self.framelab[i,:,0],self.framelab[i,:,1],
                 self.framelab[i,:,2])=self.convertolab(self.frame[i,:,0],
                 self.frame[i,:,1],self.frame[i,:,2])
        
    def getlablims(self,xmin,xmax,ymin,ymax,zmin,zmax):
        "get extrema of laboratory coordinates of frame"
        self.xlmin=min(min(self.framelab[0,:,0]),min(self.framelab[1,:,0]),
                       min(self.framelab[2,:,0]),min(self.framelab[3,:,0]),
                       min(self.framelab[4,:,0]),min(self.framelab[5,:,0]))
        self.ylmin=min(min(self.framelab[0,:,1]),min(self.framelab[1,:,1]),
                       min(self.framelab[2,:,1]),min(self.framelab[3,:,1]),
                       min(self.framelab[4,:,1]),min(self.framelab[5,:,1]))
        self.zlmin=min(min(self.framelab[0,:,2]),min(self.framelab[1,:,2]),
                       min(self.framelab[2,:,2]),min(self.framelab[3,:,2]),
                       min(self.framelab[4,:,2]),min(self.framelab[5,:,2]))
        self.xlmax=max(max(self.framelab[0,:,0]),max(self.framelab[1,:,0]),
                       max(self.framelab[2,:,0]),max(self.framelab[3,:,0]),
                       max(self.framelab[4,:,0]),max(self.framelab[5,:,0]))
        self.ylmax=max(max(self.framelab[0,:,1]),max(self.framelab[1,:,1]),
                       max(self.framelab[2,:,1]),max(self.framelab[3,:,1]),
                       max(self.framelab[4,:,1]),max(self.framelab[5,:,1]))
        self.zlmax=max(max(self.framelab[0,:,2]),max(self.framelab[1,:,2]),
                       max(self.framelab[2,:,2]),max(self.framelab[3,:,2]),
                       max(self.framelab[4,:,2]),max(self.framelab[5,:,2]))


    def projectlab(self,x,y,z,xoffset,yoffset):
        """
Project x,y,z lab arrays onto the window.  Assumes eye is aligned with
middle of object, offset by percentages xoffset, yoffset of
maximum x, y extents
        """
        xmid=.5*(self.xlmin+self.xlmax)
        ymid=.5*(self.ylmin+self.ylmax)
        xoff = xoffset*(self.xlmax-self.xlmin)
        yoff = yoffset*(self.ylmax-self.ylmin)
        xproj=(x-xmid-xoff)*self.winz/z
        yproj=(y-ymid-yoff)*self.winz/z
        return xproj,yproj

    def project(self,x,y,z,xoffset,yoffset):
        """
Project x,y,z local arrays onto the window. Assumes eye is aligned with
middle of window
        """
        (xlab,ylab,zlab)=self.convertolab(x,y,z)
        (xproj,yproj)=self.projectlab(xlab,ylab,zlab,xoffset,yoffset)
        return (xproj,yproj)

    def plot3d(self,x,y,z,color="black",autoframe=true,linetype=1,marker="\1",
               xoffset=defaultxoffset,yoffset=defaultyoffset):
        """
The plotting function.  See doc(Plarr3d) for explanation of arguments
        """
        #  if autoframe=false, call makeframe manually
        #  to make a frame.  Will set autoscale=true if no frame has
        #  ever been made for this instance; otherwise will use
        #  the last frame calculated.

        if not self.calledmakeframe:
            print "WARNING: no frame pre-calculated; setting autoframe = true"
        if self.calledmakeframe and not autoframe :
            print "WARNING: autoframe = false, working with last frame"
        if  autoframe or not self.calledmakeframe:
            print "about to make frame"
            xmin=min(x)
            xmax=max(x)
            ymin=min(y)
            ymax=max(y)
            zmin=min(z)
            zmax=max(z)
            # make the frame
            self.makeframe(xmin,xmax,ymin,ymax,zmin,zmax)

        # plot the frame
        pldefault(marks=0)
        for i in range(12):
#            (xfp,yfp)=self.project(self.framelab[i,:,0],
#                  self.framelab[i,:,1],self.framelab[i,:,2])
            (xfp,yfp)=self.project(self.frame[i,:,0],
                  self.frame[i,:,1],self.frame[i,:,2],xoffset,yoffset)
            plg(yfp,xfp,color=self.framecolor)
        (xp,yp) = self.project(x,y,z,xoffset,yoffset)
        plg(yp,xp,color=color,type=linetype,marker=marker)
#      make scales have unit ratio
        limits(square=1)
#      turn off frame and ticks
        gridxy(0x200)
        
    def plotstereo(self,x,y,z,lcolor="cyan",rcolor="red",autoframe=true,
                   linetype=1,marker="\1", stereooffset=defaultstereooffset):
       """
makes stereo (2-color separation) plots of arrays.  Arguments
as in plot3d, except lcolor and rcolor (defaults "cyan" and "red")
are colors for left and right eye images, and xoffset (defaulted
to module's stereooffset value) should be set to percentage horizontal
shift of right eye viewpoint relative to object's lab frame x extent
Note, lcolor and rcolor can be strings of standard colors, or
three-tuples of rgb values (range 0-255).
       """
       # if specified rgb values for lcolor or rcolor, create
       # a custom palette with the rgb values
       lcoloruse=lcolor;rcoloruse=rcolor
       if type(lcolor) == type((1,2)) or type(lcolor) == type([1,2]) \
         or type(lcolor) == type(dumarr) \
         or type(rcolor) == type((1,2)) or type(rcolor) == type([1,2]) \
         or type(rcolor) == type(dumarr):
          (lcoloruse,rcoloruse) = makepalette(lcolor,rcolor)
       origframecolor=self.framecolor
       #plot right frame
       self.framecolor=rcoloruse


       self.plot3d(x,y,z,color=rcoloruse,autoframe=autoframe,
          linetype=linetype,marker=marker,xoffset=stereooffset,yoffset=.01)
       #plot left frame
       self.framecolor=lcoloruse
       self.plot3d(x,y,z,color=lcoloruse,autoframe=autoframe,
          linetype=linetype,marker=marker,xoffset=0.)
       self.framecolor=origframecolor

def makepalette(lcolor,rcolor):
# make a custom palette if lcolor or rcolor are tuples
# The color variable needs to be set to lcolor or rcolor if
# they are strings, but set to 0 or 1 if lcolor or rcolor
# are 3-tuples.  (And 1 is only used for rcolor if both
# lcolor and rcolor are specified as 3-tuples; otherwise
# which ever of lcolor and rcolor is set to a 3-tuple, the
# corresponding color variable is set to 0)
   haveleftrgb = 0
   redarr=[];greenarr=[];bluearr=[]
   if type(lcolor) is type("a"):
      lcoloruse=lcolor
   else:
      haveleftrgb=1
      lcoloruse=0
      redarr.append(lcolor[0])
      greenarr.append(lcolor[1])
      bluearr.append(lcolor[2])
   if type(rcolor) is type("a"):
      rcoloruse=rcolor
   else:
      rcoloruse=haveleftrgb
#     so now rcoloruse will be 1 if both lcolor and rcolor are
#     3-tuples; 0 otherwise.
      redarr.append(rcolor[0])
      greenarr.append(rcolor[1])
      bluearr.append(rcolor[2])
   palette(redarr,greenarr,bluearr)
   return (lcoloruse,rcoloruse)
