plarr3d_version = "$Id: plarr3d.py,v 1.1 2001/12/06 21:16:21 rcohen Exp $"
def plarr3ddoc():
   print """
This module contains a class Plarr3d with methods to plot a 3-D array
of points conneted or not by lines.  For more info, import Plarr3d and type
doc(Plarr3d).
   """
from Numeric import *
from gist import *
true=1; false=0

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
   auxiliary function p.makeframe(xmin,xmax,ymin,ymax,zmin,zmax)

 p.plot3d(x,y,z,color,autoframe,type,marker)

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
    type (default 1): connect the points if 1; display markers at
      points if 0.
    marker (default "\1"): marker to display at points if they are
      not being connected.  "\1" plots dots at points; see gist
      documentation for other possibilties

 p.makeframe(xmin,xmax,ymin,ymax,zmin,zmax):
    Generates a frame with minimum and maximum values as specified.
    If it has never been called, autoframe in plotframe will be set
    to "true".

    """

    def __init__(self,winz=5.,objz=15.,theta=10.,phi=0.):
        """
plots to window 1 (scaled) unit by 1 unit.
winz is z position of window, relative to eye
objz is z position of origin of local coordinates
theta is inclination of local z axis relative to plotting z
    (rotation w.r.t. lab x axis)
phi is rotation of local coord. about local y axis, relative
    to plotting window's x,y.
        """
        self.winz=winz
        self.objz=objz
        thetarad=pi*theta/180.
        phirad=phi*pi/180.
        self.costh = cos(thetarad)
        self.sinth = sin(thetarad)
        self.cosphi=cos(phirad)
        self.sinphi = sin(phirad)
        self.framelab=zeros([12,2,3],"d")
        self.frame=zeros([12,2,3],"d")
        self.frame0=zeros([12,2,3],"d")
        self.scaleframe=zeros([12,2,3],"d")
        self.framecolor="black"
        self.calledmakeframe = false


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
        self.getlablims(xmin,xmax,ymin,ymax,zmin,zmax)
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


    def projectlab(self,x,y,z):
        """
Project x,y,z lab arrays onto the window.  Assumes eye is aligned with
middle of object
        """
        xmid=.5*(self.xlmin+self.xlmax)
        self.ymid=.5*(self.ylmin+self.ylmax)
        xproj=(x-xmid)*self.winz/z
        yproj=(y-self.ymid)*self.winz/z
        return xproj,yproj

    def project(self,x,y,z):
        """
Project x,y,z local arrays onto the window Assumes eye is aligned with
middle of window
        """
        (xlab,ylab,zlab)=self.convertolab(x,y,z)
        (xproj,yproj)=self.projectlab(xlab,ylab,zlab)
        return (xproj,yproj)

    def plot3d(self,x,y,z,color="black",autoframe=true,type=1,marker="\1"):
        """
The plotting function.  See doc(Plarr3d) for explanation of arguments
        """
        # if autoframe=false, call getscalefac manually
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
                  self.frame[i,:,1],self.frame[i,:,2])
            plg(yfp,xfp,color=self.framecolor)
        (xp,yp) = self.project(x,y,z)
        plg(yp,xp,color=color,type=type,marker=marker)
#      make scales have unit ratio
        limits(square=1)
#      turn off frame and ticks
        gridxy(0x200)
        
