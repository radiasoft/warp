"""This module contains a class Plarr3d with methods to plot a 3-D array of points connected or not by lines.  Also included is a method to make a color-separated stereoscopic plot. For more info, import Plarr3d and type doc(Plarr3d)."""
def plarr3ddoc():
    import plarr3d
    print plarr3d.__doc__
from numpy import *
from gist import *
import collections
true=1; false=0
defaultxoffset=0.
defaultyoffset=0.
defaultstereooffset=.1
# label offsets as fraction of total frame size
xlabeloffset = -0.05
yzlabeloffset = -0.03
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
     stereooffset,framelabels,labelprecision):

 p.plot3d(x,y,z,color,autoframe,linetype,marker,xoffset,yoffset,framelabels,labelprecision)

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
    xoffset, yoffset (default 0) are horizontal and vertical offsets of the eye relative to
      the object's origin.   This is provided mainly for use in plotstereo (3D stereo plotting).
    framelabels (default true): if true, labels "x", "y" and "z" axes and the plotted limits
    labelprecision (default 2): number of significant digits in the printed axis limits

 p.makeframe(xmin,xmax,ymin,ymax,zmin,zmax):
    Generates a frame with minimum and maximum values as specified.
    If it has never been called, autoframe in plotframe will be set
    to "true".

 p.plotstereo(x,y,z,lcolor,rcolor,autoframe,linetype,marker,
     stereooffset,framelabels,labelprecision):
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

     Note these plotting routines turn off the default plotting frame.  To turn it back on
     (to make conventional plots) you can invoke the convenience method
     restoreframe() after fma() and before your next plot.

    """

    def __init__(self,winz=5.,objz=30.,theta=-10.,phi=0.):
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
        x1=x*self.cosphi+z*self.sinphi
        z1=-x*self.sinphi+z*self.cosphi
        y1=y
        # now rotate by theta about x axis:
        ylab = y1*self.costh+z1*self.sinth
        #zlab = z1*self.costh - y1*self.sinth + self.objz
        zlab = z1*self.costh - y1*self.sinth
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
        self.frame[3,1,:]=array([xmax,ymax,zmax])
        self.frame[3,0,:]=array([xmin,ymax,zmax])
        self.frame[4,1,:]=array([xmax,ymax,zmax])
        self.frame[4,0,:]=array([xmax,ymin,zmax])
        self.frame[5,1,:]=array([xmax,ymax,zmax])
        self.frame[5,0,:]=array([xmax,ymax,zmin])
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


    def projectlab(self,x,y,z):
        """
        Project x,y,z lab arrays onto the window.  Assumes eye is aligned with
        middle of object, offset by percentages self.xoffset, self.yoffset of
        maximum x, y extents
        """
        xmid=.5*(self.xlmin+self.xlmax)
        ymid=.5*(self.ylmin+self.ylmax)
        xoff = self.xoffset*(self.xlmax-self.xlmin)
        yoff = self.yoffset*(self.ylmax-self.ylmin)
        ztoeye = self.objz-z
        xproj=(x-xmid-xoff)*self.winz/ztoeye
        yproj=(y-ymid-yoff)*self.winz/ztoeye
        return xproj,yproj

    def project(self,x,y,z):
        """
        Project x,y,z local arrays onto the window. Assumes eye is aligned with
        middle of window
        """
        (xlab,ylab,zlab)=self.convertolab(x,y,z)
        (xproj,yproj)=self.projectlab(xlab,ylab,zlab)
        return (xproj,yproj)

    def outsidelines(self):
        """
        Determine which frame edges should get labels, i.e. ones that appear on the
        outside, if possible, in the projection onto the window,  to the bottom
        or left of the box for x and y axes.   To avoid possible label collision,
        pick z boundary on max x side, top or bottom.
        Depending on orienation, the outside x boundary is frame edge 8 or 0;
        the outside y boundary is 1 or 7 if z increases toward viewer, otherwise
        4 or 6; and the outside z boundary is 5 or 9 or neither
        (e.g. if we are looking directly into the high-z
        end of the frame), in which case we arbitrarily pick side 9) for z increasing
        toward the viewer; otherwise z boundary is 2 or 10.
        """
        # Test to see if side 0 appears below side 8
        x0 = self.frame[0,0,0];y0=self.frame[0,0,1];z0 = self.frame[0,0,2]
        x8 = self.frame[8,0,0];y8=self.frame[8,0,1];z8=self.frame[8,0,2]
        if (self.project(x0,y0,z0)[1] < self.project(x8,y8,z8)[1]):
            self.xlabelside = 0
        else:
            self.xlabelside = 8

        # Test to see if z is increasing toward the observer or not
        self.zincr_tord_obs = (sign(self.cosphi)+1)/2.
        if self.zincr_tord_obs:
            # Test to see if side 1 appears to left of  side 7
            x1 = self.frame[1,0,0];y1=self.frame[1,0,1];z1 = self.frame[1,0,2]
            x7 = self.frame[7,0,0];y7=self.frame[7,0,1];z7=self.frame[7,0,2]
            if (self.project(x1,y1,z1)[0] < self.project(x7,y7,z7)[0]):
                self.ylabelside = 1
            else:
                self.ylabelside = 7
        else:
            # Test to see if side 4 appears to left of  side 6
            x4 = self.frame[4,0,0];y4=self.frame[4,0,1];z4 = self.frame[4,0,2]
            x6 = self.frame[6,0,0];y6=self.frame[6,0,1];z6=self.frame[6,0,2]
            if (self.project(x4,y4,z4)[0] < self.project(x6,y6,z6)[0]):
                self.ylabelside = 4
            else:
                self.ylabelside = 6

        if self.zincr_tord_obs:
            # Test to see if upper end of side 10 appears above and to right, or below and
            # to left, of lower end.
            x10u = self.frame[10,1,0];y10u =self.frame[10,1,1];z10u = self.frame[10,1,2]
            x10l = self.frame[10,0,0];y10l =self.frame[10,0,1];z10l = self.frame[10,0,2]
            if ((self.project(x10u,y10u,z10u)[0] < self.project(x10l,y10l,z10l)[0]) and \
                (self.project(x10u,y10u,z10u)[1] < self.project(x10l,y10l,z10l)[1])) or \
                ((self.project(x10u,y10u,z10u)[0] > self.project(x10l,y10l,z10l)[0]) and \
                 (self.project(x10u,y10u,z10u)[1] > self.project(x10l,y10l,z10l)[1])):
                self.zlabelside = 10
            else:
                self.zlabelside = 2
        else:
            # Test to see if upper end of side 5 appears above and to right, or below and
            # to left, of lower end
            x5u = self.frame[5,1,0];y5u =self.frame[5,1,1];z5u = self.frame[5,1,2]
            x5l = self.frame[5,0,0];y5l =self.frame[5,0,1];z5l = self.frame[5,0,2]
            if ((self.project(x5u,y5u,z5u)[0] < self.project(x5l,y5l,z5l)[0]) and \
                (self.project(x5u,y5u,z5u)[1] < self.project(x5l,y5l,z5l)[1])) or \
                ((self.project(x5u,y5u,z5u)[0] > self.project(x5l,y5l,z5l)[0]) and \
                 (self.project(x5u,y5u,z5u)[1] > self.project(x5l,y5l,z5l)[1])):
                self.zlabelside = 5
            else:
                self.zlabelside = 9


    def printframelabel(self,labelside):
        """
        prints extrema and axis label on a side of a frame
        """
        labelprecision = self.labelprecision
        extraspacebeg=0.0;extraspaceend=0.0;extraspacelabel=0.0
        if (labelside == 0) or (labelside == 8):
            # x axis
            axislabel = "x"
            offsets = array([0,xlabeloffset])
            # create text for label of starting and ending x positions
            begpostxt = "%.*g" %(labelprecision,self.frame[labelside,0,0])
            endpostxt = "%.*g" %(labelprecision,self.frame[labelside,1,0])
            # Note x axis has labels offset vertically
            labeljustify = "CA"
            if self.zincr_tord_obs:
                begjustify = "LA"
                endjustify = "RA"
            else:
                endjustify = "LA"
                begjustify = "RA"
        if (labelside == 1) or (labelside == 7) or (labelside == 4) or (labelside == 6):
            axislabel = "y"
            offsets = array([yzlabeloffset,0])
            # create text for label of starting and ending y positions
            begpostxt = "%.*g" %(labelprecision,self.frame[labelside,0,1])
            endpostxt = "%.*g" %(labelprecision,self.frame[labelside,1,1])
            # Note y axis has labels offset horizontally
            labeljustify = "RH"
            begjustify = "RA"
            endjustify = "RT"
        if (labelside == 5) or (labelside == 9) or (labelside == 2) or (labelside == 10):
            axislabel = "z"
            offsets = array([yzlabeloffset,0])
            # print "z labelside, offsets = ", labelside,offsets
            # if labelside = 5 need a bit of extra space to left of begpostxt and label
            # if labelside = 10 need a bit of extra space to left of endpostxt and label
            if (labelside == 9):
                extraspacelabel = -.01
                extraspacebeg = -.01
            if (labelside == 5):
                extraspacelabel = -.01
                extraspacebeg = -.01
            # create text for label of starting and ending z positions
            begpostxt = "%.*g" %(labelprecision,self.frame[labelside,0,2])
            endpostxt = "%.*g" %(labelprecision,self.frame[labelside,1,2])
            # Note z axis has labels offset horizontally
            labeljustify = "RH"
            if (self.theta < 0 and self.zincr_tord_obs) or (self.theta > 0 and not self.zincr_tord_obs):
                begjustify = "RT"
                endjustify = "RA"
            else:
                endjustify = "RT"
                begjustify = "RA"
        # calculate absolute values of offsets by multiplying by x and y ranges of plotted frame
        offsets=offsets*self.winrange

        # get position, projected onto window,  of middle of the side, to label the axis
        midposx = 0.5*(self.frame[labelside,0,0]+self.frame[labelside,1,0])
        midposy = 0.5*(self.frame[labelside,0,1]+self.frame[labelside,1,1])
        midposz = 0.5*(self.frame[labelside,0,2]+self.frame[labelside,1,2])
        axislabelposition = array(self.project(midposx,midposy,midposz))+offsets
        plt(axislabel,axislabelposition[0]+extraspacelabel,axislabelposition[1],justify=labeljustify,tosys=1,color=self.framecolor)
        #
        # get position, projected onto window, of start and end of the labelled sides
        (xfp,yfp)=self.project(self.frame[labelside,:,0],self.frame[labelside,:,1],  \
                               self.frame[labelside,:,2])
        xfp = xfp + offsets[0]
        yfp = yfp + offsets[1]
        plt(begpostxt,xfp[0]+extraspacebeg,yfp[0],justify=begjustify,tosys=1,color=self.framecolor)
        plt(endpostxt,xfp[1]+extraspaceend,yfp[1],justify=endjustify,tosys=1,color=self.framecolor)
        #
        # make sure the plot has room for the labels
        winrangebar = sum(self.winrange)/2.
        winsafety =0.03*winrangebar
        rawlimits = limits(square=1)
        limits(rawlimits[0]-winsafety,rawlimits[1],rawlimits[2]-winsafety,rawlimits[3])

    def printframelabels(self,labelprecision=2):
        """
        Prints labels on the frame.
        labelprecision is number of digits in printed frame limit labels
        """
        self.labelprecision = labelprecision
        # Determine edges on which to place labels
        self.outsidelines()
        self.printframelabel(self.xlabelside)
        self.printframelabel(self.ylabelside)
        self.printframelabel(self.zlabelside)

    def plot3d(self,x,y,z,color="black",autoframe=true,linetype=1,marker="\1",
        xoffset=defaultxoffset,yoffset=defaultyoffset,framelabels=true,labelprecision=2):
        """
        The plotting function.  See doc(Plarr3d) for explanation of arguments
        """
        #  if autoframe=false, call makeframe manually
        #  to make a frame.  Will set autoscale=true if no frame has
        #  ever been made for this instance; otherwise will use
        #  the last frame calculated.

        self.xoffset = xoffset
        self.yoffset = yoffset
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

        #calculate extrema of frame in window for scaling of plotted frame label offsets
        xfmin=0;xfmax=0;yfmin=0;yfmax=0
        for i in range(12):
            #   (xfp,yfp)=self.project(self.framelab[i,:,0],
            #         self.framelab[i,:,1],self.framelab[i,:,2])
            (xfp,yfp)=self.project(self.frame[i,:,0],
                  self.frame[i,:,1],self.frame[i,:,2])
            xfmin=min(xfmin,xfp[0])
            yfmin=min(yfmin,yfp[0])
            xfmax=max(xfmax,xfp[1])
            yfmax=max(yfmax,yfp[1])
            plg(yfp,xfp,color=self.framecolor)
        #store the range of the frame projected onto the window
        self.winrange=array([xfmax-xfmin,yfmax-yfmin])
        #
        (xp,yp) = self.project(x,y,z)
        plg(yp,xp,color=color,type=linetype,marker=marker)
        #  make scales have unit ratio
        limits(square=1)
        # turn off frame and ticks
        gridxy(0x200)
        # turn on axis labels if requested
        if (framelabels == true):
            self.printframelabels(labelprecision)

    def plotstereo(self,x,y,z,lcolor="cyan",rcolor="red",autoframe=true,
                   linetype=1,marker="\1", stereooffset=defaultstereooffset,
                   framelabels=true,labelprecision=2):
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
        if (isinstance(lcolor,collections.Sequence) or isinstance(lcolor,ndarray)
          or isinstance(rcolor,collections.Sequence) or isinstance(rcolor,ndarray)):
            (lcoloruse,rcoloruse) = makepalette(lcolor,rcolor)
        origframecolor=self.framecolor
        #plot right frame
        self.framecolor=rcoloruse


        self.plot3d(x,y,z,color=rcoloruse,autoframe=autoframe,
           linetype=linetype,marker=marker,xoffset=stereooffset,yoffset=.01,
           framelabels=framelabels,labelprecision=labelprecision)
        #plot left frame
        self.framecolor=lcoloruse
        self.plot3d(x,y,z,color=lcoloruse,autoframe=autoframe,
           linetype=linetype,marker=marker,xoffset=0.,
           framelabels=framelabels,labelprecision=labelprecision)
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
    if isinstance(lcolor,basestring):
        lcoloruse=lcolor
    else:
        haveleftrgb=1
        lcoloruse=0
        redarr.append(lcolor[0])
        greenarr.append(lcolor[1])
        bluearr.append(lcolor[2])
    if isinstance(rcolor,basestring):
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

def restoreframe():
    # restore the default WARP boxed plotting frame
    gridxy(0x62b)
