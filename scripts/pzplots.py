from warp import *
pzplots_version = "$Id: pzplots.py,v 1.8 2002/01/24 23:07:30 dave Exp $"

def pzplotsdoc():
  print """
pzpnum: Plot no. of simulation particles versus Z
pzppcell: Plot no. of simulation particles per cell versus Z
pzxbar: Plot mean X coordinate versus Z
pzybar: Plot mean Y coordinate versus Z
pzzbar: Plot mean axial location versus Z
pzxpbar: Plot mean X' versus Z
pzypbar: Plot mean Y' versus Z
pzvxbar: Plot mean Vx versus Z
pzvybar: Plot mean Vy versus Z
pzvzbar: Plot mean Vz versus Z
pzxybar: Plot mean product of X  and Y  versus Z
pzxypbar: Plot mean product of X  and Y' versus Z
pzyxpbar: Plot mean product of Y  and X' versus Z
pzxpypbar: Plot mean product of X' and Y' versus Z
pzxsqbar: Plot mean X-squared versus Z
pzysqbar: Plot mean Y-squared versus Z
pzzsqbar: Plot mean Z-squared versus Z
pzxpsqbar: Plot mean X' squared versus Z
pzypsqbar: Plot mean Y' squared versus Z
pzvxsqbar: Plot mean Vx squared versus Z
pzvysqbar: Plot mean Vy squared versus Z
pzvzsqbar: Plot mean Vz squared versus Z
pzxxpbar: Plot mean product of X and X' versus Z
pzyypbar: Plot mean product of Y and Y' versus Z
pzzvzbar: Plot mean product of Z and Vz versus Z
pzxvzbar: Plot mean product of X and Vz versus Z
pzyvzbar: Plot mean product of Y and Vz versus Z
pzvxvzbar: Plot mean product of Vx and Vz versus Z
pzvyvzbar: Plot mean product of Vy and Vz versus Z
pzxrms: Plot rMS X versus Z
pzyrms: Plot rMS Y versus Z
pzzrms: Plot rMS Z versus Z
pzxprms: Plot rMS X' versus Z
pzyprms: Plot rMS Y' versus Z
pzepsx: Plot x-X' emittance versus Z
pzepsy: Plot y-Y' emittance versus Z
pzepsz: Plot z-Z' emittance versus Z
pzepsnx: Plot x-X' normalized emittance versus Z
pzepsny: Plot y-Y' normalized emittance versus Z
pzepsnz: Plot z-Z' normalized emittance versus Z
pzepsg: Plot generalized emittance versus Z
pzepsh: Plot generalized emittance versus Z
pzepsng: Plot generalized normalized emittance versus Z
pzepsnh: Plot generalized normalized emittance versus Z
pzvxrms: Plot true RMS Vx versus Z
pzvyrms: Plot true RMS Vy versus Z
pzvzrms: Plot true RMS Vz versus Z
pzxxpslope: Plot slope of x-x' phase space versus Z
pzyypslope: Plot slope of y-y' phase space versus Z
pzrhomid: Plot charge dens. on axis versus Z
pzrhomax: Plot charge dens. max-over-X,Y versus Z
pzcurr: Plot beam current versus Z
pzegap: Plot gap electric field versus Z
pzlchg: Plot line charge versus Z
pzvzofz: Plot mean axial velocity versus Z
pzezax: Plot Z electric field on axis versus Z
pzphiax: Plot electrostatic potential on axis versus Z
pzrhoax: Plot charge density on axis versus Z
pzenvx: Plot beam X envelope (twice Xrms) versus Z
pzenvy: Plot beam Y envelope (twice Yrms) versus Z
pzxedge: Plot beam X envelope (twice Xrms) versus Z
pzxpedge: Plot beam X' envelope versus Z
pzyedge: Plot beam Y envelope (twice Yrms) versus Z
pzypedge: Plot beam Y' envelope versus Z
pzxedges: Plot beam X edges (centroid +- twice Xrms) versus Z
pzyedges: Plot beam Y edges (centroid +- twice Yrms) versus Z
pzenvxp: Plot beam X' envelope (2*xxpbar/xrms) versus Z
pzenvyp: Plot beam Y' envelope (2*yypbar/yrms) versus Z
  """

##########################################################################
def pzpnum(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots pnumz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.pnumz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("No. of simulation particles versus Z",titleb,"(number)")

##########################################################################
def pzppcell(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
             marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots number of particles per cell versus z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  ppcell = top.pnumz/(pi*top.xrmsz*top.yrmsz/(w3d.dx*w3d.dy))
  warpplg(ppcell,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("No. of particles per cell versus Z",titleb,"(number)")

##########################################################################
def pzxbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean X coordinate versus Z",titleb,"(m)")

##########################################################################
def pzybar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots ybarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.ybarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Y coordinate versus Z",titleb,"(m)")

##########################################################################
def pzzbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots zbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.zbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean axial location versus Z",titleb,"(m)")

##########################################################################
def pzxpbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xpbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xpbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean X' versus Z",titleb,"(rad)")

##########################################################################
def pzypbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots ypbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.ypbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Y' versus Z",titleb,"(rad)")

##########################################################################
def pzvxbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vxbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vxbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Vx versus Z",titleb,"(m/s)")

##########################################################################
def pzvybar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vybarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vybarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Vy versus Z",titleb,"(m/s)")

##########################################################################
def pzvzbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vzbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vzbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Vz versus Z",titleb,"(m/s)")

##########################################################################
def pzxybar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xybarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xybarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of X  and Y  versus Z",titleb,"(m^2)")

##########################################################################
def pzxypbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xypbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xypbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of X  and Y' versus Z",titleb,"(m-rad)")

##########################################################################
def pzyxpbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots yxpbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.yxpbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of Y  and X' versus Z",titleb,"(m-rad)")

##########################################################################
def pzxpypbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xpypbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xpypbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of X' and Y' versus Z",titleb,"(rad^2)")

##########################################################################
def pzxsqbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xsqbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xsqbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean X-squared versus Z",titleb,"(m^2)")

##########################################################################
def pzysqbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots ysqbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.ysqbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Y-squared versus Z",titleb,"(m^2)")

##########################################################################
def pzzsqbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots zsqbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.zsqbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Z-squared versus Z",titleb,"(m^2)")

##########################################################################
def pzxpsqbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xpsqbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xpsqbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean X' squared versus Z",titleb,"(rad^2)")

##########################################################################
def pzypsqbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots ypsqbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.ypsqbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Y' squared versus Z",titleb,"(rad^2)")

##########################################################################
def pzvxsqbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vxsqbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vxsqbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Vx squared versus Z",titleb,"((m/s)^2)")

##########################################################################
def pzvysqbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vysqbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vysqbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Vy squared versus Z",titleb,"((m/s)^2)")

##########################################################################
def pzvzsqbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vzsqbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vzsqbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Vz squared versus Z",titleb,"((m/s)^2)")

##########################################################################
def pzxxpbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xxpbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xxpbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of X and X' versus Z",titleb,"(m-rad)")

##########################################################################
def pzyypbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots yypbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.yypbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of Y and Y' versus Z",titleb,"(m-rad)")

##########################################################################
def pzzvzbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots zvzbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.zvzbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of Z and Vz versus Z",titleb,"(m^2/s)")

##########################################################################
def pzxvzbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xvzbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xvzbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of X and Vz versus Z",titleb,"(m^2/s)")

##########################################################################
def pzyvzbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots yvzbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.yvzbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of Y and Vz versus Z",titleb,"(m^2/s)")

##########################################################################
def pzvxvzbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vxvzbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vxvzbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of Vx and Vz versus Z",titleb,"((m/s)^2)")

##########################################################################
def pzvyvzbar(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vyvzbarz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vyvzbarz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean product of Vy and Vz versus Z",titleb,"((m/s)^2)")

##########################################################################
def pzxrms(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xrmsz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xrmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("RMS X versus Z",titleb,"(m)")

##########################################################################
def pzyrms(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots yrmsz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.yrmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("RMS Y versus Z",titleb,"(m)")

##########################################################################
def pzzrms(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots zrmsz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.zrmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("RMS Z versus Z",titleb,"(m)")

##########################################################################
def pzxprms(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots xprmsz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xprmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("RMS X' versus Z",titleb,"(rad)")

##########################################################################
def pzyprms(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots yprmsz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.yprmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("RMS Y' versus Z",titleb,"(rad)")

##########################################################################
def pzepsx(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epsxz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epsxz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("X-X' emittance versus Z",titleb,"(pi-m-rad)")

##########################################################################
def pzepsy(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epsyz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epsyz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Y-Y' emittance versus Z",titleb,"(pi-m-rad)")

##########################################################################
def pzepsz(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epszz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epszz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Z-Z' emittance versus Z",titleb,"(pi-m-rad)")

##########################################################################
def pzepsnx(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epsnxz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epsnxz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("X-X' normalized emittance versus Z",titleb,"(pi-mm-mrad)")

##########################################################################
def pzepsny(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epsnyz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epsnyz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Y-Y' normalized emittance versus Z",titleb,"(pi-mm-mrad)")

##########################################################################
def pzepsnz(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epsnzz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epsnzz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Z-Z' normalized emittance versus Z",titleb,"(pi-mm-mrad)")

##########################################################################
def pzepsg(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epsgz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epsgz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Generalized emittance versus Z",titleb,"(pi-m-rad)")

##########################################################################
def pzepsh(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epshz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epshz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Generalized emittance versus Z",titleb,"(pi-m-rad)")

##########################################################################
def pzepsng(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epsngz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epsngz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Generalized normalized emittance versus Z",titleb,"(pi-mm-mrad)")

##########################################################################
def pzepsnh(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots epsnhz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.epsnhz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Generalized normalized emittance versus Z",titleb,"(pi-mm-mrad)")

##########################################################################
def pzvxrms(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vxrmsz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vxrmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("True RMS Vx versus Z",titleb,"(m/s)")

##########################################################################
def pzvyrms(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vyrmsz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vyrmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("True RMS Vy versus Z",titleb,"(m/s)")

##########################################################################
def pzvzrms(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots vzrmsz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vzrmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("True RMS Vz versus Z",titleb,"(m/s)")

##########################################################################
def pzxxpslope(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
               marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plot slope of x-x' phase space versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  sxz = (top.xxpbarz - top.xbarz*top.xpbarz)/ \
        where(greater(top.xrmsz,0.),top.xrmsz**2,1.)
  warpplg(sxz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Slope of x-x' phase space",titleb,"(1)")

##########################################################################
def pzyypslope(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
               marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plot slope of y-y' phase space versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  syz = (top.yypbarz - top.ybarz*top.ypbarz)/ \
        where(greater(top.yrmsz,0.),top.yrmsz**2,1.)
  warpplg(syz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Slope of y-y' phase space",titleb,"(1)")

##########################################################################
def pzrhomid(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots rhomidz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.rhomidz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Charge dens. on axis versus Z",titleb,"(C/m^3)")

##########################################################################
def pzrhomax(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots rhomaxz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.rhomaxz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Charge dens. max-over-X,Y versus Z",titleb,"(C/m^3)")

##########################################################################
def pzcurr(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots current along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zplmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.curr,zoffset+top.zplmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam Current",titleb,"(Amps)")
ppcurr = pzcurr

##########################################################################
def pzegap(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots smeared Ez along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zplmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.egap,zoffset+top.zplmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Gap Electric Field",titleb,"(V/m)")

##########################################################################
def pzlchg(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots linecharge along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zplmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.linechg,zoffset+top.zplmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Line Charge",titleb,"(C/m^2)")
pplchg = pzlchg

##########################################################################
def pzvzofz(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots Vz along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zplmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.vzofz,zoffset+top.zplmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Axial Velocity",titleb,"(m/s)")
ppvzofz = pzvzofz

##########################################################################
def pzezax(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots Self Ez along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zplmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.ezax,zoffset+top.zplmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Z Electric Field on Axis",titleb,"(V/m)")
ppezax = pzezax

##########################################################################
def pzphiax(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots electrostatic potential along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zplmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.phiax,zoffset+top.zplmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Electrostatic Potential on Axis",titleb,"(V)")
  if lframe:
    if ((top.phiplmin != 0.0)&(top.phiplmax == 0.0)):
      limits(top.zzmin,top.zzmax,top.phiplmin)
    elif ((top.phiplmin == 0.0)&(top.phiplmax != 0.0)):
      limits(top.zzmin,top.zzmax,max(top.phiax),top.phiplmax)
    elif ((top.phiplmin != 0.0)&(top.phiplmax != 0.0)):
      limits(top.zzmin,top.zzmax,top.phiplmin,top.phiplmax)
ppphiax = pzphiax

##########################################################################
def pzrhoax(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots space-charge density along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zplmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.rhoax,zoffset+top.zplmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Charge Density on Axis",titleb,"(C)")
pprhoax = pzrhoax

##########################################################################
def pzenvx(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
           marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots beam X envelope (twice X rms) versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(2.*top.xrmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam X envelope (2*rms)",titleb,"(m)")
pzxedge = pzenvx

##########################################################################
def pzxpedge(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
             marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots beam X' envelope versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xpedgez = (top.xxpbarz-top.xbarz*top.xpbarz)/ \
            where(greater(top.xrmsz,0.),top.xrmsz,1.)
  warpplg(xpedgez,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam X' envelope",titleb,"(m)")

##########################################################################
def pzenvy(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
           marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots beam Y envelope (twice Y rms) versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(2.*top.yrmsz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam Y envelope (2*rms)",titleb,"(m)")
pzyedge = pzenvy

##########################################################################
def pzypedge(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
             marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots beam Y' envelope versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  ypedgez = (top.yypbarz-top.ybarz*top.ypbarz)/ \
            where(greater(top.yrmsz,0.),top.yrmsz,1.)
  warpplg(ypedgez,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam Y' envelope",titleb,"(m)")

##########################################################################
def pzxedges(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
             marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots beam X edges (centroid +- twice X rms) versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.xbarz+2.*top.xrmsz,zoffset+top.zmntmesh/zscale,color=color,
          linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  warpplg(top.xbarz-2.*top.xrmsz,zoffset+top.zmntmesh/zscale,color=color,
          linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam X edges (xbar+-2*rms)",titleb,"(m)")

##########################################################################
def pzyedges(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
             marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plots beam Y edges (centroid +- twice Y rms) versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  warpplg(top.ybarz+2.*top.yrmsz,zoffset+top.zmntmesh/zscale,color=color,
          linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  warpplg(top.ybarz-2.*top.yrmsz,zoffset+top.zmntmesh/zscale,color=color,
          linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam Y edges (ybar+-2*rms)",titleb,"(m)")

##########################################################################
def pzenvxp(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plot beam X' envelope (2*xxpbar/xrms) versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  sxz = 2.*(top.xxpbarz - top.xbarz*top.xpbarz)/ \
        where(greater(top.xrmsz,0.),top.xrmsz,1.)
  warpplg(sxz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam X' envelope",titleb,"(rad)")

##########################################################################
def pzenvyp(zoffset=0.,zscale=1.,color="fg",linetype="solid",marks=0,
            marker=None,msize=1.,width=1.,lframe=0,titleb=None,titles=1):
  """Plot beam Y' envelope (2*yypbar/yrms) versus Z
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + top.zmntmesh/zscale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles"""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  syz = 2.*(top.yypbarz - top.ybarz*top.ypbarz)/ \
        where(greater(top.yrmsz,0.),top.yrmsz,1.)
  warpplg(syz,zoffset+top.zmntmesh/zscale,color=color,linetype=linetype,
          marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Beam Y' envelope",titleb,"(rad)")

