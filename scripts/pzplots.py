"""
RMS:
  pzxrms: Plots RMS X versus Z
  pzyrms: Plots RMS Y versus Z
  pzzrms: Plots RMS Z versus Z
  pzrrms: Plots RMS R versus Z
  pzxprms: Plots RMS X' versus Z
  pzyprms: Plots RMS Y' versus Z
  pzvxrms: Plots true RMS Vx versus Z
  pzvyrms: Plots true RMS Vy versus Z
  pzvzrms: Plots true RMS Vz versus Z
  pzepsx: Plots x-X' emittance versus Z
  pzepsy: Plots y-Y' emittance versus Z
  pzepsz: Plots z-Z' emittance versus Z
  pzepsnx: Plots x-X' normalized emittance versus Z
  pzepsny: Plots y-Y' normalized emittance versus Z
  pzepsnz: Plots z-Z' normalized emittance versus Z
  pzepsg: Plots generalized emittance versus Z
  pzepsh: Plots generalized emittance versus Z
  pzepsng: Plots generalized normalized emittance versus Z
  pzepsnh: Plots generalized normalized emittance versus Z
  pzxxpslope: Plots slope of x-x' phase space versus Z
  pzyypslope: Plots slope of y-y' phase space versus Z
Envelope:
  pzenvx: Plots beam X envelope (twice Xrms) versus Z
  pzenvy: Plots beam Y envelope (twice Yrms) versus Z
  pzxedge: Plots beam X envelope (twice Xrms) versus Z
  pzxpedge: Plots beam X' envelope versus Z
  pzyedge: Plots beam Y envelope (twice Yrms) versus Z
  pzypedge: Plots beam Y' envelope versus Z
  pzredge: Plots beam R envelope (root 2 Rrms) versus Z
  pzxedges: Plots beam X edges (centroid +- twice Xrms) versus Z
  pzyedges: Plots beam Y edges (centroid +- twice Yrms) versus Z
  pzredges: Plots beam R edges (+- root 2 Rrms) versus Z
  pzenvxp: Plots beam X' envelope (2*xxpbar/xrms) versus Z
  pzenvyp: Plots beam Y' envelope (2*yypbar/yrms) versus Z
Means:
  pzxbar: Plots mean X coordinate versus Z
  pzybar: Plots mean Y coordinate versus Z
  pzzbar: Plots mean axial location versus Z
  pzxpbar: Plots mean X' versus Z
  pzypbar: Plots mean Y' versus Z
  pzvxbar: Plots mean Vx versus Z
  pzvybar: Plots mean Vy versus Z
  pzvzbar: Plots mean Vz versus Z
  pzxybar: Plots mean product of X  and Y  versus Z
  pzxypbar: Plots mean product of X  and Y' versus Z
  pzyxpbar: Plots mean product of Y  and X' versus Z
  pzxpypbar: Plots mean product of X' and Y' versus Z
  pzxsqbar: Plots mean X-squared versus Z
  pzysqbar: Plots mean Y-squared versus Z
  pzzsqbar: Plots mean Z-squared versus Z
  pzxpsqbar: Plots mean X' squared versus Z
  pzypsqbar: Plots mean Y' squared versus Z
  pzvxsqbar: Plots mean Vx squared versus Z
  pzvysqbar: Plots mean Vy squared versus Z
  pzvzsqbar: Plots mean Vz squared versus Z
  pzxxpbar: Plots mean product of X and X' versus Z
  pzyypbar: Plots mean product of Y and Y' versus Z
  pzzvzbar: Plots mean product of Z and Vz versus Z
  pzxvzbar: Plots mean product of X and Vz versus Z
  pzyvzbar: Plots mean product of Y and Vz versus Z
  pzvxvzbar: Plots mean product of Vx and Vz versus Z
  pzvyvzbar: Plots mean product of Vy and Vz versus Z
Miscellaneous:
  pzcurr: Plots beam current versus Z
  pzlchg: Plots line charge versus Z
  pzvzofz: Plots mean axial velocity versus Z
  pzrhomid: Plots charge dens. on axis versus Z
  pzrhomax: Plots charge dens. max-over-X,Y versus Z
  pzrhoax: Plots charge density on axis versus Z
  pzphiax: Plots electrostatic potential on axis versus Z
  pzegap: Plots gap electric field versus Z
  pzezax: Plots Z electric field on axis versus Z
  pzpnum: Plots no. of simulation particles versus Z
  pzppcell: Plots no. of simulation particles per cell versus Z
"""

from warp import *
import __main__
pzplots_version = "$Id: pzplots.py,v 1.21 2008/05/01 23:56:17 dave Exp $"

def pzplotsdoc():
  import pzplots
  print pzplots.__doc__

def setzdiagsflag(flag):
  "Turns on or off the various z diagnostics"
  w3d.lsrhoax3d = flag
  w3d.lgtlchg3d = flag
  w3d.lgetvzofz = flag
  w3d.lgetese3d = flag
  w3d.lsphiax3d = flag
  w3d.lsezax3d  = flag
  w3d.lsetcurr  = flag

###########################################################################
def _extractvar(name,varsuffix=None,pkg='top',ff=None):
  """
Helper function which, given a name, returns the appropriate data. Note that
name could actually be the variable itself, in which case, it is just
returned.
  """
  if type(name) == StringType:
    # --- if varsuffix is specified, try to evaluate the name with the
    # --- suffix. If ok, return the result, otherwise, default to the
    # --- fortran variable in the specified package.
    if varsuffix is not None:
      vname = name + str(varsuffix)
      try:    result = ff.read(vname)
      except: result = None
      if result is not None: return result
      try:    result = __main__.__dict__[vname]
      except: result = None
      if result is not None: return result
    try:    result = ff.read(name+'@'+pkg)
    except: result = None
    if result is not None: return result
    return getattr(packageobject(pkg),name)
  else:
    return name

def _extractvarkw(name,kw,pkg='top'):
  return _extractvar(name,kw.get('varsuffix',None),pkg=pkg)

def _gettitler(js):
  if js == -1: return "All species"
  else:        return "Species %d"%js

##########################################################################
def pzpnum(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots pnumz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  pnumz = _extractvar('pnumz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(pnumz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("No. of simulation particles versus Z",titleb,"(number)",
            _gettitler(js))

##########################################################################
def pzppcell(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots number of particles per cell versus z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  pnumz = _extractvar('pnumz',varsuffix,'top',ff)[...,js]
  xrmsz = _extractvar('xrmsz',varsuffix,'top',ff)[...,js]
  yrmsz = _extractvar('yrmsz',varsuffix,'top',ff)[...,js]
  rrmsz = _extractvar('rrmsz',varsuffix,'top',ff)[...,js]
  dx = _extractvar('dx',varsuffix,'w3d',ff)
  dy = _extractvar('dy',varsuffix,'w3d',ff)
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  if w3d.ny > 0:
    beamarea = 4.*pi*xrmsz*yrmsz
    beamarea = where(beamarea==0.,1.,beamarea)
    ppcell = pnumz/(beamarea/(dx*dy))*scale
    if w3d.l2symtry: ppcell = 2.*ppcell
    if w3d.l4symtry: ppcell = 4.*ppcell
  else:
    beamradius = sqrt(2.)*rrmsz
    beamradius = where(beamradius==0.,1.,beamradius)
    ppcell = pnumz/(beamradius/dx)*scale
  plg(ppcell,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("No. of particles per cell versus Z",titleb,"(number)",
            _gettitler(js))

##########################################################################
def pzxbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xbarz = _extractvar('xbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean X coordinate versus Z",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzybar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots ybarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  ybarz = _extractvar('ybarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(ybarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Y coordinate versus Z",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzzbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots zbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  zbarz = _extractvar('zbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(zbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean axial location versus Z",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzxpbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xpbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xpbarz = _extractvar('xpbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xpbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean X' versus Z",titleb,"(rad)",
            _gettitler(js))

##########################################################################
def pzypbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots ypbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  ypbarz = _extractvar('ypbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(ypbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Y' versus Z",titleb,"(rad)",
            _gettitler(js))

##########################################################################
def pzvxbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vxbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vxbarz = _extractvar('vxbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vxbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Vx versus Z",titleb,"(m/s)",
            _gettitler(js))

##########################################################################
def pzvybar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vybarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vybarz = _extractvar('vybarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vybarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Vy versus Z",titleb,"(m/s)",
            _gettitler(js))

##########################################################################
def pzvzbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vzbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vzbarz = _extractvar('vzbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vzbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Vz versus Z",titleb,"(m/s)",
            _gettitler(js))

##########################################################################
def pzxybar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xybarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xybarz = _extractvar('xybarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xybarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of X  and Y  versus Z",titleb,"(m^2)",
            _gettitler(js))

##########################################################################
def pzxypbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xypbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xypbarz = _extractvar('xypbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xypbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of X  and Y' versus Z",titleb,"(m-rad)",
            _gettitler(js))

##########################################################################
def pzyxpbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots yxpbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yxpbarz = _extractvar('yxpbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(yxpbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of Y  and X' versus Z",titleb,"(m-rad)",
            _gettitler(js))

##########################################################################
def pzxpypbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xpypbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xpypbarz = _extractvar('xpypbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xpypbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of X' and Y' versus Z",titleb,"(rad^2)",
            _gettitler(js))

##########################################################################
def pzxsqbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xsqbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xsqbarz = _extractvar('xsqbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xsqbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean X-squared versus Z",titleb,"(m^2)",
            _gettitler(js))

##########################################################################
def pzysqbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots ysqbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  ysqbarz = _extractvar('ysqbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(ysqbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Y-squared versus Z",titleb,"(m^2)",
            _gettitler(js))

##########################################################################
def pzzsqbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots zsqbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  zsqbarz = _extractvar('zsqbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(zsqbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Z-squared versus Z",titleb,"(m^2)",
            _gettitler(js))

##########################################################################
def pzxpsqbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xpsqbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xpsqbarz = _extractvar('xpsqbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xpsqbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean X' squared versus Z",titleb,"(rad^2)",
            _gettitler(js))

##########################################################################
def pzypsqbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots ypsqbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  ypsqbarz = _extractvar('ypsqbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(ypsqbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Y' squared versus Z",titleb,"(rad^2)",
            _gettitler(js))

##########################################################################
def pzvxsqbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vxsqbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vxsqbarz = _extractvar('vxsqbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vxsqbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Vx squared versus Z",titleb,"((m/s)^2)",
            _gettitler(js))

##########################################################################
def pzvysqbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vysqbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vysqbarz = _extractvar('vysqbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vysqbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Vy squared versus Z",titleb,"((m/s)^2)",
            _gettitler(js))

##########################################################################
def pzvzsqbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vzsqbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vzsqbarz = _extractvar('vzsqbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vzsqbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean Vz squared versus Z",titleb,"((m/s)^2)",
            _gettitler(js))

##########################################################################
def pzxxpbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xxpbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xxpbarz = _extractvar('xxpbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xxpbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of X and X' versus Z",titleb,"(m-rad)",
            _gettitler(js))

##########################################################################
def pzyypbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots yypbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yypbarz = _extractvar('yypbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(yypbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of Y and Y' versus Z",titleb,"(m-rad)",
            _gettitler(js))

##########################################################################
def pzzvzbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots zvzbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  zvzbarz = _extractvar('zvzbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(zvzbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of Z and Vz versus Z",titleb,"(m^2/s)",
            _gettitler(js))

##########################################################################
def pzxvzbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xvzbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xvzbarz = _extractvar('xvzbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xvzbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of X and Vz versus Z",titleb,"(m^2/s)",
            _gettitler(js))

##########################################################################
def pzyvzbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots yvzbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yvzbarz = _extractvar('yvzbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(yvzbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of Y and Vz versus Z",titleb,"(m^2/s)",
            _gettitler(js))

##########################################################################
def pzvxvzbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vxvzbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vxvzbarz = _extractvar('vxvzbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vxvzbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of Vx and Vz versus Z",titleb,"((m/s)^2)",
            _gettitler(js))

##########################################################################
def pzvyvzbar(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vyvzbarz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vyvzbarz = _extractvar('vyvzbarz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vyvzbarz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Mean product of Vy and Vz versus Z",titleb,"((m/s)^2)",
            _gettitler(js))

##########################################################################
def pzxrms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xrmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xrmsz = _extractvar('xrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("RMS X versus Z",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzyrms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots yrmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yrmsz = _extractvar('yrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(yrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("RMS Y versus Z",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzzrms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots zrmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  zrmsz = _extractvar('zrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(zrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("RMS Z versus Z",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzrrms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots rrmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  rrmsz = _extractvar('rrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(rrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("RMS R versus Z",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzxprms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots xprmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xprmsz = _extractvar('xprmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xprmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("RMS X' versus Z",titleb,"(rad)",
            _gettitler(js))

##########################################################################
def pzyprms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots yprmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yprmsz = _extractvar('yprmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(yprmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("RMS Y' versus Z",titleb,"(rad)",
            _gettitler(js))

##########################################################################
def pzepsx(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epsxz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epsxz = _extractvar('epsxz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epsxz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("X-X' emittance versus Z",titleb,"(!p-m-rad)",
            _gettitler(js))

##########################################################################
def pzepsy(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epsyz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epsyz = _extractvar('epsyz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epsyz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Y-Y' emittance versus Z",titleb,"(!p-m-rad)",
            _gettitler(js))

##########################################################################
def pzepsz(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epszz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epszz = _extractvar('epszz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epszz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Z-Z' emittance versus Z",titleb,"(!p-m-rad)",
            _gettitler(js))

##########################################################################
def pzepsnx(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epsnxz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epsnxz = _extractvar('epsnxz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epsnxz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("X-X' normalized emittance versus Z",titleb,"(!p-mm-mrad)",
            _gettitler(js))

##########################################################################
def pzepsny(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epsnyz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epsnyz = _extractvar('epsnyz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epsnyz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Y-Y' normalized emittance versus Z",titleb,"(!p-mm-mrad)",
            _gettitler(js))

##########################################################################
def pzepsnz(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epsnzz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epsnzz = _extractvar('epsnzz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epsnzz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Z-Z' normalized emittance versus Z",titleb,"(!p-mm-mrad)",
            _gettitler(js))

##########################################################################
def pzepsg(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epsgz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epsgz = _extractvar('epsgz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epsgz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Generalized emittance versus Z",titleb,"(!p-m-rad)",
            _gettitler(js))

##########################################################################
def pzepsh(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epshz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epshz = _extractvar('epshz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epshz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Generalized emittance versus Z",titleb,"(!p-m-rad)",
            _gettitler(js))

##########################################################################
def pzepsng(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epsngz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epsngz = _extractvar('epsngz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epsngz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
     ptitles("Generalized normalized emittance versus Z",titleb,"(!p-mm-mrad)",
             _gettitler(js))

##########################################################################
def pzepsnh(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots epsnhz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  epsnhz = _extractvar('epsnhz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(epsnhz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
     ptitles("Generalized normalized emittance versus Z",titleb,"(!p-mm-mrad)",
             _gettitler(js))

##########################################################################
def pzvxrms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vxrmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vxrmsz = _extractvar('vxrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vxrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("True RMS Vx versus Z",titleb,"(m/s)",
            _gettitler(js))

##########################################################################
def pzvyrms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vyrmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vyrmsz = _extractvar('vyrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vyrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("True RMS Vy versus Z",titleb,"(m/s)",
            _gettitler(js))

##########################################################################
def pzvzrms(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots vzrmsz along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vzrmsz = _extractvar('vzrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(vzrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("True RMS Vz versus Z",titleb,"(m/s)",
            _gettitler(js))

##########################################################################
def pzxxpslope(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots slope of x-x' phase space versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xxpbarz = _extractvar('xxpbarz',varsuffix,'top',ff)[...,js]
  xbarz = _extractvar('xbarz',varsuffix,'top',ff)[...,js]
  xpbarz = _extractvar('xpbarz',varsuffix,'top',ff)[...,js]
  xrmsz = _extractvar('xrmsz',varsuffix,'top',ff)[...,js]
  sxz = (xxpbarz - xbarz*xpbarz)/ \
        where(greater(xrmsz,0.),xrmsz**2,1.)*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(sxz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Slope of x-x' phase space",titleb,"(1)",
            _gettitler(js))

##########################################################################
def pzyypslope(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots slope of y-y' phase space versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yypbarz = _extractvar('yypbarz',varsuffix,'top',ff)[...,js]
  ybarz = _extractvar('ybarz',varsuffix,'top',ff)[...,js]
  ypbarz = _extractvar('ypbarz',varsuffix,'top',ff)[...,js]
  yrmsz = _extractvar('yrmsz',varsuffix,'top',ff)[...,js]
  syz = (yypbarz - ybarz*ypbarz)/ \
        where(greater(yrmsz,0.),yrmsz**2,1.)*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(syz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Slope of y-y' phase space",titleb,"(1)",
            _gettitler(js))

##########################################################################
def pzrhomid(zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots rhomidz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  rhomidz = _extractvar('rhomidz',varsuffix,'top',ff)*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(rhomidz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Charge dens. on axis versus Z",titleb,"(C/m^3)")

##########################################################################
def pzrhomax(zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots rhomaxz along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  rhomaxz = _extractvar('rhomaxz',varsuffix,'top',ff)*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(rhomaxz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Charge dens. max-over-X,Y versus Z",titleb,"(C/m^3)")

##########################################################################
def pzcurr(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots current along z-axis
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zplmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  curr = _extractvar('curr',varsuffix,'top',ff)[...,js]*scale
  zplmesh = _extractvar('zplmesh',varsuffix,'top',ff)
  plg(curr,zoffset+zplmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam Current",titleb,"(Amps)",
            _gettitler(js))
ppcurr = pzcurr

##########################################################################
def pzegap(zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots smeared Ez along z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zplmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  egap = _extractvar('egap',varsuffix,'top',ff)*scale
  zplmesh = _extractvar('zplmesh',varsuffix,'top',ff)
  plg(egap,zoffset+zplmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Gap Electric Field",titleb,"(V/m)")

##########################################################################
def pzlchg(zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots linecharge along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zplmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  linechg = _extractvar('linechg',varsuffix,'top',ff)*scale
  zplmesh = _extractvar('zplmesh',varsuffix,'top',ff)
  plg(linechg,zoffset+zplmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Line Charge",titleb,"(C/m)")
pplchg = pzlchg
pzlinechg = pzlchg

##########################################################################
def pzvzofz(zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots Vz along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zplmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  vzofz = _extractvar('vzofz',varsuffix,'top',ff)*scale
  zplmesh = _extractvar('zplmesh',varsuffix,'top',ff)
  plg(vzofz,zoffset+zplmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Mean Axial Velocity",titleb,"(m/s)")
ppvzofz = pzvzofz

##########################################################################
def pzezax(zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots Self Ez along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zplmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  ezax = _extractvar('ezax',varsuffix,'top',ff)*scale
  zplmesh = _extractvar('zplmesh',varsuffix,'top',ff)
  plg(ezax,zoffset+zplmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Z Electric Field on Axis",titleb,"(V/m)")
ppezax = pzezax

##########################################################################
def pzphiax(zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots electrostatic potential along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zplmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  phiax = _extractvar('phiax',varsuffix,'top',ff)*scale
  zplmesh = _extractvar('zplmesh',varsuffix,'top',ff)
  plg(phiax,zoffset+zplmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Electrostatic Potential on Axis",titleb,"(V)")
  if lframe:
    phiplmin = _extractvar('phiplmin',varsuffix,'top',ff)*scale
    phiplmax = _extractvar('phiplmax',varsuffix,'top',ff)*scale
    zzmin = _extractvar('zzmin',varsuffix,'top',ff)
    zzmax = _extractvar('zzmax',varsuffix,'top',ff)
    if ((phiplmin != 0.0)&(phiplmax == 0.0)):
      limits(zzmin,zzmax,phiplmin)
    elif ((phiplmin == 0.0)&(phiplmax != 0.0)):
      limits(zzmin,zzmax,max(phiax),phiplmax)
    elif ((phiplmin != 0.0)&(phiplmax != 0.0)):
      limits(zzmin,zzmax,phiplmin,phiplmax)
ppphiax = pzphiax

##########################################################################
def pzrhoax(zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots space-charge density along the z-axis
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zplmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  rhoax = _extractvar('rhoax',varsuffix,'top',ff)*scale
  zplmesh = _extractvar('zplmesh',varsuffix,'top',ff)
  plg(rhoax,zoffset+zplmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles: ptitles("Charge Density on Axis",titleb,"(C)")
pprhoax = pzrhoax

##########################################################################
def pzenvx(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam X envelope (twice X rms) versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xrmsz = _extractvar('xrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(2.*xrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam X envelope (2*rms)",titleb,"(m)",
            _gettitler(js))
pzxedge = pzenvx

##########################################################################
def pzxpedge(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam X' envelope versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xxpbarz = _extractvar('xxpbarz',varsuffix,'top',ff)[...,js]
  xbarz = _extractvar('xbarz',varsuffix,'top',ff)[...,js]
  xpbarz = _extractvar('xpbarz',varsuffix,'top',ff)[...,js]
  xrmsz = _extractvar('xrmsz',varsuffix,'top',ff)[...,js]
  xpedgez = (xxpbarz-xbarz*xpbarz)/ \
            where(greater(xrmsz,0.),xrmsz,1.)*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xpedgez,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam X' envelope",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzenvy(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam Y envelope (twice Y rms) versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yrmsz = _extractvar('yrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(2.*yrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam Y envelope (2*rms)",titleb,"(m)",
            _gettitler(js))
pzyedge = pzenvy

##########################################################################
def pzypedge(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam Y' envelope versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yypbarz = _extractvar('yypbarz',varsuffix,'top',ff)[...,js]
  ybarz = _extractvar('ybarz',varsuffix,'top',ff)[...,js]
  ypbarz = _extractvar('ypbarz',varsuffix,'top',ff)[...,js]
  yrmsz = _extractvar('yrmsz',varsuffix,'top',ff)[...,js]
  ypedgez = (yypbarz-ybarz*ypbarz)/ \
            where(greater(yrmsz,0.),yrmsz,1.)*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(ypedgez,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam Y' envelope",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzenvr(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
           marks=0,marker=None,msize=1.,width=1.,lframe=0,
           titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam R envelope (root 2 R rms) versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  rrmsz = _extractvar('rrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(sqrt(2.)*rrmsz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam R envelope (sqrt(2)*rms)",titleb,"(m)",
            _gettitler(js))
pzredge = pzenvr

##########################################################################
def pzxedges(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam X edges (centroid +- twice X rms) versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xbarz = _extractvar('xbarz',varsuffix,'top',ff)[...,js]*scale
  xrmsz = _extractvar('xrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(xbarz+2.*xrmsz,zoffset+zmntmesh/zscale,color=color,
      linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  plg(xbarz-2.*xrmsz,zoffset+zmntmesh/zscale,color=color,
      linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam X edges (xbar+-2*rms)",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzyedges(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam Y edges (centroid +- twice Y rms) versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  ybarz = _extractvar('ybarz',varsuffix,'top',ff)[...,js]*scale
  yrmsz = _extractvar('yrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(ybarz+2.*yrmsz,zoffset+zmntmesh/zscale,color=color,
      linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  plg(ybarz-2.*yrmsz,zoffset+zmntmesh/zscale,color=color,
      linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam Y edges (ybar+-2*rms)",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzredges(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam R edges (+- root 2 R rms) versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  rrmsz = _extractvar('rrmsz',varsuffix,'top',ff)[...,js]*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(+2.*rrmsz,zoffset+zmntmesh/zscale,color=color,
      linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  plg(-2.*rrmsz,zoffset+zmntmesh/zscale,color=color,
      linetype=linetype,marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam R edges (+-sqrt(2)*rms)",titleb,"(m)",
            _gettitler(js))

##########################################################################
def pzenvxp(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam X' envelope (2*xxpbar/xrms) versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  xxpbarz = _extractvar('xxpbarz',varsuffix,'top',ff)[...,js]
  xbarz = _extractvar('xbarz',varsuffix,'top',ff)[...,js]
  xpbarz = _extractvar('xpbarz',varsuffix,'top',ff)[...,js]
  xrmsz = _extractvar('xrmsz',varsuffix,'top',ff)[...,js]
  sxz = 2.*(xxpbarz - xbarz*xpbarz)/ \
        where(greater(xrmsz,0.),xrmsz,1.)*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(sxz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam X' envelope",titleb,"(rad)",
            _gettitler(js))

##########################################################################
def pzenvyp(js=-1,zoffset=0.,zscale=1.,scale=1.,color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
  """Plots beam Y' envelope (2*yypbar/yrms) versus Z
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - zoffset=0: offset added to axis
  - zscale=1: scale of axis
    plots versus zoffset + zmntmesh/zscale
  - scale=1.: factor to scale data by
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="Z": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
  if zscale == 0.: raise "zscale must be nonzero"
  if titleb is None:
    if zscale == 1.: titleb = "Z (m)"
    else: titleb = "Z"
  yypbarz = _extractvar('yypbarz',varsuffix,'top',ff)[...,js]
  ybarz = _extractvar('ybarz',varsuffix,'top',ff)[...,js]
  ypbarz = _extractvar('ypbarz',varsuffix,'top',ff)[...,js]
  yrmsz = _extractvar('yrmsz',varsuffix,'top',ff)[...,js]
  syz = 2.*(yypbarz - ybarz*ypbarz)/ \
        where(greater(yrmsz,0.),yrmsz,1.)*scale
  zmntmesh = _extractvar('zmntmesh',varsuffix,'top',ff)
  plg(syz,zoffset+zmntmesh/zscale,color=color,linetype=linetype,
      marks=marks,marker=marker,msize=msize,width=width)
  if titles:
    ptitles("Beam Y' envelope",titleb,"(rad)",
            _gettitler(js))


##########################################################################
def pzplotstest(**kw):
  """
Test all pzplots routines.
  """
  apply(pzpnum,(),kw);fma()
  apply(pzppcell,(),kw);fma()
  apply(pzxbar,(),kw);fma()
  apply(pzybar,(),kw);fma()
  apply(pzzbar,(),kw);fma()
  apply(pzxpbar,(),kw);fma()
  apply(pzypbar,(),kw);fma()
  apply(pzvxbar,(),kw);fma()
  apply(pzvybar,(),kw);fma()
  apply(pzvzbar,(),kw);fma()
  apply(pzxybar,(),kw);fma()
  apply(pzxypbar,(),kw);fma()
  apply(pzyxpbar,(),kw);fma()
  apply(pzxpypbar,(),kw);fma()
  apply(pzxsqbar,(),kw);fma()
  apply(pzysqbar,(),kw);fma()
  apply(pzzsqbar,(),kw);fma()
  apply(pzxpsqbar,(),kw);fma()
  apply(pzypsqbar,(),kw);fma()
  apply(pzvxsqbar,(),kw);fma()
  apply(pzvysqbar,(),kw);fma()
  apply(pzvzsqbar,(),kw);fma()
  apply(pzxxpbar,(),kw);fma()
  apply(pzyypbar,(),kw);fma()
  apply(pzzvzbar,(),kw);fma()
  apply(pzxvzbar,(),kw);fma()
  apply(pzyvzbar,(),kw);fma()
  apply(pzvxvzbar,(),kw);fma()
  apply(pzvyvzbar,(),kw);fma()
  apply(pzxrms,(),kw);fma()
  apply(pzyrms,(),kw);fma()
  apply(pzzrms,(),kw);fma()
  apply(pzxprms,(),kw);fma()
  apply(pzyprms,(),kw);fma()
  apply(pzepsx,(),kw);fma()
  apply(pzepsy,(),kw);fma()
  apply(pzepsz,(),kw);fma()
  apply(pzepsnx,(),kw);fma()
  apply(pzepsny,(),kw);fma()
  apply(pzepsnz,(),kw);fma()
  apply(pzepsg,(),kw);fma()
  apply(pzepsh,(),kw);fma()
  apply(pzepsng,(),kw);fma()
  apply(pzepsnh,(),kw);fma()
  apply(pzvxrms,(),kw);fma()
  apply(pzvyrms,(),kw);fma()
  apply(pzvzrms,(),kw);fma()
  apply(pzxxpslope,(),kw);fma()
  apply(pzyypslope,(),kw);fma()
  apply(pzrhomid,(),kw);fma()
  apply(pzrhomax,(),kw);fma()
  apply(pzcurr,(),kw);fma()
  apply(pzegap,(),kw);fma()
  apply(pzlchg,(),kw);fma()
  apply(pzvzofz,(),kw);fma()
  apply(pzezax,(),kw);fma()
  apply(pzphiax,(),kw);fma()
  apply(pzrhoax,(),kw);fma()
  apply(pzenvx,(),kw);fma()
  apply(pzenvy,(),kw);fma()
  apply(pzxedge,(),kw);fma()
  apply(pzxpedge,(),kw);fma()
  apply(pzyedge,(),kw);fma()
  apply(pzypedge,(),kw);fma()
  apply(pzxedges,(),kw);fma()
  apply(pzyedges,(),kw);fma()
  apply(pzenvxp,(),kw);fma()
  apply(pzenvyp,(),kw);fma()
