"""Functions to plot lab window moments data

ppnumlw: Plots pnumlw as a function of time
pxbarlw: Plots xbarlw as a function of time
pybarlw: Plots ybarlw as a function of time
pvzbarlw: Plots vzbarlw as a function of time
pepsxlw: Plots epsxlw as a function of time
pepsylw: Plots epsylw as a function of time
pepzslw: Plots epzslw as a function of time
pvxrmslw: Plots vxrmslw as a function of time
pvyrmslw: Plots vyrmslw as a function of time
pvzrmslw: Plots vzrmslw as a function of time
pxrmslw: Plots xrmslw as a function of time
pyrmslw: Plots yrmslw as a function of time
prrmslw: Plots rrmslw as a function of time
pxxpbarlw: Plots xxpbarlw as a function of time
pyypbarlw: Plots yypbarlw as a function of time
pcurrlw: Plots currlw as a function of time
plostparslw: Plots lostparslw as a function of time
plinechglw: Plots linechglw as a function of time

"""

from warp import *
import __main__
lwplots_version = "$Id: lwplots.py,v 1.2 2009/04/03 21:58:17 dave Exp $"

def lwplotsdoc():
    import lwplots
    print lwplots.__doc__

###########################################################################
def addlabwindow(zlw):
    """Adds a new lab window moments calculation point at the given location."""
    # --- Find first non-zero value of top.zlw
    iz = argmin(abs(top.zlw))
    if top.zlw[iz] == 0.:
        top.zlw[iz] = zlw
    else:
        # --- More space is needed
        top.labwn += 1
        gchange('Lab_Moments')
        top.zlw[-1] = zlw

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

def _gettitler(ilw,js):
    if js == -1: return "All species in lab window %d"%ilw
    else:        return "Species %d in lab window %d"%(js,ilw)

##########################################################################

def ppnumlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots pnumlw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    pnumlw = _extractvar('pnumlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(pnumlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("No. of simulation particles versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pxbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xbarlw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xbarlw = _extractvar('xbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("X bar versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pybarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots ybarlw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    ybarlw = _extractvar('ybarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(ybarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Y bar versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pvzbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vzbarlw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vzbarlw = _extractvar('vzbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vzbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Vz bar versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pepsxlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsxlw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsxlw = _extractvar('epsxlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsxlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("X emittance versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pepsylw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsylw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsylw = _extractvar('epsylw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsylw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Y emittance versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pepzslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epzslw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epzslw = _extractvar('epzslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epzslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Z emittance versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pvxrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vxrmslw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vxrmslw = _extractvar('vxrmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vxrmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Vx rms versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pvyrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vyrmslw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vyrmslw = _extractvar('vyrmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vyrmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Vy rms versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pvzrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vzrmslw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vzrmslw = _extractvar('vzrmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vzrmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Vz rms versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pxrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xrmslw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xrmslw = _extractvar('xrmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xrmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("X rms versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pyrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots yrmslw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    yrmslw = _extractvar('yrmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(yrmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Y rms versus time",titleb,"(number)",
                _gettitler(ilw,js))

def prrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots rrmslw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    rrmslw = _extractvar('rrmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(rrmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("R rms versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pxxpbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xxpbarlw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xxpbarlw = _extractvar('xxpbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xxpbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("X-X' bar versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pyypbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots yypbarlw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    yypbarlw = _extractvar('yypbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(yypbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Y-Y' bar versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pcurrlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots currlw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    currlw = _extractvar('currlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(currlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Current versus time",titleb,"(number)",
                _gettitler(ilw,js))

def plostparslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
                color="fg",linetype="solid",
                marks=0,marker=None,msize=1.,width=1.,lframe=0,
                titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots lostparslw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    lostparslw = _extractvar('lostparslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(lostparslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Number of lost particles versus time",titleb,"(number)",
                _gettitler(ilw,js))

def plinechglw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots linechglw as a function of time
  - ilw: lab window number to plot, must be specified
  - js=-1: species number, zero based. When -1, plots data combined from all
           species
  - toffset=0: offset added to time axis
  - tscale=1: scale of time axis - plots versus (timelw+toffset)/tscale
  - scale=1.: factor to scale data by - plots data/scale
  - color='fg': curve color
  - linetype='solid': line type
  - marks=0: turns on identifying marks on the curve
  - marker=None: marker type (see gist manual for the list)
  - msize=1: marker size
  - width=1: line width
  - lframe=0: specifies whether or not to set plot limits
  - titleb="t": bottom title
  - titles=1: specifies whether or not to plot titles
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data to plot."""
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    linechglw = _extractvar('linechglw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(linechglw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Line charge versus time",titleb,"(number)",
                _gettitler(ilw,js))

##########################################################################
def lwplotstest(ilw,**kw):
    """
Test all lwplots routines.
    """
    apply(ppnumlw,(ilw),kw);fma()
    apply(pxbarlw,(ilw),kw);fma()
    apply(pybarlw,(ilw),kw);fma()
    apply(pvzbarlw,(ilw),kw);fma()
    apply(pepsxlw,(ilw),kw);fma()
    apply(pepsylw,(ilw),kw);fma()
    apply(pepzslw,(ilw),kw);fma()
    apply(pvxrmslw,(ilw),kw);fma()
    apply(pvyrmslw,(ilw),kw);fma()
    apply(pvzrmslw,(ilw),kw);fma()
    apply(pxrmslw,(ilw),kw);fma()
    apply(pyrmslw,(ilw),kw);fma()
    apply(prrmslw,(ilw),kw);fma()
    apply(pxxpbarlw,(ilw),kw);fma()
    apply(pyypbarlw,(ilw),kw);fma()
    apply(pcurrlw,(ilw),kw);fma()
    apply(plostparslw,(ilw),kw);fma()
    apply(plinechglw,(ilw),kw);fma()

