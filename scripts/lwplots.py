"""Functions to plot lab window moments data

getlw: returns the data for a given moment

ppnumlw: Plots number of particles as a function of time
pxbarlw: Plots X bar as a function of time
pybarlw: Plots Y bar as a function of time
pzbarlw: Plots Z bar as a function of time
pxpbarlw: Plots X' bar as a function of time
pypbarlw: Plots Y' bar as a function of time
pvxbarlw: Plots Vx bar as a function of time
pvybarlw: Plots Vy bar as a function of time
pvzbarlw: Plots Vz bar as a function of time
pxybarlw: Plots XY bar as a function of time
pxypbarlw: Plots XY' bar as a function of time
pyxpbarlw: Plots YX' bar as a function of time
pxpypbarlw: Plots X'Y' bar as a function of time
pxsqbarlw: Plots XX bar as a function of time
pysqbarlw: Plots YY bar as a function of time
pzsqbarlw: Plots ZZ bar as a function of time
pxpsqbarlw: Plots X'X' bar as a function of time
pypsqbarlw: Plots Y'Y' bar as a function of time
pvxsqbarlw: Plots VxVx bar as a function of time
pvysqbarlw: Plots VyVy bar as a function of time
pvzsqbarlw: Plots Vz*Vz bar as a function of time
pxxpbarlw: Plots XX' bar as a function of time
pyypbarlw: Plots YY' bar as a function of time
pzvzbarlw: Plots ZVz bar as a function of time
pxvzbarlw: Plots XVz bar as a function of time
pyvzbarlw: Plots YVz bar as a function of time
pvxvzbarlw: Plots VxVz bar as a function of time
pvyvzbarlw: Plots VyVz bar as a function of time
pxrmslw: Plots X RMS as a function of time
pyrmslw: Plots Y RMS as a function of time
pzrmslw: Plots Z RMS as a function of time
prrmslw: Plots R RMS as a function of time
pxprmslw: Plots X' RMS as a function of time
pyprmslw: Plots Y' RMS as a function of time
pepsxlw: Plots X emittance as a function of time
pepsylw: Plots Y emittance as a function of time
pepszlw: Plots Z emittance as a function of time
pepsnxlw: Plots normalized X emittance as a function of time
pepsnylw: Plots normalized Y emittance as a function of time
pepsnzlw: Plots normalized Z emittance as a function of time
pepsrlw: Plots R emittance as a function of time
pepsglw: Plots generalized emittance as a function of time
pepshlw: Plots generalized emittance as a function of time
pepsnrlw: Plots R emittance as a function of time
pepsnglw: Plots normalized generalized emittance as a function of time
pepsnhlw: Plots normalized generalized emittance as a function of time
pvxrmslw: Plots Vx RMS as a function of time
pvyrmslw: Plots Vy RMS as a function of time
pvzrmslw: Plots Vz RMS as a function of time
pcurrlw: Plots current as a function of time
plostparslw: Plots number of lost particles as a function of time
plinechglw: Plots line-charge as a function of time

"""

from warp import *
import __main__
lwplots_version = "$Id: lwplots.py,v 1.6 2010/08/19 18:16:24 dave Exp $"

def lwplotsdoc():
    import lwplots
    print lwplots.__doc__

###########################################################################
def addlabwindow(zlw):
    """Adds a new lab window moments calculation point at the given location."""
    # --- Find first non-set value of top.zlw
    if top.nlabwn > 0:
        iz = argmax(top.zlw)
        if top.zlw[iz] == largepos:
            top.zlw[iz] = zlw
            return iz
    # --- More space is needed
    top.nlabwn += 1
    gchange('Lab_Moments')
    top.zlw[-1] = zlw
    return top.nlabwn-1

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

def getlw(name,ilw,js=-1,varsuffix=None,ff=None):
    """Returns the specified lab window data
  - name: name of the moment to return, for example 'xrms'
  - ilw: lab window number to return, must be specified
  - js=-1: species number, zero based. When -1, returns the data combined from
           all species
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables
  - ff=None: An opened file object can be specified as the place from which to
             get the data.
    """
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    if name[-2:] != 'lw': name = name + 'lw'
    return _extractvar(name,varsuffix,'top',ff)[s,ilw,js]

##########################################################################

def ppnumlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots pnumlw, number of particles, as a function of time
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
             get the data to plot.
    """
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
        ptitles("number of particles versus time",titleb,"(number)",
                _gettitler(ilw,js))

def pxbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xbarlw, X bar, as a function of time
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
             get the data to plot.
    """
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
        ptitles("X bar versus time",titleb,"(m)",
                _gettitler(ilw,js))

def pybarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots ybarlw, Y bar, as a function of time
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
             get the data to plot.
    """
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
        ptitles("Y bar versus time",titleb,"(m)",
                _gettitler(ilw,js))

def pzbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots zbarlw, Z bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    zbarlw = _extractvar('zbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(zbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Z bar versus time",titleb,"(m)",
                _gettitler(ilw,js))

def pxpbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xpbarlw, X' bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xpbarlw = _extractvar('xpbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xpbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("X' bar versus time",titleb,"(rad)",
                _gettitler(ilw,js))

def pypbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots ypbarlw, Y' bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    ypbarlw = _extractvar('ypbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(ypbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Y' bar versus time",titleb,"(rad)",
                _gettitler(ilw,js))

def pvxbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vxbarlw, Vx bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vxbarlw = _extractvar('vxbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vxbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Vx bar versus time",titleb,"(m/s)",
                _gettitler(ilw,js))

def pvybarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vybarlw, Vy bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vybarlw = _extractvar('vybarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vybarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Vy bar versus time",titleb,"(m/s)",
                _gettitler(ilw,js))

def pvzbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vzbarlw, Vz bar, as a function of time
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
             get the data to plot.
    """
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
        ptitles("Vz bar versus time",titleb,"(m/s)",
                _gettitler(ilw,js))

def pxybarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xybarlw, XY bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xybarlw = _extractvar('xybarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xybarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("XY bar versus time",titleb,"(m**2)",
                _gettitler(ilw,js))

def pxypbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xypbarlw, XY' bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xypbarlw = _extractvar('xypbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xypbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("XY' bar versus time",titleb,"(m-rad)",
                _gettitler(ilw,js))

def pyxpbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots yxpbarlw, YX' bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    yxpbarlw = _extractvar('yxpbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(yxpbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("YX' bar versus time",titleb,"(m-rad)",
                _gettitler(ilw,js))

def pxpypbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xpypbarlw, X'Y' bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xpypbarlw = _extractvar('xpypbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xpypbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("X'Y' bar versus time",titleb,"(rad-rad)",
                _gettitler(ilw,js))

def pxsqbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xsqbarlw, XX bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xsqbarlw = _extractvar('xsqbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xsqbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("XX bar versus time",titleb,"(m**2)",
                _gettitler(ilw,js))

def pysqbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots ysqbarlw, YY bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    ysqbarlw = _extractvar('ysqbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(ysqbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("YY bar versus time",titleb,"(m**2)",
                _gettitler(ilw,js))

def pzsqbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots zsqbarlw, ZZ bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    zsqbarlw = _extractvar('zsqbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(zsqbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("ZZ bar versus time",titleb,"(m**2)",
                _gettitler(ilw,js))

def pxpsqbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xpsqbarlw, X'X' bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xpsqbarlw = _extractvar('xpsqbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xpsqbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("X'X' bar versus time",titleb,"(rad-rad)",
                _gettitler(ilw,js))

def pypsqbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots ypsqbarlw, Y'Y' bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    ypsqbarlw = _extractvar('ypsqbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(ypsqbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Y'Y' bar versus time",titleb,"(rad-rad)",
                _gettitler(ilw,js))

def pvxsqbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vxsqbarlw, VxVx bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vxsqbarlw = _extractvar('vxsqbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vxsqbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("VxVx bar versus time",titleb,"(m**2/s**2)",
                _gettitler(ilw,js))

def pvysqbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vysqbarlw, VyVy bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vysqbarlw = _extractvar('vysqbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vysqbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("VyVy bar versus time",titleb,"(m**2/s**2)",
                _gettitler(ilw,js))

def pvzsqbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vzsqbarlw, VzVz bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vzsqbarlw = _extractvar('vzsqbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vzsqbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("VzVz bar versus time",titleb,"(m**2/s**2)",
                _gettitler(ilw,js))

def pxxpbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xxpbarlw, XX' bar, as a function of time
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
             get the data to plot.
    """
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
        ptitles("XX' bar versus time",titleb,"(m-rad)",
                _gettitler(ilw,js))

def pyypbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots yypbarlw, YY' bar, as a function of time
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
             get the data to plot.
    """
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
        ptitles("YY' bar versus time",titleb,"(m-rad)",
                _gettitler(ilw,js))

def pzvzbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots zvzbarlw, ZVz bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    zvzbarlw = _extractvar('zvzbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(zvzbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("ZVz bar versus time",titleb,"(m-m/s)",
                _gettitler(ilw,js))

def pxvzbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xvzbarlw, XVz bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xvzbarlw = _extractvar('xvzbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xvzbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("XVz bar versus time",titleb,"(m-m/s)",
                _gettitler(ilw,js))

def pyvzbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
              color="fg",linetype="solid",
              marks=0,marker=None,msize=1.,width=1.,lframe=0,
              titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots yvzbarlw, YVz bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    yvzbarlw = _extractvar('yvzbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(yvzbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("YVz bar versus time",titleb,"(m-m/s)",
                _gettitler(ilw,js))

def pvxvzbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vxvzbarlw, VxVz bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vxvzbarlw = _extractvar('vxvzbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vxvzbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("VxVz bar versus time",titleb,"(m**2/s**2)",
                _gettitler(ilw,js))

def pvyvzbarlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vyvzbarlw, VyVz bar, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    vyvzbarlw = _extractvar('vyvzbarlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(vyvzbarlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("VyVz bar versus time",titleb,"(m**2/s**2)",
                _gettitler(ilw,js))

def pxrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xrmslw, X RMS, as a function of time
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
             get the data to plot.
    """
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
        ptitles("X RMS versus time",titleb,"(m)",
                _gettitler(ilw,js))

def pyrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots yrmslw, Y RMS, as a function of time
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
             get the data to plot.
    """
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
        ptitles("Y RMS versus time",titleb,"(m)",
                _gettitler(ilw,js))

def pzrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots zrmslw, Z RMS, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    zrmslw = _extractvar('zrmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(zrmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Z RMS versus time",titleb,"(m)",
                _gettitler(ilw,js))

def prrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots rrmslw, R RMS, as a function of time
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
             get the data to plot.
    """
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
        ptitles("R RMS versus time",titleb,"(m)",
                _gettitler(ilw,js))

def pxprmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots xprmslw, X' RMS, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    xprmslw = _extractvar('xprmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(xprmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("X' RMS versus time",titleb,"(rad)",
                _gettitler(ilw,js))

def pyprmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots yprmslw, Y' RMS, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    yprmslw = _extractvar('yprmslw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(yprmslw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Y' RMS versus time",titleb,"(rad)",
                _gettitler(ilw,js))

def pepsxlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsxlw, X emittance, as a function of time
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
             get the data to plot.
    """
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
        ptitles("X emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsylw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsylw, Y emittance, as a function of time
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
             get the data to plot.
    """
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
        ptitles("Y emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepszlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epszlw, Z emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epszlw = _extractvar('epszlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epszlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("Z emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsnxlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsnxlw, normalized X emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsnxlw = _extractvar('epsnxlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsnxlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("normalized X emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsnylw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsnylw, normalized Y emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsnylw = _extractvar('epsnylw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsnylw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("normalized Y emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsnzlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsnzlw, normalized Z emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsnzlw = _extractvar('epsnzlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsnzlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("normalized Z emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsrlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsrlw, R emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsrlw = _extractvar('epsrlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsrlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("R emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsglw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsglw, generalized emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsglw = _extractvar('epsglw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsglw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("generalized emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepshlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epshlw, generalized emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epshlw = _extractvar('epshlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epshlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("generalized emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsnrlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsnrlw, R emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsnrlw = _extractvar('epsnrlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsnrlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("R emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsnglw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsnglw, normalized generalized emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsnglw = _extractvar('epsnglw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsnglw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("normalized generalized emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pepsnhlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots epsnhlw, normalized generalized emittance, as a function of time
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
             get the data to plot.
    """
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    epsnhlw = _extractvar('epsnhlw',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(epsnhlw/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("normalized generalized emittance versus time",titleb,"(pi-m-rad)",
                _gettitler(ilw,js))

def pvxrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vxrmslw, Vx RMS, as a function of time
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
             get the data to plot.
    """
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
        ptitles("Vx RMS versus time",titleb,"(m/s)",
                _gettitler(ilw,js))

def pvyrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vyrmslw, Vy RMS, as a function of time
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
             get the data to plot.
    """
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
        ptitles("Vy RMS versus time",titleb,"(m/s)",
                _gettitler(ilw,js))

def pvzrmslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
             color="fg",linetype="solid",
             marks=0,marker=None,msize=1.,width=1.,lframe=0,
             titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots vzrmslw, Vz RMS, as a function of time
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
             get the data to plot.
    """
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
        ptitles("Vz RMS versus time",titleb,"(m/s)",
                _gettitler(ilw,js))

def pcurrlw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
            color="fg",linetype="solid",
            marks=0,marker=None,msize=1.,width=1.,lframe=0,
            titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots currlw, current, as a function of time
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
             get the data to plot.
    """
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
        ptitles("current versus time",titleb,"(A)",
                _gettitler(ilw,js))

def plostparslw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
                color="fg",linetype="solid",
                marks=0,marker=None,msize=1.,width=1.,lframe=0,
                titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots lostparslw, number of lost particles, as a function of time
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
             get the data to plot.
    """
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
        ptitles("number of lost particles versus time",titleb,"(number)",
                _gettitler(ilw,js))

def plinechglw(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
               color="fg",linetype="solid",
               marks=0,marker=None,msize=1.,width=1.,lframe=0,
               titleb=None,titles=1,varsuffix=None,ff=None):
    """Plots linechglw, line-charge, as a function of time
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
             get the data to plot.
    """
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
        ptitles("line-charge versus time",titleb,"(C)",
                _gettitler(ilw,js))

##########################################################################
def lwplotstest(ilw,**kw):
    """
Test all lwplots routines.
    """
    apply(ppnumlw,(ilw),kw);fma()
    apply(pxbarlw,(ilw),kw);fma()
    apply(pybarlw,(ilw),kw);fma()
    apply(pzbarlw,(ilw),kw);fma()
    apply(pxpbarlw,(ilw),kw);fma()
    apply(pypbarlw,(ilw),kw);fma()
    apply(pvxbarlw,(ilw),kw);fma()
    apply(pvybarlw,(ilw),kw);fma()
    apply(pvzbarlw,(ilw),kw);fma()
    apply(pxybarlw,(ilw),kw);fma()
    apply(pxypbarlw,(ilw),kw);fma()
    apply(pyxpbarlw,(ilw),kw);fma()
    apply(pxpypbarlw,(ilw),kw);fma()
    apply(pxsqbarlw,(ilw),kw);fma()
    apply(pysqbarlw,(ilw),kw);fma()
    apply(pzsqbarlw,(ilw),kw);fma()
    apply(pxpsqbarlw,(ilw),kw);fma()
    apply(pypsqbarlw,(ilw),kw);fma()
    apply(pvxsqbarlw,(ilw),kw);fma()
    apply(pvysqbarlw,(ilw),kw);fma()
    apply(pvzsqbarlw,(ilw),kw);fma()
    apply(pxxpbarlw,(ilw),kw);fma()
    apply(pyypbarlw,(ilw),kw);fma()
    apply(pzvzbarlw,(ilw),kw);fma()
    apply(pxvzbarlw,(ilw),kw);fma()
    apply(pyvzbarlw,(ilw),kw);fma()
    apply(pvxvzbarlw,(ilw),kw);fma()
    apply(pvyvzbarlw,(ilw),kw);fma()
    apply(pxrmslw,(ilw),kw);fma()
    apply(pyrmslw,(ilw),kw);fma()
    apply(pzrmslw,(ilw),kw);fma()
    apply(prrmslw,(ilw),kw);fma()
    apply(pxprmslw,(ilw),kw);fma()
    apply(pyprmslw,(ilw),kw);fma()
    apply(pepsxlw,(ilw),kw);fma()
    apply(pepsylw,(ilw),kw);fma()
    apply(pepszlw,(ilw),kw);fma()
    apply(pepsnxlw,(ilw),kw);fma()
    apply(pepsnylw,(ilw),kw);fma()
    apply(pepsnzlw,(ilw),kw);fma()
    apply(pepsrlw,(ilw),kw);fma()
    apply(pepsglw,(ilw),kw);fma()
    apply(pepshlw,(ilw),kw);fma()
    apply(pepsnrlw,(ilw),kw);fma()
    apply(pepsnglw,(ilw),kw);fma()
    apply(pepsnhlw,(ilw),kw);fma()
    apply(pvxrmslw,(ilw),kw);fma()
    apply(pvyrmslw,(ilw),kw);fma()
    apply(pvzrmslw,(ilw),kw);fma()
    apply(pcurrlw,(ilw),kw);fma()
    apply(plostparslw,(ilw),kw);fma()
    apply(plinechglw,(ilw),kw);fma()




###########################################################################
# --- Below is the code that was used to generate the functions above. This
# --- was saved in case the function needs to be changed. This can be used
# --- again instead of modifying all of the functions above.
func =\
"""def p%(lwm)s(ilw,js=-1,toffset=0.,tscale=1.,scale=1.,
%(space)scolor="fg",linetype="solid",
%(space)smarks=0,marker=None,msize=1.,width=1.,lframe=0,
%(space)stitleb=None,titles=1,varsuffix=None,ff=None):
    \"\"\"Plots %(lwm)s, %(comment)s, as a function of time
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
             get the data to plot.
    \"\"\"
    if tscale == 0.: raise "tscale must be nonzero"
    if titleb is None:
        if tscale == 1.: titleb = "t (s)"
        else: titleb = "t"
    ilabwn = _extractvar('ilabwn',varsuffix,'top',ff)[ilw,js]
    s = s_[:ilabwn]
    %(lwm)s = _extractvar('%(lwm)s',varsuffix,'top',ff)[s,ilw,js]
    timelw = _extractvar('timelw',varsuffix,'top',ff)[s,ilw,js]
    plg(%(lwm)s/scale,(toffset+timelw)/tscale,color=color,linetype=linetype,
        marks=marks,marker=marker,msize=msize,width=width)
    if titles:
        ptitles("%(comment)s versus time",titleb,"(%(units)s)",
                _gettitler(ilw,js))
"""

lwlist = [
["pnumlw","number of particles","number"],
["xbarlw","X bar","m"],
["ybarlw","Y bar","m"],
["zbarlw","Z bar","m"],
["xpbarlw","X' bar","rad"],
["ypbarlw","Y' bar","rad"],
["vxbarlw","Vx bar","m/s"],
["vybarlw","Vy bar","m/s"],
["vzbarlw","Vz bar","m/s"],
["xybarlw","XY bar","m**2"],
["xypbarlw","XY' bar","m-rad"],
["yxpbarlw","YX' bar","m-rad"],
["xpypbarlw","X'Y' bar","rad-rad"],
["xsqbarlw","XX bar","m**2"],
["ysqbarlw","YY bar","m**2"],
["zsqbarlw","ZZ bar","m**2"],
["xpsqbarlw","X'X' bar","rad-rad"],
["ypsqbarlw","Y'Y' bar","rad-rad"],
["vxsqbarlw","VxVx bar","m**2/s**2"],
["vysqbarlw","VyVy bar","m**2/s**2"],
["vzsqbarlw","VzVz bar","m**2/s**2"],
["xxpbarlw","XX' bar","m-rad"],
["yypbarlw","YY' bar","m-rad"],
["zvzbarlw","ZVz bar","m-m/s"],
["xvzbarlw","XVz bar","m-m/s"],
["yvzbarlw","YVz bar","m-m/s"],
["vxvzbarlw","VxVz bar","m**2/s**2"],
["vyvzbarlw","VyVz bar","m**2/s**2"],
["xrmslw","X RMS","m"],
["yrmslw","Y RMS","m"],
["zrmslw","Z RMS","m"],
["rrmslw","R RMS","m"],
["xprmslw","X' RMS","rad"],
["yprmslw","Y' RMS","rad"],
["epsxlw","X emittance","pi-m-rad"],
["epsylw","Y emittance","pi-m-rad"],
["epszlw","Z emittance","pi-m-rad"],
["epsnxlw","normalized X emittance","pi-m-rad"],
["epsnylw","normalized Y emittance","pi-m-rad"],
["epsnzlw","normalized Z emittance","pi-m-rad"],
["epsrlw","R emittance","pi-m-rad"],
["epsglw","generalized emittance","pi-m-rad"],
["epshlw","generalized emittance","pi-m-rad"],
["epsnrlw","R emittance","pi-m-rad"],
["epsnglw","normalized generalized emittance","pi-m-rad"],
["epsnhlw","normalized generalized emittance","pi-m-rad"],
["vxrmslw","Vx RMS","m/s"],
["vyrmslw","Vy RMS","m/s"],
["vzrmslw","Vz RMS","m/s"],
["currlw","current","A"],
["lostparslw","number of lost particles","number"],
["linechglw","line-charge","C"]]

"""
for lwm,comment,units in lwlist:
  space = (6+len(lwm))*" "
  print func%locals()
"""
