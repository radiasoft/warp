from warp import *
from mplot import *
import __main__
histplots_version = "$Id: histplots.py,v 1.17 2002/10/09 18:19:40 dave Exp $"

hpbasictext = """
  - absc: Data for the abscissa. Defaults to either thist or hzbeam
  - xmin, xmax, ymin, ymax: Plot limits, defaults to extrema of data
  - titlet, titleb, titlel, titler: Plot titles, defaults to blank except
    titleb  which is set to name of bottom axis when absc not input
  - xscale=1.0: Scale factor for abscissa
  - xoffset=0.0: Offset for abscissa
  - yscale=1.0: Scale factor for oordinate
  - yoffset=0.0: Offset for oordinate
  - lnormalized=0: When true, scales data by initial value
  - istart=0: Time index at which to start the plots
  - iend=jhist: Time index at which to end the plots
  - istep=1: Step size of time data plotted
  - lhzbeam=0: When true, plots data versus hzbeam instead of thist
  - lvsz=0: Old name for lhzbeam
  - lzshift=0: specifies whether the z-axis is shifted by the window
              location
  - logplot=0: When true, make a log plot
  - color='fg': Color of the curve plotted
  - marks=0: Marks to place on curve
  - marker=None: Marker to place on curve
  - msize=1.0: Marker size
  - width=1.0: Line width
  - linetype='solid': Line type
  - titles=1: When false, no titles are printed
  - plsysval=1: Plot system to make plot in (quadrant plots for example)
    see work.gs for the plot system values.
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables"""
hpbasicwintext = (
"""  - iw=0: Window to chose from""" + hpbasictext +
"""  - lzshift=0: specifies whether the z-axis is shifted by the window
              location""")
hpbasicconttext = (
"""
  - absc: Data for the abscissa. Defaults to either thist or hzbeam
  - xmin, xmax, ymin, ymax: Plot limits, defaults to extrema of data
  - titlet, titleb, titlel, titler: Plot titles, defaults to blank except
    titleb  which is set to name of bottom axis when absc not input
  - xscale=1.0: Scale factor for abscissa
  - xoffset=0.0: Offset for abscissa
  - yscale=1.0: Scale factor for oordinate
  - yoffset=0.0: Offset for oordinate
  - istart=0: Time index at which to start the plots
  - iend=jhist: Time index at which to end the plots
  - istep=max(iend/20,1): Step size of time data plotted
  - jstart=0: Z index at which to start the plots
  - jend=nzzarr: Z index at which to end the plots
  - jstep=max(jend/32,1): Z step size of time data plotted
  - lhzbeam=0: When true, plots data versus hzbeam instead of thist
  - logplot=0: When true, make a log plot
  - color='fg': Color of the curve plotted
  - marks=0: Marks to place on curve
  - marker=None: Marker to place on curve
  - msize=1.0: Marker size
  - titles=1: Line width
  - width=1.0: Line width
  - linetype='solid': Line type
  - levs=None: number of contours levels
  - titles=true: When false, no titles are printed
  - filled=0: When true, plot filled contours
  - varsuffix=None: When specified, variables with that suffix are used
                    instead of the fortran variables""")
hpzarraytext = (
  """
  - contour=0 when 1 plots contours
  - overlay=0 when 1 overlays values from different times
  - iz when specified, plots history at that location""" + hpbasicconttext +
  """ - For mountain range plot options, run doc(mountainplot1)""")

###########################################################################
def hpdoc():
  print """
What follows is a list of all possible arguments to any of the history
plotting routines, along with their default values.
  """ + hpbasictext

###########################################################################
def _extractvar(name,varsuffix=None,pkg='top'):
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
      try:
        result = eval(vname,__main__.__dict__)
      except:
        result = None
      if result is not None: return result
    return eval(pkg+'.'+name,globals())
  else:
    return name

def _extractvarkw(name,kw,pkg='top'):
  return _extractvar(name,kw.get('varsuffix',None),pkg=pkg)

###########################################################################
def hpbasic(oord,kwdict={},**kw):
  """
This is the basic history plot with all of its possible arguments. The
parsing of the arguments is done in one place here and each of the history
plot functions just passes most of the arguments into this function. The
only required argument of course is the data to be plotted.
  """
  kwdefaults = {'absc':None,'xmin':'e','xmax':'e','ymin':'e','ymax':'e',
                'titlet':'','titleb':'','titlel':'','titler':'',
                'xscale':1.0,'xoffset':0.0,'yscale':array([1.0]),'yoffset':0.0,
                'lnormalized':0,'istart':0,'iend':'jhist','istep':1,
                'lhzbeam':0,'lvsz':0,'logplot':0,
                'color':'fg','marks':0,'marker':None,'msize':1.0,'titles':1,
                'plsysval':1,'width':1.,'linetype':'solid',
                'varsuffix':None}
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  badargs = checkarguments(kwvalues,kwdefaults)
  if badargs: raise "bad argument ",string.join(badargs.keys())
  for arg in kwvalues.keys(): exec(arg+" = kwvalues['"+arg+"']")

  iend = _extractvar(iend,varsuffix)
  oord = _extractvar(oord,varsuffix)

  # --- Now complete the setup
  if lnormalized:
    yscale = 1./oord[...,istart]
    if not titlel: titlel = titlet + " over initial value"
  if type(yscale) != type(array([])): yscale = array([yscale])
  if not absc:
    if lhzbeam or lvsz:
      hzbeam = _extractvar("hzbeam",varsuffix)
      absc = hzbeam[istart:iend+1:istep]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "Z (m)"
        else:
          titleb = "Z"
    else:
      thist = _extractvar("thist",varsuffix)
      absc = thist[istart:iend+1:istep]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "time (s)"
        else:
          titleb = "time"
  if logplot:
    logxy(1,0)
    oord = log10(maximum(10e-12,
                 oord[...,istart:iend+1:istep]*yscale[...,NewAxis]+yoffset))
    if not titler:
      titler = "logarithmic scale"
  else:
    oord = oord[...,istart:iend+1:istep]*yscale[...,NewAxis]+yoffset

  # --- Now actually make the plot after all of that ado.   
  pla(transpose(oord),absc,color=color,msize=msize,marks=marks,marker=marker,
      width=width,linetype=linetype,decomposed=0)
  if titles: ptitles(titlet,titleb,titlel,titler,plsysval)
  limits(xmin,xmax,ymin,ymax)
  if logplot: logxy(0,0)


#############################################################################
# Basic history plot for data in windows
def hpbasicwin(oord,iw=0,kwdict={},**kw):
  kw.update(kwdict)
  oord = _extractvarkw(oord,kw)
  zwindows = _extractvarkw('zwindows',kw)
  if 'lzshift' in kw.keys():
    lzshift = kw['lzshift']
    del kw['lzshift']
  else:
    lzshift = 0
  if lzshift: kw['xoffset'] = 0.5*(zwindows[0,iw]+zwindows[1,iw])
  if iw == 0:
    kw['titlet'] = kw['titlet'] + " for whole beam"
    hpbasic(oord[iw,:],kw)
  elif iw > 0:
    kw['titler']="z%1d = %6.3f"%(iw,.5*(zwindows[0,iw]+zwindows[1,iw]))
    hpbasic(oord[iw,:],kw)
  else:
    nzwind = _extractvarkw('nzwind',kw)
    for i in range(1,nzwind):
      plsys(3+(i-1)%4)
      kw['titler'] = "z%1d = %6.3f"%(i,.5*(zwindows[0,i]+zwindows[1,i]))
      kw['plsysval'] = 3+(i-1)%4
      hpbasic(oord[i,:],kw)
      if (i-1)%4 == 3: fma()
    fma()

###########################################################################
# This is the basic history contour plot with all of its possible arguments.
# The parsing of the arguments is done in one place here and each of the
# history plot functions (declared further down) just pass most of the
# arguments into this function. The only required argument of course of the
# data to be plotted.
def hpbasiccont(oord,oordmesh,kwdict={},**kw):
  kwdefaults = {'absc':None,'xmin':'e','xmax':'e','ymin':'e','ymax':'e',
                'titlet':'','titleb':'','titlel':'','titler':'',
                'xscale':1.0,'xoffset':0.0,'yscale':1.0,'yoffset':0.0,
                'istart':0,'iend':'jhist','istep':None,'jstart':0,
                'jend':'nzzarr','jstep':None,'lhzbeam':0,'logplot':0,
                'color':'fg','marks':0,'marker':None,'msize':1.0,
                'titles':1,'levs':10,'filled':0,'width':1.,'linetype':'solid',
                'varsuffix':None}
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")
  badargs = checkarguments(kwvalues,kwdefaults)
  if badargs: raise "bad argument ",string.join(badargs.keys())

  iend = _extractvar(iend,varsuffix)
  jend = _extractvar(jend,varsuffix)
  oord = _extractvar(oord,varsuffix)
  oordmesh = _extractvar(oordmesh,varsuffix)

  # --- Some special arguments
  if istep is None: istep = max(iend/20,1)
  if jstep is None: jstep = max(jend/32,1)

  # --- Check through the possible keyword arguments
  for arg in kw.keys():
    if arg in keywords:
      exec(arg+" = kw['"+arg+"']")
      del kw[arg]
  # --- Raise an error if there are keywords left over
  if len(kw) > 0:
    badkwlist = 'unexpected keyword argument:'
    if len(kw) > 1: badkwlist = badkwlist[:-1] + 's:'
    for arg in kw.keys():
      badkwlist = badkwlist + ' ' + arg
    raise TypeError, badkwlist


  # --- Now complete the setup
  if not absc:
    if lhzbeam:
      hzbeam = _extractvar('hzbeam',varsuffix)
      absc = hzbeam[istart:iend+1:istep]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "Z (m)"
        else:
          titleb = "Z"
    else:
      thist = _extractvar('thist',varsuffix)
      absc = thist[istart:iend+1:istep]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "time (s)"
        else:
          titleb = "time"
  if logplot:
    oord=log10(maximum(10e-12,
               transpose(oord[jstart:jend+1:jstep,istart:iend+1:istep])*
               yscale+yoffset))
    if not titler:
      titler = "logarithmic scale"
  else:
    oord=transpose(oord[jstart:jend+1:jstep,istart:iend+1:istep]*yscale+yoffset)

  # --- Now actually make the plot after all of that ado.   
  if filled:
    plotfc(oord,absc,oordmesh[jstart:jend+1:jstep],contours=levs)
  else:
    plotc(oord,absc,oordmesh[jstart:jend+1:jstep],color=color,levs=levs,
          width=width,linetype=linetype)
  if titles: ptitles(titlet,titleb,titlel,titler)
  if xmin == 'e': xmin = min(oordmesh[jstart:jend+1:jstep])
  if xmax == 'e': xmax = max(oordmesh[jstart:jend+1:jstep])
  if ymin == 'e': ymin = min(absc)
  if ymax == 'e': ymax = max(absc)
  limits(xmin,xmax,ymin,ymax)

###########################################################################
# --- plot hzarray in various ways
def hpzarray(hzarray,contour=0,overlay=0,iz=None,kwdict={},**kw):
  """
Plots data in various ways. By default, makes a mountain range plot.
  - contour=0 when 1 plots contours
  - overlay=0 when 1 overlays values from different times
  - iz when specified, plots history at that location
  - For mountain range plot options, run doc(mountainplot1)
  """
  kw.update(kwdict)
  #kw.setdefault('titleb','time')
  hzarray = _extractvarkw(hzarray,kw)
  zmntmesh = _extractvarkw('zmntmesh',kw)
  jhist = _extractvarkw('jhist',kw)
  if iz:
    nz = _extractvarkw('nz',kw,'w3d')
    if iz == nz/2:
      kw['titlet']=kw['titlet']+" at Center"
    else:
      kw['titlet']=kw['titlet']+" at iz = %d"%iz
    kw['titler'] = 'z = %6.3f'%zmntmesh[iz]
    hpbasic(hzarray[iz,:],kw)
    return
  elif contour:
    hpbasiccont(hzarray,zmntmesh,kw)
  elif overlay:
    zmesh = _extractvarkw('zmesh',kw,'w3d')
    kw.setdefault('istep',max(jhist/10,1))
    kw.setdefault('ord',zmesh)
    kw.setdefault('iend',jhist)
    mountainplot1(kw['titlet'],hzarray,kw)
  else:
    kw.setdefault('iend',jhist)
    mountainplot1(kw['titlet'],hzarray,kw)

###########################################################################
###########################################################################
###########################################################################

def hpzbeam(kwdict={},**kw):
  "Beam frame location."
  kw.update(kwdict)
  kw['titlet']="Beam frame location"
  kw['titlel']="(m)"
  hpbasic('hzbeam',kw)
if sys.version[:5] != "1.5.1":
  hpzbeam.__doc__ = hpzbeam.__doc__ + hpbasictext


def hpvbeam(kwdict={},**kw):
  "Beam frame velocity."
  kw.update(kwdict)
  kw['titlet']="Beam frame velocity"
  kw['titlel']="(m/s)"
  hpbasic('hvbeam',kw)
if sys.version[:5] != "1.5.1":
  hpvbeam.__doc__ = hpvbeam.__doc__ + hpbasictext


def hpbmlen(kwdict={},**kw):
  "RMS beam length."
  kw.update(kwdict)
  kw['titlet']="RMS beam length"
  kw['titlel']="(m)"
  hpbasic('hbmlen',kw)
if sys.version[:5] != "1.5.1":
  hpbmlen.__doc__ = hpbmlen.__doc__ + hpbasictext


def hpefld(kwdict={},**kw):
  "Field energy."
  kw.update(kwdict)
  kw['titlet']="Field energy"
  kw['titlel']="(J)"
  hpbasic('hefld',kw)
if sys.version[:5] != "1.5.1":
  hpefld.__doc__ = hpefld.__doc__ + hpbasictext


def hpekzmbe(kwdict={},**kw):
  "Total Z Kinetic energy minus beam energy."
  kw.update(kwdict)
  kw['titlet']="Total Z Kinetic energy minus beam energy"
  kw['titlel']="(J)"
  hpbasic('hekzmbe',kw)
if sys.version[:5] != "1.5.1":
  hpekzmbe.__doc__ = hpekzmbe.__doc__ + hpbasictext


def hpekzbeam(kwdict={},**kw):
  "Z Kinetic energy in the beam frame."
  kw.update(kwdict)
  kw['titlet']="Z Kinetic energy in the beam frame"
  kw['titlel']="(J)"
  hpbasic('hekzbeam',kw)
if sys.version[:5] != "1.5.1":
  hpekzbeam.__doc__=hpekzbeam.__doc__+hpbasictext


def hpekperp(kwdict={},**kw):
  "Perp Kinetic energy."
  kw.update(kwdict)
  kw['titlet']="Perp Kinetic energy"
  kw['titlel']="(J)"
  hpbasic('hekperp',kw)
if sys.version[:5] != "1.5.1":
  hpekperp.__doc__ = hpekperp.__doc__ + hpbasictext


def hpekinz(iw=0,kwdict={},**kw):
  "Z Kinetic energy (units of MV)"
  kw.update(kwdict)
  kw['titlet']="Z Kinetic energy"
  kw['titlel']="(MV)"
  lrelativ = _extractvarkw('lrelativ',kw)
  hvzbar = _extractvarkw('hvzbar',kw)
  aion = _extractvarkw('aion',kw)
  if lrelativ:
    u = (hvzbar/clight)**2
    e = (aion*amu*clight**2)*(u/(sqrt(1.-u)+1.-u))/jperev              
  else:
    e = 0.5*aion*amu*hvzbar**2/jperev
  hpbasicwin(e*1.e-6,iw,kw)
if sys.version[:5] != "1.5.1":
  hpekinz.__doc__ = hpekinz.__doc__ + hpbasicwintext


def hpekin(iw=0,kwdict={},**kw):
  "Total Kinetic energy (units of MV)"
  kw.update(kwdict)
  kw['titlet']="Total Kinetic energy"
  kw['titlel']="(MV)"
  lrelativ = _extractvarkw('lrelativ',kw)
  hvxbar = _extractvarkw('hvxbar',kw)
  hvybar = _extractvarkw('hvybar',kw)
  hvzbar = _extractvarkw('hvzbar',kw)
  aion = _extractvarkw('aion',kw)
  if lrelativ:
    u = (hvxbar**2 + hvybar**2 + hvzbar**2)/clight**2
    e = (aion*amu*clight**2)*(u/(sqrt(1.-u)+1.-u))/jperev              
  else:
    u = hvxbar**2 + hvybar**2 + hvzbar**2
    e = 0.5*aion*amu*u/jperev
  hpbasicwin(e*1.e-6,iw,kw)
if sys.version[:5] != "1.5.1":
  hpekin.__doc__ = hpekin.__doc__ + hpbasicwintext


def hpepsx(iw=0,kwdict={},**kw):
  "X emittance."
  kw.update(kwdict)
  kw['titlet']="X emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin('hepsx',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsx.__doc__ = hpepsx.__doc__ + hpbasicwintext


def hpepsy(iw=0,kwdict={},**kw):
  "Y emittance."
  kw.update(kwdict)
  kw['titlet']="Y emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin('hepsy',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsy.__doc__ = hpepsy.__doc__ + hpbasicwintext


def hpepsz(iw=0,kwdict={},**kw):
  "Z emittance."
  kw.update(kwdict)
  kw['titlet']="Z emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin('hepsz',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsz.__doc__ = hpepsz.__doc__ + hpbasicwintext


def hpepsnx(iw=0,kwdict={},**kw):
  "X normalized emittance."
  kw.update(kwdict)
  kw['titlet']="X normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin('hepsnx',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsnx.__doc__ = hpepsnx.__doc__ + hpbasicwintext


def hpepsny(iw=0,kwdict={},**kw):
  "Y normalized emittance."
  kw.update(kwdict)
  kw['titlet']="Y normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin('hepsny',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsny.__doc__ = hpepsny.__doc__ + hpbasicwintext


def hpepsnz(iw=0,kwdict={},**kw):
  "Z normalized emittance."
  kw.update(kwdict)
  kw['titlet']="Z normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin('hepsnz',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsnz.__doc__ =  hpepsnz.__doc__ + hpbasicwintext


def hpepsg(iw=0,kwdict={},**kw):
  "Generalized emittance."
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin('hepsg',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsg.__doc__ = hpepsg.__doc__ + hpbasicwintext


def hpepsh(iw=0,kwdict={},**kw):
  "Generalized emittance."
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin('hepsh',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsh.__doc__ = hpepsh.__doc__ + hpbasicwintext


def hpepsng(iw=0,kwdict={},**kw):
  "Generalized normalized emittance."
  kw.update(kwdict)
  kw['titlet']="Generalized normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin('hepsng',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsng.__doc__ = hpepsng.__doc__ + hpbasicwintext


def hpepsnh(iw=0,kwdict={},**kw):
  "Generalized normalized emittance."
  kw.update(kwdict)
  kw['titlet']="Generalized normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin('hepsnh',iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsnh.__doc__ = hpepsnh.__doc__ + hpbasicwintext


def hppnum(iw=0,kwdict={},**kw):
  "Number of particles."
  kw.update(kwdict)
  kw['titlet']="Number of particles"
  kw['titlel']=""
  hpbasicwin('hpnum',iw,kw)
if sys.version[:5] != "1.5.1":
  hppnum.__doc__ = hppnum.__doc__ + hpbasicwintext


def hprhomid(iw=0,kwdict={},**kw):
  "Charge density on axis."
  kw.update(kwdict)
  kw['titlet']="Charge density on axis"
  kw['titlel']="(C/m**3)"
  hpbasicwin('hrhomid',iw,kw)
if sys.version[:5] != "1.5.1":
  hprhomid.__doc__ = hprhomid.__doc__ + hpbasicwintext


def hprhomax(iw=0,kwdict={},**kw):
  "Charge density max."
  kw.update(kwdict)
  kw['titlet']="Charge density max"
  kw['titlel']="(C/m**3)"
  hpbasicwin('hrhomax',iw,kw)
if sys.version[:5] != "1.5.1":
  hprhomax.__doc__ = hprhomax.__doc__ + hpbasicwintext


def hpxbar(iw=0,kwdict={},**kw):
  "True mean x."
  kw.update(kwdict)
  kw['titlet']="True mean x"
  kw['titlel']="(m)"
  hpbasicwin('hxbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxbar.__doc__ = hpxbar.__doc__ + hpbasicwintext


def hpybar(iw=0,kwdict={},**kw):
  "True mean y."
  kw.update(kwdict)
  kw['titlet']="True mean y"
  kw['titlel']="(m)"
  hpbasicwin('hybar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpybar.__doc__ = hpybar.__doc__ + hpbasicwintext


def hpxybar(iw=0,kwdict={},**kw):
  "True mean xy."
  kw.update(kwdict)
  kw['titlet']="True mean xy"
  kw['titlel']="(m**2)"
  hpbasicwin('hxybar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxybar.__doc__ = hpxybar.__doc__ + hpbasicwintext


def hpxrms(iw=0,kwdict={},**kw):
  "True RMS x."
  kw.update(kwdict)
  kw['titlet']="True RMS x"
  kw['titlel']="(m)"
  hpbasicwin('hxrms',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxrms.__doc__ = hpxrms.__doc__ + hpbasicwintext


def hpyrms(iw=0,kwdict={},**kw):
  "True RMS y."
  kw.update(kwdict)
  kw['titlet']="True RMS y"
  kw['titlel']="(m)"
  hpbasicwin('hyrms',iw,kw)
if sys.version[:5] != "1.5.1":
  hpyrms.__doc__ = hpyrms.__doc__ + hpbasicwintext


def hpxprms(iw=0,kwdict={},**kw):
  "True RMS x'."
  kw.update(kwdict)
  kw['titlet']="True RMS x'"
  kw['titlel']="(rad)"
  hpbasicwin('hxprms',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxprms.__doc__ = hpxprms.__doc__ + hpbasicwintext


def hpyprms(iw=0,kwdict={},**kw):
  "True RMS y'."
  kw.update(kwdict)
  kw['titlet']="True RMS y'"
  kw['titlel']="(rad)"
  hpbasicwin('hyprms',iw,kw)
if sys.version[:5] != "1.5.1":
  hpyprms.__doc__ = hpyprms.__doc__ + hpbasicwintext


def hpxsqbar(iw=0,kwdict={},**kw):
  "Mean x squared."
  kw.update(kwdict)
  kw['titlet']="Mean x squared"
  kw['titlel']="(m**2)"
  hpbasicwin('hxsqbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxsqbar.__doc__ = hpxsqbar.__doc__ + hpbasicwintext


def hpysqbar(iw=0,kwdict={},**kw):
  "Mean y squared."
  kw.update(kwdict)
  kw['titlet']="Mean y squared"
  kw['titlel']="(m**2)"
  hpbasicwin('hysqbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpysqbar.__doc__ = hpysqbar.__doc__ + hpbasicwintext


def hpvxbar(iw=0,kwdict={},**kw):
  "Mean vx."
  kw.update(kwdict)
  kw['titlet']="Mean vx"
  kw['titlel']="(m/s)"
  hpbasicwin('hvxbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpvxbar.__doc__ = hpvxbar.__doc__ + hpbasicwintext


def hpvybar(iw=0,kwdict={},**kw):
  "Mean vy."
  kw.update(kwdict)
  kw['titlet']="Mean vy"
  kw['titlel']="(m/s)"
  hpbasicwin('hvybar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpvybar.__doc__ = hpvybar.__doc__ + hpbasicwintext


# --- plot mean z velocity in beam or lab frame vs time
def hpvzbar(iw=0,beamframe=1,kwdict={},**kw):
  """Mean Z Velocity (beam frame or lab frame)
  - beamframe=1: when true, plot Vz relative to beam frame (vzbar - vbeam)"""
  kw.update(kwdict)
  kw['titlel']="(m/s)"
  vbeam = _extractvarkw('vbeam',kw)
  hvzbar = _extractvarkw('hvzbar',kw)
  if beamframe:
    kw['titlet']="Mean Z Velocity (beam frame)"
    kw['titler']="vbeam = %6.3e"%vbeam
    hpbasicwin(hvzbar-vbeam,iw,kw)
  else:
    kw['titlet']="Mean Z Velocity"
    hpbasicwin('hvzbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpvzbar.__doc__ = hpvzbar.__doc__ + hpbasicwintext


def hpxpbar(iw=0,kwdict={},**kw):
  "Mean x'."
  kw.update(kwdict)
  kw['titlet']="Mean x'"
  kw['titlel']="(rad)"
  hpbasicwin('hxpbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxpbar.__doc__ = hpxpbar.__doc__ + hpbasicwintext


def hpypbar(iw=0,kwdict={},**kw):
  "Mean y'."
  kw.update(kwdict)
  kw['titlet']="Mean y'"
  kw['titlel']="(rad)"
  hpbasicwin('hypbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpypbar.__doc__ = hpypbar.__doc__ + hpbasicwintext


def hpvxrms(iw=0,kwdict={},**kw):
  "True RMS vx."
  kw.update(kwdict)
  kw['titlet']="True RMS vx"
  kw['titlel']="(m/s)"
  hpbasicwin('hvxrms',iw,kw)
if sys.version[:5] != "1.5.1":
  hpvxrms.__doc__ = hpvxrms.__doc__ + hpbasicwintext


def hpvyrms(iw=0,kwdict={},**kw):
  "True RMS vy."
  kw.update(kwdict)
  kw['titlet']="True RMS vy"
  kw['titlel']="(m/s)"
  hpbasicwin('hvyrms',iw,kw)
if sys.version[:5] != "1.5.1":
  hpvyrms.__doc__ = hpvyrms.__doc__ + hpbasicwintext


def hpvzrms(iw=0,kwdict={},**kw):
  "True RMS vz."
  kw.update(kwdict)
  kw['titlet']="True RMS vz"
  kw['titlel']="(m/s)"
  hpbasicwin('hvzrms',iw,kw)
if sys.version[:5] != "1.5.1":
  hpvzrms.__doc__ = hpvzrms.__doc__ + hpbasicwintext


def hpxpsqbar(iw=0,kwdict={},**kw):
  "Mean x' squared."
  kw.update(kwdict)
  kw['titlet']="Mean x' squared"
  kw['titlel']="(rad**2)"
  hpbasicwin('hxpsqbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxpsqbar.__doc__ = hpxpsqbar.__doc__ + hpbasicwintext


def hpypsqbar(iw=0,kwdict={},**kw):
  "Mean y' squared."
  kw.update(kwdict)
  kw['titlet']="Mean y' squared"
  kw['titlel']="(rad**2)"
  hpbasicwin('hypsqbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpypsqbar.__doc__ = hpypsqbar.__doc__ + hpbasicwintext


def hpxxpbar(iw=0,kwdict={},**kw):
  "Mean x*x'."
  kw.update(kwdict)
  kw['titlet']="Mean x*x'"
  kw['titlel']="(m-rad)"
  hpbasicwin('hxxpbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxxpbar.__doc__ = hpxxpbar.__doc__ + hpbasicwintext


def hpyypbar(iw=0,kwdict={},**kw):
  "Mean y*y'."
  kw.update(kwdict)
  kw['titlet']="Mean y*y'"
  kw['titlel']="(m-rad)"
  hpbasicwin('hyypbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpyypbar.__doc__ = hpyypbar.__doc__ + hpbasicwintext


def hpxypbar(iw=0,kwdict={},**kw):
  "Mean x*y'."
  kw.update(kwdict)
  kw['titlet']="Mean x*y'"
  kw['titlel']="(m-rad)"
  hpbasicwin('hxypbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxypbar.__doc__ = hpxypbar.__doc__ + hpbasicwintext


def hpyxpbar(iw=0,kwdict={},**kw):
  "Mean y*x'."
  kw.update(kwdict)
  kw['titlet']="Mean y*x'"
  kw['titlel']="(m**2)"
  hpbasicwin('hyxpbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpyxpbar.__doc__ = hpyxpbar.__doc__ + hpbasicwintext


def hpxpypbar(iw=0,kwdict={},**kw):
  "Mean x'*y'."
  kw.update(kwdict)
  kw['titlet']="Mean x'*y'"
  kw['titlel']="(rad**2)"
  hpbasicwin('hxpypbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxpypbar.__doc__ = hpxpypbar.__doc__ + hpbasicwintext


def hpxvzbar(iw=0,kwdict={},**kw):
  "Mean x*vz."
  kw.update(kwdict)
  kw['titlet']="Mean x*vz"
  kw['titlel']="(m-m/s)"
  hpbasicwin('hxvzbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpxvzbar.__doc__ = hpxvzbar.__doc__ + hpbasicwintext


def hpyvzbar(iw=0,kwdict={},**kw):
  "Mean y*vz."
  kw.update(kwdict)
  kw['titlet']="Mean y*vz"
  kw['titlel']="(m-m/s)"
  hpbasicwin('hyvzbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpyvzbar.__doc__ = hpyvzbar.__doc__ + hpbasicwintext


def hpvxvzbar(iw=0,kwdict={},**kw):
  "Mean vx*vz."
  kw.update(kwdict)
  kw['titlet']="Mean vx*vz"
  kw['titlel']="((m/s)**2)"
  hpbasicwin('hvxvzbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpvxvzbar.__doc__ = hpvxvzbar.__doc__ + hpbasicwintext


def hpvyvzbar(iw=0,kwdict={},**kw):
  "Mean vy*vz."
  kw.update(kwdict)
  kw['titlet']="Mean vy*vz"
  kw['titlel']="((m/s)**2)"
  hpbasicwin('hvyvzbar',iw,kw)
if sys.version[:5] != "1.5.1":
  hpvyvzbar.__doc__ = hpvyvzbar.__doc__ + hpbasicwintext


def hplinechg(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Line charge density."
  lhlinechg = _extractvarkw('lhlinechg',kw)
  if not lhlinechg: return
  kw.update(kwdict)
  kw['titlet']="Line charge density"
  hpzarray('hlinechg',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hplinechg.__doc__ = hplinechg.__doc__ + hpzarraytext


def hpvzofz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz versus space and time."
  lhvzofz = _extractvarkw('lhvzofz',kw)
  if not lhvzofz: return
  kw.update(kwdict)
  kw['titlet']="Vz versus space and time"
  hpzarray('hvzofz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvzofz.__doc__ = hpvzofz.__doc__ + hpzarraytext


def hpepsxz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X emittance."
  lhepsxz = _extractvarkw('lhepsxz',kw)
  if not lhepsxz: return
  kw.update(kwdict)
  kw['titlet']="X emittance"
  hpzarray('hepsxz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsxz.__doc__ = hpepsxz.__doc__ + hpzarraytext


def hpepsyz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y emittance."
  lhepsyz = _extractvarkw('lhepsyz',kw)
  if not lhepsyz: return
  kw.update(kwdict)
  kw['titlet']="Y emittance"
  hpzarray('hepsyz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsyz.__doc__ = hpepsyz.__doc__ + hpzarraytext


def hpepsnxz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X normalized emittance."
  lhepsnxz = _extractvarkw('lhepsnxz',kw)
  if not lhepsnxz: return
  kw.update(kwdict)
  kw['titlet']="X normalized emittance"
  hpzarray('hepsnxz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsnxz.__doc__ = hpepsnxz.__doc__ + hpzarraytext


def hpepsnyz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y normalized emittance."
  lhepsnyz = _extractvarkw('lhepsnyz',kw)
  if not lhepsnyz: return
  kw.update(kwdict)
  kw['titlet']="Y normalized emittance"
  hpzarray('hepsnyz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsnyz.__doc__ = hpepsnyz.__doc__ + hpzarraytext


def hpepsgz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized emittance."
  lhepsgz = _extractvarkw('lhepsgz',kw)
  if not lhepsgz: return
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  hpzarray('hepsgz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsgz.__doc__ = hpepsgz.__doc__ + hpzarraytext


def hpepshz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized emittance."
  lhepshz = _extractvarkw('lhepshz',kw)
  if not lhepshz: return
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  hpzarray('hepshz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepshz.__doc__ = hpepshz.__doc__ + hpzarraytext


def hpepsngz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized nrmlzd emittance."
  lhepsngz = _extractvarkw('lhepsngz',kw)
  if not lhepsngz: return
  kw.update(kwdict)
  kw['titlet']="Generalized nrmlzd emittance"
  hpzarray('hepsngz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsngz.__doc__ = hpepsngz.__doc__ + hpzarraytext


def hpepsnhz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized nrmlzd emittance."
  lhepsnhz = _extractvarkw('lhepsnhz',kw)
  if not lhepsnhz: return
  kw.update(kwdict)
  kw['titlet']="Generalized nrmlzd emittance"
  hpzarray('hepsnhz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsnhz.__doc__ = hpepsnhz.__doc__ + hpzarraytext


def hpxbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X bar."
  lhxbarz = _extractvarkw('lhxbarz',kw)
  if not lhxbarz: return
  kw.update(kwdict)
  kw['titlet']="X bar"
  hpzarray('hxbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxbarz.__doc__ = hpxbarz.__doc__ + hpzarraytext


def hpybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y bar."
  lhybarz = _extractvarkw('lhybarz',kw)
  if not lhybarz: return
  kw.update(kwdict)
  kw['titlet']="Y bar"
  hpzarray('hybarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpybarz.__doc__ = hpybarz.__doc__ + hpzarraytext


def hpxybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XY bar."
  lhxybarz = _extractvarkw('lhxybarz',kw)
  if not lhxybarz: return
  kw.update(kwdict)
  kw['titlet']="XY bar"
  hpzarray('hxybarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxybarz.__doc__ = hpxybarz.__doc__ + hpzarraytext


def hpxrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X rms."
  lhxrmsz = _extractvarkw('lhxrmsz',kw)
  if not lhxrmsz: return
  kw.update(kwdict)
  kw['titlet']="X rms"
  hpzarray('hxrmsz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxrmsz.__doc__ = hpxrmsz.__doc__ + hpzarraytext


def hpyrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y rms."
  lhyrmsz = _extractvarkw('lhyrmsz',kw)
  if not lhyrmsz: return
  kw.update(kwdict)
  kw['titlet']="Y rms"
  hpzarray('hyrmsz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyrmsz.__doc__ = hpyrmsz.__doc__ + hpzarraytext


def hpxprmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X' rms."
  lhxprmsz = _extractvarkw('lhxprmsz',kw)
  if not lhxprmsz: return
  kw.update(kwdict)
  kw['titlet']="X' rms"
  hpzarray('hxprmsz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxprmsz.__doc__ = hpxprmsz.__doc__ + hpzarraytext


def hpyprmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y' rms."
  lhyprmsz = _extractvarkw('lhyprmsz',kw)
  if not lhyprmsz: return
  kw.update(kwdict)
  kw['titlet']="Y' rms"
  hpzarray('hyprmsz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyprmsz.__doc__ = hpyprmsz.__doc__ + hpzarraytext


def hpxsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X**2 bar."
  lhxsqbarz = _extractvarkw('lhxsqbarz',kw)
  if not lhxsqbarz: return
  kw.update(kwdict)
  kw['titlet']="X**2 bar"
  hpzarray('hxsqbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxsqbarz.__doc__ = hpxsqbarz.__doc__ + hpzarraytext


def hpysqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y**2 bar."
  lhysqbarz = _extractvarkw('lhysqbarz',kw)
  if not lhysqbarz: return
  kw.update(kwdict)
  kw['titlet']="Y**2 bar"
  hpzarray('hysqbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpysqbarz.__doc__ = hpysqbarz.__doc__ + hpzarraytext


def hpvxbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vx bar."
  lhvxbarz = _extractvarkw('lhvxbarz',kw)
  if not lhvxbarz: return
  kw.update(kwdict)
  kw['titlet']="Vx bar"
  hpzarray('hvxbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvxbarz.__doc__ = hpvxbarz.__doc__ + hpzarraytext


def hpvybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vy bar."
  lhvybarz = _extractvarkw('lhvybarz',kw)
  if not lhvybarz: return
  kw.update(kwdict)
  kw['titlet']="Vy bar"
  hpzarray('hvybarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvybarz.__doc__ = hpvybarz.__doc__ + hpzarraytext


def hpvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz bar."
  lhvzbarz = _extractvarkw('lhvzbarz',kw)
  if not lhvzbarz: return
  kw.update(kwdict)
  kw['titlet']="Vz bar"
  hpzarray('hvzbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvzbarz.__doc__ = hpvzbarz.__doc__ + hpzarraytext


def hpxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X' bar."
  lhxpbarz = _extractvarkw('lhxpbarz',kw)
  if not lhxpbarz: return
  kw.update(kwdict)
  kw['titlet']="X' bar"
  hpzarray('hxpbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxpbarz.__doc__ = hpxpbarz.__doc__ + hpzarraytext


def hpypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y' bar."
  lhypbarz = _extractvarkw('lhypbarz',kw)
  if not lhypbarz: return
  kw.update(kwdict)
  kw['titlet']="Y' bar"
  hpzarray('hypbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpypbarz.__doc__ = hpypbarz.__doc__ + hpzarraytext


def hpvxrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vx rms."
  lhvxrmsz = _extractvarkw('lhvxrmsz',kw)
  if not lhvxrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vx rms"
  hpzarray('hvxrmsz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvxrmsz.__doc__ = hpvxrmsz.__doc__ + hpzarraytext


def hpvyrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vy rms."
  lhvyrmsz = _extractvarkw('lhvyrmsz',kw)
  if not lhvyrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vy rms"
  hpzarray('hvyrmsz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvyrmsz.__doc__ = hpvyrmsz.__doc__ + hpzarraytext


def hpvzrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz rms."
  lhvzrmsz = _extractvarkw('lhvzrmsz',kw)
  if not lhvzrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vz rms"
  hpzarray('hvzrmsz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvzrmsz.__doc__ = hpvzrmsz.__doc__ + hpzarraytext


def hpxpsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X'**2 bar."
  lhxpsqbarz = _extractvarkw('lhxpsqbarz',kw)
  if not lhxpsqbarz: return
  kw.update(kwdict)
  kw['titlet']="X'**2 bar"
  hpzarray('hxpsqbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxpsqbarz.__doc__ = hpxpsqbarz.__doc__ + hpzarraytext


def hpypsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y'**2 bar."
  lhypsqbarz = _extractvarkw('lhypsqbarz',kw)
  if not lhypsqbarz: return
  kw.update(kwdict)
  kw['titlet']="Y'**2 bar"
  hpzarray('hypsqbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpypsqbarz.__doc__ = hpypsqbarz.__doc__ + hpzarraytext


def hpxxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XX' bar."
  lhxxpbarz = _extractvarkw('lhxxpbarz',kw)
  if not lhxxpbarz: return
  kw.update(kwdict)
  kw['titlet']="XX' bar"
  hpzarray('hxxpbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxxpbarz.__doc__ = hpxxpbarz.__doc__ + hpzarraytext


def hpyypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YY' bar."
  lhyypbarz = _extractvarkw('lhyypbarz',kw)
  if not lhyypbarz: return
  kw.update(kwdict)
  kw['titlet']="YY' bar"
  hpzarray('hyypbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyypbarz.__doc__ = hpyypbarz.__doc__ + hpzarraytext


def hpxypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XY' bar."
  lhxypbarz = _extractvarkw('lhxypbarz',kw)
  if not lhxypbarz: return
  kw.update(kwdict)
  kw['titlet']="XY' bar"
  hpzarray('hxypbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxypbarz.__doc__ = hpxypbarz.__doc__ + hpzarraytext


def hpyxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YX' bar."
  lhyxpbarz = _extractvarkw('lhyxpbarz',kw)
  if not lhyxpbarz: return
  kw.update(kwdict)
  kw['titlet']="YX' bar"
  hpzarray('hyxpbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyxpbarz.__doc__ = hpyxpbarz.__doc__ + hpzarraytext


def hpxpypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X'Y' bar."
  lhxpypbarz = _extractvarkw('lhxpypbarz',kw)
  if not lhxpypbarz: return
  kw.update(kwdict)
  kw['titlet']="X'Y' bar"
  hpzarray('hxpypbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxpypbarz.__doc__ = hpxpypbarz.__doc__ + hpzarraytext


def hpxvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XVz bar."
  lhxvzbarz = _extractvarkw('lhxvzbarz',kw)
  if not lhxvzbarz: return
  kw.update(kwdict)
  kw['titlet']="XVz bar"
  hpzarray('hxvzbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxvzbarz.__doc__ = hpxvzbarz.__doc__ + hpzarraytext


def hpyvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YVz bar."
  lhyvzbarz = _extractvarkw('lhyvzbarz',kw)
  if not lhyvzbarz: return
  kw.update(kwdict)
  kw['titlet']="YVz bar"
  hpzarray('hyvzbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyvzbarz.__doc__ = hpyvzbarz.__doc__ + hpzarraytext


def hpvxvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "VxVz bar."
  lhvxvzbarz = _extractvarkw('lhvxvzbarz',kw)
  if not lhvxvzbarz: return
  kw.update(kwdict)
  kw['titlet']="VxVz bar"
  hpzarray('hvxvzbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvxvzbarz.__doc__ = hpvxvzbarz.__doc__ + hpzarraytext


def hpvyvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "VyVz bar."
  lhvyvzbarz = _extractvarkw('lhvyvzbarz',kw)
  if not lhvyvzbarz: return
  kw.update(kwdict)
  kw['titlet']="VyVz bar"
  hpzarray('hvyvzbarz',contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvyvzbarz.__doc__ = hpvyvzbarz.__doc__ + hpzarraytext


def hptotalke(kwdict={},**kw):
  "Total Kinetic Energy."
  kw.update(kwdict)
  kw['titlet']="Total Kinetic Energy"
  kw['titlel']="(J)"
  hekzbeam = _extractvarkw('hekzbeam',kw)
  hekperp = _extractvarkw('hekperp',kw)
  hpbasic(hekzbeam+hekperp,kw)


def hptotale(kwdict={},**kw):
  "Total Energy."
  kw.update(kwdict)
  kw['titlet']="Total Energy"
  kw['titlel']="(J)"
  hefld = _extractvarkw('hefld',kw)
  hekzbeam = _extractvarkw('hekzbeam',kw)
  hekperp = _extractvarkw('hekperp',kw)
  hpbasic(hefld+hekzbeam+hekperp,kw)


def hpthermale(iw=0,kwdict={},**kw):
  "Z Thermal Energy."
  kw.update(kwdict)
  kw['titlet']="Z Thermal Energy"
  kw['titlel']="(J)"
  sm = _extractvarkw('sm',kw)
  sw = _extractvarkw('sw',kw)
  sp_fract = _extractvarkw('sp_fract',kw)
  hpnum = _extractvarkw('hpnum',kw)
  hvzrms = _extractvarkw('hvzrms',kw)
  hpbasicwin(0.5*sum(sm*sw*sp_fract)*hpnum*hvzrms**2,iw,kw)
if sys.version[:5] != "1.5.1":
  hpthermale.__doc__ = hpthermale.__doc__ + hpbasicwintext


def hpeps6d(iw=0,kwdict={},**kw):
  "6-D Emittance."
  kw.update(kwdict)
  kw['titlet']="6-D Emittance"
  kw['titlel']="((pi-m-rad)**3)"
  hepsx = _extractvarkw('hepsx',kw)
  hepsy = _extractvarkw('hepsy',kw)
  hepsz = _extractvarkw('hepsz',kw)
  hpbasicwin(hepsx*hepsy*hepsz,iw,kw)
if sys.version[:5] != "1.5.1":
  hpeps6d.__doc__ = hpeps6d.__doc__ + hpbasicwintext


def hpepst(iw=0,kwdict={},**kw):
  "Transverse Emittance."
  kw.update(kwdict)
  kw['titlet']="Transverse Emittance"
  kw['titlel']="(pi-m-rad)"
  hepsx = _extractvarkw('hepsx',kw)
  hepsy = _extractvarkw('hepsy',kw)
  hpbasicwin(sqrt(hepsx*hepsy),iw,kw)
if sys.version[:5] != "1.5.1":
  hpepst.__doc__ = hpepst.__doc__ + hpbasicwintext


def hpepsnt(iw=0,kwdict={},**kw):
  "Normalized Transverse Emittance."
  kw.update(kwdict)
  kw['titlet']="Normalized Transverse Emittance"
  kw['titlel']="(pi-mm-mrad)"
  hepsnx = _extractvarkw('hepsnx',kw)
  hepsny = _extractvarkw('hepsny',kw)
  hpbasicwin(sqrt(hepsnx*hepsny),iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsnt.__doc__ = hpepsnt.__doc__ + hpbasicwintext


def hpxedge(iw=0,kwdict={},**kw):
  "X Beam Edge."
  kw.update(kwdict)
  kw['titlet']="X Beam Edge"
  kw['titlel']="(m)"
  hxrms = _extractvarkw('hxrms',kw)
  hpbasicwin(2.*hxrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxedge.__doc__ = hpxedge.__doc__ + hpbasicwintext


def hpxpedge(iw=0,kwdict={},**kw):
  "X' at Beam Edge."
  kw.update(kwdict)
  kw['titlet']="X' at Beam Edge"
  kw['titlel']="(m)"
  hxxpbar = _extractvarkw('hxxpbar',kw)
  hxbar = _extractvarkw('hxbar',kw)
  hxpbar = _extractvarkw('hxpbar',kw)
  hxrms = _extractvarkw('hxrms',kw)
  xpedge = (hxxpbar-hxbar*hxpbar)/where(greater(hxrms,0.),hxrms,1.)
  hpbasicwin(xpedge,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxpedge.__doc__ = hpxpedge.__doc__ + hpbasicwintext


def hpyedge(iw=0,kwdict={},**kw):
  "Y Beam Edge."
  kw.update(kwdict)
  kw['titlet']="Y Beam Edge"
  kw['titlel']="(m)"
  hyrms = _extractvarkw('hyrms',kw)
  hpbasicwin(2.*hyrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyedge.__doc__ = hpyedge.__doc__ + hpbasicwintext


def hpypedge(iw=0,kwdict={},**kw):
  "Y' at Beam Edge."
  kw.update(kwdict)
  kw['titlet']="Y' at Beam Edge"
  kw['titlel']="(m)"
  hyypbar = _extractvarkw('hyypbar',kw)
  hybar = _extractvarkw('hybar',kw)
  hypbar = _extractvarkw('hypbar',kw)
  hyrms = _extractvarkw('hyrms',kw)
  ypedge = (hyypbar-hybar*hypbar)/where(greater(hyrms,0.),hyrms,1.)
  hpbasicwin(ypedge,iw,kw)
if sys.version[:5] != "1.5.1":
  hpypedge.__doc__ = hpypedge.__doc__ + hpbasicwintext


def hpxedges(iw=0,kwdict={},**kw):
  "X Beam Edges plus centroid."
  kw.update(kwdict)
  kw['titlet']="X Beam Edges plus centroid"
  kw['titlel']="(m)"
  hxrms = _extractvarkw('hxrms',kw)
  hxbar = _extractvarkw('hxbar',kw)
  hpbasicwin(+2.*hxrms+hxbar,iw,kw)
  hpbasicwin(-2.*hxrms+hxbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxedge.__doc__ = hpxedge.__doc__ + hpbasicwintext


def hpyedges(iw=0,kwdict={},**kw):
  "Y Beam Edges plus centroid."
  kw.update(kwdict)
  kw['titlet']="Y Beam Edges plus centroid"
  kw['titlel']="(m)"
  hyrms = _extractvarkw('hyrms',kw)
  hybar = _extractvarkw('hybar',kw)
  hpbasicwin(+2.*hyrms+hybar,iw,kw)
  hpbasicwin(-2.*hyrms+hybar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyedge.__doc__ = hpyedge.__doc__ + hpbasicwintext


def hpenvx(iw=0,kwdict={},**kw):
  "X Beam Edge."
  kw.update(kwdict)
  kw['titlet']="X Beam Edge"
  kw['titlel']="(m)"
  hxrms = _extractvarkw('hxrms',kw)
  hpbasicwin(2.*hxrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpenvx.__doc__ = hpenvx.__doc__ + hpbasicwintext


def hpenvy(iw=0,kwdict={},**kw):
  "Y Beam Edge."
  kw.update(kwdict)
  kw['titlet']="Y Beam Edge"
  kw['titlel']="(m)"
  hyrms = _extractvarkw('hyrms',kw)
  hpbasicwin(2.*hyrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpenvy.__doc__ = hpenvy.__doc__ + hpbasicwintext

# --- Plots of current
def hpcurr(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Current."
  lhcurrz = _extractvarkw('lhcurrz',kw)
  lhlinechg = _extractvarkw('lhlinechg',kw)
  lhvzofz = _extractvarkw('lhvzofz',kw)
  if lhcurrz:
    kw.update(kwdict)
    kw['titlet']="Current"
    hcurrz = _extractvarkw('hcurrz',kw)
    hpzarray(hcurrz,contour,overlay,iz,kw)
  elif lhlinechg and lhvzofz:
    kw.update(kwdict)
    kw['titlet']="Current"
    hlinechg = _extractvarkw('hlinechg',kw)
    hvzofz = _extractvarkw('hvzofz',kw)
    hpzarray(hlinechg*hvzofz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpcurr.__doc__ = hpcurr.__doc__ + hpzarraytext

def histplotsdoc():
  """
hpdoc(): Prints complete list of arguments for histplot routines
hptotalke(): Total Kinetic Energy
hptotale(): Total Energy
hpthermale(): Z Thermal Energy
hpeps6d(): 6-D Emittance
hpepst(): Transverse Emittance
hpepsnt(): Normalized Transverse Emittance
hpxedge(): X Beam Edge (twice rms)
hpxpedge(): X' Beam Edge
hpyedge(): Y Beam Edge (twice rms)
hpypedge(): Y' Beam Edge
hpxedges(): X Beam Edges plus centroid
hpyedges(): Y Beam Edges plus centroid
hpenvx = hpxedge
hpenvy = hpyedge
hpzbeam(): Beam frame location
hpvbeam(): Beam frame velocity
hpbmlen(): RMS beam length
hpefld(): Field energy
hpekzmbe(): Total Z Kinetic energy minus beam energy
hpekzbeam(): Z Kinetic energy in the beam frame
hpekperp(): Perp Kinetic energy
hpekinz(): Z Kinetic energy (units of MV)
hpekin(): Total kinetic energy (units of MV)
hpepsx(): X emittance
hpepsy(): Y emittance
hpepsz(): Z emittance
hpepsnx(): X normalized emittance
hpepsny(): Y normalized emittance
hpepsnz(): Z normalized emittance
hpepsg(): Generalized emittance
hpepsh(): Generalized emittance
hpepsng(): Generalized normalized emittance
hpepsnh(): Generalized normalized emittance
hppnum(): Number of particles
hprhomid(): Charge density on axis
hprhomax(): Charge density max
hpxbar(): True mean x
hpybar(): True mean y
hpxybar(): True mean xy
hpxrms(): True RMS x
hpyrms(): True RMS y
hpxprms(): True RMS x'
hpyprms(): True RMS y'
hpxsqbar(): Mean x squared
hpysqbar(): Mean y squared
hpvxbar(): Mean vx
hpvybar(): Mean vy
hpvzbar(): Mean vz
hpxpbar(): Mean x'
hpypbar(): Mean y'
hpvxrms(): True RMS vx
hpvyrms(): True RMS vy
hpvzrms(): True RMS vz
hpxpsqbar(): Mean x' squared
hpypsqbar(): Mean y' squared
hpxxpbar(): Mean x*x'
hpyypbar(): Mean y*y'
hpxypbar(): Mean x*y'
hpyxpbar(): Mean y*x'
hpxpypbar(): Mean x'*y'
hpxvzbar(): Mean x*vz
hpyvzbar(): Mean y*vz
hpvxvzbar(): Mean vx*vz
hpvyvzbar(): Mean vy*vz
hplinechg(): Line charge density
hpvzofz(): Vz versus space and time
hpepsxz(): X emittance
hpepsyz(): Y emittance
hpepsnxz(): X normalized emittance
hpepsnyz(): Y normalized emittance
hpepsgz(): Generalized emittance
hpepshz(): Generalized emittance
hpepsngz(): Generalized nrmlzd emittance
hpepsnhz(): Generalized nrmlzd emittance
hpxbarz(): X bar
hpybarz(): Y bar
hpxybarz(): XY bar
hpxrmsz(): X rms
hpyrmsz(): Y rms
hpxprmsz(): X' rms
hpyprmsz(): Y' rms
hpxsqbarz(): X**2 bar
hpysqbarz(): Y**2 bar
hpvxbarz(): Vx bar
hpvybarz(): Vy bar
hpvzbarz(): Vz bar
hpxpbarz(): X' bar
hpypbarz(): Y' bar
hpvxrmsz(): Vx rms
hpvyrmsz(): Vy rms
hpvzrmsz(): Vz rms
hpxpsqbarz(): X'**2 bar
hpypsqbarz(): Y'**2 bar
hpxxpbarz(): XX' bar
hpyypbarz(): YY' bar
hpxypbarz(): XY' bar
hpyxpbarz(): YX' bar
hpxpypbarz(): X'Y' bar
hpxvzbarz(): XVz bar
hpyvzbarz(): YVz bar
hpvxvzbarz(): VxVz bar
hpvyvzbarz(): VyVz bar
  """
  print histplotsdoc.__doc__

def histplotstest(**kw):
  """
Test all histplots
  """
  apply(hptotalke,(),kw);fma()
  apply(hptotalke,(),kw);fma()
  apply(hptotale,(),kw);fma()
  apply(hpthermale,(),kw);fma()
  apply(hpeps6d,(),kw);fma()
  apply(hpepst,(),kw);fma()
  apply(hpepsnt,(),kw);fma()
  apply(hpxedge,(),kw);fma()
  apply(hpxpedge,(),kw);fma()
  apply(hpyedge,(),kw);fma()
  apply(hpypedge,(),kw);fma()
  apply(hpxedges,(),kw);fma()
  apply(hpyedges,(),kw);fma()
  apply(hpenvx,(),kw);fma()
  apply(hpenvy,(),kw);fma()
  apply(hpzbeam,(),kw);fma()
  apply(hpvbeam,(),kw);fma()
  apply(hpbmlen,(),kw);fma()
  apply(hpefld,(),kw);fma()
  apply(hpekzmbe,(),kw);fma()
  apply(hpekzbeam,(),kw);fma()
  apply(hpekperp,(),kw);fma()
  apply(hpekinz,(),kw);fma()
  apply(hpekin,(),kw);fma()
  apply(hpepsx,(),kw);fma()
  apply(hpepsy,(),kw);fma()
  apply(hpepsz,(),kw);fma()
  apply(hpepsnx,(),kw);fma()
  apply(hpepsny,(),kw);fma()
  apply(hpepsnz,(),kw);fma()
  apply(hpepsg,(),kw);fma()
  apply(hpepsh,(),kw);fma()
  apply(hpepsng,(),kw);fma()
  apply(hpepsnh,(),kw);fma()
  apply(hppnum,(),kw);fma()
  apply(hprhomid,(),kw);fma()
  apply(hprhomax,(),kw);fma()
  apply(hpxbar,(),kw);fma()
  apply(hpybar,(),kw);fma()
  apply(hpxybar,(),kw);fma()
  apply(hpxrms,(),kw);fma()
  apply(hpyrms,(),kw);fma()
  apply(hpxprms,(),kw);fma()
  apply(hpyprms,(),kw);fma()
  apply(hpxsqbar,(),kw);fma()
  apply(hpysqbar,(),kw);fma()
  apply(hpvxbar,(),kw);fma()
  apply(hpvybar,(),kw);fma()
  apply(hpvzbar,(),kw);fma()
  apply(hpxpbar,(),kw);fma()
  apply(hpypbar,(),kw);fma()
  apply(hpvxrms,(),kw);fma()
  apply(hpvyrms,(),kw);fma()
  apply(hpvzrms,(),kw);fma()
  apply(hpxpsqbar,(),kw);fma()
  apply(hpypsqbar,(),kw);fma()
  apply(hpxxpbar,(),kw);fma()
  apply(hpyypbar,(),kw);fma()
  apply(hpxypbar,(),kw);fma()
  apply(hpyxpbar,(),kw);fma()
  apply(hpxpypbar,(),kw);fma()
  apply(hpxvzbar,(),kw);fma()
  apply(hpyvzbar,(),kw);fma()
  apply(hpvxvzbar,(),kw);fma()
  apply(hpvyvzbar,(),kw);fma()
  apply(hplinechg,(),kw);fma()
  apply(hpvzofz,(),kw);fma()
  apply(hpepsxz,(),kw);fma()
  apply(hpepsyz,(),kw);fma()
  apply(hpepsnxz,(),kw);fma()
  apply(hpepsnyz,(),kw);fma()
  apply(hpepsgz,(),kw);fma()
  apply(hpepshz,(),kw);fma()
  apply(hpepsngz,(),kw);fma()
  apply(hpepsnhz,(),kw);fma()
  apply(hpxbarz,(),kw);fma()
  apply(hpybarz,(),kw);fma()
  apply(hpxybarz,(),kw);fma()
  apply(hpxrmsz,(),kw);fma()
  apply(hpyrmsz,(),kw);fma()
  apply(hpxprmsz,(),kw);fma()
  apply(hpyprmsz,(),kw);fma()
  apply(hpxsqbarz,(),kw);fma()
  apply(hpysqbarz,(),kw);fma()
  apply(hpvxbarz,(),kw);fma()
  apply(hpvybarz,(),kw);fma()
  apply(hpvzbarz,(),kw);fma()
  apply(hpxpbarz,(),kw);fma()
  apply(hpypbarz,(),kw);fma()
  apply(hpvxrmsz,(),kw);fma()
  apply(hpvyrmsz,(),kw);fma()
  apply(hpvzrmsz,(),kw);fma()
  apply(hpxpsqbarz,(),kw);fma()
  apply(hpypsqbarz,(),kw);fma()
  apply(hpxxpbarz,(),kw);fma()
  apply(hpyypbarz,(),kw);fma()
  apply(hpxypbarz,(),kw);fma()
  apply(hpyxpbarz,(),kw);fma()
  apply(hpxpypbarz,(),kw);fma()
  apply(hpxvzbarz,(),kw);fma()
  apply(hpyvzbarz,(),kw);fma()
  apply(hpvxvzbarz,(),kw);fma()
  apply(hpvyvzbarz,(),kw);fma()

