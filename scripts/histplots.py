from warp import *
from mplot import *
histplots_version = "$Id: histplots.py,v 1.11 2002/01/24 23:08:01 dave Exp $"

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
  - iend=top.jhist: Time index at which to end the plots
  - istep=1: Step size of time data plotted
  - lhzbeam=1: When true, plots data versus hzbeam instead of thist
  - lvsz=1: Old name for lhzbeam
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
    see work.gs for the plot system values."""
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
  - iend=top.jhist: Time index at which to end the plots
  - istep=max(iend/20,1): Step size of time data plotted
  - jstart=0: Z index at which to start the plots
  - jend=top.nzzarr: Z index at which to end the plots
  - jstep=max(jend/32,1): Z step size of time data plotted
  - lhzbeam=1: When true, plots data versus hzbeam instead of thist
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
  - filled=0: When true, plot filled contours""")
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
                'lnormalized':0,'istart':0,'iend':top.jhist,'istep':1,
                'lhzbeam':1,'lvsz':1,'logplot':0,
                'color':'fg','marks':0,'marker':None,'msize':1.0,'titles':1,
                'plsysval':1,'width':1.,'linetype':'solid'}
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  badargs = checkarguments(kwvalues,kwdefaults)
  if badargs: raise "bad argument ",string.join(badargs.keys())
  for arg in kwvalues.keys(): exec(arg+" = kwvalues['"+arg+"']")

  # --- Now complete the setup
  if lnormalized:
    yscale = 1./oord[...,istart]
    if not titlel: titlel = titlet + " over initial value"
  if type(yscale) != type(array([])): yscale = array([yscale])
  if not absc:
    if (not lhzbeam or not lvsz):
      absc = top.thist[istart:iend+1:istep]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "time (s)"
        else:
          titleb = "time"
    else:
      absc = top.hzbeam[istart:iend+1:istep]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "Z (m)"
        else:
          titleb = "Z"
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
  if 'lzshift' in kw.keys():
    lzshift = kw['lzshift']
    del kw['lzshift']
  else:
    lzshift = 0
  if lzshift: kw['xoffset'] = 0.5*(top.zwindows[0,iw]+top.zwindows[1,iw])
  if iw == 0:
    kw['titlet'] = kw['titlet'] + " for whole beam"
    hpbasic(oord[iw,:],kw)
  elif iw > 0:
    kw['titler']="z%1d = %6.3f"%(iw,.5*(top.zwindows[0,iw]+top.zwindows[1,iw]))
    hpbasic(oord[iw,:],kw)
  else:
    for i in range(1,top.nzwind):
      plsys(3+(i-1)%4)
      kw['titler'] = "z%1d = %6.3f"%(i,.5*(top.zwindows[0,i]+top.zwindows[1,i]))
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
                'istart':0,'iend':top.jhist,'istep':None,'jstart':0,
                'jend':top.nzzarr,'jstep':None,'lhzbeam':1,'logplot':0,
                'color':'fg','marks':0,'marker':None,'msize':1.0,
                'titles':1,'levs':10,'filled':0,'width':1.,'linetype':'solid'}
  kwvalues = kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  for arg in kwdefaults.keys(): exec(arg+" = kwvalues['"+arg+"']")
  badargs = checkarguments(kwvalues,kwdefaults)
  if badargs: raise "bad argument ",string.join(badargs.keys())

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
    if (lhzbeam):
      absc = top.hzbeam[istart:iend+1:istep]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "Z (m)"
        else:
          titleb = "Z"
    else:
      absc = top.thist[istart:iend+1:istep]*xscale + xoffset
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
  kw['titleb']="Z"
  if iz:
    if iz == w3d.nz/2:
      kw['titlet']=kw['titlet']+" at Center"
    else:
      kw['titlet']=kw['titlet']+" at iz = %d"%iz
    kw['titler'] = 'z = %6.3f'%top.zmntmesh[iz]
    hpbasic(hzarray[iz,:],kw)
    return
  elif contour:
    hpbasiccont(hzarray,top.zmntmesh,kw)
  elif overlay:
    kw['istep'] = max(top.jhist/10,1)
    kw['ord'] = w3d.zmesh
    kw['jhist'] = top.jhist
    mountainplot1(kw['titlet'],hzarray,kw)
  else:
    kw['jhist'] = top.jhist
    mountainplot1(kw['titlet'],hzarray,kw)

###########################################################################
###########################################################################
###########################################################################

def hpzbeam(kwdict={},**kw):
  "Beam frame location."
  kw.update(kwdict)
  kw['titlet']="Beam frame location"
  kw['titlel']="(m)"
  hpbasic(top.hzbeam,kw)
if sys.version[:5] != "1.5.1":
  hpzbeam.__doc__ = hpzbeam.__doc__ + hpbasictext


def hpvbeam(kwdict={},**kw):
  "Beam frame velocity."
  kw.update(kwdict)
  kw['titlet']="Beam frame velocity"
  kw['titlel']="(m/s)"
  hpbasic(top.hvbeam,kw)
if sys.version[:5] != "1.5.1":
  hpvbeam.__doc__ = hpvbeam.__doc__ + hpbasictext


def hpbmlen(kwdict={},**kw):
  "RMS beam length."
  kw.update(kwdict)
  kw['titlet']="RMS beam length"
  kw['titlel']="(m)"
  hpbasic(top.hbmlen,kw)
if sys.version[:5] != "1.5.1":
  hpbmlen.__doc__ = hpbmlen.__doc__ + hpbasictext


def hpefld(kwdict={},**kw):
  "Field energy."
  kw.update(kwdict)
  kw['titlet']="Field energy"
  kw['titlel']="(J)"
  hpbasic(top.hefld,kw)
if sys.version[:5] != "1.5.1":
  hpefld.__doc__ = hpefld.__doc__ + hpbasictext


def hpekzmbe(kwdict={},**kw):
  "Total Z Kinetic energy minus beam energy."
  kw.update(kwdict)
  kw['titlet']="Total Z Kinetic energy minus beam energy"
  kw['titlel']="(J)"
  hpbasic(top.hekzmbe,kw)
if sys.version[:5] != "1.5.1":
  hpekzmbe.__doc__ = hpekzmbe.__doc__ + hpbasictext


def hpekzbeam(kwdict={},**kw):
  "Z Kinetic energy in the beam frame."
  kw.update(kwdict)
  kw['titlet']="Z Kinetic energy in the beam frame"
  kw['titlel']="(J)"
  hpbasic(top.hekzbeam,kw)
if sys.version[:5] != "1.5.1":
  hpekzbeam.__doc__=hpekzbeam.__doc__+hpbasictext


def hpekperp(kwdict={},**kw):
  "Perp Kinetic energy."
  kw.update(kwdict)
  kw['titlet']="Perp Kinetic energy"
  kw['titlel']="(J)"
  hpbasic(top.hekperp,kw)
if sys.version[:5] != "1.5.1":
  hpekperp.__doc__ = hpekperp.__doc__ + hpbasictext


def hpekinz(iw=0,kwdict={},**kw):
  "Z Kinetic energy (units of MV)"
  kw.update(kwdict)
  kw['titlet']="Z Kinetic energy"
  kw['titlel']="(MV)"
  if top.lrelativ:
    u = (top.hvzbar/clight)**2
    e = (top.aion*amu*clight**2)*(u/(sqrt(1.-u)+1.-u))/jperev              
  else:
    e = 0.5*top.aion*amu*top.hvzbar**2/jperev
  hpbasicwin(e*1.e-6,iw,kw)
if sys.version[:5] != "1.5.1":
  hpekinz.__doc__ = hpekinz.__doc__ + hpbasicwintext


def hpekin(iw=0,kwdict={},**kw):
  "Total Kinetic energy (units of MV)"
  kw.update(kwdict)
  kw['titlet']="Total Kinetic energy"
  kw['titlel']="(MV)"
  if top.lrelativ:
    u = (top.hvxbar**2 + top.hvybar**2 + top.hvzbar**2)/clight**2
    e = (top.aion*amu*clight**2)*(u/(sqrt(1.-u)+1.-u))/jperev              
  else:
    u = top.hvxbar**2 + top.hvybar**2 + top.hvzbar**2
    e = 0.5*top.aion*amu*u/jperev
  hpbasicwin(e*1.e-6,iw,kw)
if sys.version[:5] != "1.5.1":
  hpekin.__doc__ = hpekin.__doc__ + hpbasicwintext


def hpepsx(iw=0,kwdict={},**kw):
  "X emittance."
  kw.update(kwdict)
  kw['titlet']="X emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsx,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsx.__doc__ = hpepsx.__doc__ + hpbasicwintext


def hpepsy(iw=0,kwdict={},**kw):
  "Y emittance."
  kw.update(kwdict)
  kw['titlet']="Y emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsy,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsy.__doc__ = hpepsy.__doc__ + hpbasicwintext


def hpepsz(iw=0,kwdict={},**kw):
  "Z emittance."
  kw.update(kwdict)
  kw['titlet']="Z emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsz,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsz.__doc__ = hpepsz.__doc__ + hpbasicwintext


def hpepsnx(iw=0,kwdict={},**kw):
  "X normalized emittance."
  kw.update(kwdict)
  kw['titlet']="X normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsnx,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsnx.__doc__ = hpepsnx.__doc__ + hpbasicwintext


def hpepsny(iw=0,kwdict={},**kw):
  "Y normalized emittance."
  kw.update(kwdict)
  kw['titlet']="Y normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsny,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsny.__doc__ = hpepsny.__doc__ + hpbasicwintext


def hpepsnz(iw=0,kwdict={},**kw):
  "Z normalized emittance."
  kw.update(kwdict)
  kw['titlet']="Z normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsnz,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsnz.__doc__ =  hpepsnz.__doc__ + hpbasicwintext


def hpepsg(iw=0,kwdict={},**kw):
  "Generalized emittance."
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsg,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsg.__doc__ = hpepsg.__doc__ + hpbasicwintext


def hpepsh(iw=0,kwdict={},**kw):
  "Generalized emittance."
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsh,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsh.__doc__ = hpepsh.__doc__ + hpbasicwintext


def hpepsng(iw=0,kwdict={},**kw):
  "Generalized normalized emittance."
  kw.update(kwdict)
  kw['titlet']="Generalized normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsng,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsng.__doc__ = hpepsng.__doc__ + hpbasicwintext


def hpepsnh(iw=0,kwdict={},**kw):
  "Generalized normalized emittance."
  kw.update(kwdict)
  kw['titlet']="Generalized normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsnh,iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsnh.__doc__ = hpepsnh.__doc__ + hpbasicwintext


def hppnum(iw=0,kwdict={},**kw):
  "Number of particles."
  kw.update(kwdict)
  kw['titlet']="Number of particles"
  kw['titlel']=""
  hpbasicwin(top.hpnum,iw,kw)
if sys.version[:5] != "1.5.1":
  hppnum.__doc__ = hppnum.__doc__ + hpbasicwintext


def hprhomid(iw=0,kwdict={},**kw):
  "Charge density on axis."
  kw.update(kwdict)
  kw['titlet']="Charge density on axis"
  kw['titlel']="(C/m**3)"
  hpbasicwin(top.hrhomid,iw,kw)
if sys.version[:5] != "1.5.1":
  hprhomid.__doc__ = hprhomid.__doc__ + hpbasicwintext


def hprhomax(iw=0,kwdict={},**kw):
  "Charge density max."
  kw.update(kwdict)
  kw['titlet']="Charge density max"
  kw['titlel']="(C/m**3)"
  hpbasicwin(top.hrhomax,iw,kw)
if sys.version[:5] != "1.5.1":
  hprhomax.__doc__ = hprhomax.__doc__ + hpbasicwintext


def hpxbar(iw=0,kwdict={},**kw):
  "True mean x."
  kw.update(kwdict)
  kw['titlet']="True mean x"
  kw['titlel']="(m)"
  hpbasicwin(top.hxbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxbar.__doc__ = hpxbar.__doc__ + hpbasicwintext


def hpybar(iw=0,kwdict={},**kw):
  "True mean y."
  kw.update(kwdict)
  kw['titlet']="True mean y"
  kw['titlel']="(m)"
  hpbasicwin(top.hybar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpybar.__doc__ = hpybar.__doc__ + hpbasicwintext


def hpxybar(iw=0,kwdict={},**kw):
  "True mean xy."
  kw.update(kwdict)
  kw['titlet']="True mean xy"
  kw['titlel']="(m**2)"
  hpbasicwin(top.hxybar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxybar.__doc__ = hpxybar.__doc__ + hpbasicwintext


def hpxrms(iw=0,kwdict={},**kw):
  "True RMS x."
  kw.update(kwdict)
  kw['titlet']="True RMS x"
  kw['titlel']="(m)"
  hpbasicwin(top.hxrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxrms.__doc__ = hpxrms.__doc__ + hpbasicwintext


def hpyrms(iw=0,kwdict={},**kw):
  "True RMS y."
  kw.update(kwdict)
  kw['titlet']="True RMS y"
  kw['titlel']="(m)"
  hpbasicwin(top.hyrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyrms.__doc__ = hpyrms.__doc__ + hpbasicwintext


def hpxprms(iw=0,kwdict={},**kw):
  "True RMS x'."
  kw.update(kwdict)
  kw['titlet']="True RMS x'"
  kw['titlel']="(rad)"
  hpbasicwin(top.hxprms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxprms.__doc__ = hpxprms.__doc__ + hpbasicwintext


def hpyprms(iw=0,kwdict={},**kw):
  "True RMS y'."
  kw.update(kwdict)
  kw['titlet']="True RMS y'"
  kw['titlel']="(rad)"
  hpbasicwin(top.hyprms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyprms.__doc__ = hpyprms.__doc__ + hpbasicwintext


def hpxsqbar(iw=0,kwdict={},**kw):
  "Mean x squared."
  kw.update(kwdict)
  kw['titlet']="Mean x squared"
  kw['titlel']="(m**2)"
  hpbasicwin(top.hxsqbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxsqbar.__doc__ = hpxsqbar.__doc__ + hpbasicwintext


def hpysqbar(iw=0,kwdict={},**kw):
  "Mean y squared."
  kw.update(kwdict)
  kw['titlet']="Mean y squared"
  kw['titlel']="(m**2)"
  hpbasicwin(top.hysqbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpysqbar.__doc__ = hpysqbar.__doc__ + hpbasicwintext


def hpvxbar(iw=0,kwdict={},**kw):
  "Mean vx."
  kw.update(kwdict)
  kw['titlet']="Mean vx"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvxbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpvxbar.__doc__ = hpvxbar.__doc__ + hpbasicwintext


def hpvybar(iw=0,kwdict={},**kw):
  "Mean vy."
  kw.update(kwdict)
  kw['titlet']="Mean vy"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvybar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpvybar.__doc__ = hpvybar.__doc__ + hpbasicwintext


# --- plot mean z velocity in beam or lab frame vs time
def hpvzbar(iw=0,beamframe=1,kwdict={},**kw):
  """Mean Z Velocity (beam frame or lab frame)
  - beamframe=1: when true, plot Vz relative to beam frame (vzbar - vbeam)"""
  kw.update(kwdict)
  kw['titlel']="(m/s)"
  if beamframe:
    kw['titlet']="Mean Z Velocity (beam frame)"
    kw['titler']="vbeam = %6.3e"%top.vbeam
    hpbasicwin(top.hvzbar-top.vbeam,iw,kw)
  else:
    kw['titlet']="Mean Z Velocity"
    hpbasicwin(top.hvzbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpvzbar.__doc__ = hpvzbar.__doc__ + hpbasicwintext


def hpxpbar(iw=0,kwdict={},**kw):
  "Mean x'."
  kw.update(kwdict)
  kw['titlet']="Mean x'"
  kw['titlel']="(rad)"
  hpbasicwin(top.hxpbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxpbar.__doc__ = hpxpbar.__doc__ + hpbasicwintext


def hpypbar(iw=0,kwdict={},**kw):
  "Mean y'."
  kw.update(kwdict)
  kw['titlet']="Mean y'"
  kw['titlel']="(rad)"
  hpbasicwin(top.hypbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpypbar.__doc__ = hpypbar.__doc__ + hpbasicwintext


def hpvxrms(iw=0,kwdict={},**kw):
  "True RMS vx."
  kw.update(kwdict)
  kw['titlet']="True RMS vx"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvxrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpvxrms.__doc__ = hpvxrms.__doc__ + hpbasicwintext


def hpvyrms(iw=0,kwdict={},**kw):
  "True RMS vy."
  kw.update(kwdict)
  kw['titlet']="True RMS vy"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvyrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpvyrms.__doc__ = hpvyrms.__doc__ + hpbasicwintext


def hpvzrms(iw=0,kwdict={},**kw):
  "True RMS vz."
  kw.update(kwdict)
  kw['titlet']="True RMS vz"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvzrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpvzrms.__doc__ = hpvzrms.__doc__ + hpbasicwintext


def hpxpsqbar(iw=0,kwdict={},**kw):
  "Mean x' squared."
  kw.update(kwdict)
  kw['titlet']="Mean x' squared"
  kw['titlel']="(rad**2)"
  hpbasicwin(top.hxpsqbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxpsqbar.__doc__ = hpxpsqbar.__doc__ + hpbasicwintext


def hpypsqbar(iw=0,kwdict={},**kw):
  "Mean y' squared."
  kw.update(kwdict)
  kw['titlet']="Mean y' squared"
  kw['titlel']="(rad**2)"
  hpbasicwin(top.hypsqbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpypsqbar.__doc__ = hpypsqbar.__doc__ + hpbasicwintext


def hpxxpbar(iw=0,kwdict={},**kw):
  "Mean x*x'."
  kw.update(kwdict)
  kw['titlet']="Mean x*x'"
  kw['titlel']="(m-rad)"
  hpbasicwin(top.hxxpbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxxpbar.__doc__ = hpxxpbar.__doc__ + hpbasicwintext


def hpyypbar(iw=0,kwdict={},**kw):
  "Mean y*y'."
  kw.update(kwdict)
  kw['titlet']="Mean y*y'"
  kw['titlel']="(m-rad)"
  hpbasicwin(top.hyypbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyypbar.__doc__ = hpyypbar.__doc__ + hpbasicwintext


def hpxypbar(iw=0,kwdict={},**kw):
  "Mean x*y'."
  kw.update(kwdict)
  kw['titlet']="Mean x*y'"
  kw['titlel']="(m-rad)"
  hpbasicwin(top.hxypbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxypbar.__doc__ = hpxypbar.__doc__ + hpbasicwintext


def hpyxpbar(iw=0,kwdict={},**kw):
  "Mean y*x'."
  kw.update(kwdict)
  kw['titlet']="Mean y*x'"
  kw['titlel']="(m**2)"
  hpbasicwin(top.hyxpbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyxpbar.__doc__ = hpyxpbar.__doc__ + hpbasicwintext


def hpxpypbar(iw=0,kwdict={},**kw):
  "Mean x'*y'."
  kw.update(kwdict)
  kw['titlet']="Mean x'*y'"
  kw['titlel']="(rad**2)"
  hpbasicwin(top.hxpypbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxpypbar.__doc__ = hpxpypbar.__doc__ + hpbasicwintext


def hpxvzbar(iw=0,kwdict={},**kw):
  "Mean x*vz."
  kw.update(kwdict)
  kw['titlet']="Mean x*vz"
  kw['titlel']="(m-m/s)"
  hpbasicwin(top.hxvzbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxvzbar.__doc__ = hpxvzbar.__doc__ + hpbasicwintext


def hpyvzbar(iw=0,kwdict={},**kw):
  "Mean y*vz."
  kw.update(kwdict)
  kw['titlet']="Mean y*vz"
  kw['titlel']="(m-m/s)"
  hpbasicwin(top.hyvzbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyvzbar.__doc__ = hpyvzbar.__doc__ + hpbasicwintext


def hpvxvzbar(iw=0,kwdict={},**kw):
  "Mean vx*vz."
  kw.update(kwdict)
  kw['titlet']="Mean vx*vz"
  kw['titlel']="((m/s)**2)"
  hpbasicwin(top.hvxvzbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpvxvzbar.__doc__ = hpvxvzbar.__doc__ + hpbasicwintext


def hpvyvzbar(iw=0,kwdict={},**kw):
  "Mean vy*vz."
  kw.update(kwdict)
  kw['titlet']="Mean vy*vz"
  kw['titlel']="((m/s)**2)"
  hpbasicwin(top.hvyvzbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpvyvzbar.__doc__ = hpvyvzbar.__doc__ + hpbasicwintext


def hplinechg(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Line charge density."
  if not top.lhlinechg: return
  kw.update(kwdict)
  kw['titlet']="Line charge density"
  hpzarray(top.hlinechg,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hplinechg.__doc__ = hplinechg.__doc__ + hpzarraytext


def hpvzofz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz versus space and time."
  if not top.lhvzofz: return
  kw.update(kwdict)
  kw['titlet']="Vz versus space and time"
  hpzarray(top.hvzofz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvzofz.__doc__ = hpvzofz.__doc__ + hpzarraytext


def hpepsxz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X emittance."
  if not top.lhepsxz: return
  kw.update(kwdict)
  kw['titlet']="X emittance"
  hpzarray(top.hepsxz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsxz.__doc__ = hpepsxz.__doc__ + hpzarraytext


def hpepsyz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y emittance."
  if not top.lhepsyz: return
  kw.update(kwdict)
  kw['titlet']="Y emittance"
  hpzarray(top.hepsyz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsyz.__doc__ = hpepsyz.__doc__ + hpzarraytext


def hpepsnxz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X normalized emittance."
  if not top.lhepsnxz: return
  kw.update(kwdict)
  kw['titlet']="X normalized emittance"
  hpzarray(top.hepsnxz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsnxz.__doc__ = hpepsnxz.__doc__ + hpzarraytext


def hpepsnyz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y normalized emittance."
  if not top.lhepsnyz: return
  kw.update(kwdict)
  kw['titlet']="Y normalized emittance"
  hpzarray(top.hepsnyz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsnyz.__doc__ = hpepsnyz.__doc__ + hpzarraytext


def hpepsgz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized emittance."
  if not top.lhepsgz: return
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  hpzarray(top.hepsgz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsgz.__doc__ = hpepsgz.__doc__ + hpzarraytext


def hpepshz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized emittance."
  if not top.lhepshz: return
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  hpzarray(top.hepshz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepshz.__doc__ = hpepshz.__doc__ + hpzarraytext


def hpepsngz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized nrmlzd emittance."
  if not top.lhepsngz: return
  kw.update(kwdict)
  kw['titlet']="Generalized nrmlzd emittance"
  hpzarray(top.hepsngz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsngz.__doc__ = hpepsngz.__doc__ + hpzarraytext


def hpepsnhz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized nrmlzd emittance."
  if not top.lhepsnhz: return
  kw.update(kwdict)
  kw['titlet']="Generalized nrmlzd emittance"
  hpzarray(top.hepsnhz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpepsnhz.__doc__ = hpepsnhz.__doc__ + hpzarraytext


def hpxbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X bar."
  if not top.lhxbarz: return
  kw.update(kwdict)
  kw['titlet']="X bar"
  hpzarray(top.hxbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxbarz.__doc__ = hpxbarz.__doc__ + hpzarraytext


def hpybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y bar."
  if not top.lhybarz: return
  kw.update(kwdict)
  kw['titlet']="Y bar"
  hpzarray(top.hybarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpybarz.__doc__ = hpybarz.__doc__ + hpzarraytext


def hpxybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XY bar."
  if not top.lhxybarz: return
  kw.update(kwdict)
  kw['titlet']="XY bar"
  hpzarray(top.hxybarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxybarz.__doc__ = hpxybarz.__doc__ + hpzarraytext


def hpxrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X rms."
  if not top.lhxrmsz: return
  kw.update(kwdict)
  kw['titlet']="X rms"
  hpzarray(top.hxrmsz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxrmsz.__doc__ = hpxrmsz.__doc__ + hpzarraytext


def hpyrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y rms."
  if not top.lhyrmsz: return
  kw.update(kwdict)
  kw['titlet']="Y rms"
  hpzarray(top.hyrmsz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyrmsz.__doc__ = hpyrmsz.__doc__ + hpzarraytext


def hpxprmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X' rms."
  if not top.lhxprmsz: return
  kw.update(kwdict)
  kw['titlet']="X' rms"
  hpzarray(top.hxprmsz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxprmsz.__doc__ = hpxprmsz.__doc__ + hpzarraytext


def hpyprmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y' rms."
  if not top.lhyprmsz: return
  kw.update(kwdict)
  kw['titlet']="Y' rms"
  hpzarray(top.hyprmsz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyprmsz.__doc__ = hpyprmsz.__doc__ + hpzarraytext


def hpxsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X**2 bar."
  if not top.lhxsqbarz: return
  kw.update(kwdict)
  kw['titlet']="X**2 bar"
  hpzarray(top.hxsqbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxsqbarz.__doc__ = hpxsqbarz.__doc__ + hpzarraytext


def hpysqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y**2 bar."
  if not top.lhysqbarz: return
  kw.update(kwdict)
  kw['titlet']="Y**2 bar"
  hpzarray(top.hysqbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpysqbarz.__doc__ = hpysqbarz.__doc__ + hpzarraytext


def hpvxbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vx bar."
  if not top.lhvxbarz: return
  kw.update(kwdict)
  kw['titlet']="Vx bar"
  hpzarray(top.hvxbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvxbarz.__doc__ = hpvxbarz.__doc__ + hpzarraytext


def hpvybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vy bar."
  if not top.lhvybarz: return
  kw.update(kwdict)
  kw['titlet']="Vy bar"
  hpzarray(top.hvybarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvybarz.__doc__ = hpvybarz.__doc__ + hpzarraytext


def hpvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz bar."
  if not top.lhvzbarz: return
  kw.update(kwdict)
  kw['titlet']="Vz bar"
  hpzarray(top.hvzbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvzbarz.__doc__ = hpvzbarz.__doc__ + hpzarraytext


def hpxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X' bar."
  if not top.lhxpbarz: return
  kw.update(kwdict)
  kw['titlet']="X' bar"
  hpzarray(top.hxpbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxpbarz.__doc__ = hpxpbarz.__doc__ + hpzarraytext


def hpypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y' bar."
  if not top.lhypbarz: return
  kw.update(kwdict)
  kw['titlet']="Y' bar"
  hpzarray(top.hypbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpypbarz.__doc__ = hpypbarz.__doc__ + hpzarraytext


def hpvxrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vx rms."
  if not top.lhvxrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vx rms"
  hpzarray(top.hvxrmsz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvxrmsz.__doc__ = hpvxrmsz.__doc__ + hpzarraytext


def hpvyrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vy rms."
  if not top.lhvyrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vy rms"
  hpzarray(top.hvyrmsz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvyrmsz.__doc__ = hpvyrmsz.__doc__ + hpzarraytext


def hpvzrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz rms."
  if not top.lhvzrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vz rms"
  hpzarray(top.hvzrmsz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvzrmsz.__doc__ = hpvzrmsz.__doc__ + hpzarraytext


def hpxpsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X'**2 bar."
  if not top.lhxpsqbarz: return
  kw.update(kwdict)
  kw['titlet']="X'**2 bar"
  hpzarray(top.hxpsqbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxpsqbarz.__doc__ = hpxpsqbarz.__doc__ + hpzarraytext


def hpypsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y'**2 bar."
  if not top.lhypsqbarz: return
  kw.update(kwdict)
  kw['titlet']="Y'**2 bar"
  hpzarray(top.hypsqbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpypsqbarz.__doc__ = hpypsqbarz.__doc__ + hpzarraytext


def hpxxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XX' bar."
  if not top.lhxxpbarz: return
  kw.update(kwdict)
  kw['titlet']="XX' bar"
  hpzarray(top.hxxpbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxxpbarz.__doc__ = hpxxpbarz.__doc__ + hpzarraytext


def hpyypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YY' bar."
  if not top.lhyypbarz: return
  kw.update(kwdict)
  kw['titlet']="YY' bar"
  hpzarray(top.hyypbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyypbarz.__doc__ = hpyypbarz.__doc__ + hpzarraytext


def hpxypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XY' bar."
  if not top.lhxypbarz: return
  kw.update(kwdict)
  kw['titlet']="XY' bar"
  hpzarray(top.hxypbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxypbarz.__doc__ = hpxypbarz.__doc__ + hpzarraytext


def hpyxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YX' bar."
  if not top.lhyxpbarz: return
  kw.update(kwdict)
  kw['titlet']="YX' bar"
  hpzarray(top.hyxpbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyxpbarz.__doc__ = hpyxpbarz.__doc__ + hpzarraytext


def hpxpypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X'Y' bar."
  if not top.lhxpypbarz: return
  kw.update(kwdict)
  kw['titlet']="X'Y' bar"
  hpzarray(top.hxpypbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxpypbarz.__doc__ = hpxpypbarz.__doc__ + hpzarraytext


def hpxvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XVz bar."
  if not top.lhxvzbarz: return
  kw.update(kwdict)
  kw['titlet']="XVz bar"
  hpzarray(top.hxvzbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpxvzbarz.__doc__ = hpxvzbarz.__doc__ + hpzarraytext


def hpyvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YVz bar."
  if not top.lhyvzbarz: return
  kw.update(kwdict)
  kw['titlet']="YVz bar"
  hpzarray(top.hyvzbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpyvzbarz.__doc__ = hpyvzbarz.__doc__ + hpzarraytext


def hpvxvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "VxVz bar."
  if not top.lhvxvzbarz: return
  kw.update(kwdict)
  kw['titlet']="VxVz bar"
  hpzarray(top.hvxvzbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvxvzbarz.__doc__ = hpvxvzbarz.__doc__ + hpzarraytext


def hpvyvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "VyVz bar."
  if not top.lhvyvzbarz: return
  kw.update(kwdict)
  kw['titlet']="VyVz bar"
  hpzarray(top.hvyvzbarz,contour,overlay,iz,kw)
if sys.version[:5] != "1.5.1":
  hpvyvzbarz.__doc__ = hpvyvzbarz.__doc__ + hpzarraytext


def hptotalke(kwdict={},**kw):
  "Total Kinetic Energy."
  kw.update(kwdict)
  kw['titlet']="Total Kinetic Energy"
  kw['titlel']="(J)"
  hpbasic(top.hekzbeam+top.hekperp,kw)


def hptotale(kwdict={},**kw):
  "Total Energy."
  kw.update(kwdict)
  kw['titlet']="Total Energy"
  kw['titlel']="(J)"
  hpbasic(top.hefld+top.hekzbeam+top.hekperp,kw)


def hpthermale(iw=0,kwdict={},**kw):
  "Z Thermal Energy."
  kw.update(kwdict)
  kw['titlet']="Z Thermal Energy"
  kw['titlel']="(J)"
  hpbasicwin(0.5*sum(top.sm*top.sw*top.sp_fract)*top.hpnum*top.hvzrms**2,iw,kw)
if sys.version[:5] != "1.5.1":
  hpthermale.__doc__ = hpthermale.__doc__ + hpbasicwintext


def hpeps6d(iw=0,kwdict={},**kw):
  "6-D Emittance."
  kw.update(kwdict)
  kw['titlet']="6-D Emittance"
  kw['titlel']="((pi-m-rad)**3)"
  hpbasicwin(top.hepsx*top.hepsy*top.hepsz,iw,kw)
if sys.version[:5] != "1.5.1":
  hpeps6d.__doc__ = hpeps6d.__doc__ + hpbasicwintext


def hpepst(iw=0,kwdict={},**kw):
  "Transverse Emittance."
  kw.update(kwdict)
  kw['titlet']="Transverse Emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(sqrt(top.hepsx*top.hepsy),iw,kw)
if sys.version[:5] != "1.5.1":
  hpepst.__doc__ = hpepst.__doc__ + hpbasicwintext


def hpepsnt(iw=0,kwdict={},**kw):
  "Normalized Transverse Emittance."
  kw.update(kwdict)
  kw['titlet']="Normalized Transverse Emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(sqrt(top.hepsnx*top.hepsny),iw,kw)
if sys.version[:5] != "1.5.1":
  hpepsnt.__doc__ = hpepsnt.__doc__ + hpbasicwintext


def hpxedge(iw=0,kwdict={},**kw):
  "X Beam Edge."
  kw.update(kwdict)
  kw['titlet']="X Beam Edge"
  kw['titlel']="(m)"
  hpbasicwin(2.*top.hxrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxedge.__doc__ = hpxedge.__doc__ + hpbasicwintext


def hpxpedge(iw=0,kwdict={},**kw):
  "X' at Beam Edge."
  kw.update(kwdict)
  kw['titlet']="X' at Beam Edge"
  kw['titlel']="(m)"
  xpedge = (top.hxxpbar-top.hxbar*top.hxpbar)/ \
           where(greater(top.hxrms,0.),top.hxrms,1.)
  hpbasicwin(xpedge,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxpedge.__doc__ = hpxpedge.__doc__ + hpbasicwintext


def hpyedge(iw=0,kwdict={},**kw):
  "Y Beam Edge."
  kw.update(kwdict)
  kw['titlet']="Y Beam Edge"
  kw['titlel']="(m)"
  hpbasicwin(2.*top.hyrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyedge.__doc__ = hpyedge.__doc__ + hpbasicwintext


def hpypedge(iw=0,kwdict={},**kw):
  "Y' at Beam Edge."
  kw.update(kwdict)
  kw['titlet']="Y' at Beam Edge"
  kw['titlel']="(m)"
  ypedge = (top.hyypbar-top.hybar*top.hypbar)/ \
           where(greater(top.hyrms,0.),top.hyrms,1.)
  hpbasicwin(ypedge,iw,kw)
if sys.version[:5] != "1.5.1":
  hpypedge.__doc__ = hpypedge.__doc__ + hpbasicwintext


def hpxedges(iw=0,kwdict={},**kw):
  "X Beam Edges plus centroid."
  kw.update(kwdict)
  kw['titlet']="X Beam Edges plus centroid"
  kw['titlel']="(m)"
  hpbasicwin(+2.*top.hxrms+top.hxbar,iw,kw)
  hpbasicwin(-2.*top.hxrms+top.hxbar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpxedge.__doc__ = hpxedge.__doc__ + hpbasicwintext


def hpyedges(iw=0,kwdict={},**kw):
  "Y Beam Edges plus centroid."
  kw.update(kwdict)
  kw['titlet']="Y Beam Edges plus centroid"
  kw['titlel']="(m)"
  hpbasicwin(+2.*top.hyrms+top.hybar,iw,kw)
  hpbasicwin(-2.*top.hyrms+top.hybar,iw,kw)
if sys.version[:5] != "1.5.1":
  hpyedge.__doc__ = hpyedge.__doc__ + hpbasicwintext


def hpenvx(iw=0,kwdict={},**kw):
  "X Beam Edge."
  kw.update(kwdict)
  kw['titlet']="X Beam Edge"
  kw['titlel']="(m)"
  hpbasicwin(2.*top.hxrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpenvx.__doc__ = hpenvx.__doc__ + hpbasicwintext


def hpenvy(iw=0,kwdict={},**kw):
  "Y Beam Edge."
  kw.update(kwdict)
  kw['titlet']="Y Beam Edge"
  kw['titlel']="(m)"
  hpbasicwin(2.*top.hyrms,iw,kw)
if sys.version[:5] != "1.5.1":
  hpenvy.__doc__ = hpenvy.__doc__ + hpbasicwintext

# --- Plots of current
def hpcurr(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Current."
  if top.lhcurrz:
    kw.update(kwdict)
    kw['titlet']="Current"
    hpzarray(top.hcurrz,contour,overlay,iz,kw)
  elif top.lhlinechg and top.lhvzofz:
    kw.update(kwdict)
    kw['titlet']="Current"
    hpzarray(top.hlinechg*top.hvzofz,contour,overlay,iz,kw)
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

