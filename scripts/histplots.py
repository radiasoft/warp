from warp import *
from mplot import *
histplots_version = "$Id: histplots.py,v 1.2 2001/01/25 22:10:34 dave Exp $"

###########################################################################
def hpdoc():
  """
What follows is a list of all possible arguments to any of the history
plotting routines, along with their default values.
  - absc Data for the abscissa. Defaults to either thist or hzbeam
  - xmin, xmax, ymin, ymax Plot limits, defaults to extrema of data
  - titlet, titleb, titlel, titler Plot titles, defaults to blank except
    titleb which is set to name of bottom axis when absc not input
  - xscale = 1.0 Scale factor for abscissa
  - xoffset = 0.0 Offset for abscissa
  - yscale = 1.0 Scale factor for oordinate
  - yoffset = 0.0 Offset for oordinate
  - lnormalized = 0 When true, scales data by initial value
  - i1 = 0 Time index at which to start the plots
  - i2 = top.jhist Time index at which to end the plots
  - nnhp = 1 Step size of time data plotted
  - lhzbeam = 0 When true, plots data versus hzbeam instead of thist
  - logplot = 0 When true, make a log plot
  - color = 'fg' Color of the curve plotted
  - marks = 0 Marks to place on curve
  - marker = None Marker to place on curve
  - msize = 1.0 Marker size
  - width = 1.0 Line width
  - titles = 1 When false, no titles are printed
  - plsysval = 1 Plot system to make plot in (quadrant plots for example)
    see work.gs for the plot system values.
  """
  print hpdoc.__doc__

###########################################################################
def hpbasic(oord,kwdict={},**kw):
  """
This is the basic history plot with all of its possible arguments. The
parsing of the arguments is done in one place here and each of the history
plot functions just passes most of the arguments into this function. The
only required argument of course is the data to be plotted.
  """
  keywords = ['absc','xmin','xmax','ymin','ymax','titlet','titleb','titlel',
              'titler','xscale','xoffset','yscale','yoffset','lnormalized',
              'i1','i2','nnhp','lhzbeam','logplot','color','marks','marker',
              'msize','titles','plsysval','width']

  # --- Add together dictionaries
  kw.update(kwdict)

  # --- Set the default values of input
  absc = None
  xmin = 'e'
  xmax = 'e'
  ymin = 'e'
  ymax = 'e'
  titlet = ''
  titleb = ''
  titlel = ''
  titler = ''
  xscale = 1.0
  xoffset = 0.0
  yscale = array([1.0])
  yoffset = 0.0
  lnormalized = 0
  i1 = 0
  i2 = top.jhist
  nnhp = 1
  lhzbeam = 0
  logplot = 0
  color = 'fg'
  marks = 0
  marker = None
  msize = 1.0
  titles = 1
  plsysval = 1
  width = 1.
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
  if lnormalized:
    yscale = 1./oord[...,i1]
    if not titlel: titlel = titlet + " over initial value"
  if type(yscale) != type(array([])): yscale = array([yscale])
  if not absc:
    if (lhzbeam):
      absc = top.hzbeam[i1:i2+1:nnhp]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "Z (m)"
        else:
          titleb = "Z"
    else:
      absc = top.thist[i1:i2+1:nnhp]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "time (s)"
        else:
          titleb = "time"
  if logplot:
    logxy(1,0)
    oord = log10(maximum(10e-12,
                         oord[...,i1:i2+1:nnhp]*yscale[...,NewAxis]+yoffset))
    if not titler:
      titler = "logarithmic scale"
  else:
    oord = oord[...,i1:i2+1:nnhp]*yscale[...,NewAxis]+yoffset

  # --- Now actually make the plot after all of that ado.   
  pla(transpose(oord),absc,color=color,msize=msize,marks=marks,marker=marker,
      width=width)
  if titles: ptitles(titlet,titleb,titlel,titler,plsysval)
  limits(xmin,xmax,ymin,ymax)
  if logplot: logxy(0,0)


#############################################################################
# Basic history plot for data in windows
def hpbasicwin(oord,iw=0,kwdict={},**kw):
  kw.update(kwdict)
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
  kw.update(kwdict)
  keywords = ['absc','xmin','xmax','ymin','ymax','titlet','titleb','titlel',
              'titler','xscale','xoffset','yscale','yoffset','i1','i2','inct',
              'j1','j2','incz','lhzbeam','logplot','color','marks','marker',
              'msize','titles','levs','filled']
  # --- Set the default values of input
  absc = None
  xmin = 'e'
  xmax = 'e'
  ymin = 'e'
  ymax = 'e'
  titlet = ''
  titleb = ''
  titlel = ''
  titler = ''
  xscale = 1.0
  xoffset = 0.0
  yscale = 1.0
  yoffset = 0.0
  i1 = 0
  i2 = top.jhist
  inct = max(i2/20,1)
  j1 = 0
  j2 = top.nzzarr
  incz = max(j2/32,1)
  lhzbeam = 0
  logplot = 0
  color = 'fg'
  marks = 0
  marker = None
  msize = 1.0
  titles = 1
  levs = None
  filled = 0
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
      absc = top.hzbeam[i1:i2+1:inct]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "Z (m)"
        else:
          titleb = "Z"
    else:
      absc = top.thist[i1:i2+1:inct]*xscale + xoffset
      if not titleb:
        if (xscale == 1.):
          titleb = "time (s)"
        else:
          titleb = "time"
  if logplot:
    oord = log10(maximum(10e-12,transpose(oord[j1:j2+1:incz,i1:i2+1:inct])*
                                yscale+yoffset))
    if not titler:
      titler = "logarithmic scale"
  else:
    oord = transpose(oord[j1:j2+1:incz,i1:i2+1:inct]*yscale+yoffset)

  # --- Now actually make the plot after all of that ado.   
  if filled:
    plotfc(oord,absc,oordmesh[j1:j2+1:incz],color=color,contours=levs)
  else:
    plotc(oord,absc,oordmesh[j1:j2+1:incz],color=color,levs=levs)
  if titles: ptitles(titlet,titleb,titlel,titler)
  if xmin == 'e': xmin = min(oordmesh[j1:j2+1:incz])
  if xmax == 'e': xmax = max(oordmesh[j1:j2+1:incz])
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
    hpbasic(hzarray[iz,:],kw)
    return
  elif contour:
    hpbasiccont(hzarray,top.zmntmesh,kw)
  elif overlay:
    kw['nnhp'] = max(top.jhist/10,1)
    kw['absc'] = w3d.zmesh
    hpbasic(hzarray,kw)
  else:
    mountainplot1(kw['titlet'],hzarray,kw)

###########################################################################
###########################################################################
###########################################################################

def hpzbeam(kwdict={},**kw):
  "Beam frame location. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Beam frame location"
  kw['titlel']="(m)"
  hpbasic(top.hzbeam,kw)


def hpvbeam(kwdict={},**kw):
  "Beam frame velocity. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Beam frame velocity"
  kw['titlel']="(m/s)"
  hpbasic(top.hvbeam,kw)


def hpbmlen(kwdict={},**kw):
  "RMS beam length. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="RMS beam length"
  kw['titlel']="(m)"
  hpbasic(top.hbmlen,kw)


def hpefld(kwdict={},**kw):
  "Field energy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Field energy"
  kw['titlel']="(J)"
  hpbasic(top.hefld,kw)


def hpekzmbe(kwdict={},**kw):
  "Total Z Kinetic energy minus beam energy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Total Z Kinetic energy minus beam energy"
  kw['titlel']="(J)"
  hpbasic(top.hekzmbe,kw)


def hpekzbeam(kwdict={},**kw):
  "Z Kinetic energy in the beam frame. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Z Kinetic energy in the beam frame"
  kw['titlel']="(J)"
  hpbasic(top.hekzbeam,kw)


def hpekperp(kwdict={},**kw):
  "Perp Kinetic energy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Perp Kinetic energy"
  kw['titlel']="(J)"
  hpbasic(top.hekperp,kw)


def hpepsx(iw=0,kwdict={},**kw):
  "X emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="X emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsx,iw,kw)


def hpepsy(iw=0,kwdict={},**kw):
  "Y emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Y emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsy,iw,kw)


def hpepsz(iw=0,kwdict={},**kw):
  "Z emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Z emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsz,iw,kw)


def hpepsnx(iw=0,kwdict={},**kw):
  "X normalized emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="X normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsnx,iw,kw)


def hpepsny(iw=0,kwdict={},**kw):
  "Y normalized emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Y normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsny,iw,kw)


def hpepsnz(iw=0,kwdict={},**kw):
  "Z normalized emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Z normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsnz,iw,kw)


def hpepsg(iw=0,kwdict={},**kw):
  "Generalized emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsg,iw,kw)


def hpepsh(iw=0,kwdict={},**kw):
  "Generalized emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(top.hepsh,iw,kw)


def hpepsng(iw=0,kwdict={},**kw):
  "Generalized normalized emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Generalized normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsng,iw,kw)


def hpepsnh(iw=0,kwdict={},**kw):
  "Generalized normalized emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Generalized normalized emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(top.hepsnh,iw,kw)


def hppnum(iw=0,kwdict={},**kw):
  "Number of particles. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Number of particles"
  kw['titlel']=""
  hpbasicwin(top.hpnum,iw,kw)


def hprhomid(iw=0,kwdict={},**kw):
  "Charge density on axis. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Charge density on axis"
  kw['titlel']="(C/m**3)"
  hpbasicwin(top.hrhomid,iw,kw)


def hprhomax(iw=0,kwdict={},**kw):
  "Charge density max. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Charge density max"
  kw['titlel']="(C/m**3)"
  hpbasicwin(top.hrhomax,iw,kw)


def hpxbar(iw=0,kwdict={},**kw):
  "True mean x. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True mean x"
  kw['titlel']="(m)"
  hpbasicwin(top.hxbar,iw,kw)


def hpybar(iw=0,kwdict={},**kw):
  "True mean y. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True mean y"
  kw['titlel']="(m)"
  hpbasicwin(top.hybar,iw,kw)


def hpxybar(iw=0,kwdict={},**kw):
  "True mean xy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True mean xy"
  kw['titlel']="(m**2)"
  hpbasicwin(top.hxybar,iw,kw)


def hpxrms(iw=0,kwdict={},**kw):
  "True RMS x. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True RMS x"
  kw['titlel']="(m)"
  hpbasicwin(top.hxrms,iw,kw)


def hpyrms(iw=0,kwdict={},**kw):
  "True RMS y. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True RMS y"
  kw['titlel']="(m)"
  hpbasicwin(top.hyrms,iw,kw)


def hpxprms(iw=0,kwdict={},**kw):
  "True RMS x'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True RMS x'"
  kw['titlel']="(rad)"
  hpbasicwin(top.hxprms,iw,kw)


def hpyprms(iw=0,kwdict={},**kw):
  "True RMS y'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True RMS y'"
  kw['titlel']="(rad)"
  hpbasicwin(top.hyprms,iw,kw)


def hpxsqbar(iw=0,kwdict={},**kw):
  "Mean x squared. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean x squared"
  kw['titlel']="(m**2)"
  hpbasicwin(top.hxsqbar,iw,kw)


def hpysqbar(iw=0,kwdict={},**kw):
  "Mean y squared. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean y squared"
  kw['titlel']="(m**2)"
  hpbasicwin(top.hysqbar,iw,kw)


def hpvxbar(iw=0,kwdict={},**kw):
  "Mean vx. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean vx"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvxbar,iw,kw)


def hpvybar(iw=0,kwdict={},**kw):
  "Mean vy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean vy"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvybar,iw,kw)


def hpvzbar(iw=0,kwdict={},**kw):
  "Mean vz. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean vz"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvzbar,iw,kw)


def hpxpbar(iw=0,kwdict={},**kw):
  "Mean x'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean x'"
  kw['titlel']="(rad)"
  hpbasicwin(top.hxpbar,iw,kw)


def hpypbar(iw=0,kwdict={},**kw):
  "Mean y'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean y'"
  kw['titlel']="(rad)"
  hpbasicwin(top.hypbar,iw,kw)


def hpvxrms(iw=0,kwdict={},**kw):
  "True RMS vx. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True RMS vx"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvxrms,iw,kw)


def hpvyrms(iw=0,kwdict={},**kw):
  "True RMS vy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True RMS vy"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvyrms,iw,kw)


def hpvzrms(iw=0,kwdict={},**kw):
  "True RMS vz. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="True RMS vz"
  kw['titlel']="(m/s)"
  hpbasicwin(top.hvzrms,iw,kw)


def hpxpsqbar(iw=0,kwdict={},**kw):
  "Mean x' squared. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean x' squared"
  kw['titlel']="(rad**2)"
  hpbasicwin(top.hxpsqbar,iw,kw)


def hpypsqbar(iw=0,kwdict={},**kw):
  "Mean y' squared. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean y' squared"
  kw['titlel']="(rad**2)"
  hpbasicwin(top.hypsqbar,iw,kw)


def hpxxpbar(iw=0,kwdict={},**kw):
  "Mean x*x'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean x*x'"
  kw['titlel']="(m-rad)"
  hpbasicwin(top.hxxpbar,iw,kw)


def hpyypbar(iw=0,kwdict={},**kw):
  "Mean y*y'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean y*y'"
  kw['titlel']="(m-rad)"
  hpbasicwin(top.hyypbar,iw,kw)


def hpxypbar(iw=0,kwdict={},**kw):
  "Mean x*y'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean x*y'"
  kw['titlel']="(m-rad)"
  hpbasicwin(top.hxypbar,iw,kw)


def hpyxpbar(iw=0,kwdict={},**kw):
  "Mean y*x'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean y*x'"
  kw['titlel']="(m**2)"
  hpbasicwin(top.hyxpbar,iw,kw)


def hpxpypbar(iw=0,kwdict={},**kw):
  "Mean x'*y'. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean x'*y'"
  kw['titlel']="(rad**2)"
  hpbasicwin(top.hxpypbar,iw,kw)


def hpxvzbar(iw=0,kwdict={},**kw):
  "Mean x*vz. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean x*vz"
  kw['titlel']="(m-m/s)"
  hpbasicwin(top.hxvzbar,iw,kw)


def hpyvzbar(iw=0,kwdict={},**kw):
  "Mean y*vz. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean y*vz"
  kw['titlel']="(m-m/s)"
  hpbasicwin(top.hyvzbar,iw,kw)


def hpvxvzbar(iw=0,kwdict={},**kw):
  "Mean vx*vz. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean vx*vz"
  kw['titlel']="((m/s)**2)"
  hpbasicwin(top.hvxvzbar,iw,kw)


def hpvyvzbar(iw=0,kwdict={},**kw):
  "Mean vy*vz. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Mean vy*vz"
  kw['titlel']="((m/s)**2)"
  hpbasicwin(top.hvyvzbar,iw,kw)


def hplinechg(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Line charge density. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhlinechg: return
  kw.update(kwdict)
  kw['titlet']="Line charge density"
  hpzarray(top.hlinechg,contour,overlay,iz,kw)


def hpvzofz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz versus space and time. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvzofz: return
  kw.update(kwdict)
  kw['titlet']="Vz versus space and time"
  hpzarray(top.hvzofz,contour,overlay,iz,kw)


def hpepsxz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X emittance. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhepsxz: return
  kw.update(kwdict)
  kw['titlet']="X emittance"
  hpzarray(top.hepsxz,contour,overlay,iz,kw)


def hpepsyz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y emittance. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhepsyz: return
  kw.update(kwdict)
  kw['titlet']="Y emittance"
  hpzarray(top.hepsyz,contour,overlay,iz,kw)


def hpepsnxz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X normalized emittance. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhepsnxz: return
  kw.update(kwdict)
  kw['titlet']="X normalized emittance"
  hpzarray(top.hepsnxz,contour,overlay,iz,kw)


def hpepsnyz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y normalized emittance. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhepsnyz: return
  kw.update(kwdict)
  kw['titlet']="Y normalized emittance"
  hpzarray(top.hepsnyz,contour,overlay,iz,kw)


def hpepsgz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized emittance. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhepsgz: return
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  hpzarray(top.hepsgz,contour,overlay,iz,kw)


def hpepshz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized emittance. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhepshz: return
  kw.update(kwdict)
  kw['titlet']="Generalized emittance"
  hpzarray(top.hepshz,contour,overlay,iz,kw)


def hpepsngz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized nrmlzd emittance. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhepsngz: return
  kw.update(kwdict)
  kw['titlet']="Generalized nrmlzd emittance"
  hpzarray(top.hepsngz,contour,overlay,iz,kw)


def hpepsnhz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Generalized nrmlzd emittance. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhepsnhz: return
  kw.update(kwdict)
  kw['titlet']="Generalized nrmlzd emittance"
  hpzarray(top.hepsnhz,contour,overlay,iz,kw)


def hpxbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxbarz: return
  kw.update(kwdict)
  kw['titlet']="X bar"
  hpzarray(top.hxbarz,contour,overlay,iz,kw)


def hpybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhybarz: return
  kw.update(kwdict)
  kw['titlet']="Y bar"
  hpzarray(top.hybarz,contour,overlay,iz,kw)


def hpxybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XY bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxybarz: return
  kw.update(kwdict)
  kw['titlet']="XY bar"
  hpzarray(top.hxybarz,contour,overlay,iz,kw)


def hpxrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X rms. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxrmsz: return
  kw.update(kwdict)
  kw['titlet']="X rms"
  hpzarray(top.hxrmsz,contour,overlay,iz,kw)


def hpyrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y rms. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhyrmsz: return
  kw.update(kwdict)
  kw['titlet']="Y rms"
  hpzarray(top.hyrmsz,contour,overlay,iz,kw)


def hpxprmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X' rms. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxprmsz: return
  kw.update(kwdict)
  kw['titlet']="X' rms"
  hpzarray(top.hxprmsz,contour,overlay,iz,kw)


def hpyprmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y' rms. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhyprmsz: return
  kw.update(kwdict)
  kw['titlet']="Y' rms"
  hpzarray(top.hyprmsz,contour,overlay,iz,kw)


def hpxsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X**2 bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxsqbarz: return
  kw.update(kwdict)
  kw['titlet']="X**2 bar"
  hpzarray(top.hxsqbarz,contour,overlay,iz,kw)


def hpysqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y**2 bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhysqbarz: return
  kw.update(kwdict)
  kw['titlet']="Y**2 bar"
  hpzarray(top.hysqbarz,contour,overlay,iz,kw)


def hpvxbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vx bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvxbarz: return
  kw.update(kwdict)
  kw['titlet']="Vx bar"
  hpzarray(top.hvxbarz,contour,overlay,iz,kw)


def hpvybarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vy bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvybarz: return
  kw.update(kwdict)
  kw['titlet']="Vy bar"
  hpzarray(top.hvybarz,contour,overlay,iz,kw)


def hpvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvzbarz: return
  kw.update(kwdict)
  kw['titlet']="Vz bar"
  hpzarray(top.hvzbarz,contour,overlay,iz,kw)


def hpxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X' bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxpbarz: return
  kw.update(kwdict)
  kw['titlet']="X' bar"
  hpzarray(top.hxpbarz,contour,overlay,iz,kw)


def hpypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y' bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhypbarz: return
  kw.update(kwdict)
  kw['titlet']="Y' bar"
  hpzarray(top.hypbarz,contour,overlay,iz,kw)


def hpvxrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vx rms. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvxrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vx rms"
  hpzarray(top.hvxrmsz,contour,overlay,iz,kw)


def hpvyrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vy rms. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvyrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vy rms"
  hpzarray(top.hvyrmsz,contour,overlay,iz,kw)


def hpvzrmsz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Vz rms. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvzrmsz: return
  kw.update(kwdict)
  kw['titlet']="Vz rms"
  hpzarray(top.hvzrmsz,contour,overlay,iz,kw)


def hpxpsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X'**2 bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxpsqbarz: return
  kw.update(kwdict)
  kw['titlet']="X'**2 bar"
  hpzarray(top.hxpsqbarz,contour,overlay,iz,kw)


def hpypsqbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Y'**2 bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhypsqbarz: return
  kw.update(kwdict)
  kw['titlet']="Y'**2 bar"
  hpzarray(top.hypsqbarz,contour,overlay,iz,kw)


def hpxxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XX' bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxxpbarz: return
  kw.update(kwdict)
  kw['titlet']="XX' bar"
  hpzarray(top.hxxpbarz,contour,overlay,iz,kw)


def hpyypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YY' bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhyypbarz: return
  kw.update(kwdict)
  kw['titlet']="YY' bar"
  hpzarray(top.hyypbarz,contour,overlay,iz,kw)


def hpxypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XY' bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxypbarz: return
  kw.update(kwdict)
  kw['titlet']="XY' bar"
  hpzarray(top.hxypbarz,contour,overlay,iz,kw)


def hpyxpbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YX' bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhyxpbarz: return
  kw.update(kwdict)
  kw['titlet']="YX' bar"
  hpzarray(top.hyxpbarz,contour,overlay,iz,kw)


def hpxpypbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "X'Y' bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxpypbarz: return
  kw.update(kwdict)
  kw['titlet']="X'Y' bar"
  hpzarray(top.hxpypbarz,contour,overlay,iz,kw)


def hpxvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "XVz bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhxvzbarz: return
  kw.update(kwdict)
  kw['titlet']="XVz bar"
  hpzarray(top.hxvzbarz,contour,overlay,iz,kw)


def hpyvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "YVz bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhyvzbarz: return
  kw.update(kwdict)
  kw['titlet']="YVz bar"
  hpzarray(top.hyvzbarz,contour,overlay,iz,kw)


def hpvxvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "VxVz bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvxvzbarz: return
  kw.update(kwdict)
  kw['titlet']="VxVz bar"
  hpzarray(top.hvxvzbarz,contour,overlay,iz,kw)


def hpvyvzbarz(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "VyVz bar. For more help, run hpdoc() or doc(hpzarray)."
  if not top.lhvyvzbarz: return
  kw.update(kwdict)
  kw['titlet']="VyVz bar"
  hpzarray(top.hvyvzbarz,contour,overlay,iz,kw)


def hptotalke(kwdict={},**kw):
  "Total Kinetic Energy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Total Kinetic Energy"
  kw['titlel']="(J)"
  hpbasic(top.hekzbeam+top.hekperp,kw)


def hptotale(kwdict={},**kw):
  "Total Energy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Total Energy"
  kw['titlel']="(J)"
  hpbasic(top.hefld+top.hekzbeam+top.hekperp,kw)


def hpthermale(iw=0,kwdict={},**kw):
  "Z Thermal Energy. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Z Thermal Energy"
  kw['titlel']="(J)"
  hpbasicwin(0.5*sum(top.sm*top.sw*top.sp_fract)*top.hpnum*top.hvzrms**2,iw,kw)


def hpeps6d(iw=0,kwdict={},**kw):
  "6-D Emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="6-D Emittance"
  kw['titlel']="((pi-m-rad)**3)"
  hpbasicwin(top.hepsx*top.hepsy*top.hepsz,iw,kw)


def hpepst(iw=0,kwdict={},**kw):
  "Transverse Emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Transverse Emittance"
  kw['titlel']="(pi-m-rad)"
  hpbasicwin(sqrt(top.hepsx*top.hepsy),iw,kw)


def hpepsnt(iw=0,kwdict={},**kw):
  "Normalized Transverse Emittance. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Normalized Transverse Emittance"
  kw['titlel']="(pi-mm-mrad)"
  hpbasicwin(sqrt(top.hepsnx*top.hepsny),iw,kw)


def hpxedge(iw=0,kwdict={},**kw):
  "X Beam Edge. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="X Beam Edge"
  kw['titlel']="(m)"
  hpbasicwin(2.*top.hxrms,iw,kw)


def hpyedge(iw=0,kwdict={},**kw):
  "Y Beam Edge. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Y Beam Edge"
  kw['titlel']="(m)"
  hpbasicwin(2.*top.hyrms,iw,kw)


def hpenvx(iw=0,kwdict={},**kw):
  "X Beam Edge. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="X Beam Edge"
  kw['titlel']="(m)"
  hpbasicwin(2.*top.hxrms,iw,kw)


def hpenvy(iw=0,kwdict={},**kw):
  "Y Beam Edge. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="Y Beam Edge"
  kw['titlel']="(m)"
  hpbasicwin(2.*top.hyrms,iw,kw)

# --- plot mean z velocity in beam frame vs time
def hpvzbar(iw=0,beamframe=1,kwdict={},**kw):
  kw.update(kwdict)
  if beamframe:
    t1="Mean Z Velocity (beam frame)"
    kw['titler']="vbeam = %6.3e"%top.vbeam
    hpbasicwin(top.hvzbar-top.vbeam,iw,t1,t1,t1,kw)
  else:
    t1="Mean Z Velocity (lab frame)"
    hpbasicwin(top.hvzbar,iw,t1,t1,t1,kw)

# --- Plots of current
def hpcurr(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "Current. For a complete list of arguments, run hpdoc()."
  if not top.lhlinechg and not top.lhvzofz: return
  kw.update(kwdict)
  kw['titlet']="Current"
  hpzarray(top.hlinechg*top.hvzofz,contour,overlay,iz,kw)

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
hpyedge(): Y Beam Edge (twice rms)
hpenvx = hpxedge
hpenvy = hpyedge
hpzbeam(): Beam frame location
hpvbeam(): Beam frame velocity
hpbmlen(): RMS beam length
hpefld(): Field energy
hpekzmbe(): Total Z Kinetic energy minus beam energy
hpekzbeam(): Z Kinetic energy in the beam frame
hpekperp(): Perp Kinetic energy
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
hpptotalke(): Total Kinetic Energy
hpptotale(): Total Energy
hppthermale(): Z Thermal Energy
hppeps6d(): 6-D Emittance
hppepst(): Transverse Emittance
hppepsnt(): Normalized Transverse Emittance
hppxedge(): X Beam Edge
hppyedge(): Y Beam Edge
hppenvx(): X Beam Edge
hppenvy(): Y Beam Edge
  """
  print histplotsdoc.__doc__

