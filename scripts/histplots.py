from warp import *
from mplot import *
histplots_version = "$Id: histplots.py,v 1.1 2000/10/16 18:34:19 dave Exp $"

###########################################################################
def hpdoc():
  """
What follows is a list of all possible arguments to any of the history
plotting routines, along with their default values.
  - oord Data to be plotted (Only required argument)
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

#############################################################################
# Now define the many plot functions.

# --- List of all of the history variables and their docs
# --- It would be better to generate these list directly from
# --- top.varlist("Hist"). The list is done here to clean up the documentation.
hist1list = [["hzbeam","Beam frame location","(m)"],
             ["hvbeam","Beam frame velocity","(m/s)"],
             ["hbmlen","RMS beam length","(m)"],
             ["hefld","Field energy","(J)"],
             ["hekzmbe","Total Z Kinetic energy minus beam energy","(J)"],
             ["hekzbeam","Z Kinetic energy in the beam frame","(J)"],
             ["hekperp","Perp Kinetic energy","(J)"]]
hist2list = [["hepsx","X emittance","(pi-m-rad)"],
             ["hepsy","Y emittance","(pi-m-rad)"],
             ["hepsz","Z emittance","(pi-m-rad)"],
             ["hepsnx","X normalized emittance","(pi-mm-mrad)"],
             ["hepsny","Y normalized emittance","(pi-mm-mrad)"],
             ["hepsnz","Z normalized emittance","(pi-mm-mrad)"],
             ["hepsg","Generalized emittance","(pi-m-rad)"],
             ["hepsh","Generalized emittance","(pi-m-rad)"],
             ["hepsng","Generalized normalized emittance","(pi-mm-mrad)"],
             ["hepsnh","Generalized normalized emittance","(pi-mm-mrad)"],
             ["hpnum","Number of particles",""],
             ["hrhomid","Charge density on axis","(C/m**3)"],
             ["hrhomax","Charge density max","(C/m**3)"],
             ["hxbar","True mean x","(m)"],
             ["hybar","True mean y","(m)"],
             ["hxybar","True mean xy","(m**2)"],
             ["hxrms","True RMS x","(m)"],
             ["hyrms","True RMS y","(m)"],
             ["hxprms","True RMS x'","(rad)"],
             ["hyprms","True RMS y'","(rad)"],
             ["hxsqbar","Mean x squared","(m**2)"],
             ["hysqbar","Mean y squared","(m**2)"],
             ["hvxbar","Mean vx","(m/s)"],
             ["hvybar","Mean vy","(m/s)"],
             ["hvzbar","Mean vz","(m/s)"],
             ["hxpbar","Mean x'","(rad)"],
             ["hypbar","Mean y'","(rad)"],
             ["hvxrms","True RMS vx","(m/s)"],
             ["hvyrms","True RMS vy","(m/s)"],
             ["hvzrms","True RMS vz","(m/s)"],
             ["hxpsqbar","Mean x' squared","(rad**2)"],
             ["hypsqbar","Mean y' squared","(rad**2)"],
             ["hxxpbar","Mean x*x'","(m-rad)"],
             ["hyypbar","Mean y*y'","(m-rad)"],
             ["hxypbar","Mean x*y'","(m-rad)"],
             ["hyxpbar","Mean y*x'","(m**2)"],
             ["hxpypbar","Mean x'*y'","(rad**2)"],
             ["hxvzbar","Mean x*vz","(m-m/s)"],
             ["hyvzbar","Mean y*vz","(m-m/s)"],
             ["hvxvzbar","Mean vx*vz","((m/s)**2)"],
             ["hvyvzbar","Mean vy*vz","((m/s)**2)"]]
hist3list = [["hlinechg","Line charge density"],
             ["hvzofz","Vz versus space and time"],
             ["hepsxz","X emittance"],
             ["hepsyz","Y emittance"],
             ["hepsnxz","X normalized emittance"],
             ["hepsnyz","Y normalized emittance"],
             ["hepsgz","Generalized emittance"],
             ["hepshz","Generalized emittance"],
             ["hepsngz","Generalized nrmlzd emittance"],
             ["hepsnhz","Generalized nrmlzd emittance"],
             ["hxbarz","X bar"],
             ["hybarz","Y bar"],
             ["hxybarz","XY bar"],
             ["hxrmsz","X rms"],
             ["hyrmsz","Y rms"],
             ["hxprmsz","X' rms"],
             ["hyprmsz","Y' rms"],
             ["hxsqbarz","X**2 bar"],
             ["hysqbarz","Y**2 bar"],
             ["hvxbarz","Vx bar"],
             ["hvybarz","Vy bar"],
             ["hvzbarz","Vz bar"],
             ["hxpbarz","X' bar"],
             ["hypbarz","Y' bar"],
             ["hvxrmsz","Vx rms"],
             ["hvyrmsz","Vy rms"],
             ["hvzrmsz","Vz rms"],
             ["hxpsqbarz","X'**2 bar"],
             ["hypsqbarz","Y'**2 bar"],
             ["hxxpbarz","XX' bar"],
             ["hyypbarz","YY' bar"],
             ["hxypbarz","XY' bar"],
             ["hyxpbarz","YX' bar"],
             ["hxpypbarz","X'Y' bar"],
             ["hxvzbarz","XVz bar"],
             ["hyvzbarz","YVz bar"],
             ["hvxvzbarz","VxVz bar"],
             ["hvyvzbarz","VyVz bar"]]

hist4list=[["hptotalke","Total Kinetic Energy","top.hekzbeam+top.hekperp",
            "(J)"],
           ["hptotale","Total Energy","top.hefld+top.hekzbeam+top.hekperp",
                       "(J)"]]

hist5list=[["hpthermale","Z Thermal Energy",
            "0.5*sum(top.sm*top.sw*top.sp_fract)*top.hpnum*top.hvzrms**2",
            "(J)"],
           ["hpeps6d","6-D Emittance","top.hepsx*top.hepsy*top.hepsz",
            "((pi-m-rad)**3)"],
           ["hpepst","Transverse Emittance","sqrt(top.hepsx*top.hepsy)",
            "(pi-m-rad)"],
           ["hpepsnt","Normalized Transverse Emittance",
            "sqrt(top.hepsnx*top.hepsny)","(pi-mm-mrad)"],
           ["hpxedge","X Beam Edge","2.*top.hxrms","(m)"],
           ["hpyedge","Y Beam Edge","2.*top.hyrms","(m)"],
           ["hpenvx","X Beam Edge","2.*top.hxrms","(m)"],
           ["hpenvy","Y Beam Edge","2.*top.hyrms","(m)"]]


###########################################################################
# --- Template for global history plots
dd1 = """
def %s(kwdict={},**kw):
  "%s. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="%s"
  kw['titlel']="%s"
  hpbasic(%s,kw)
"""

# --- Template for window history plots
dd2 = """
def %s(iw=0,kwdict={},**kw):
  "%s. For a complete list of arguments, run hpdoc()."
  kw.update(kwdict)
  kw['titlet']="%s"
  kw['titlel']="%s"
  hpbasicwin(%s,iw,kw)
"""

# --- Template for zarrays history plots
dd3 = """
def %s(contour=0,overlay=0,iz=None,kwdict={},**kw):
  "%s. For more help, run hpdoc() or doc(hpzarray)."
  if not top.%s: return
  kw.update(kwdict)
  kw['titlet']="%s"
  hpzarray(%s,contour,overlay,iz,kw)
"""

###########################################################################
# --- Create functions for all of the hist variables
for (vname,vdoc,vunits) in hist1list:
  exec(dd1%('hp'+vname[1:],vdoc,vdoc,vunits,'top.'+vname))
for (vname,vdoc,vunits) in hist2list:
  exec(dd2%('hp'+vname[1:],vdoc,vdoc,vunits,'top.'+vname))
for (vname,vdoc) in hist3list:
  exec(dd3%('hp'+vname[1:],vdoc,'l'+vname,vdoc,'top.'+vname))
for (vname,vdoc,vexp,vunits) in hist4list:
  exec(dd1%(vname,vdoc,vdoc,vunits,vexp))
for (vname,vdoc,vexp,vunits) in hist5list:
  exec(dd2%(vname,vdoc,vdoc,vunits,vexp))

###########################################################################
# --- Create some extra specific functions for special cases.

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


###########################################################################
# --- Documentation that list all callable routines
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
hpenvy = hpyedge"""
  print histplotsdoc.__doc__

if sys.version[:5] != "1.5.1":
  for (vname,vdoc,vunits) in hist1list:
    histplotsdoc.__doc__ = histplotsdoc.__doc__+'hp'+vname[1:]+"(): "+vdoc+"\n"
  for (vname,vdoc,vunits) in hist2list:
    histplotsdoc.__doc__ = histplotsdoc.__doc__+'hp'+vname[1:]+"(): "+vdoc+"\n"
  for (vname,vdoc) in hist3list:
    histplotsdoc.__doc__ = histplotsdoc.__doc__+'hp'+vname[1:]+"(): "+vdoc+"\n"
  for (vname,vdoc,vexp,vunits) in hist4list:
    histplotsdoc.__doc__ = histplotsdoc.__doc__+'hp'+vname[1:]+"(): "+vdoc+"\n"
  for (vname,vdoc,vexp,vunits) in hist5list:
    histplotsdoc.__doc__ = histplotsdoc.__doc__+'hp'+vname[1:]+"(): "+vdoc+"\n"
