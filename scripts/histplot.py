from histplots import *
histplot_version = "$Id: histplot.py,v 1.1 2000/10/16 18:34:19 dave Exp $"

def histplotdoc():
  print """
Contains the command histplot() which plots a standard set of history plots.
"""

def histplot(**kw):
  hpzbeam(kwdict=kw);fma()
  hpvbeam(kwdict=kw);fma()
  hpefld(kwdict=kw);fma()
  hpekzmbe(kwdict=kw);fma()
  hpekzbeam(kwdict=kw);fma()
  hpekzbeam(logplot=1,kwdict=kw);fma()
  hpekperp(kwdict=kw);fma()
  hpekperp(logplot=1,kwdict=kw);fma()
  hptotalke(kwdict=kw);fma()
  hptotale(kwdict=kw);fma()
  hpthermale(kwdict=kw);fma()
  hpthermale(logplot=1,kwdict=kw);fma()
  hplinechg(iz=w3d.iz_axis,kwdict=kw);fma()
  hpbmlen(kwdict=kw);fma()
  hppnum(iw=0,kwdict=kw);fma()
  hppnum(iw=-1,kwdict=kw)
  hpepsz(kwdict=kw);fma()
  hpeps6d(kwdict=kw);fma()
  for n in range(top.nzwind):
    hpepsx(iw=n,lims=1,kwdict=kw);fma()
    hpepsx(iw=n,lims=0,kwdict=kw);fma()
  for n in range(top.nzwind):
    hpepsy(iw=n,lims=1,kwdict=kw);fma()
    hpepsy(iw=n,lims=0,kwdict=kw);fma()
  for n in range(top.nzwind):
    hpepst(iw=n,lims=1,kwdict=kw);fma()
    hpepst(iw=n,lims=0,kwdict=kw);fma()
  hpepsx(iw=-1,lims=0,kwdict=kw)
  hpepsy(iw=-1,lims=0,kwdict=kw)
  hpepst(iw=-1,lims=0,kwdict=kw)
  for n in range(top.nzwind):
    hpepsnx(iw=n,lims=1,kwdict=kw);fma()
    hpepsnx(iw=n,lims=0,kwdict=kw);fma()
  for n in range(top.nzwind):
    hpepsny(iw=n,lims=1,kwdict=kw);fma()
    hpepsny(iw=n,lims=0,kwdict=kw);fma()
  for n in range(top.nzwind):
    hpepsnt(iw=n,lims=1,kwdict=kw);fma()
    hpepsnt(iw=n,lims=0,kwdict=kw);fma()
  hpepsnx(iw=-1,lims=0,kwdict=kw)
  hpepsny(iw=-1,lims=0,kwdict=kw)
  hpepsnt(iw=-1,lims=0,kwdict=kw)
  hpxbar(iw=0,kwdict=kw);fma()
  hpxbar(iw=-1,kwdict=kw);fma()
  hpybar(iw=0,kwdict=kw);fma()
  hpybar(iw=-1,kwdict=kw)
  hpxedge(iw=0 ,kwdict=kw);fma()
  hpxedge(iw=-1,kwdict=kw)
  hpyedge(iw=0 ,kwdict=kw);fma()
  hpyedge(iw=-1,kwdict=kw)
  hpvzbar(iw=0 ,beamframe=1,kwdict=kw);fma()
  hpvzbar(iw=-1,beamframe=1,kwdict=kw)
  hpvxrms(iw=0 ,kwdict=kw);fma()
  hpvxrms(iw=-1,kwdict=kw)
  hpvyrms(iw=0 ,kwdict=kw);fma()
  hpvyrms(iw=-1,kwdict=kw)
  hpvzrms(iw=0 ,kwdict=kw);fma()
  hpvzrms(iw=0 ,logplot=1,kwdict=kw);fma()
  hpvzrms(iw=-1,kwdict=kw)
  hprhomid(iw=0 ,kwdict=kw);fma()
  hprhomid(iw=-1,kwdict=kw)
  hprhomax(iw=0 ,kwdict=kw);fma()
  hprhomax(iw=-1,kwdict=kw)
  hpvzofz(contour=1,kwdict=kw);fma()
  hplinechg(contour=1,kwdict=kw);fma()
  hpcurr(contour=1,kwdict=kw);fma()
  hpvzofz(overlay=1,kwdict=kw);fma()
  hplinechg(overlay=1,kwdict=kw);fma()
  hpcurr(overlay=1,kwdict=kw);fma()
  hpvzofz(kwdict=kw);fma()
  hplinechg(kwdict=kw);fma()
  hpcurr(kwdict=kw);fma()

if sys.version[:5] != "1.5.1":
  histplot.__doc__ = hpdoc
