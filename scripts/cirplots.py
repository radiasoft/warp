from warp import *
cirplots_version = "$Id: cirplots.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"

##############################################################################
# Handy plots for the circe package

def hcplot1(ii,istep,lvss,data,
            titlet=None,titleb=None,titlel=None,titler=None):
  if ii == -1:
    i1 = 0
    i2 = cir.nit
  else:
    i1 = ii
    i2 = ii
  for i in xrange(i1,i2,istep):
    if lvss:
      plg(data[i,:cir.jhcir+1],cir.hscir[:cir.jhcir+1])
    else:
      plg(data[i,:cir.jhcir+1],cir.htcir[cir.nit/2,:cir.jhcir+1])
  if titlet or titleb or titlel or titler:
    if not titlet: titlet = "circe"
    if not titleb:
      if lvss:
        titleb = "s (m)"
      else:
        titleb = "time (s)"
    ptitles(titlet,titleb,titlel,titler)

def hcplot2(istep,jstep,hoffset,voffset,lvss,data,
            titlet=None,titleb=None,titlel=None,titler=None):
  if lvss:
    scale = hvzcir[cir.nit/2,:cir.jhcir+1]*top.clight
  else:
    scale = ones(cir.jhcir+1)
  for j in xrange(0,cir.jhcir+1,jstep):
    if hoffset:
      plg(voffset*j + data[::istep,j],
          scale[j]*(hoffset*j + cir.htcir[::istep,j]))
    else:
      plg(voffset*j + data[::istep,j],
          scale[j]*(cir.htcir[::istep,j]-cir.htcir[cir.nit/2,j]))
  if titlet or titleb or titlel or titler:
    if not titlet: titlet = "circe"
    if not titleb:
      if lvss:
        titleb = "s (m)"
      else:
        titleb = "time (s)"
    ptitles(titlet,titleb,titlel,titler)

def hca(ii=cir.nit/2,istep=1,lvss=1,
        titlet=None,titleb=None,titlel=None,titler=None):
  if not titlel: titlel = "a (m)"
  hcplot1(ii,istep,lvss,cir.hacir,titlet,titleb,titlel,titler)

def hcb(ii=cir.nit/2,istep=1,lvss=1,
        titlet=None,titleb=None,titlel=None,titler=None):
  if not titlel: titlel = "b (m)"
  hcplot1(ii,istep,lvss,cir.hbcir,titlet,titleb,titlel,titler)

def hcap(ii=cir.nit/2,istep=1,lvss=1,
         titlet=None,titleb=None,titlel=None,titler=None):
  if not titlel: titlel = "a' (rad)"
  hcplot1(ii,istep,lvss,cir.hapcir,titlet,titleb,titlel,titler)

def hcbp(ii=cir.nit/2,istep=1,lvss=1,
         titlet=None,titleb=None,titlel=None,titler=None):
  if not titlel: titlel = "a' (rad)"
  hcplot1(ii,istep,lvss,cir.hbpcir,titlet,titleb,titlel,titler)

def hccur(istep=1,jstep=1,hoffset=0,voffset=0,lvss=0,
          titlet=None,titleb=None,titlel=None,titler=None):
  if not titlel: titlel = "current (Amps)"
  hcplot2(istep,jstep,hoffset,voffset,lvss,cir.hcur,titlet,titleb,titlel,titler)

def hcvz(istep=1,jstep=1,hoffset=0,voffset=0,lvss=0,
         titlet=None,titleb=None,titlel=None,titler=None):
  if not titlel: titlel = "Vz/Clight (1)"
  hcplot2(istep,jstep,hoffset,voffset,lvss,cir.hvzcir,titlet,titleb,titlel,titler)

def hcden(istep=1,jstep=1,hoffset=0,voffset=0,lvss=0,
          titlet=None,titleb=None,titlel=None,titler=None):
  if not titlel: titlel = "Line charge density/Clight"
  hcplot2(istep,jstep,hoffset,voffset,lvss,cir.hden,titlet,titleb,titlel,titler)

def hcdq(istep=1,jstep=1,hoffset=0,voffset=0,lvss=0,
         titlet=None,titleb=None,titlel=None,titler=None):
  if not titlel: titlel = "Charge per slice (Coulombs)"
  hcplot2(istep,jstep,hoffset,voffset,lvss,cir.hdq,titlet,titleb,titlel,titler)
