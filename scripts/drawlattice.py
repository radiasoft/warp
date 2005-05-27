"""Creates the function drawlattice which plots the lattice.
"""
from warp import *
drawlattice_version = "$Id: drawlattice.py,v 1.6 2005/05/27 22:38:27 dave Exp $"
def drawlatticedoc():
  import drawlattice
  print drawlattice.__doc__

#############################################################################
def _getelem(elem,zlatmin,zlatmax):
  ne = getattr(top,'n'+elem)
  nn = 0
  if ne >= 0:
    # --- Note that the elemzs and elemze will be unallocated if ne < 0.
    zs = getattr(top,elem+'zs')
    ze = getattr(top,elem+'ze')
    ii = compress(logical_and(ze >= zlatmin,zs <= zlatmax), iota(0,len(zs)))
    if len(ii) > 0: nn = ii[0]
  else:
    ii = []
  return ii,nn

def _plotele(x,z,c,n1,n2):
  if c is None: return
  if n1 > 0:
    cc = (c*ones(n1)).astype('b')
    nn = n2*ones(n1)
    plfp(cc,x,z,nn)
    z = array(z)
    x = array(x)
    z.shape = (n1,n2)
    x.shape = (n1,n2)
    pla(transpose(x),transpose(z))

def _addelement(zs,ze,zend,dh,zele,xele,zaxis,xaxis,zz,xx,zl,xl):
  dl = ze - zs
  z0 = zaxis[-1] + (zs - zend)
  x0 = xaxis[-1]
  z1 = z0 + dl
  x1 = x0
  zz = zz + list(z0 + zele*dl)
  xx = xx + list(x0 + xele*dh)
  zaxis = zaxis + [z0,z1]
  xaxis = xaxis + [x0,x1]
  zl = zl + [0.5*(z0+z1)]
  xl = xl + [0.5*(x0+x1)]
  return zaxis,xaxis,zz,xx,zl,xl

def _makearcs(x,z,narc):
  xx,zz = [],[]
  for i in xrange(len(x)-1):
    xx = xx + list(span(x[i],x[i+1],narc))
    zz = zz + list(span(z[i],z[i+1],narc))
  return array(xx),array(zz)
  

#############################################################################
def drawlattice(zlatmin=0,zlatmax=None,ilinflg=0,ilabflg=1,ratio=None,narc=10,
                height=None,
                zquad=None,xquad=None,quadlab=None,quadcolor=10,
                zhele=None,xhele=None,helelab=None,helecolor=35,
                zemlt=None,xemlt=None,emltlab=None,emltcolor=60,
                zmmlt=None,xmmlt=None,mmltlab=None,mmltcolor=85,
                zbgrd=None,xbgrd=None,bgrdlab=None,bgrdcolor=110,
                zpgrd=None,xpgrd=None,pgrdlab=None,pgrdcolor=135,
                zaccl=None,xaccl=None,accllab=None,acclcolor=160,
                zdipo=None,xdipo=None,dipolab=None,dipocolor=185,
                zbend=None,xbend=None,bendlab=None,bendcolor=20):
  """
Draws the lattice.  All lattice elements starting in the interval 
 - zlatmin=0,zlatmax: lattice elements within the range are plotted.
                      zlatmax defaults to zlatperi when it is nonzero,
                      or maximum z of elements.
 - ilinflg=0: when true, bends are straightened out. The default is to plot
              elements in laboratory frame, showing bends.
 - ilabflg=1: when true, labels on each element are plotted
 - ratio=1.5: ratio of height to length of elements
 - narc=10: number of points in bend arcs, any curved line in a bend
 - height=(zlatmax-zlatmin)/10.: height of elements blocks, in meters

For each element type, the following inputs can be optionally given to affect
how the element is displayed. Here 'elem' refers is replaced by one of the
element names, such as 'quad' or 'bgrd'.:
 - zelem, xelem: coordinates of image to draw, normalized to
                 0 <= z <= 1, with 0 at elemzs, and 1 at elemze
               -.5 <= x <= 0.5
 - elemlab: label for each element, must be either a string or a function
            that takes one arguement, the integer id of the element being
            plotted.
 - elemcolor: color index for each element - index is relative to the active
              palette.

Some lattice element symbols and labels are as follows:

 element       symbol       label/comments

 drift         -----        none, if bent, represents a bent beam pipe 
                                  with no lattice element.
                ______      
               |      |     f   , d   ; focus, defocus electric quad 
 quadrupole    |------|     F   , D   ; focus, defocus magnetic quad
               |      |    
                ------

 emlt          /-----\ 
               |-----|
               \-----/
                ______
               |      |     p    ; electric dipole
 dipole        |------|     P    ; magnetic dipole
                \    /
                  \/
                _
               |  -         
 accelerator   |-----       A    ; acceleration module 
               |_ - 

Dipoles are drawn consistent with the existence of any bent beam pipes. 
The routine assumes no particular ordering between different element
type, and should draw any general lattice.
  """

  # --- Set default value of zlatmax
  if zlatmax is None:
    if top.zlatperi > 0:
      zlatmax = top.zlatperi
    else:
      zlatmax = zlatmin
      if top.quads: zlatmax = max(max(top.quadze),zlatmax)
      if top.heles: zlatmax = max(max(top.heleze),zlatmax)
      if top.emlts: zlatmax = max(max(top.emltze),zlatmax)
      if top.mmlts: zlatmax = max(max(top.mmltze),zlatmax)
      if top.bgrds: zlatmax = max(max(top.bgrdze),zlatmax)
      if top.pgrds: zlatmax = max(max(top.pgrdze),zlatmax)
      if top.dipos: zlatmax = max(max(top.dipoze),zlatmax)
      if top.bends: zlatmax = max(max(top.bendze),zlatmax)
      if top.accls: zlatmax = max(max(top.acclze),zlatmax)

  # --- Shapes for the elements
  # --- Note that 0 <= z <= 1 and -0.5 <= x <= +0.5
  if zquad is None:
    zquad = array([0., 0. , 1. ,  1. ,  0. , 0.])
    xquad = array([0., 0.5, 0.5, -0.5, -0.5, 0.])
    if ilinflg == 0: xquad,zquad = _makearcs(xquad,zquad,narc)
  if zhele is None:
    zhele = array([0., 0. , 1. ,  1. ,  0. , 0.])
    xhele = array([0., 0.5, 0.5, -0.5, -0.5, 0.])
    if ilinflg == 0: xhele,zhele = _makearcs(xhele,zhele,narc)
  if zemlt is None:
    zemlt = array([0., 0.  , 0.1, 0.9, 1.  ,  1.  ,  0.9,  0.1,  0.  , 0.])
    xemlt = array([0., 0.4 , 0.5, 0.5, 0.4 , -0.4 , -0.5, -0.5, -0.4 , 0.])
    if ilinflg == 0: xemlt,zemlt = _makearcs(xemlt,zemlt,narc)
  if zmmlt is None:
    zmmlt = array([0., 0.  , 0.2, 0.8, 1.  ,  1.  ,  0.8,  0.2,  0.  , 0.])
    xmmlt = array([0., 0.25, 0.5, 0.5, 0.25, -0.25, -0.5, -0.5, -0.25, 0.])
    if ilinflg == 0: xmmlt,zmmlt = _makearcs(xmmlt,zmmlt,narc)
  if zbgrd is None:
    zbgrd = array([0.,0.,0.333,0.333,0.667,0.667,1.,1.,0.667,0.667,0.333,0.333,0.,0.,1.,1.,0.,0.])
    xbgrd = array([0.,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.167,0.167,-0.167,-0.167,0.])
    if ilinflg == 0: xbgrd,zbgrd = _makearcs(xbgrd,zbgrd,narc)
  if zpgrd is None:
    zpgrd = array([0.,0.,0.333,0.333,0.667,0.667,1.,1.,0.667,0.667,0.333,0.333,0.,0.])
    xpgrd = array([0.,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.])
    if ilinflg == 0: xpgrd,zpgrd = _makearcs(xpgrd,zpgrd,narc)
  if zaccl is None:
    zaccl = array([0., 0. , 1.,  0. , 0.])
    xaccl = array([0., 0.5, 0., -0.5, 0.])
    if ilinflg == 0: xaccl,zaccl = _makearcs(xaccl,zaccl,narc)
  if zdipo is None:
    zdipo = array([0., 0. , 1. , 1.,  0.5, 0.])
    xdipo = array([0., 0.5, 0.5, 0., -0.5, 0.])
    if ilinflg == 0: xdipo,zdipo = _makearcs(xdipo,zdipo,narc)
  if zbend is None:
    zbend = array([0.] + list(1.*iota(0,narc)/narc) +
                         list(1.*iota(narc,0,-1)/narc) + [0.])
    xbend = array([0.] + (narc+1)*[0.5] + (narc+1)*[-0.5] + [0.])

  # --- Labels for elements
  if quadlab is None:
    def quadlab(nq):
      if top.quadde[nq] != 0.:
        if top.quadde[nq] > 0: return "f%d"%nq
        else:                  return "d%d"%nq
      else:
        if top.quaddb[nq] > 0: return "F%d"%nq
        else:                  return "D%d"%nq
  elif type(quadlab) is not FunctionType:
    def quadlab(nq,quadlab=quadlab):
      return quadlab

  if helelab is None:
    def helelab(nh):
        return "H%d"%nh
  elif type(helelab) is not FunctionType:
    def helelab(nh,helelab=helelab):
      return helelab

  if dipolab is None:
    def dipolab(nd):
      if top.dipoex[nd] != 0.: return "p%d"%nd
      else:                    return "P%d"%nd
  elif type(dipolab) is not FunctionType:
    def dipolab(nd,dipolab=dipolab):
      return dipolab

  if accllab is None:
    def accllab(na):
        return "A%d"%na
  elif type(accllab) is not FunctionType:
    def accllab(na,accllab=accllab):
      return accllab

  if emltlab is None:
    def emltlab(ne):
        return "E%d"%ne
  elif type(emltlab) is not FunctionType:
    def emltlab(ne,emltlab=emltlab):
      return emltlab

  if mmltlab is None:
    def mmltlab(nm):
        return "M%d"%nm
  elif type(mmltlab) is not FunctionType:
    def mmltlab(nm,mmltlab=mmltlab):
      return mmltlab

  if bgrdlab is None:
    def bgrdlab(nb):
        return "B%d"%nb
  elif type(bgrdlab) is not FunctionType:
    def bgrdlab(nb,bgrdlab=bgrdlab):
      return bgrdlab

  if pgrdlab is None:
    def pgrdlab(np):
        return "P%d"%np
  elif type(pgrdlab) is not FunctionType:
    def pgrdlab(np,pgrdlab=pgrdlab):
      return pgrdlab

  if bendlab is not None and type(bendlab) is not FunctionType:
    def bendlab(np,bendlab=bendlab):
      return bendlab

  # --- Create work arrays
  zaxis = []
  xaxis = []
  zq = []
  xq = []
  zh = []
  xh = []
  ze = []
  xe = []
  zm = []
  xm = []
  zb = []
  xb = []
  zp = []
  xp = []
  za = []
  xa = []
  zd = []
  xd = []
  zc = []
  xc = []
  zl = []
  xl = []
  cl = []

  # --- find element indices in plot region
  iq,nq = _getelem('quad',zlatmin,zlatmax)
  ih,nh = _getelem('hele',zlatmin,zlatmax)
  ie,ne = _getelem('emlt',zlatmin,zlatmax)
  im,nm = _getelem('mmlt',zlatmin,zlatmax)
  ib,nb = _getelem('bgrd',zlatmin,zlatmax)
  ip,np = _getelem('pgrd',zlatmin,zlatmax)
  id,nd = _getelem('dipo',zlatmin,zlatmax)
  ia,na = _getelem('accl',zlatmin,zlatmax)
  ic,nc = _getelem('bend',zlatmin,zlatmax)

  # --- Get maximum element length, and set height proportional to that
  # --- so all elements are the same height.
  if ratio is not None:
    hmax = 0.
    if len(iq) > 0: hmax = max(hmax,max(take(top.quadze-top.quadzs,iq)))
    if len(ih) > 0: hmax = max(hmax,max(take(top.heleze-top.helezs,ih)))
    if len(ie) > 0: hmax = max(hmax,max(take(top.emltze-top.emltzs,ie)))
    if len(im) > 0: hmax = max(hmax,max(take(top.mmltze-top.mmltzs,im)))
    if len(ib) > 0: hmax = max(hmax,max(take(top.bgrdze-top.bgrdzs,ib)))
    if len(ip) > 0: hmax = max(hmax,max(take(top.pgrdze-top.pgrdzs,ip)))
    if len(id) > 0: hmax = max(hmax,max(take(top.dipoze-top.dipozs,id)))
    if len(ia) > 0: hmax = max(hmax,max(take(top.acclze-top.acclzs,ia)))
    if len(ic) > 0: hmax = max(hmax,max(take(top.bendze-top.bendzs,ic)))
    dh = ratio*hmax
  elif height is not None:
    dh = height
  else:
    dh = (zlatmax - zlatmin)/10.

  # --- First point on axis
  zaxis.append(zlatmin)
  xaxis.append(0.)

  # --- Initialize "counters"
  zend = zlatmin

  # --- loop over all elements in range
  for n in xrange(len(iq)+len(ih)+len(id)+len(ic)+len(ia)+len(ie)+len(im)+len(ib)+len(ip)):

    # --- find element with minimum z and flag with ilatspc  
    zqmin=zhmin=zemin=zmmin=zbmin=zpmin=zdmin=zamin=zcmin=zlatmax
    if len(iq) > 0 and nq in iq: zqmin = top.quadzs[nq]
    if len(ih) > 0 and nh in ih: zhmin = top.helezs[nh]
    if len(ie) > 0 and ne in ie: zemin = top.emltzs[ne]
    if len(im) > 0 and nm in im: zmmin = top.mmltzs[nm]
    if len(ib) > 0 and nb in ib: zbmin = top.bgrdzs[nb]
    if len(ip) > 0 and np in ip: zpmin = top.pgrdzs[np]
    if len(id) > 0 and nd in id: zdmin = top.dipozs[nd]
    if len(ia) > 0 and na in ia: zamin = top.acclzs[na]
    if len(ic) > 0 and nc in ic: zcmin = top.bendzs[nc]
    ii = argmin([zqmin,zhmin,zemin,zmmin,zbmin,zpmin,zamin,zcmin,zdmin])
    ilatspc =   ['q'  ,'h'  ,'e'  ,'m'  ,'b'  ,'p'  ,'a'  ,'c'  ,'d'][ii]

    if ilatspc == 'q':
      # --- load quadrapole element (focussing or defocussing)
      zaxis,xaxis,zq,xq,zl,xl = _addelement(top.quadzs[nq],top.quadze[nq],
                                            zend,dh,zquad,xquad,
                                            zaxis,xaxis,zq,xq,zl,xl)
      cl.append(quadlab(nq))
      zend = top.quadze[nq]
      nq = nq + 1
    if ilatspc == 'h':
      # --- load hele element
      zaxis,xaxis,zh,xh,zl,xl = _addelement(top.helezs[nh],top.heleze[nh],
                                            zend,dh,zhele,xhele,
                                            zaxis,xaxis,zh,xh,zl,xl)
      cl.append(helelab(nh))
      zend = top.heleze[nh]
      nh = nh + 1
    elif ilatspc == 'e':
      # --- load emlt element 
      zaxis,xaxis,ze,xe,zl,xl = _addelement(top.emltzs[ne],top.emltze[ne],
                                            zend,dh,zemlt,xemlt,
                                            zaxis,xaxis,ze,xe,zl,xl)
      zend = top.emltze[ne]
      cl.append(emltlab(ne))
      ne = ne + 1
    elif ilatspc == 'm':
      # --- load mmlt element 
      zaxis,xaxis,zm,xm,zl,xl = _addelement(top.mmltzs[nm],top.mmltze[nm],
                                            zend,dh,zmmlt,xmmlt,
                                            zaxis,xaxis,zm,xm,zl,xl)
      zend = top.mmltze[nm]
      cl.append(mmltlab(nm))
      nm = nm + 1
    elif ilatspc == 'b':
      # --- load bgrd element 
      zaxis,xaxis,zb,xb,zl,xl = _addelement(top.bgrdzs[nb],top.bgrdze[nb],
                                            zend,dh,zbgrd,xbgrd,
                                            zaxis,xaxis,zb,xb,zl,xl)
      zend = top.bgrdze[nb]
      cl.append(bgrdlab(nb))
      nb = nb + 1
    elif ilatspc == 'p':
      # --- load pgrd element 
      zaxis,xaxis,zp,xp,zl,xl = _addelement(top.pgrdzs[np],top.pgrdze[np],
                                            zend,dh,zpgrd,xpgrd,
                                            zaxis,xaxis,zp,xp,zl,xl)
      zend = top.pgrdze[np]
      cl.append(pgrdlab(np))
      np = np + 1
    elif ilatspc == 'a':
      # --- load accelerator element 
      zaxis,xaxis,za,xa,zl,xl = _addelement(top.acclzs[na],top.acclze[na],
                                            zend,dh,zaccl,xaccl,
                                            zaxis,xaxis,za,xa,zl,xl)
      zend = top.acclze[na]
      cl.append(accllab(na))
      na = na + 1
    elif ilatspc == 'd':
      # --- load dipole element
      dl = top.dipoze[nd] - top.dipozs[nd]
      if nc-1 in ic and top.bendzs[nc-1] == top.dipozs[nd]: sf = 0.5
      else:                                                 sf = 1.0
      z0 = zaxis[-1] + (top.dipozs[nd] - zend)
      x0 = xaxis[-1]
      zd = zd + list(z0 + zdipo*dl)
      xd = xd + list(x0 + sf*xdipo*dh)
      z1 = z0 + dl
      x1 = x0
      zaxis = zaxis + [z0,z1]
      xaxis = xaxis + [x0,x1]
      zend = top.dipoze[nd] 
      if dipolab is not None:
        zl = zl + [0.5*(z0+z1)]
        xl = xl + [0.5*(x0+x1)]
        cl = cl + [dipolab(nd)]
      nd = nd + 1
    elif ilatspc == 'c':
      # --- load bent beam pipe with no element 
      dl = top.bendze[nc] - top.bendzs[nc]
      z0 = zaxis[-1] + (top.bendzs[nc] - zend)
      x0 = xaxis[-1]
      z1 = z0 + dl
      x1 = x0
      zc = zc + list(z0 + zbend*dl)
      xc = xc + list(x0 + xbend*dh)
      zaxis = zaxis + [z0,z1]
      xaxis = xaxis + [x0,x1]
      zend = top.bendze[nc]
      if bendlab is not None:
        zl = zl + [0.5*(z0+z1)]
        xl = xl + [0.5*(x0+x1)]
        cl = cl + [bendlab(nc)]
      nc = nc + 1
    
  # --- Add the last point
  zaxis.append(zaxis[-1] + (zlatmax - zend))
  xaxis.append(xaxis[-1])

  if ilinflg == 0:
    xq,zq = array(xq),array(zq)
    xh,zh = array(xh),array(zh)
    xe,ze = array(xe),array(ze)
    xm,zm = array(xm),array(zm)
    xb,zb = array(xb),array(zb)
    xp,zp = array(xp),array(zp)
    xa,za = array(xa),array(za)
    xc,zc = array(xc),array(zc)
    xd,zd = array(xd),array(zd)
    xaxis,zaxis = _makearcs(xaxis,zaxis,narc)
    xl,zl = array(xl),array(zl)

    tolabfrm(zlatmin,len(xq),xq,zq)
    tolabfrm(zlatmin,len(xh),xh,zh)
    tolabfrm(zlatmin,len(xe),xe,ze)
    tolabfrm(zlatmin,len(xm),xm,zm)
    tolabfrm(zlatmin,len(xb),xb,zb)
    tolabfrm(zlatmin,len(xp),xp,zp)
    tolabfrm(zlatmin,len(xa),xa,za)
    tolabfrm(zlatmin,len(xc),xc,zc)
    tolabfrm(zlatmin,len(xd),xd,zd)
    tolabfrm(zlatmin,len(xaxis),xaxis,zaxis)
    tolabfrm(zlatmin,len(xl),xl,zl)

  # --- plot the lattice
  _plotele(xc,zc,bendcolor,len(ic),len(zbend))
  _plotele(xq,zq,quadcolor,len(iq),len(zquad))
  _plotele(xh,zh,helecolor,len(ih),len(zhele))
  _plotele(xe,ze,emltcolor,len(ie),len(zemlt))
  _plotele(xm,zm,mmltcolor,len(im),len(zmmlt))
  _plotele(xb,zb,bgrdcolor,len(ib),len(zbgrd))
  _plotele(xp,zp,pgrdcolor,len(ip),len(zpgrd))
  _plotele(xa,za,acclcolor,len(ia),len(zaccl))
  _plotele(xd,zd,dipocolor,len(id),len(zdipo))

  # --- Plot the axis
  plg(xaxis,zaxis)
  # --- Plot the labels
  if ilabflg:
    for z,x,c in map(None,zl,xl,cl):
      plt(c,z,x,tosys=1,justify="CA")
  ptitles(titleb="meters",titlel="meters")

