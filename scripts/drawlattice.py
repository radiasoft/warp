"""
The function drawlattice which plots the lattice.
"""
from warp import *
drawlattice_version = "$Id: drawlattice.py,v 1.1 2002/05/31 23:17:24 dave Exp $"
def drawlatticedoc():
  import drawlattice
  print drawlattice.__doc__


#############################################################################
def _getelem(ll,zs,ze,zlatmin,zlatmax):
  nn = 0
  if ll:
    ii = compress(logical_and(ze >= zlatmin,zs <= zlatmax), iota(0,len(zs)))
    if len(ii) > 0: nn = ii[0]
  else:
    ii = []
  return ii,nn

def _plotele(x,z,c,n1,n2):
  if n1 > 0:
    cc = (c*ones(n1)).astype('b')
    nn = n2*ones(n1)
    plfp(cc,x,z,nn)
    z = array(z)
    x = array(x)
    z.shape = (n1,n2)
    x.shape = (n1,n2)
    pla(transpose(x),transpose(z))

def _addelement(zs,ze,zend,ang,dh,zele,xele,zaxis,xaxis,zz,xx,zl,xl):
  dl = ze - zs
  z0 = zaxis[-1] + (zs - zend)*cos(ang)
  x0 = xaxis[-1] + (zs - zend)*sin(ang)
  z1 = z0 + dl*cos(ang)
  x1 = x0 + dl*sin(ang)
  zz = zz + list(z0 + zele*dl*cos(ang) - xele*dh*sin(ang))
  xx = xx + list(x0 + xele*dh*cos(ang) + zele*dl*sin(ang))
  zaxis = zaxis + [z0,z1]
  xaxis = xaxis + [x0,x1]
  zl = zl + [0.5*(z0+z1)]
  xl = xl + [0.5*(x0+x1)]
  return zaxis,xaxis,zz,xx,zl,xl

#############################################################################
def drawlattice(zlatmin=0,zlatmax=None,ilinflg=0,ilabflg=1,ratio=1.5,narc=10,
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
                      zlatmax defaults to zlatperi, or max z of elements.
 - ilinflg=0: when true, bends are straightened out.
 - ratio=1.5: ratio of height to lenght of elements
 - narc=10: number of points in bend arcs
For each element type are the following inputs:
 - zelem, xelem: coordinates of image to draw, normalized to
                 0 <= z <= 1, with 0 at elemzs, and 1 at elemze
               -.5 <= x <= 0.5
 - elemlab: label for each element
 - elemcolor: color index for each element

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
  if zquad is None:
    zquad = array([0., 0. , 1. ,  1. ,  0. , 0.])
    xquad = array([0., 0.5, 0.5, -0.5, -0.5, 0.])
  if zhele is None:
    zhele = array([0., 0. , 1. ,  1. ,  0. , 0.])
    xhele = array([0., 0.5, 0.5, -0.5, -0.5, 0.])
  if zemlt is None:
    zemlt = array([0., 0.  , 0.1, 0.9, 1.  ,  1.  ,  0.9,  0.1,  0.  , 0.])
    xemlt = array([0., 0.4 , 0.5, 0.5, 0.4 , -0.4 , -0.5, -0.5, -0.4 , 0.])
  if zmmlt is None:
    zmmlt = array([0., 0.  , 0.2, 0.8, 1.  ,  1.  ,  0.8,  0.2,  0.  , 0.])
    xmmlt = array([0., 0.25, 0.5, 0.5, 0.25, -0.25, -0.5, -0.5, -0.25, 0.])
  if zbgrd is None:
    zbgrd = array([0.,0.,0.333,0.333,0.667,0.667,1.,1.,0.667,0.667,0.333,0.333,0.,0.,1.,1.,0.,0.])
    xbgrd = array([0.,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.25,0.25,-0.25,-0.25,0.])
  if zpgrd is None:
    zpgrd = array([0.,0.,0.333,0.333,0.667,0.667,1.,1.,0.667,0.667,0.333,0.333,0.,0.])
    xpgrd = array([0.,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.])
  if zaccl is None:
    zaccl = array([0., 0. , 1.,  0. , 0.])
    xaccl = array([0., 0.5, 0., -0.5, 0.])
  if zdipo is None:
    zdipo = array([0., 0. , 1. , 1.,  0.5, 0.])
    xdipo = array([0., 0.5, 0.5, 0., -0.5, 0.])
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
  iq,nq = _getelem(top.quads,top.quadzs,top.quadze,zlatmin,zlatmax)
  ih,nh = _getelem(top.heles,top.helezs,top.heleze,zlatmin,zlatmax)
  ie,ne = _getelem(top.emlts,top.emltzs,top.emltze,zlatmin,zlatmax)
  im,nm = _getelem(top.mmlts,top.mmltzs,top.mmltze,zlatmin,zlatmax)
  ib,nb = _getelem(top.bgrds,top.bgrdzs,top.bgrdze,zlatmin,zlatmax)
  ip,np = _getelem(top.pgrds,top.pgrdzs,top.pgrdze,zlatmin,zlatmax)
  id,nd = _getelem(top.dipos,top.dipozs,top.dipoze,zlatmin,zlatmax)
  ia,na = _getelem(top.accls,top.acclzs,top.acclze,zlatmin,zlatmax)
  ic,nc = _getelem(top.bends,top.bendzs,top.bendze,zlatmin,zlatmax)

  # --- Get maximum element length, and set height proportional to that
  # --- so all elements are the same height.
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

  # --- First point on axis
  zaxis.append(zlatmin)
  xaxis.append(0.)

  # --- Initialize "counters"
  ang  = 0.
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
                                            zend,ang,dh,zquad,xquad,
                                            zaxis,xaxis,zq,xq,zl,xl)
      cl.append(quadlab(nq))
      zend = top.quadze[nq]
      nq = nq + 1
    if ilatspc == 'h':
      # --- load hele element
      zaxis,xaxis,zh,xh,zl,xl = _addelement(top.helezs[nh],top.heleze[nh],
                                            zend,ang,dh,zhele,xhele,
                                            zaxis,xaxis,zh,xh,zl,xl)
      cl.append(helelab(nh))
      zend = top.heleze[nh]
      nh = nh + 1
    elif ilatspc == 'e':
      # --- load emlt element 
      zaxis,xaxis,ze,xe,zl,xl = _addelement(top.emltzs[ne],top.emltze[ne],
                                            zend,ang,dh,zemlt,xemlt,
                                            zaxis,xaxis,ze,xe,zl,xl)
      zend = top.emltze[ne]
      cl.append(emltlab(ne))
      ne = ne + 1
    elif ilatspc == 'm':
      # --- load mmlt element 
      zaxis,xaxis,zm,xm,zl,xl = _addelement(top.mmltzs[nm],top.mmltze[nm],
                                            zend,ang,dh,zmmlt,xmmlt,
                                            zaxis,xaxis,zm,xm,zl,xl)
      zend = top.mmltze[nm]
      cl.append(mmltlab(nm))
      nm = nm + 1
    elif ilatspc == 'b':
      # --- load bgrd element 
      zaxis,xaxis,zb,xb,zl,xl = _addelement(top.bgrdzs[nb],top.bgrdze[nb],
                                            zend,ang,dh,zbgrd,xbgrd,
                                            zaxis,xaxis,zb,xb,zl,xl)
      zend = top.bgrdze[nb]
      cl.append(bgrdlab(nb))
      nb = nb + 1
    elif ilatspc == 'p':
      # --- load pgrd element 
      zaxis,xaxis,zp,xp,zl,xl = _addelement(top.pgrdzs[np],top.pgrdze[np],
                                            zend,ang,dh,zpgrd,xpgrd,
                                            zaxis,xaxis,zp,xp,zl,xl)
      zend = top.pgrdze[np]
      cl.append(pgrdlab(np))
      np = np + 1
    elif ilatspc == 'a':
      # --- load accelerator element 
      zaxis,xaxis,za,xa,zl,xl = _addelement(top.acclzs[na],top.acclze[na],
                                            zend,ang,dh,zaccl,xaccl,
                                            zaxis,xaxis,za,xa,zl,xl)
      zend = top.acclze[na]
      cl.append(accllab(na))
      na = na + 1
    elif ilatspc == 'd':
      # --- load dipole element
      dl = top.dipoze[nd] - top.dipozs[nd]
      if not ilinflg and nc-1 in ic and top.bendzs[nc-1] == top.dipozs[nd]:
        angb = (top.bendze[nc-1]-top.bendzs[nc-1])/top.bendrc[nc-1]
        dlt  = dl*(2./angb)*sin(angb/2.)
        z0 = zaxis[-1] - dlt*cos(ang+angb-angb/2.)
        x0 = xaxis[-1] - dlt*sin(ang+angb-angb/2.)
        zcent = z0 + top.bendrc[nc-1]*sin(ang+angb)
        xcent = x0 - top.bendrc[nc-1]*cos(ang+angb)
        aa = ang+angb - angb*zdipo
        rc = top.bendrc[nc-1]
        zd = zd + list(zcent - (rc + 0.5*sign(dh,rc)*xdipo)*sin(aa))
        xd = xd + list(xcent + (rc + 0.5*sign(dh,rc)*xdipo)*cos(aa))
        z1 = z0 + dlt*cos(ang+angb-angb/2.)
        x1 = x0 + dlt*sin(ang+angb-angb/2.)
      else:
        if nc-1 in ic and top.bendzs[nc-1] == top.dipozs[nd]: sf = 0.5
        else:                                                 sf = 1.0
        z0 = zaxis[-1] + (top.dipozs[nd] - zend)*cos(ang)
        x0 = xaxis[-1] + (top.dipozs[nd] - zend)*sin(ang)
        zd = zd + list(z0 + zdipo*dl*cos(ang) - sf*xdipo*dh*sin(ang))
        xd = xd + list(x0 + sf*xdipo*dh*cos(ang) + zdipo*dl*sin(ang))
        z1 = z0 + dl*cos(ang)
        x1 = x0 + dl*sin(ang)
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
      if ilinflg == 1:
        z0 = zaxis[-1] + (top.bendzs[nc] - zend)
        x0 = xaxis[-1]
        z1 = z0 + dl
        x1 = x0
        zc = zc + list(z0 + zbend*dl)
        xc = xc + list(x0 + xbend*dh)
        zaxis = zaxis + [z0,z1]
        xaxis = xaxis + [x0,x1]
      else:
        angb = dl/top.bendrc[nc]
        dlt  = dl*(2./angb)*sin(angb/2.)
        z0 = zaxis[-1] + (top.bendzs[nc] - zend)*cos(ang)
        x0 = xaxis[-1] + (top.bendzs[nc] - zend)*sin(ang)
        z1 = z0 + dlt*cos(ang-angb/2.)
        x1 = x0 + dlt*sin(ang-angb/2.)
        zcent = z0 + top.bendrc[nc]*sin(ang)
        xcent = x0 - top.bendrc[nc]*cos(ang)
        zc = zc + list(zcent - (top.bendrc[nc]+dh*xbend)*sin(ang - angb*zbend))
        xc = xc + list(xcent + (top.bendrc[nc]+dh*xbend)*cos(ang - angb*zbend))
        aa = ang - angb*iota(0,narc)/narc
        zaxis = zaxis + list(zcent - top.bendrc[nc]*sin(aa))
        xaxis = xaxis + list(xcent + top.bendrc[nc]*cos(aa))
        ang = ang - angb
      zend = top.bendze[nc]
      if bendlab is not None:
        zl = zl + [0.5*(z0+z1)]
        xl = xl + [0.5*(x0+x1)]
        cl = cl + [bendlab(nc)]
      nc = nc + 1
    
  # --- Add the last point
  zaxis.append(zaxis[-1] + (zlatmax - zend)*cos(ang))
  xaxis.append(xaxis[-1] + (zlatmax - zend)*sin(ang))

  # --- plot the lattice
  _plotele(xq,zq,quadcolor,len(iq),len(zquad))
  _plotele(xh,zh,helecolor,len(ih),len(zhele))
  _plotele(xe,ze,emltcolor,len(ie),len(zemlt))
  _plotele(xm,zm,mmltcolor,len(im),len(zmmlt))
  _plotele(xb,zb,bgrdcolor,len(ib),len(zbgrd))
  _plotele(xp,zp,pgrdcolor,len(ip),len(zpgrd))
  _plotele(xa,za,acclcolor,len(ia),len(zaccl))
  _plotele(xc,zc,bendcolor,len(ic),len(zbend))
  _plotele(xd,zd,dipocolor,len(id),len(zdipo))

  # --- Plot the axis
  plg(xaxis,zaxis)
  # --- Plot the labels
  for z,x,c in map(None,zl,xl,cl): plt(c,z,x,tosys=1,justify="CA")
  ptitles(titleb="meters",titlel="meters")

