from warp import *
fringedquads_version = "$Id: fringedquads.py,v 1.11 2008/11/21 20:49:47 dave Exp $"
# --- Set up quadrupoles with fringes.
# --- Currently uses form proportional to tanh(cot(z)), which is essentially
# --- a linear falloff with rounded corners to match derivatives.
# --- Extracted from IRE decks and generalized.
# --- DPG 7/6/1999

def cot(x):
  return cos(x)/sin(x)

def fringedquads(fringelen=4.5e-2,fringescale=0.5,lpseudooct=true,
                 usequads=false,firstfringe=0,npoints=1000,lclear=true,
                 lscale=true,fringe=None,fringep=None,fringepp=None):
  """
Set up quadrupoles with fringes based on hard edged quadrupoles, either
quad or hele.  Currently uses form proportional to tanh(cot(z)), which is
essentially a linear falloff with rounded corners to match derivatives.
  - fringelen = 4.5e-2, total length of the fringe 
  - fringescale = 0.5, fraction of the fringe which extends beyond the end of
                       the hard edged quadrupole
  - lpseudooct = true, turns pseudo-octopole term on
  - usequads = false, when true, forces the use of quad elements
  - firstfringe = 0, first quad which should have fringes added
  - npoints = 1000, total number of points in axial profile
  - lclear = true, forces the hard edged element to be zerod out
  - lscale = true, scales field so that the integral of the quad strength
                  over the element is the same as the hard edged element.
                  Otherwise field at center is set directly to field of hard
                  edged element.
  - fringe = function returning the fringe. The function must have a
             single argument n. The function must return a 1-D array of
             length n+1 with the fringe at evenly spaced points along
             the fringed, including the end points.  (The value of n will
             depend on npoints and the fringelen.) The profile should satisfy
             f(zmin)=0, f(zmax)=1, f'(zmin)=f'(zmax)=0. The default form
             of the fringe is tanh(cot(z)) (where 0<=z<=pi).
  - fringep = first derivative of fringe, scaled so f'(i) ~ (f(i+1)-f(i-1))/2
  - fringepp = second derivative of fringe
Note that fringelen,fringescale, and lscale can all be arrays the same shape
as the hard edged element arrays so that the parameters can be varied from
element to element. If fringe is supplied but either fringep or fringepp are
not, then the derivatives will be done with a finite difference of fringe.
  """

  # --- Check for quads or heles. Copy either the quad or hele data into
  # --- temporary arrays (using same name as quad elements). This avoids
  # --- having to check which type of element to use in the work below.
  quads = false
  heles = false
  if top.quads or usequads:
    quads = true
    quadzs = top.quadzs[firstfringe:]
    quadze = top.quadze[firstfringe:]
    quadap = top.quadap[firstfringe:]
    quadde = top.quadde[firstfringe:]
    quaddb = top.quaddb[firstfringe:]
    quadox = top.qoffx[firstfringe:]
    quadoy = top.qoffy[firstfringe:]
    quadrr = top.quadrr[firstfringe:]
    quadrl = top.quadrl[firstfringe:]
    quadgl = top.quadgl[firstfringe:]
    quadgp = top.quadgp[firstfringe:]
    quadpw = top.quadpw[firstfringe:]
    quadpa = top.quadpa[firstfringe:]
    quadpe = zeros(top.nquad+1-firstfringe,'d')
    quadpm = zeros(top.nquad+1-firstfringe,'d')
  elif top.heles:
    heles = true
    quadzs = []
    quadze = []
    quadap = []
    quadde = []
    quaddb = []
    quadox = []
    quadoy = []
    quadrr = []
    quadrl = []
    quadgl = []
    quadgp = []
    quadpw = []
    quadpa = []
    quadpe = []
    quadpm = []
    for ih in xrange(firstfringe,top.nhele+1):
      for im in xrange(top.nhmlt):
        if top.hele_n[im,ih]==2 and top.hele_v[im,ih]==0:
          quadzs.append(top.helezs[ih])
          quadze.append(top.heleze[ih])
          quadap.append(top.heleap[ih])
          quadde.append(top.heleae[im,ih])
          quaddb.append(top.heleam[im,ih])
          quadox.append(top.heleox[ih])
          quadoy.append(top.heleoy[ih])
          quadrr.append(top.helerr[ih])
          quadrl.append(top.helerl[ih])
          quadgl.append(top.helegl[ih])
          quadgp.append(top.helegp[ih])
          quadpw.append(top.helepw[ih])
          quadpa.append(top.helepa[ih])
          quadpe.append(top.helepe[im,ih])
          quadpm.append(top.helepm[im,ih])
          break # --- Just to prevent overlapping elements.
    quadzs = array(quadzs)
    quadze = array(quadze)
    quadap = array(quadap)
    quadde = array(quadde)
    quaddb = array(quaddb)
    quadox = array(quadox)
    quadoy = array(quadoy)
    quadrr = array(quadrr)
    quadrl = array(quadrl)
    quadgl = array(quadgl)
    quadgp = array(quadgp)
    quadpw = array(quadpw)
    quadpa = array(quadpa)
    quadpe = array(quadpe)
    quadpm = array(quadpm)
  else:
    raise "No quad elements were specified, maybe need to pass usequads=true"

  # --- Make sure that the fringelen and scale are arrays.
  fringelen = fringelen*ones(len(quadzs),'d')
  fringescale = fringescale*ones(len(quadzs),'d')
  lscale = lscale*ones(len(quadzs),'d')

  # --- Get the number of each kind of quad, electric and magnetic and
  # --- allocate the space appropriately.
  top.nemlt = len(oldnonzero(quadde)) - 1
  top.nmmlt = len(oldnonzero(quaddb)) - 1
  top.neerr = top.nemlt
  top.nmerr = top.nmmlt

  gchange("Lattice",0)

  # --- Setup the multipole components space
  if len(oldnonzero(quadde)) > 0:
    top.nemltsets = top.nemlt + 1
    top.nesmult = 1 # Number of multipole components for the electric quads
    if lpseudooct: top.nesmult = 2
    top.nzemltmax = npoints

  if len(oldnonzero(quaddb)) > 0:
    top.nmmltsets = top.nmmlt + 1
    top.nmsmult = 1 # Number of multipole components for the magnetic quad
    if lpseudooct: top.nmsmult = 2
    top.nzmmltmax = npoints

  gchange("Mult_data",0)

  # --- Set parameters constant for all elements
  if len(oldnonzero(quadde)) > 0:
    top.nzemlt = npoints
    top.emlt_n[0] = 2.
    top.emlt_v[0] = 0.
    if (top.nesmult == 2):
      top.emlt_n[1] = 2.
      top.emlt_v[1] = 1.
  if len(oldnonzero(quaddb)) > 0:
    top.nzmmlt = npoints
    top.mmlt_n[0] = 2.
    top.mmlt_v[0] = 0.
    if (top.nmsmult == 2):
      top.mmlt_n[1] = 2.
      top.mmlt_v[1] = 1.

  # --- Function that defines the profile of the fringe fields
  # --- This function should have the following properties
  # --- f(0) = 0, f(1) = 1.
  # --- f'(0) = 0, f'(1) = 0.
  if not fringe:
    # --- use the default form
    def fringe(n):
      zz = pi*iota(1,n-1)/n
      ff = zeros(n+1,Float)
      ff[n] = 1.
      ff[1:n] = 0.5*(1. - tanh(cot(zz)))
      return ff
    def fringep(n):
      zz = pi*iota(1,n-1)/n
      ffp = zeros(n+1,Float)
      ffp[1:n] = pi/n*0.5*(1.-tanh(cot(zz))**2)*(1.+cot(zz)**2)
      return ffp
    def fringepp(n):
      zz = pi*iota(1,n-1)/n
      ffpp = zeros(n+1,Float)
      ffpp[1:n]=(pi/n)**2*0.5*(2.*tanh(cot(zz))*(1.-tanh(cot(zz))**2)*
        (1.+cot(zz)**2)**2-(1.-tanh(cot(zz))**2)*2.*cot(zz)*(1.+cot(zz)**2))
      return ffpp
  if not fringep:
    # --- Use finite difference of fringe
    def fringep(n,fringe=fringe):
      ffp = zeros(n+1,Float)
      ff = fringe(n)
      ffp[1:n] = (ff[2:] - ff[:-2])/2.
      return ffp
    if not fringepp:
      # --- If fringepp is also not defined, use finite difference of fringe
      def fringepp(n,fringe=fringe):
        ffpp = zeros(n+1,Float)
        ff = fringe(n)
        ffpp[1:n] = (ff[2:] + ff[:-2] - 2.*ff[1:-1])/2.
        return ffpp
  if not fringepp:
    # --- If fringep was defined but not fringepp, then use finite difference
    # --- on fringep directly.
    def fringepp(n,fringep=fringep):
      ffpp = zeros(n+1,Float)
      ffp = fringep(n)
      ffpp[1:n] = (ffp[2:] - ffp[:-2])/2.
      return ffpp

  # --- Loop over quads
  neq = -1
  nmq = -1
  for iq in xrange(len(quadzs)):

    # --- Set up the lattice elements
    if quadde[iq] != 0.:
      # --- Electric quads
      neq = neq + 1
      top.emltzs[neq] = quadzs[iq] - fringelen[iq]*fringescale[iq]
      top.emltze[neq] = quadze[iq] + fringelen[iq]*fringescale[iq]
      top.emltid[neq] = neq + 1
      top.emltox[neq] = quadox[iq]
      top.emltoy[neq] = quadoy[iq]
      top.emltap[neq] = quadap[iq]
      top.emltrr[neq] = quadrr[iq]
      top.emltrl[neq] = quadrl[iq]
      top.emltgl[neq] = quadgl[iq]
      top.emltgp[neq] = quadgp[iq]
      top.emltpw[neq] = quadpw[iq]
      top.emltpa[neq] = quadpa[iq]
      top.emltph[neq] = quadpe[iq]
      top.dzemlt[neq] = (top.emltze[neq] - top.emltzs[neq])/top.nzemltmax
      top.esemlt[:,0,neq] = 1.
      nend = int(fringelen[iq]/top.dzemlt[neq])
      nmax = top.nzemltmax
      top.esemlt[0:nend+1,0,neq] = fringe(nend)
      top.esemlt[nmax:nmax-nend-1:-1,0,neq] = fringe(nend)
      top.esemltp[0:nend+1,0,neq] = fringep(nend)
      top.esemltp[nmax:nmax-nend-1:-1,0,neq] = -fringep(nend)
      if lscale[iq]:
        scalefactor = ((quadze[iq] - quadzs[iq])*quadde[iq]/
                      (sum(top.esemlt[:,0,neq])*top.dzemlt[neq]))
      else:
        scalefactor = quadde[iq]
      top.esemlt[:,0,neq] = top.esemlt[:,0,neq]*scalefactor
      top.esemltp[:,0,neq] = top.esemltp[:,0,neq]*scalefactor/top.dzemlt[neq]
      if (top.nesmult == 2):
        top.esemlt[:,1,neq] = 0.
        top.esemlt[0:nend+1,1,neq] = fringepp(nend)
        top.esemlt[nmax:nmax-nend-1:-1,1,neq]=fringepp(nend)
        top.esemlt[:,1,neq] = -(top.esemlt[:,1,neq]*scalefactor/
         (top.dzemlt[neq]**2*(4.*top.emlt_v[1]*(top.emlt_n[1]+top.emlt_v[1]))))
    if quaddb[iq] != 0.:
      # --- Magnetic quads
      nmq = nmq + 1
      top.mmltzs[nmq] = quadzs[iq] - fringelen[iq]*fringescale[iq]
      top.mmltze[nmq] = quadze[iq] + fringelen[iq]*fringescale[iq]
      top.mmltid[nmq] = nmq + 1
      top.mmltox[nmq] = quadox[iq]
      top.mmltoy[nmq] = quadoy[iq]
      top.mmltap[nmq] = quadap[iq]
      top.mmltph[nmq] = quadpm[iq]
      top.dzmmlt[nmq] = (top.mmltze[nmq] - top.mmltzs[nmq])/top.nzmmltmax
      top.msmmlt[:,0,nmq] = 1.
      nend = int(fringelen[iq]/top.dzmmlt[nmq])
      nmax = top.nzmmltmax
      top.msmmlt[0:nend+1,0,nmq] = fringe(nend)
      top.msmmlt[nmax:nmax-nend-1:-1,0,nmq] = fringe(nend)
      top.msmmltp[0:nend+1,0,nmq] = fringep(nend)
      top.msmmltp[nmax:nmax-nend-1:-1,0,nmq] = -fringep(nend)
      if lscale[iq]:
        scalefactor = ((quadze[iq] - quadzs[iq])*quaddb[iq]/
                      (sum(top.msmmlt[:,0,nmq])*top.dzmmlt[nmq]))
      else:
        scalefactor = quaddb[iq]
      top.msmmlt[:,0,nmq] = top.msmmlt[:,0,nmq]*scalefactor
      top.msmmltp[:,0,nmq] = top.msmmltp[:,0,nmq]*scalefactor/top.dzmmlt[nmq]
      if (top.nmsmult == 2):
        top.msmmlt[:,1,nmq] = 0.
        top.msmmlt[0:nend+1,1,nmq] = fringepp(nend)
        top.msmmlt[nmax:nmax-nend-1:-1,1,nmq]=fringepp(nend)
        top.msmmlt[:,1,nmq] = -(top.msmmlt[:,1,nmq]*scalefactor/
         (top.dzmmlt[nmq]**2*(4.*top.mmlt_v[1]*(top.mmlt_n[1]+top.mmlt_v[1]))))

  # --- Turn off the hard edged elements by switching off the flag.
  # --- Only do this if all of the hard edged elements have fringes added.
  if quads and firstfringe == 0: top.quads = false
  if heles and firstfringe == 0: top.heles = false

  # --- If requested, delete the hard edged data.
  if lclear:
    if quads:
      top.nquad = firstfringe - 1
      top.nqerr = firstfringe - 1
      gchange("Lattice")
    if heles:
      top.nhele = firstfringe - 1
      gchange("Lattice")
  else:
    if quads:
      top.quadde[firstfringe:] = 0.
      top.quaddb[firstfringe:] = 0.
    if heles:
      top.heleae[:,firstfringe:] = 0.
      top.heleam[:,firstfringe:] = 0.
      top.heleep[:,firstfringe:] = 0.
      top.helemp[:,firstfringe:] = 0.

  # --- Turn on the fringed elements
  if top.nemlt > 0: top.emlts = true
  if top.nmmlt > 0: top.mmlts = true

  # --- Make sure internal lattice arrays set properly
  resetlat()
  setlatt()

def testfringedequads():
  ppp=zeros((w3d.nx+1,w3d.ny+1,101),'d')
  pp1=zeros((w3d.nx+1,w3d.ny+1,101),'d')
  xx=w3d.xmesh*ones(w3d.nx+1,'d')[:,newaxis]-0.5*(w3d.xmmax+w3d.xmmin)
  yy=w3d.ymesh[:,newaxis]*ones(w3d.ny+1,'d')-0.5*(w3d.ymmax+w3d.ymmin)
  rr=sqrt(xx**2+yy**2)
  tt=arctan(yy/xx)
  tt[w3d.nx/2,w3d.ny/2]=0.
  ireg=ones((w3d.nx+1,w3d.ny+1),'l')
  for i in xrange(101):
    ppp[:,:,i] = top.esemlt[i,0,0]*rr**2*cos(2*(tt+top.emltpe[i]))
    pp1[:,:,i] = top.esemlt[i,1,0]*rr**4*cos(2*(tt+top.emltpe[i]))

