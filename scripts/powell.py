from warp import sign
powell_version = "$Id: powell.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"
# This is a conversion of routines from Numerical Recipes

########################################################################
def mnbrak(ax,bx,func):
  gold = 1.618034
  glimit = 100.
  tiny = 1.e-20
  fa = func(ax)
  fb = func(bx)
  if(fb>fa):
    ax,bx = bx,ax
    fa,fb = fb,fa
  cx = bx+gold*(bx-ax)
  fc = func(cx)
  while 1:
    if(fb>=fc):
      r = (bx-ax)*(fb-fc)
      q = (bx-cx)*(fb-fa)
      u = bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),q-r))
      ulim = bx+glimit*(cx-bx)
      if((bx-u)*(u-cx)>0.):
        fu = func(u)
        if(fu<fc):
          ax = bx
          fa = fb
          bx = u
          fb = fu
          continue
        elif(fu>fb):
          cx = u
          fc = fu
          continue
        u = cx+gold*(cx-bx)
        fu = func(u)
      elif((cx-u)*(u-ulim)>0.):
        fu = func(u)
        if(fu<fc):
          bx = cx
          cx = u
          u = cx+gold*(cx-bx)
          fb = fc
          fc = fu
          fu = func(u)
      elif((u-ulim)*(ulim-cx)>=0.):
        u = ulim
        fu = func(u)
      else:
        u = cx+gold*(cx-bx)
        fu = func(u)
      ax = bx
      bx = cx
      cx = u
      fa = fb
      fb = fc
      fc = fu
    return (ax,bx,cx,fa,fb,fc)


########################################################################
def brent(ax,bx,cx,f,tol=1.e-10):
  cgold = .3819660
  zeps = 1.0e-10
  a = min(ax,cx)
  b = max(ax,cx)
  v = bx
  w = v
  x = v
  e = 0.
  fx = f(x)
  fv = fx
  fw = fx
  for iter in xrange(100):
    xm = 0.5*(a+b)
    tol1 = tol*abs(x)+zeps
    tol2 = 2.*tol1
    if(abs(x-xm)<=(tol2-.5*(b-a))): return (x,fx)
    dontgoto1 = 1
    if(abs(e)>tol1):
      r = (x-w)*(fx-fv)
      q = (x-v)*(fx-fw)
      p = (x-v)*q-(x-w)*r
      q = 2.*(q-r)
      if(q>0.): p = -p
      q = abs(q)
      etemp = e
      e = d
      if not (abs(p)>=abs(.5*q*etemp) or p<=q*(a-x) or p>=q*(b-x)):
        d = p/q
        u = x+d
        if(u-a<tol2 or b-u<tol2): d = sign(tol1,xm-x)
	dontgoto1 = 0
    if dontgoto1:
      if(x>=xm):
        e = a-x
      else:
        e = b-x
      d = cgold*e
    if(abs(d)>=tol1):
      u = x+d
    else:
      u = x+sign(tol1,d)
    fu = f(u)
    if(fu<=fx):
      if(u>=x):
        a = x
      else:
        b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
    else:
      if(u<x):
        a = u
      else:
        b = u
      if(fu<=fw or w==x):
        v = w
        fv = fw
        w = u
        fw = fu
      elif(fu<=fv or v==x or v==w):
        v = u
        fv = fu
  print 'brent exceed maximum iterations.'

########################################################################
def linmin(p,xi,tol=1.e-10):
  global pcom,xicom
  pcom = p
  xicom = xi
  ax = 0.
  xx = 1.
  bx = 2.
  (ax,xx,bx,fa,fx,fb) =  mnbrak(ax,xx,f1dim)
  (xmin,fret) = brent(ax,xx,bx,f1dim,tol)
  xi = xmin*xi
  p = p + xi
  return (p,xi,fret)


 
########################################################################
def f1dim(x):
  global pcom,xicom
  return func(pcom + x*xicom)

def getcom():
  global pcom,xicom
  return (pcom,xicom)

########################################################################
def powell(p,xi,localfunc,ftol=1.e-10,itmax=100):
  global func
  func = localfunc

  fret=func(p)
  pt = p + 0.
  iter=0
  while 1:
    iter=iter+1
    fp=fret
    ibig=0
    delta=0.
    for i in xrange(len(p)):
      xit = xi[:,i] + 0.
      (p,xit,fret) = linmin(p,xit,ftol)
      if(abs(fp-fret)>delta):
        delta=abs(fp-fret)
        ibig=i
    if(2.*abs(fp-fret)<=ftol*(abs(fp)+abs(fret)) or iter >= itmax):
      return (p,xi,fret,iter)
    print iter,fret
    if(iter==itmax): print 'powell exceeding maximum iterations.'
    ptt = 2.*p-pt
    xit = p-pt
    pt = p + 0.
    fptt=func(ptt)
    if(fptt>=fp): continue
    t = 2.*(fp-2.*fret+fptt)*(fp-fret-delta)**2-delta*(fp-fptt)**2
    if(t>=0.): continue
    (p,xit,fret) = linmin(p,xit,ftol)
    xi[:,ibig]=xit
