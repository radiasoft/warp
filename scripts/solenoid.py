"""Functions to generate solenoid lattice elements.
The solenoid field is obtained from the analytic field profile of a cylindrical
current sheet. Thanks to Ed Lee for the formula. The following functions are
provided:
 - addsolenoid(zi,zf,ri,ro,current,nzpoints=10000,fringelen=10.): adds a
     single solenoid
 - B0, B0p, B0pp, B0ppp: return B0 and its derivatives
 - Bz, Br: return the B field at a given r,z coordinate (using same truncation)

The field on axis is given by
  B0(z) = (k*mu0*((l - 2*z)/sqrt(4*R**2 + (l - 2*z)**2) +
                  (l + 2*z)/sqrt(4*R**2 + (l + 2*z)**2)))/2.
where k is the current in units Ampere-turns/meter, mu0 has the standard
meaning, and R is the radius of the current sheet. Starting with this
expression, the field off axis is given by the multipole expansion.

Bz(r,z) = B0 - B0''*r**2/4 + ...
Br(r,z) = -B0'*r/2 + B0'''*r**3/16 - ...

As currently written, the series is truncated and only the terms shown are
included, up B0'''.
"""
from warp import *
from lattice import addnewmmlt
solenoid_version = "$Id: solenoid.py,v 1.6 2007/10/09 23:28:41 dave Exp $"

def solenoiddoc():
  import solenoid
  print solenoid.__doc__

def B0(z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return (k*mu0*((l - 2*z)/sqrt(c3) + (l + 2*z)/sqrt(c4)))/2.

def B0p(z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 4*k*mu0*R**2*(-(c3)**-1.5 + (c4)**-1.5)

def B0pp(z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return (24*k*mu0*R**2*(-l/(c3)**2.5 + (2*z)/(c3)**2.5 - 
                          l/(c4)**2.5 - (2*z)/(c4)**2.5))

def B0ppp(z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return (192*k*mu0*R**2*(-(c3)**-2.5 + (c4)**-2.5 + 
                   R**2*(5/(c3)**3.5 - 5/(c4)**3.5)))

def B0p4(z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 1920*k*mu0*R**2*(l*(-c3**-3.5 - c4**-3.5 +
                             7*(c3**-4.5 + c4**-4.5)*R**2) +
                          2*(c3**-3.5 - c4**-3.5 +
                             (-7/c3**4.5 + 7/c4**4.5)*R**2)*z)

def B0p5(z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 23040*k*mu0*R**2*(-c3**-3.5 + c4**-3.5 +
                           14*(c3**-4.5 - c4**-4.5)*R**2 +
                           42*(-c3**-5.5 + c4**-5.5)*R**4)

def B0p6(z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 322560*k*mu0*R**2*(l*(-c3**-4.5 - c4**-4.5 +
                               18*(c3**-5.5 + c4**-5.5)*R**2 +
                               (-66/c3**6.5 - 66/c4**6.5)*R**4) + 
                            2*(c3**-4.5 - c4**-4.5 +
                               18*(-c3**-5.5 + c4**-5.5)*R**2 +
                               66*(c3**-6.5 - c4**-6.5)*R**4)*z)

def B0p7(z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 5160960*k*mu0*R**2*(-c3**-4.5 + c4**-4.5 +
                             27*(c3**-5.5 - c4**-5.5)*R**2 +
                             198*(-c3**-6.5 + c4**-6.5)*R**4 +
                             429*(c3**-7.5 - c4**-7.5)*R**6)

def Bz(r,z,zcent,k,R,l):
  z = z - zcent
  c1 = l - 2*z
  c2 = l + 2*z
  c3 = 4*R**2 + c1**2
  c4 = 4*R**2 + c2**2
  return ((k*mu0*(c1/sqrt(c3) + c2/sqrt(c4)))/2.
         + 6*k*mu0*r**2*R**2*(2*z*(-(c3)**-2.5 + (c4)**-2.5) + 
                              l*((c3)**-2.5 + (c4)**-2.5)))
# return ((k*mu0*((l - 2*z)/sqrt(4*R**2 + (l - 2*z)**2) + 
#            (l + 2*z)/sqrt(4*R**2 + (l + 2*z)**2)))/2.\
#        + 6*k*mu0*r**2*R**2*
#        (2*z*(-(4*R**2 + (l - 2*z)**2)**-2.5 + 
#             (4*R**2 + (l + 2*z)**2)**-2.5) + 
#          l*((4*R**2 + (l - 2*z)**2)**-2.5 + 
#             (4*R**2 + (l + 2*z)**2)**-2.5)))
# return ((k*mu0*((l/2. - z)/sqrt(R**2 + (l/2. - z)**2) + 
#            (l/2. + z)/sqrt(R**2 + (l/2. + z)**2)))/2.\
#        - (k*mu0*r**2*
#          ((-2*(l/2. - z))/
#             (R**2 + (l/2. - z)**2)**1.5 + 
#            (-(R**2 + (l/2. - z)**2)**-1.5 + 
#               (3*(l/2. - z)**2)/
#                (R**2 + (l/2. - z)**2)**2.5)*(l/2. - z)\
#             - (2*(l/2. + z))/
#             (R**2 + (l/2. + z)**2)**1.5 + 
#            (l/2. + z)*
#             ((3*(l/2. + z)**2)/
#                (R**2 + (l/2. + z)**2)**2.5 - 
#               (R**2 + (l/2. + z)**2)**-1.5)))/8.)

def Br(r,z,zcent,k,R,l):
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return (2*k*mu0*r*R**2*((c3)**-1.5 - (c4)**-1.5 + 
                          6*r**2*(-(c3)**-2.5 + (c4)**-2.5 + 
                          R**2*(5/(c3)**3.5 - 5/(c4)**3.5))))
# return (-(k*mu0*r*(-(1/sqrt(R**2 + (l/2. - z)**2)) + 
#             (l/2. - z)**2/
#              (R**2 + (l/2. - z)**2)**1.5 - 
#             (l/2. + z)**2/
#              (R**2 + (l/2. + z)**2)**1.5 + 
#             1/sqrt(R**2 + (l/2. + z)**2)))/4. + 
#       (k*mu0*r**3*(3*
#             ((R**2 + (l/2. - z)**2)**-1.5 - 
#               (3*(l/2. - z)**2)/
#                (R**2 + (l/2. - z)**2)**2.5) + 
#            ((-9*(l/2. - z))/
#                (R**2 + (l/2. - z)**2)**2.5 + 
#               (15*(l/2. - z)**3)/
#                (R**2 + (l/2. - z)**2)**3.5)*(l/2. - z)\
#             + (l/2. + z)*
#             ((-15*(l/2. + z)**3)/
#                (R**2 + (l/2. + z)**2)**3.5 + 
#               (9*(l/2. + z))/
#                (R**2 + (l/2. + z)**2)**2.5) + 
#            3*((3*(l/2. + z)**2)/
#                (R**2 + (l/2. + z)**2)**2.5 - 
#               (R**2 + (l/2. + z)**2)**-1.5)))/32.)


def addsolenoid(zi,zf,ri,ro=None,maxbz=None,current=0.,
                nzpoints=10000,fringelen=10.,
                nsheets=1,
                **kw):
  """
Adds a solenoid element as an mmlt lattice element.
 - zi: z start of the current sheet
 - zf: z end of the current sheet
 - ri: inner radius of the sheet
 - ro=ri: outer radius of the sheet (note that only (ri+ro)/2 is actually used)
 - nsheets=1: number of current sheets
 - maxbz: maximum Bz field in T; used to calculate current if specified
 - current: current in the sheet, in units of Ampere-turns/meter; ignored for nonzero maxbz
 - nzpoints=10000: number of points in the table generated
 - fringelen=10.: length of region before and after the current sheet to
                  include the field fringe, in units of the sheet radius
Note that the actual sheet radius is given be (ri+ro)/2. The aperture is given
by ri. The fringelen uses the actual sheet radius.
  """
  if ro is None: ro = ri
  if nsheets == 1:
    rsheets = array([(ri + ro)/2.])
  else:
    rsheets = ri + arange(nsheets)*(ro-ri)/(nsheets-1)

  zcent = (zi + zf)/2.
  l = (zf - zi)
  
  if maxbz is not None:
    current = maxbz/(mu0*sum(1./sqrt(4.0*(rsheets/l)**2 + 1.0)))
  if current == 0.0:
    print 'warning: solenoid with zero current at zcent = %-6.3f' % zcent

  zs = zi - fringelen*max(rsheets)
  ze = zf + fringelen*max(rsheets)
  ap = kw.get('ap',ri)
  z = span(zs,ze,nzpoints+1)
  ms  = zeros((nzpoints+1,2),'d')
  msp = zeros((nzpoints+1,2),'d')

  for R in rsheets:
    ms[:,0]  += B0(z,zcent,current,R,l)
    msp[:,0] += B0p(z,zcent,current,R,l)
    ms[:,1]  += B0pp(z,zcent,current,R,l)
    msp[:,1] += B0ppp(z,zcent,current,R,l)
    # --- This is a slowly converging series for radius approaching R
    # --- so having extra terms doesn't help, and can be worse since the
    # --- terms get larger at first.
    #ms[:,2]  += B0p4(z,zcent,current,R,l)
    #msp[:,2] += B0p5(z,zcent,current,R,l)
    #ms[:,3]  += B0p6(z,zcent,current,R,l)
    #msp[:,3] += B0p7(z,zcent,current,R,l)

  nn = [0,0]
  vv = [0,1]
  return addnewmmlt(zs,ze,ap,ms=ms,msp=msp,nn=nn,vv=vv,**kw)

addnewsolenoid = addsolenoid

