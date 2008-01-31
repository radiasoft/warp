"""Functions to generate solenoid lattice elements.

The following functions are available:
 - addsolenoid: multipole expansion
 - addgriddedsolenoid: gridded data via solve of del**2 A = mu0*J

"""
from warp import *
from lattice import addnewmmlt,addnewbgrd
solenoid_version = "$Id: solenoid.py,v 1.8 2008/01/31 01:09:42 dave Exp $"

def solenoiddoc():
  import solenoid
  print solenoid.__doc__

# --- Functions for the multipole representation
def B0(z,zcent,k,R,l):
  "Bz on axis"
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return (k*mu0*((l - 2*z)/sqrt(c3) + (l + 2*z)/sqrt(c4)))/2.

def B0p(z,zcent,k,R,l):
  "First z derivative of Bz on axis"
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 4*k*mu0*R**2*(-(c3)**-1.5 + (c4)**-1.5)

def B0pp(z,zcent,k,R,l):
  "Second z derivative of Bz on axis"
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return (24*k*mu0*R**2*(-l/(c3)**2.5 + (2*z)/(c3)**2.5 - 
                          l/(c4)**2.5 - (2*z)/(c4)**2.5))

def B0ppp(z,zcent,k,R,l):
  "Third z derivative of Bz on axis"
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return (192*k*mu0*R**2*(-(c3)**-2.5 + (c4)**-2.5 + 
                   R**2*(5/(c3)**3.5 - 5/(c4)**3.5)))

def B0p4(z,zcent,k,R,l):
  "Fourth z derivative of Bz on axis"
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 1920*k*mu0*R**2*(l*(-c3**-3.5 - c4**-3.5 +
                             7*(c3**-4.5 + c4**-4.5)*R**2) +
                          2*(c3**-3.5 - c4**-3.5 +
                             (-7/c3**4.5 + 7/c4**4.5)*R**2)*z)

def B0p5(z,zcent,k,R,l):
  "Fifth z derivative of Bz on axis"
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 23040*k*mu0*R**2*(-c3**-3.5 + c4**-3.5 +
                           14*(c3**-4.5 - c4**-4.5)*R**2 +
                           42*(-c3**-5.5 + c4**-5.5)*R**4)

def B0p6(z,zcent,k,R,l):
  "Sixth z derivative of Bz on axis"
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
  "Seventh z derivative of Bz on axis"
  z = z - zcent
  c3 = 4*R**2 + (l - 2*z)**2
  c4 = 4*R**2 + (l + 2*z)**2
  return 5160960*k*mu0*R**2*(-c3**-4.5 + c4**-4.5 +
                             27*(c3**-5.5 - c4**-5.5)*R**2 +
                             198*(-c3**-6.5 + c4**-6.5)*R**4 +
                             429*(c3**-7.5 - c4**-7.5)*R**6)

def addsolenoid(zi,zf,ri,ro=None,maxbz=None,current=0.,
                nzpoints=10000,fringelen=10.,
                nsheets=1,v=1,
                **kw):
  """
Adds a solenoid element represented as a multipole expansion of the field on
axis. This creates a mmlt lattice element.
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
 - v=1: number of non-linear terms to include. max value is 3, though
        it is not recommended to use v>1.
Note that the actual sheet radius is given be (ri+ro)/2. The aperture is given
by ri. The fringelen uses the actual sheet radius.

The solenoid field is obtained from the analytic field profile of a cylindrical
current sheet. The field on axis is given by
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
  z = zs + (ze - zs)*arange(nzpoints+1)/nzpoints
  ms  = zeros((nzpoints+1,v+1),'d')
  msp = zeros((nzpoints+1,v+1),'d')

  for R in rsheets:
    ms[:,0]  += B0(z,zcent,current,R,l)
    msp[:,0] += B0p(z,zcent,current,R,l)
    if v >= 1:
      ms[:,1]  += B0pp(z,zcent,current,R,l)
      msp[:,1] += B0ppp(z,zcent,current,R,l)
    # --- This is a slowly converging series for radius approaching R
    # --- so having extra terms doesn't help, and can be worse since the
    # --- terms get larger at first.
    if v >= 2:
      ms[:,2]  += B0p4(z,zcent,current,R,l)
      msp[:,2] += B0p5(z,zcent,current,R,l)
    if v >= 3:
      ms[:,3]  += B0p6(z,zcent,current,R,l)
      msp[:,3] += B0p7(z,zcent,current,R,l)

  nn = zeros(v+1,'l')
  vv = arange(v+1,dtype='l')
  return addnewmmlt(zs,ze,ap,ms=ms,msp=msp,nn=nn,vv=vv,**kw)

addnewsolenoid = addsolenoid


#============================================================================
#============================================================================
#============================================================================
def addgriddedsolenoid(zcenter=None,length=None,rinner=None,router=None,
                       bzmax=None,
                       lcylindrical=true,
                       saveBsolver=false,scalebz=true,tol=1.e-5,
                       nx=None,ny=None,nz=None,
                       dx=None,dy=None,dz=None,
                       xmmin=None,xmmax=None,
                       ymmin=None,ymmax=None,
                       zmmin=None,zmmax=None,
                       l4symtry=None,l2symtry=None,
                       **kw):
  """
Adds a solenoid using gridded data, via a bgrd element.
The B field is calculated by carrying out a solve of del**2 A = mu0*J,
with the solenoid windings uniformly filled with current.
Input arguments:
 - zcenter: z center of the solenoid windings
 - length: z length of the solenoid windings
 - rinner,router: inner and outer radius of the windings
 - bzmax: Bz field on axis at the z center
 - lcylindrical=true:
 - nx,ny,nz,dx,dy,dz,xmmin,xmmax,ymmin,ymmax,zmmin,zmmax: all default from w3d
 - saveBsolver=false: when true, the field solver used to calculated
                      the fields is save as an attribute of the function.
 - scalebz=true: when true, B is scale to exactly match bzmax
 - tol=1.e-5: convergence tolerence for the field solve relative to bzmax
  """

  assert zcenter is not None,ValueError('zcenter must be given')
  assert length is not None,ValueError('length must be given')
  assert rinner is not None,ValueError('rinner must be given')
  assert router is not None,ValueError('router must be given')
  assert bzmax is not None,ValueError('bzmax must be given')

  # --- Note that the grid parameters are passed in using a separate dict
  # --- instead of using the generic kw, since kw is expected to contain
  # --- options that addnewbgrd will know about, but the field solver not know.
  # --- The values are only put in the dictionary if they are not None
  # --- since the solver only checks for the presence of the options and not
  # --- their value.
  solverdict={}
  if nx is not None: solverdict['nx'] = nx
  if ny is not None: solverdict['ny'] = ny
  if nz is not None: solverdict['nz'] = nz
  if dx is not None: solverdict['dx'] = dx
  if dy is not None: solverdict['dy'] = dy
  if dz is not None: solverdict['dz'] = dz
  if xmmin is not None: solverdict['xmmin'] = xmmin
  if xmmax is not None: solverdict['xmmax'] = xmmax
  if ymmin is not None: solverdict['ymmin'] = ymmin
  if ymmax is not None: solverdict['ymmax'] = ymmax
  if zmmin is not None: solverdict['zmmin'] = zmmin
  if zmmax is not None: solverdict['zmmax'] = zmmax
  solverdict['nslaves'] = 0
  if lcylindrical:
    solverdict['ny'] = 0
    Bsolver = MagnetostaticMG(lcylindrical=1,**solverdict)
    xx,zz = getmesh2d(0.,Bsolver.dx,Bsolver.nx,
                      Bsolver.zmmin,Bsolver.dz,Bsolver.nz)
    yy = zeros(xx.shape,'d')
  else:
    if l4symtry is not None: solverdict['l4symtry'] = l4symtry
    if l2symtry is not None: solverdict['l2symtry'] = l2symtry
    Bsolver = MagnetostaticMG(**solverdict)
    xx,yy,zz = getmesh3d(Bsolver.xmmin,Bsolver.dx,Bsolver.nx,
                         Bsolver.ymmin,Bsolver.dy,Bsolver.ny,
                         Bsolver.zmmin,Bsolver.dz,Bsolver.nz)

  if saveBsolver: addgriddedsolenoid.Bsolver = Bsolver

  # --- Get the min and max in Quadrant I
  xmin = sqrt(maximum(rinner**2 - yy**2,0.))
  xmax = sqrt(maximum(router**2 - yy**2,0.))

  ymin = sqrt(maximum(rinner**2 - xx**2,0.))
  ymax = sqrt(maximum(router**2 - xx**2,0.))

  zmin = zcenter - length/2.
  zmax = zcenter + length/2.

  if not lcylindrical:
    # --- Fix the mins on axis
    xmin = where(logical_and(rinner<yy,yy<router),-xmax,xmin)
    ymin = where(logical_and(rinner<xx,xx<router),-ymax,ymin)

    # --- Copy the min and max to the neighboring quadrante where the values
    # --- are the same.
    xmin = where(yy<0.,xmin[:,::-1,:],xmin)
    xmax = where(yy<0.,xmax[:,::-1,:],xmax)
    ymin = where(xx<0.,ymin[::-1,...],ymin)
    ymax = where(xx<0.,ymax[::-1,...],ymax)

    # --- For the inverted quadrants, the negative max becomes the min, etc
    xmin = where(xx<0.,-xmax[::-1,...],xmin)
    xmax = where(xx<0.,-xmin[::-1,...],xmax)
    ymin = where(yy<0.,-ymax[:,::-1,:],ymin)
    ymax = where(yy<0.,-ymin[:,::-1,:],ymax)

  dx = Bsolver.dx
  dy = Bsolver.dy
  dz = Bsolver.dz

  # --- Integrate over a uniform current density, using a linear weight
  # --- function for each current infinitesimal. This is equivalent to doing
  # --- linear weighting with an infinite number of particles uniformly
  # --- filling the area of the windings.
  # --- Write the weight function as
  # --- Sx = 1/dx * | 1 - (xi - x)/dx  for x <= xi <= x+dx
  # ---             | (xi - x)/dx + 1  for x-dx <= xi <= x
  # ---             | 0                otherwise
  # --- Then for each xi, the contribution is multiplied by
  # --- wx = Integral from xl to xu of Sx over x
  # --- where xl is the max of xi-dx and the lower bound of the windings
  # ---       xu is the min of xi+dx and the upper bound of the windings
  # --- In the code below, this integral is divided into two parts,
  # --- plus and minus, so wx = dxm + dxp.
  # --- The final weighting is then wx*wy*wz.

  # --- Note that in the 3D case, the loading of the current is only
  # --- approximate. The integrals over the current density use a stair-step
  # --- approximation of the windings surface.

  z2 = minimum(zmax,zz+dz)
  z1 = maximum(zmin,zz)
  dzp = (+(z2*zz/dz - 0.5*z2**2/dz + z2)
         -(z1*zz/dz - 0.5*z1**2/dz + z1))/dz
  dzp = where(z2 >= z1,dzp,0.)

  z2 = minimum(zmax,zz)
  z1 = maximum(zmin,zz-dz)
  dzm = (+(z2 - z2*zz/dz + 0.5*z2**2/dz)
         -(z1 - z1*zz/dz + 0.5*z1**2/dz))/dz
  dzm = where(z2 >= z1,dzm,0.)

  x2 = minimum(xmax,xx+dx)
  x1 = maximum(xmin,xx)
  dxp = (+(x2*xx/dx - 0.5*x2**2/dx + x2)
         -(x1*xx/dx - 0.5*x1**2/dx + x1))/dx
  dxp = where(x2 >= x1,dxp,0.)

  x2 = minimum(xmax,xx)
  x1 = maximum(xmin,xx-dx)
  dxm = (+(x2 - x2*xx/dx + 0.5*x2**2/dx)
         -(x1 - x1*xx/dx + 0.5*x1**2/dx))/dx
  dxm = where(x2 >= x1,dxm,0.)

  if lcylindrical:
    # --- In this case, y is really needed, and wy = 1
    dym = 0.5
    dyp = 0.5
  else:
    y2 = minimum(ymax,yy+dy)
    y1 = maximum(ymin,yy)
    dyp = (+(y2*yy/dy - 0.5*y2**2/dy + y2)
           -(y1*yy/dy - 0.5*y1**2/dy + y1))/dy
    dyp = where(y2 >= y1,dyp,0.)

    y2 = minimum(ymax,yy)
    y1 = maximum(ymin,yy-dy)
    dym = (+(y2 - y2*yy/dy + 0.5*y2**2/dy)
           -(y1 - y1*yy/dy + 0.5*y1**2/dy))/dy
    dym = where(y2 >= y1,dym,0.)

  ww = (dxp+dxm)*(dym+dyp)*(dzp+dzm)

  # --- Setup the source arrays for the solver
  Bsolver.setparticledomains()
  Bsolver.allocatedataarrays()
  Bsolver.zerosourcep()
  Bsolver.setsourcepforparticles(0,0,0)

  # --- Find the current given the bzmax and the area of the windings.
  # --- This assumes a large number of windings all with the same current.
  # --- The current in each winding is calculated, multiplied by the number
  # --- of windings to get the total current per meter, and divied by the
  # --- radial extent of the windings to get a current density.
  rwinds = rinner + arange(1001)/1000.*(router-rinner)
  current = (bzmax/(mu0*sum(1./sqrt(4.*(rwinds/length)**2+1.)))*
             len(rwinds)/(router-rinner))

  # --- Put the current into the solver
  if lcylindrical:
    Bsolver.sourcep[1,:,0,:] = +ww*current
  else:
    theta = arctan2(yy,xx)
    Bsolver.sourcep[0,:,:,:] = -ww*current*sin(theta)
    Bsolver.sourcep[1,:,:,:] = +ww*current*cos(theta)

  # --- Set convergence tolerence (in Tesla)
  Bsolver.mgtol = bzmax*tol*ones(3)

  # --- Force some parameters for the RZ field solver, but save their
  # --- value so they can be restored after the solve.
  solvergeom = w3d.solvergeom
  ncmax = frz.mgridrz_ncmax
  w3d.solvergeom = w3d.RZgeom
  frz.mgridrz_ncmax = 10000

  Bsolver.solve()

  w3d.solvergeom = solvergeom
  frz.mgridrz_ncmax = ncmax

  '''
  # --- Do some error checking for the RZ case
  bb = Bsolver.field
  divb = fzeros(rr.shape,'d')
  divb[1:-1,1:-1] = (1./rr[1:-1,1:-1]*bb[0,1:-1,0,1:-1] +
                       (bb[0,2:,0,1:-1] - bb[0,:-2,0,1:-1])/(2.*dx) + 
                     (bb[2,1:-1,0,2:]-bb[2,1:-1,0,:-2])/(2.*dz))
  
  curlbt = fzeros(rr.shape,'d')
  curlbt[1:-1,1:-1] = (+(bb[0,1:-1,0,2:]-bb[0,1:-1,0,:-2])/(2.*dz)
                       -(bb[2,2:,0,1:-1] - bb[2,:-2,0,1:-1])/(2.*dx))
  '''

  # --- Scale the B field to get exactly bzmax.
  if scalebz:
    bzmax_actual = max(Bsolver.field[2,0,0,:])
    Bsolver.source[...] = bzmax_actual/bzmax*Bsolver.source
    Bsolver.field[...] = bzmax_actual/bzmax*Bsolver.field

  # --- If ap is not passed in, then use the inner radius of the
  # --- solenoid windinds (which isn't necessarily a good value).
  kw.setdefault('ap',rinner)

  # --- Now add in the solenoid element
  return addnewbgrd(Bsolver.zmmin,Bsolver.zmmax,
                    dx=Bsolver.dx,dy=Bsolver.dy,
                    bx=Bsolver.field[0,...],
                    by=Bsolver.field[1,...],
                    bz=Bsolver.field[2,...],
                    rz=lcylindrical,**kw)

