"""
Contains routines for setting up plates with multiple apertures.
createapertureplate: routine to generate conductor data for field solver
ParticleScraper: class for creating particle scraping
"""
from warp import *

createapertureplate_version = "$Id: createapertureplate.py,v 1.4 2007/12/20 00:36:40 dave Exp $"
def createapertureplatedoc():
  import createapertureplate
  print createapertureplate.__doc__

def createapertureplate(xcent,ycent,zcent,rad,platewid,volt=0.,id=0,
                        platemax=None): 
  """
Generate conductor data for field solver
  - xcent: list of x centers of apertures
  - ycent: list of y centers of apertures
  - zcent: z center of plate
  - rad: list of radii of apertures
  - platewid: width of plate
  - volt=0.: voltage on the plate
  - id=0: conductor id of the plate
  - platemax=None: max radius of the plate, defaults to be big enought to
                   cover full transverse area of the mesh
  """
  f3d.gridmode = 1

  # --- Make sure there is enough space
  izplanesperplate = nint(platewid/w3d.dz)
  f3d.ncondmax = f3d.ncond + w3d.nx*w3d.ny*(izplanesperplate+1)
  f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy) + w3d.nx*3*4
  gchange('Conductor3d')

  if platemax is None:
    platemax = sqrt(w3d.xmmax**2 + w3d.ymmax**2) + max(w3d.dx,w3d.dy)

  iix = iota(0,w3d.nx)[:,NewAxis]*ones(1+w3d.ny,'l')
  iiy = iota(0,w3d.ny)*ones(1+w3d.nx,'l')[:,NewAxis]
  iix.shape = ((w3d.nx+1)*(w3d.ny+1),)
  iiy.shape = ((w3d.nx+1)*(w3d.ny+1),)
  xxx = iix*w3d.dx + w3d.xmmin
  yyy = iiy*w3d.dy + w3d.ymmin
  rrrsq = xxx**2 + yyy**2
  ix = compress(rrrsq<platemax**2,iix)
  iy = compress(rrrsq<platemax**2,iiy)

  is_start = f3d.ncond
  io_start = f3d.nocndbdy
  ie_start = f3d.necndbdy

  fuzz = 1.e-5
  # --- Set limits on points to include around each aperture
  ixmin = nint((xcent - rad - w3d.dx - w3d.xmmin)/w3d.dx)
  ixmax = nint((xcent + rad + w3d.dx - w3d.xmmin)/w3d.dx)
  iymin = nint((ycent - rad - w3d.dy - w3d.ymmin)/w3d.dy)
  iymax = nint((ycent + rad + w3d.dy - w3d.ymmin)/w3d.dy)
  ixmin = where(less(ixmin,0),0,ixmin)
  iymin = where(less(iymin,0),0,iymin)
  # --- Get location of plate center
  izl = int((zcent - platewid/2. - w3d.zmmin)/w3d.dz - fuzz)
  izr = w3d.nz - int((w3d.zmmax - zcent - platewid/2.)/w3d.dz - fuzz)
  # --- Loop over apertures
  f3d.lsrlinr = true
  nn = 1000
  f3d.npnts_sr = nn+1
  gchange("Surface_of_Rev")
  f3d.z_sr[:] = iota(-nn/2,nn/2)*platewid/nn + zcent
  for ia in xrange(len(xcent)):
    rmax = rad[ia] + 1.1*max(w3d.dx,w3d.dy)
    # --- Remove points that are near the aperture
    rr = ((ix*w3d.dx + w3d.xmmin - xcent[ia])**2 +
          (iy*w3d.dy + w3d.ymmin - ycent[ia])**2)
    ix = compress(greater(rr,rmax**2),ix)
    iy = compress(greater(rr,rmax**2),iy)
    # --- Make aperture
    f3d.r_sr[0] = rmax
    f3d.r_sr[1:-1] = rad[ia]
    f3d.r_sr[-1] = rmax
    srfrvout(" ",volt,zcent-platewid/2.,zcent+platewid/2.,xcent[ia],ycent[ia],
             rmax,true,w3d.xmmin,w3d.xmmax,w3d.ymmin,w3d.ymmax,false,
             w3d.zmmin,w3d.zmmax,top.zbeam,w3d.dx,w3d.dy,w3d.dz,
             w3d.nx,w3d.ny,w3d.nz,w3d.ix_axis,w3d.iy_axis,w3d.xmesh,w3d.ymesh,
             w3d.l2symtry,w3d.l4symtry)
  # --- Now, fill in the rest of the plate.
  izl = int((zcent - platewid/2. - w3d.zmmin)/w3d.dz + 1. - fuzz)
  wzl =    ((zcent - platewid/2. - w3d.zmmin)/w3d.dz + 1. - fuzz) - izl
  izr = int((zcent + platewid/2. - w3d.zmmin)/w3d.dz      + fuzz)
  wzr =    ((zcent + platewid/2. - w3d.zmmin)/w3d.dz      + fuzz) - izr
  if 0 <= izl <= w3d.nz:
    f3d.ixcond[f3d.ncond:f3d.ncond+len(ix)] = ix
    f3d.iycond[f3d.ncond:f3d.ncond+len(ix)] = iy
    f3d.izcond[f3d.ncond:f3d.ncond+len(ix)] = izl
    f3d.condvolt[f3d.ncond:f3d.ncond+len(ix)] = volt
    f3d.ncond = f3d.ncond + len(ix)
    if wzl > fuzz: addzsubgridpoints(ix,iy,izl-1,volt,2.,wzl)
  if 0 <= izr <= w3d.nz:
    f3d.ixcond[f3d.ncond:f3d.ncond+len(ix)] = ix
    f3d.iycond[f3d.ncond:f3d.ncond+len(ix)] = iy
    f3d.izcond[f3d.ncond:f3d.ncond+len(ix)] = izr
    f3d.condvolt[f3d.ncond:f3d.ncond+len(ix)] = volt
    f3d.ncond = f3d.ncond + len(ix)
    if wzr > fuzz: addzsubgridpoints(ix,iy,izr+1,volt,1.-wzl,2.)

  is_end = f3d.ncond
  io_end = f3d.nocndbdy
  ie_end = f3d.necndbdy
  if is_end > is_start: f3d.condnumb[is_start:is_end] = id
  if ie_end > ie_start: f3d.ecnumb[ie_start:ie_end] = id
  if io_end > io_start: f3d.ocnumb[io_start:io_end] = id

########################################################
########################################################
########################################################
def addzsubgridpoints(ix,iy,iz,volt,mz,pz):
  ieven = compress((ix+iy+iz)%2 == 0,arange(len(ix)))
  iodd  = compress((ix+iy+iz)%2 == 1,arange(len(ix)))
  if max(f3d.necndbdy+len(ieven),f3d.nocndbdy+len(iodd)) > f3d.ncndmax:
    f3d.ncndmax = max(f3d.necndbdy+len(ieven),f3d.nocndbdy+len(iodd))
    gchange("Conductor3d")
  f3d.iecndx[f3d.necndbdy:f3d.necndbdy+len(ieven)] = take(ix,ieven)
  f3d.iecndy[f3d.necndbdy:f3d.necndbdy+len(ieven)] = take(iy,ieven)
  f3d.iecndz[f3d.necndbdy:f3d.necndbdy+len(ieven)] = iz
  f3d.ecdelmx[f3d.necndbdy:f3d.necndbdy+len(ieven)] = 2.
  f3d.ecdelmy[f3d.necndbdy:f3d.necndbdy+len(ieven)] = 2.
  f3d.ecdelmz[f3d.necndbdy:f3d.necndbdy+len(ieven)] = mz
  f3d.ecdelpx[f3d.necndbdy:f3d.necndbdy+len(ieven)] = 2.
  f3d.ecdelpy[f3d.necndbdy:f3d.necndbdy+len(ieven)] = 2.
  f3d.ecdelpz[f3d.necndbdy:f3d.necndbdy+len(ieven)] = pz
  f3d.ecvolt[f3d.necndbdy:f3d.necndbdy+len(ieven)] = volt
  f3d.necndbdy = f3d.necndbdy + len(ieven)
  f3d.iocndx[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = take(ix,iodd)
  f3d.iocndy[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = take(iy,iodd)
  f3d.iocndz[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = iz
  f3d.ocdelmx[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = 2.
  f3d.ocdelmy[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = 2.
  f3d.ocdelmz[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = mz
  f3d.ocdelpx[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = 2.
  f3d.ocdelpy[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = 2.
  f3d.ocdelpz[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = pz
  f3d.ocvolt[f3d.nocndbdy:f3d.nocndbdy+len(iodd)] = volt
  f3d.nocndbdy = f3d.nocndbdy + len(iodd)

########################################################
########################################################
########################################################
# --- Create the particle scraper using griddedparticlescraper.
class ParticleScraperOld:
  """
Class for creating particle scraper for a plate with multiple apertures
  - xcent: list of x centers of apertures
  - ycent: list of y centers of apertures
  - zcent: z center of plate
  - rad: list of radii of apertures
  - platewid: width of plate
  - platemax=None: max radius of the plate, defaults to be big enought to
                   cover full transverse area of the mesh
  """
  def __init__(s,xcent,ycent,zcent,rad,platewid,platemax=None): 
    # --- Find z extent
    s.zmin = int((zcent - platewid/2. - w3d.zmmin)/w3d.dz-1)*w3d.dz + w3d.zmmin
    s.izmin = nint((s.zmin - w3d.zmmin)/w3d.dz)
    s.zmax = int((zcent + platewid/2. - w3d.zmmin)/w3d.dz+2)*w3d.dz + w3d.zmmin
    s.izmax = nint((s.zmax - w3d.zmmin)/w3d.dz)

    # --- Set default value of platemax
    if platemax is None:
      platemax = sqrt(w3d.xmmax**2 + w3d.ymmax**2) + max(w3d.dx,w3d.dy)

    # --- Set mesh size
    s.nx = w3d.nx
    s.ny = w3d.ny
    s.nz = s.izmax - s.izmin

    # --- Allocate array to hold the distances to the conductor.
    # --- Give a default value as if all points are far from a conductor
    maxdd = 2.*min(w3d.dx,w3d.dy,w3d.dz)
    s.distances = fzeros((s.nx+1,s.ny+1,s.nz+1),'d') + maxdd
    ddtemp = fzeros((s.nx+1,s.ny+1),'d') + maxdd

    # --- Get locations of the grid points
    iix = iota(0,w3d.nx)[:,NewAxis]*ones(1+w3d.ny,'l')
    iiy = iota(0,w3d.ny)*ones(1+w3d.nx,'l')[:,NewAxis]
    iix.shape = ((w3d.nx+1)*(w3d.ny+1),)
    iiy.shape = ((w3d.nx+1)*(w3d.ny+1),)
    xxx = iix*w3d.dx + w3d.xmmin
    yyy = iiy*w3d.dy + w3d.ymmin
    rrrsq = xxx**2 + yyy**2
    xx = compress(rrrsq<platemax**2,xxx)
    yy = compress(rrrsq<platemax**2,yyy)
    ix = compress(rrrsq<platemax**2,iix)
    iy = compress(rrrsq<platemax**2,iiy)
    ddlist = zeros(len(ix),'d')

    # --- For single transverse plane, for each point, find greatest distance
    # --- transversely from a conductor. Transfer the data to ddtemp.
    for ii in xrange(len(ix)):
      dd = rad - sqrt((xx[ii] - xcent)**2 + (yy[ii] - ycent)**2)
      ddlist[ii] = max(dd)
    for (iix,iiy,dd) in map(None,ix,iy,ddlist): ddtemp[iix,iiy] = dd

    # --- For each z plane, copy in the data, compared to the distances
    # --- to the conductor in z.
    for iz in iota(s.izmin,s.izmax):
       zz = iz*w3d.dz + w3d.zmmin
       dd = max(zcent - platewid/2. - zz,zz - zcent - platewid/2.)
       if dd >= 0.:
         s.distances[:,:,iz-s.izmin] = maximum(ddtemp,dd)
       elif dd < 0.:
         s.distances[:,:,iz-s.izmin] = where(greater(ddtemp,-w3d.dx),ddtemp,dd)

    # --- Install the scraper to be called before a time step
    installbeforestep(s.doscrape)

  def ppzx(s,iy=None):
    if iy is None: iy = w3d.iy_axis
    gg = gatherarray(transpose(s.distances[:,iy,:]))
    ppgeneric(grid=gg,cellarray=1,
              xmin=s.zmin,ymin=w3d.xmmin,xmax=s.zmax,ymax=w3d.xmmax)
  def ppzy(s,ix=None):
    if ix is None: ix = w3d.ix_axis
    gg = gatherarray(transpose(s.distances[ix,:,:]))
    ppgeneric(grid=gg,cellarray=1,
              xmin=s.zmin,ymin=w3d.ymmin,xmax=s.zmax,ymax=w3d.ymmax)
  def ppxy(s,iz=None):
    if iz is None: iz = w3d.iz_axis
    gg = gatherarray(s.distances[:,:,iz])
    ppgeneric(grid=gg,cellarray=1,
              xmin=w3d.xmmin,ymin=w3d.ymmin,xmax=w3d.xmmax,ymax=w3d.ymmax)

  # --- Create the particle scraper and install it
  def doscrape(s):
    for js in xrange(top.ns):
      griddedparticlescraper(js+1,s.distances,s.nx,s.ny,s.nz,
                             w3d.dx,w3d.dy,w3d.dz,w3d.xmmin,w3d.ymmin,s.zmin,
                             top.zbeam,w3d.l2symtry,w3d.l4symtry)

