from generateconductors import *
from particlescraper import *

class GeneratePlates:
  """
The following need to be defined:

w3d.dx,dy,dz,nx,ny,nz

top.xinject,top.yinject
top.vinject

These have default values if undefined:
platewid = 1.e-3
rround = 0.
pierceangle = 67.5
piercezmax = (beamletsep/2. - top.a0)*tan((90.-pierceangle)*pi/180.)
a0full = 1.
ap0full = 0.
b0full = 1.
bp0full = 0.
platesize = max(max(top.xinject),max(top.yinject)) + 0.004
rplate = 2.e-3
lsurround = false
lcurvaturemockup = false
lverbose = 1

These must be defined and be a sequence:
vplate,zplate
  """
  def __init__(self,beamletsep = 0.,
                    accllen = 0.,
                    platewid = 1.e-3,
                    rround = 0.,
                    pierceangle = 67.5,
                    piercezmax = None,
                    a0full = 1.,
                    ap0full = 0.,
                    b0full = 1.,
                    bp0full = 0.,
                    platesize = None,
                    lsurround = false,
                    lcurvaturemockup = false,
                    lverbose = 1,
                    vplate = None,
                    zplate = None,
                    rplate = None,
                    columnlen=None,
                    reset=0,
                    lscraper = 0):

    assert vplate is not None,"vplate must be specified"
    assert zplate is not None,"zplate must be specified"

    if piercezmax is None:
      piercezmax = (beamletsep/2. - top.a0)*tan((90.-pierceangle)*pi/180.)
    if platesize is None:
      platesize = max(max(top.xinject),max(top.yinject)) + 0.004
    if rplate is None:
      rplate = len(zplate)*[2.e-3]
    if columnlen is None:
      columnlen = zplate[-1]

    self.beamletsep = beamletsep
    self.accllen = accllen

    self.platewid = platewid
    self.rround = rround
    self.pierceangle = pierceangle
    self.piercezmax = piercezmax
    self.a0full = a0full
    self.ap0full = ap0full
    self.b0full = b0full
    self.bp0full = bp0full
    self.platesize = platesize
    self.lsurround = lsurround
    self.lcurvaturemockup = lcurvaturemockup
    self.lverbose = lverbose
    self.vplate = vplate
    self.zplate = zplate
    self.rplate = rplate
    self.lscraper = lscraper

    self.inittemps()

    # --- set so conductor is not recalculated
    f3d.gridmode = 1

    if reset:
      # --- delete any previous data
      f3d.ncond = 0
      f3d.necndbdy = 0
      f3d.nocndbdy = 0
      if w3d.solvergeom == w3d.RZgeom:
        frz.del_base()
        frz.init_base(w3d.nx,w3d.nz,w3d.dx,w3d.dz,w3d.xmmin,w3d.zmmin,false)

    self.createconductors()

    self.createparticlescraper()

    self.setupcurvaturemockup()

  def inittemps(s):
    # --- Preallocate the conductor arrays, making very rough guesses
    s.zleft  = array([0.]         + list(s.zplate - s.platewid/2.))
    s.zright = array([s.piercezmax] + list(s.zplate + s.platewid/2.))
    ii = compress((w3d.zmmin <= s.zleft ) & (s.zleft  <= w3d.zmmax) |
                  (w3d.zmmin <= s.zright) & (s.zright <= w3d.zmmax),
                  arange(len(s.zleft)))
    s.izleft = nint(floor(take(s.zleft,ii)/w3d.dz))
    s.izright = nint(ceil(take(s.zright,ii)/w3d.dz))

    zz = 1.*sum(s.izright - s.izleft)/w3d.nz
    f3d.ncondmax = max(f3d.ncondmax,zz*(w3d.nx+1)*(w3d.ny+1)*w3d.nz)
    f3d.ncndmax = max(f3d.ncndmax,zz*(w3d.nx+1)*(w3d.ny+1)*w3d.nz)
    gchange("Conductor3d")

    # --- Define usefull temporaries
    s.amaxinject = s.a0full - s.ap0full*s.zplate[-1]
    s.bmaxinject = s.b0full - s.bp0full*s.zplate[-1]
    if s.ap0full == 0.: s.za = largepos
    else:               s.za = s.zplate[-1] + s.a0full/(-s.ap0full)
    if s.bp0full == 0.: s.zb = largepos
    else:               s.zb = s.zplate[-1] + s.b0full/(-s.bp0full)

  # -----------------------------------------------------------------------
  # --- These routines are used to define the plate shape.
  def getplatez(self,zz,xx,yy,za,zb):
    if za > 1.e10 or zb > 1.e10:
      zl1 = zz
      zl2 = zz
    else:
      zl1 = zb - sqrt((zb-za+sqrt((za-zz)**2-xx**2))**2 - yy**2)
      zl2 = za - sqrt((za-zb+sqrt((zb-zz)**2-yy**2))**2 - xx**2)
    zlave = 0.5*(zl1 + zl2)
    return zlave

  def findplateintersect(self,z0,x0,xp0,y0,yp0,zz,za,zb):
    xx = x0 + xp0*(zz - z0)
    yy = y0 + yp0*(zz - z0)
    z1 = zz
    z2 = self.getplatez(zz,xx,yy,za,zb)
    while abs(z1-z2) > 1.e-8:
      xx = x0 + xp0*(z2 - z0)
      yy = y0 + yp0*(z2 - z0)
      z1 = z2
      z2 = self.getplatez(zz,xx,yy,za,zb)
    return (xx,yy,z2)
    
# -----------------------------------------------------------------------
  def createconductors(s):
    condgenstarttime = wtime()
    s.createsource()
    s.createapertureplates()
    s.createouterwall()

    if(w3d.solvergeom == w3d.RZgeom or w3d.solvergeom == w3d.XZgeom):
      frz.install_conductors_rz()

    f3d.ncondmax = f3d.ncond
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
    gchange("Conductor3d")

    if 1:
      # --- Set starting guess of phi using izvolt
      for iz in xrange(len(s.izvolt)):
        if top.izfsslave[me] <= iz <= top.izfsslave[me]+w3d.nz:
          w3d.phi[:,:,iz-top.izfsslave[me]+1] = s.izvolt[iz]
          if w3d.solvergeom == w3d.RZgeom:
            ng = frz.basegrid.nguardz
            frz.basegrid.phi[:,iz+ng] = s.izvolt[iz]

    if s.lverbose:
      print "setupconductors done"

    condgenendtime = wtime()
    s.condgentime = condgenendtime - condgenstarttime
    if s.lverbose:
      print "Conductor generate time = ",condgenendtime - condgenstarttime

    #if w3d.solvergeom == w3d.RZgeom: get_cond_rz(1,1)

# -----------------------------------------------------------------------
  def createsource(s):
    if s.lverbose:
      print "Create source plate"

    # --- Find max z of source plate
    xmax = s.amaxinject + s.beamletsep
    ymax = s.bmaxinject + s.beamletsep
    xpmax = xmax/s.amaxinject*s.ap0full
    ypmax = ymax/s.bmaxinject*s.bp0full
    xx,yx,zx = s.findplateintersect(s.piercezmax,xmax,xpmax,0.,0.,
                                    s.piercezmax,s.za,s.zb)
    xy,yy,zy = s.findplateintersect(s.piercezmax,0.,0.,ymax,ypmax,
                                    s.piercezmax,s.za,s.zb)
    s.sourceplatezmax = max(zx,zy)

    # --- Get list of source locations. The location obtained is the centers of
    # --- the Pierce cones.
    conezmax = s.sourceplatezmax + 1.e-2
    xs = zeros(top.ninject,'d')
    ys = zeros(top.ninject,'d')
    zs = zeros(top.ninject,'d')
    z0 = 0.
    zz = conezmax/2.
    for ij in xrange(top.ninject):
      x0 = top.xinject[ij]
      y0 = top.yinject[ij]
      xp0 = top.xpinject[ij]
      yp0 = top.xpinject[ij]
      xs[ij],ys[ij],zs[ij] = s.findplateintersect(z0,x0,xp0,y0,yp0,zz,s.za,s.zb)

    # --- Create a cylinder which covers all sources
    # --- Make plate very wide and locate it so that the front is at piercezmax
    sourcewidth = 1.e-1
    z0 = -sourcewidth/2. + s.piercezmax
    s.sourceplate = Beamletplate(s.za,s.zb,z0,sourcewidth,top.vinject[0],
                                 condid=1)

    # --- Create a cone for each emitter
    r1 = top.ainject
    r2 = top.ainject + conezmax*tan(s.pierceangle*pi/180.)
    cz = conezmax*ones(top.ninject,'d')
    theta = s.ap0full*top.xinject/s.amaxinject
    phi = s.bp0full*top.yinject/s.bmaxinject
    ss = Cones(r1,r2,cz,theta,phi,top.vinject[0],xs,ys,zs,condid=1)
    s.sourceplate = s.sourceplate - ss
    
    # --- Make a grid just big enough for the pierce cones
    # --- and generate conductor data
    sourcegrid = Grid(zmax=s.sourceplatezmax)
    sourcegrid.getdata(s.sourceplate,dfill=200.)
    sourcegrid.installdata(installrz=0)

    # --- Reset the locations of the sources
    for ij in xrange(top.ninject):
      x0 = top.xinject[ij]
      y0 = top.yinject[ij]
      xp0 = top.xpinject[ij]
      yp0 = top.xpinject[ij]
      xs[ij],ys[ij],zs[ij] = s.findplateintersect(0.,x0,xp0,y0,yp0,0.,s.za,s.zb)
    top.xinject[:] = xs
    top.yinject[:] = ys
    top.zinject[:] = zs

###########################################################################
  def createapertureplates(s):
    if s.lverbose:
      print "Setting up Pierce column"

    s.plates = []
    for ip in xrange(len(s.zplate)):
      s.createapertureplate(ip)

  def createapertureplate(s,ip):
    # --- Follow the same procedure for each aperture plate
    if s.za > 1.e10 or s.zb > 1.e10:
      plate = ZRoundedCylinder(s.platesize,s.platewid,s.rround,s.vplate[ip],
                               0.,0.,s.zplate[ip],condid=ip+2)
    else:
      z0 = s.zplate[ip]
      plate = Beamletplate(s.za,s.zb,z0,s.platewid,s.vplate[ip],condid=ip+2)
  
    # --- Get location of apertures.
    xx = top.xinject + (s.zplate[ip]-top.zinject)*top.xpinject
    yy = top.yinject + (s.zplate[ip]-top.zinject)*top.ypinject

    # --- Subtract out aperture holes
    if s.rround > 0.:
      for ij in xrange(top.ninject):
        ap = ZRoundedCylinderOut(s.rplate[ip],s.platewid,s.rround,s.vplate[ip],
                                 xx[ij],yy[ij],s.zplate[ip],condid=ip+2)
        plate = plate * ap
    else:
      rr = ones(top.ninject)*s.rplate[ip]
      ll = ones(top.ninject)*largepos
      zz = ones(top.ninject)*s.zplate[ip]
      ap = Cylinders(rr,ll,top.xpinject,top.ypinject,
                     s.vplate[ip],xx,yy,zz,condid=ip+2)
      plate = plate - ap

    s.plates.append(plate)

    # --- Make a grid just big enough for the plate
    # --- and generate conductor data
    plategrid = Grid(zmin=s.zplate[ip]-1.e-2,zmax=s.zplate[ip]+1.e-2)
    plategrid.getdata(plate,dfill=200.)
    plategrid.installdata(installrz=0)
    if s.lverbose:
      print "finished plate number ",ip

  def getizvolt(s,vinject=top.vinject[0]):
    izvolt = zeros(1+w3d.nzfull,'d')
    # --- Set voltages in diode
    zl = s.piercezmax/w3d.dz
    zr = s.zleft[1]/w3d.dz
    izl = int(zl)
    izr = int(zr)
    izvolt[:izl+1] = vinject
    izvolt[izl:izr+1] = vinject+(s.vplate[0]-vinject)*(iota(izl,izr)-zl)/(zr-zl)
    # --- Set voltages between plates
    for ip in xrange(len(s.zplate)-1):
      zl = s.zright[ip+1]/w3d.dz
      zr = s.zleft[ip+2]/w3d.dz
      for iz in iota(nint(zl),nint(zr)):
        if 0 <= iz <= w3d.nzfull:
          izvolt[iz] = s.vplate[ip]+(s.vplate[ip+1]-s.vplate[ip])*(iz - zl)/(zr - zl)
    # --- Set voltages within plates
    for ip in xrange(len(s.zplate)):
      zl = s.zleft[ip+1]/w3d.dz
      zr = s.zright[ip+1]/w3d.dz
      for iz in iota(nint(zl),nint(zr)):
        if 0 <= iz <= w3d.nzfull: izvolt[iz] = s.vplate[ip]
    # --- Set voltages between last plate and first quad
    zl = s.zright[-1]/w3d.dz
    zr = (s.zplate[-1] + s.accllen)/w3d.dz
    if   top.quadgp[0] < 0: vr = top.quadvy[0]
    elif top.quadgp[0] >= 0: vr = top.quadvx[0]
    vr = 0.
    for iz in iota(nint(zl),nint(zr)):
      if 0 <= iz <= w3d.nzfull:
        izvolt[iz] = s.vplate[-1] + (vr - s.vplate[-1])*(iz - zl)/(zr - zl)
    return izvolt

  def createouterwall(s):
    # --- Set boundary potential between plates and in long gap
    if s.lverbose:
      print "Adding surrounding pipe with voltage gradient"

    s.izvolt = s.getizvolt()

    if s.lsurround:
      def outerwall():
        f3d.srfrv_r = platesize - w3d.dx
     #srfrvout(outerwall,0.,0.,s.zplate[0],0.,0.,s.platesize,false,
     #         condid=997)
     #setconductorvoltage(s.izvolt,997)
     #srfrvout(outerwall,0.,s.zplate[0],s.zplate[-1],0.,0.,s.platesize,false,
     #         condid=998)
     #setconductorvoltage(s.izvolt,998)
      srfrvout(outerwall,0.,s.zright[-1],s.zplate[-1]+s.accllen,0.,0.,
               s.platesize,false,condid=999)
      setconductorvoltage(s.izvolt,999)

  def plotapertures(s,ip):
    xx = zeros(top.ninject,'d')
    yy = zeros(top.ninject,'d')
    plg(s.platesize*sin(arange(101)/100.*2.*pi),
        s.platesize*cos(arange(101)/100.*2.*pi))
    # --- Get location of apertures.
    xx[:] = top.xinject + (s.zplate[ip]-top.zinject)*top.xpinject
    yy[:] = top.yinject + (s.zplate[ip]-top.zinject)*top.ypinject
    for ij in xrange(top.ninject):
      plg(yy[ij] + s.rplate[ip]*sin(arange(101)/100.*2.*pi),
          xx[ij] + s.rplate[ip]*cos(arange(101)/100.*2.*pi))

  ########################################################
  def createparticlescraper(s):
    if s.lscraper:
      s.scraper = ParticleScraper(conductors=s.plates)

  ########################################################
  def setupcurvaturemockup(s):
    if not s.lcurvaturemockup: return
    # --- Add a uniform focusing element to model the transverse focusing the
    # --- would occur with the correctly spherical plates.
    if s.lcurvaturemockup and s.ap0full != 0.:
      dz = (w3d.zmmax - w3d.zmmin)/w3d.nz
      top.nemltsets = 1
      top.nesmult = 1
      top.nzemltmax = int(s.columnlen/dz)+1
      gchange("Mult_data",0)
      top.nzemlt = top.nzemltmax
      top.dzemlt = dz
      top.esemltp[:,0,0] = 0. # Set proportional to ezax below
      top.emltzs[0] = 0.
      top.emltze[0] = dz*top.nzemltmax
      top.emltid[0] = 1
      installafterfs(s.curvaturemockup)
      s.curvaturemockup()

  def curvaturemockup(s):
    sezax3d()
    ezax = top.ezax
    zzz = w3d.zmesh
    if lparallel:
      ezax = gatherallzarray(ezax)
      zzz = gatherallzmarray(zzz)
    nzm = min(w3d.nzfull,top.nzemltmax)
    ezax = ezax[:nzm+1]
    zzz = zzz[:nzm+1]
    top.esemltp[:nzm+1,0,0] = ezax/(s.a0full/s.ap0full-(max(zzz) - zzz))*(-2.)

