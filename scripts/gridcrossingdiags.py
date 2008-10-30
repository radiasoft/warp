"""Diagnostics captured as particles cross grid cells.
"""
__all__ = ['GridCrossingDiags']
from warp import *

class GridCrossingDiags(object):
  """

Sets up diagnostics at z grid cells that are gathered from particles that
cross the cell.
  - js: species of particles to include. Currently can handle only a single
        species.
  - zmmin,zmmax,dz,nz: grid parameters, all taken from w3d if not supplied.
  - nzscale=1: multiplier on nz - makes it easy to use the w3d grid parameters
               but with a differing number of grid points.
  - nr,rmax: radial extent of radial profile diagnostic. If not given,
             the radial diagnostic is not done.
  - ztarget: Z location with the radial profile is calculated.

The following quantities are calculated:
count: count of the number of particles that cross the cell each time
current: current, same as count but scaled by particle charge and 1/top.dt.
rrms:
rprms:

  """

  def __init__(self,js,zmmin=None,zmmax=None,dz=None,nz=None,nzscale=1,
               nhist=None,nr=None,rmax=None,ztarget=None):
    self.js = js
    self.zmmin = zmmin
    self.zmmax = zmmax
    self.dz = dz
    self.nz = nz
    self.nzscale = nzscale
    self.nhist = nhist
    self.nr = nr
    self.rmax = rmax
    self.ztarget = ztarget

    self.zoldpid = nextpid()

    self.count = None
    self.current = None
    self.rrms = None
    self.rprms = None

    installafterstep(self.getdiagnostics)

  def getdiagnostics(self):
    # --- Initialize grid parameters if needed. This is done here
    # --- in case the grid had not been setup yet when init was called.
    if self.zmmin is None: self.zmmin = w3d.zmmin
    if self.zmmax is None: self.zmmax = w3d.zmmax
    if self.dz is None: self.dz = w3d.dz/self.nzscale
    if self.nz is None: self.nz = w3d.nz*self.nzscale

    # --- Create handy locals.
    js = self.js
    zmmin = self.zmmin
    zmmax = self.zmmax
    dz = self.dz
    nz = self.nz

    zz = self.ztarget
    rmax = self.rmax
    nr = self.nr
    ldoradialdiag = (nr is not None)
    if ldoradialdiag: dr = rmax/nr

    zbeam = top.zbeam
    zoldpid = self.zoldpid

    # --- Initialize the data lists.
    if self.count is None:
      self.time = []
      self.count = [zeros(1+nz,'d')]
      self.current = [zeros(1+nz,'d')]
      self.rrms = [zeros(1+nz,'d')]
      self.rprms = [zeros(1+nz,'d')]
      self.dummy = zeros(1+nz,'d')
      if ldoradialdiag:
        self.rprofile = [zeros(1+nr,'d')]
        self.rprofilecount = zeros(1+nr,'d')
        self.rprofilemesh = iota(0,nr)*dr

    if self.nhist is None: nhist = top.nhist
    else:                  nhist = self.nhist

    # --- The data is gathered from top.it-nhist/2 to top.it+nhist/2-1.
    # --- This resets things at the half way point.
    if top.it%nhist == int(nhist/2):
      co = self.count[-1]
      cu = self.current[-1]
      rr = self.rrms[-1]
      rp = self.rprms[-1]
      if ldoradialdiag:
        rprof = self.rprofile[-1]
        rprofco = self.rprofilecount[-1]

      # --- Finish the calculation, gathering data from all processors and
      # --- dividing out the count.
      co[...] = parallelsum(co)
      rr[...] = parallelsum(rr)
      rr[...] = sqrt(rr/where(co==0.,1.,co))
      rp[...] = parallelsum(rp)
      rp[...] = sqrt(rp/where(co==0.,1.,co))

      if ldoradialdiag:
        rprof[...] = parallelsum(rprof)
        #rprof[0] = rprof[0]*0.75/(pi*0.5*0.5*dr**2)
        #rprof[1:] = rprof[1:]/(2.*pi*self.rprofilemesh[1:]*dr)

      # --- Scale the current appropriately.
      cu[...] = co*(top.pgroup.sq[js]*top.pgroup.sw[js]/(top.dt*nhist))

      # --- Create space for the next set of data.
      if me == 0:
        self.time.append(top.time)
        self.count.append(zeros(1+nz,'d'))
        self.current.append(zeros(1+nz,'d'))
        self.rrms.append(zeros(1+nz,'d'))
        self.rprms.append(zeros(1+nz,'d'))
        if ldoradialdiag:
          self.rprofile.append(zeros(1+nr,'d'))
      else:
        # --- On other processors, just zero out the existing arrays.
        # --- There's no reason to keep the history on all processors.
        self.count[0][:] = 0.
        self.current[0][:] = 0.
        self.rrms[0][:] = 0.
        self.rprms[0][:] = 0.
        if ldoradialdiag:
          self.rprofile[0][:] = 0.

    # --- The code below is all local and can be skipped if there are no
    # --- particles locally.
    if top.pgroup.nps[js] == 0: return

    rnew = getr(js=js,gather=0)
    rpnew = getrp(js=js,gather=0)
    znew = getz(js=js,gather=0)
    zold = getpid(js=js,id=zoldpid-1,gather=0)

    iznew = int((znew - (zbeam + zmmin))/dz)
    izold = int((zold - (zbeam + zmmin))/dz)

    icrossed = (iznew > izold)

    izc = compress(icrossed,iznew)
    rc = compress(icrossed,rnew)
    rpc = compress(icrossed,rpnew)
    np = len(izc)

    if top.wpid > 0:
      weight = getpid(js=js,id=top.wpid-1,gather=0)
      ww = compress(icrossed,weight)
    else:
      ww = ones(np,'d')

    zc = izc.astype('d')
    deposgrid1d(1,np,zc,ww,nz,self.count[-1],self.dummy,0.,nz)
    deposgrid1d(1,np,zc,ww*rc**2,nz,self.rrms[-1],self.dummy,0.,nz)
    deposgrid1d(1,np,zc,ww*rpc**2,nz,self.rprms[-1],self.dummy,0.,nz)

    if ldoradialdiag:
      izc = logical_and(zold<=zz,znew>=zz)
      rc = compress(izc,rnew)
      if len(rc) > 0:
        vz = compress(izc,getvz(js=js,gather=0))
        ee = 0.5*top.pgroup.sm[js]*vz**2
        if top.wpid > 0:
          weight = getpid(js=js,id=top.wpid-1,gather=0)
          ww = compress(izc,weight)*top.pgroup.sw[js]
        else:
          ww = top.pgroup.sw[js]
        deposgrid1d(1,len(rc),rc,ee*ww,nr,self.rprofile[-1],self.rprofilecount,
                    0.,rmax)

    # --- Save particle z positions.
    i1 = top.pgroup.ins[js] - 1
    i2 = i1 + top.pgroup.nps[js]
    top.pgroup.pid[i1:i2,zoldpid-1] = top.pgroup.zp[i1:i2]

