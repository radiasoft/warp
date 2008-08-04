
"""Defines ImplicitStep, which handles implicit time stepping"""
implicitstep_version = "$Id: implicitstep.py,v 1.9 2008/08/04 23:18:23 dave Exp $"

from warp import *
import controllers

class ImplicitStep(PackageBase):
  """
Handles implicit time stepping.
 - name="implicitstep": name to pass to the package command
 - niters=1: number of iterations of the implicit advance alogrithm
 - E0method=0: Method to use for getting the E0, the approximation to the
               future E field applied during the tilde step.
               0 is use zero E field
               1 is use E field from most recent field solve
  """

  def __init__(self,name="implicitstep",niters=1,
                    E0method=0):
    self.name = name
    registerpackage(self,self.name)

    # --- Save the input quantities
    self.niters = niters
    self.E0method = E0method

    # --- Initialize the timers
    self.timestep = 0.
    self.timedostep = 0.

    # --- Create the pid indices to hold the old E field data
    # --- and the old particle data.
    self.exoldpid = nextpid() - 1
    self.eyoldpid = nextpid() - 1
    self.ezoldpid = nextpid() - 1
    self.xoldpid  = nextpid() - 1
    self.yoldpid  = nextpid() - 1
    self.zoldpid  = nextpid() - 1
    self.uxoldpid = nextpid() - 1
    self.uyoldpid = nextpid() - 1
    self.uzoldpid = nextpid() - 1
    self.gaminvoldpid = nextpid() - 1
    # set up pid indices for predicted values of blended mover
    # WATCH OUT -- these are variables that may be used in both python
    # and fortran, and do NOT have the 1 subtracted.  So if we use
    # top.XXXpid from python where XXX is vdxold, vdyold, vdzold, bxpred,pypred,
    #  bzpred, uparoBpred, be sure to subtract 1 from it in the pid array!
    if (w3d.impinterp == 1):
      w3d.oldsetup()
    setuppgroup(top.pgroup)

  def __setstate__(self,dict):
    self.__dict__.update(dict)
    registerpackage(self,self.name)

  def saveOldParticleData(self):
    self.itold = top.it
    self.timeold = top.time
    self.zbeamold = top.zbeam
    self.zgridold = top.zgrid
    self.zgridprvold = top.zgridprv
    self.zgridndtsold = top.zgridndts.copy()
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js] and top.pgroup.nps[js] > 0:
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        top.pgroup.pid[i1:i2,self.xoldpid ] = top.pgroup.xp[i1:i2]
        top.pgroup.pid[i1:i2,self.yoldpid ] = top.pgroup.yp[i1:i2]
        top.pgroup.pid[i1:i2,self.zoldpid ] = top.pgroup.zp[i1:i2]
        top.pgroup.pid[i1:i2,self.uxoldpid] = top.pgroup.uxp[i1:i2]
        top.pgroup.pid[i1:i2,self.uyoldpid] = top.pgroup.uyp[i1:i2]
        top.pgroup.pid[i1:i2,self.uzoldpid] = top.pgroup.uzp[i1:i2]
        top.pgroup.pid[i1:i2,self.gaminvoldpid] = top.pgroup.gaminv[i1:i2]

  def restoreOldParticleData(self):
    top.it = self.itold
    top.time = self.timeold
    top.zbeam = self.zbeamold
    top.zgrid = self.zgridold
    top.zgridprv = self.zgridprvold
    top.zgridndts[:] = self.zgridndtsold
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js] and top.pgroup.nps[js] > 0:
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        if min(top.pgroup.pid[i1:i2,self.gaminvoldpid]) == 0.:
          # --- If there are some particles that were created at the tilde,
          # --- free streamed time level, they will not have old position
          # --- and velocity data. In those cases, copy the new data into
          # --- the saved data locations.
          noold = (top.pgroup.pid[i1:i2,self.gaminvoldpid] == 0)
          top.pgroup.pid[i1:i2,self.xoldpid ] = where(noold,top.pgroup.xp[i1:i2],top.pgroup.pid[i1:i2,self.xoldpid])
          top.pgroup.pid[i1:i2,self.yoldpid ] = where(noold,top.pgroup.yp[i1:i2],top.pgroup.pid[i1:i2,self.yoldpid])
          top.pgroup.pid[i1:i2,self.zoldpid ] = where(noold,top.pgroup.zp[i1:i2],top.pgroup.pid[i1:i2,self.zoldpid])
          top.pgroup.pid[i1:i2,self.uxoldpid ] = where(noold,top.pgroup.uxp[i1:i2],top.pgroup.pid[i1:i2,self.uxoldpid])
          top.pgroup.pid[i1:i2,self.uyoldpid ] = where(noold,top.pgroup.uyp[i1:i2],top.pgroup.pid[i1:i2,self.uyoldpid])
          top.pgroup.pid[i1:i2,self.uzoldpid ] = where(noold,top.pgroup.uzp[i1:i2],top.pgroup.pid[i1:i2,self.uzoldpid])
          top.pgroup.pid[i1:i2,self.gaminvoldpid ] = where(noold,top.pgroup.gaminv[i1:i2],top.pgroup.pid[i1:i2,self.gaminvoldpid])
        top.pgroup.xp [i1:i2] = top.pgroup.pid[i1:i2,self.xoldpid ]
        top.pgroup.yp [i1:i2] = top.pgroup.pid[i1:i2,self.yoldpid ]
        top.pgroup.zp [i1:i2] = top.pgroup.pid[i1:i2,self.zoldpid ]
        top.pgroup.uxp[i1:i2] = top.pgroup.pid[i1:i2,self.uxoldpid]
        top.pgroup.uyp[i1:i2] = top.pgroup.pid[i1:i2,self.uyoldpid]
        top.pgroup.uzp[i1:i2] = top.pgroup.pid[i1:i2,self.uzoldpid]
        top.pgroup.gaminv[i1:i2] = top.pgroup.pid[i1:i2,self.gaminvoldpid]

  def saveOldE(self):
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js] and top.pgroup.nps[js] > 0:
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        top.pgroup.pid[i1:i2,self.exoldpid] = top.pgroup.ex[i1:i2]
        top.pgroup.pid[i1:i2,self.eyoldpid] = top.pgroup.ey[i1:i2]
        top.pgroup.pid[i1:i2,self.ezoldpid] = top.pgroup.ez[i1:i2]

  def savePredEB(self):
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js] and top.pgroup.nps[js] > 0:
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        # From comments in constructor, need to subtract one from
        #  top.XXXpid's.
        top.pgroup.pid[i1:i2,top.expredpid-1] = top.pgroup.ex[i1:i2]
        top.pgroup.pid[i1:i2,top.eypredpid-1] = top.pgroup.ey[i1:i2]
        top.pgroup.pid[i1:i2,top.ezpredpid-1] = top.pgroup.ez[i1:i2]
        top.pgroup.pid[i1:i2,top.bxpredpid-1] = top.pgroup.bx[i1:i2]
        top.pgroup.pid[i1:i2,top.bypredpid-1] = top.pgroup.by[i1:i2]
        top.pgroup.pid[i1:i2,top.bzpredpid-1] = top.pgroup.bz[i1:i2]

  def restoreOldE(self,f=1.):
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js] and top.pgroup.nps[js] > 0:
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        top.pgroup.ex[i1:i2] = top.pgroup.pid[i1:i2,self.exoldpid]
        top.pgroup.ey[i1:i2] = top.pgroup.pid[i1:i2,self.eyoldpid]
        top.pgroup.ez[i1:i2] = top.pgroup.pid[i1:i2,self.ezoldpid]

  def averageOldAndNewEintoE(self):
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js] and top.pgroup.nps[js] > 0:
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        top.pgroup.ex[i1:i2] += top.pgroup.pid[i1:i2,self.exoldpid]
        top.pgroup.ey[i1:i2] += top.pgroup.pid[i1:i2,self.eyoldpid]
        top.pgroup.ez[i1:i2] += top.pgroup.pid[i1:i2,self.ezoldpid]
        top.pgroup.ex[i1:i2] *= 0.5
        top.pgroup.ey[i1:i2] *= 0.5
        top.pgroup.ez[i1:i2] *= 0.5

  def averageOldAndNewEintoOldE(self):
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js] and top.pgroup.nps[js] > 0:
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        top.pgroup.pid[i1:i2,self.exoldpid] += top.pgroup.ex[i1:i2]
        top.pgroup.pid[i1:i2,self.eyoldpid] += top.pgroup.ey[i1:i2]
        top.pgroup.pid[i1:i2,self.ezoldpid] += top.pgroup.ez[i1:i2]
        top.pgroup.pid[i1:i2,self.exoldpid] *= 0.5
        top.pgroup.pid[i1:i2,self.eyoldpid] *= 0.5
        top.pgroup.pid[i1:i2,self.ezoldpid] *= 0.5

  def zeroE(self):
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js] and top.pgroup.nps[js] > 0:
        i1 = top.pgroup.ins[js] - 1
        i2 = top.pgroup.ins[js] + top.pgroup.nps[js] - 1
        top.pgroup.ex[i1:i2] = 0.
        top.pgroup.ey[i1:i2] = 0.
        top.pgroup.ez[i1:i2] = 0.

  def fetche(self):
    top.lresetparticlee = true
    for js in range(top.pgroup.ns):
      if top.pgroup.limplicit[js]:
        fetche3d(top.pgroup,top.pgroup.ins[js],top.pgroup.nps[js],js+1)
    top.lresetparticlee = false

  def resetparticleb(self):
    top.pgroup.bx = 0.
    top.pgroup.by = 0.
    top.pgroup.bz = 0.

  def addextfields(self):
    # adds the external E and B fields to existing particle E and B fields
    pgroup=top.pgroup
    for js in range(pgroup.ns):
      # if top.pgroup.limplicit[js]:  why??  what happens for explicit species?
      npd=pgroup.nps[js]
      ins=pgroup.ins[js]-1
      imax=ins+npd
      x=pgroup.xp[ins:imax]
      y=pgroup.yp[ins:imax]
      z=pgroup.zp[ins:imax]
      uz=pgroup.uzp[ins:imax]
      gaminv=pgroup.gaminv[ins:imax]
      ex=pgroup.ex[ins:imax]
      ey=pgroup.ey[ins:imax]
      ez=pgroup.ez[ins:imax]
      bx=pgroup.bx[ins:imax]
      by=pgroup.by[ins:imax]
      bz=pgroup.bz[ins:imax]
      bendres=zeros(npd,"d")
      bendradi=zeros(npd,"d")
      dtr = 0.5*top.dt
      othere3d(npd,x,y,z,top.zbeam,top.zimin,top.zimax,top.straight, \
               top.ifeears,top.eears,top.eearsofz,top.dzzi,top.nzzarr, \
               top.zzmin,top.dedr,top.dexdx,top.deydy,top.dbdr, \
               top.dbxdy,top.dbydx,ex,ey,ez,bx,by,bz)

      exteb3d(npd,x,y,z,uz,gaminv,-dtr,dtr,bx,by,bz,ex,ey,ez,
              pgroup.sm[js],pgroup.sq[js],bendres,bendradi,top.gammabar,top.dt)

  def advancezgrid(self):
    # --- Accelerate grid frame.
    zcorrection = zeros(1,'d')
    acclbfrm(zcorrection)

    # --- Set timestep counter, time, and advance grid frame. The grid
    # --- frame is advanced here to be with the particles after the
    # --- position advance.
    # --- The zcorrection from the accleration of the beam frame is added on
    # --- by adding it to zgrid.
    # --- zgridprv is set here (as well as in padvnc3d) in case the
    # --- user has changed zbeam.

    top.it = top.it + 1
    if (top.lbeamcom):
      # --- In this case, it is still up to the user to set zgrid if zbeam is
      # --- set. Here, zbeam is not necessarily equal to zgrid.
      top.zgridprv = top.zgrid
    else:
      # --- Set zgridprv to zbeam so that the user only has to set zbeam.
      # --- Otherwise, zbeam is the same as zgrid, as set at the
      # --- end of padvnc3d.
      top.zgridprv = top.zbeam

    top.zgrid = top.zbeam + top.dt*top.vbeamfrm + zcorrection
    for indts in range(top.nsndts):
      if (top.it-1)%top.ndts[indts] == 0:
        # --- Only update zgridndts on steps when the group of particles
        # --- will be advanced.
        top.zgridndts[indts] = (top.zbeam +
                                top.dt*top.vbeamfrm*top.ndts[indts] +
                                zcorrection)

    # --- zgrid is integer number of dz's
    if (top.lgridqnt):
      top.zgrid = int(top.zgrid/w3d.dz + .5)*w3d.dz
      top.zgridndts[:] = int(top.zgridndts/w3d.dz + .5)*w3d.dz

  def advancezbeam(self):
    # --- Now the position of the grid can be advanced.
    # --- This is done here in case padvnc3d is called for multiple
    # --- pgroups
    if (top.lbeamcom):
      # --- Set zbeam so that it follows the center of mass of the beam.
      top.zbeam = getbeamcom(top.pgroup) - top.zbeamcomoffset
    else:
      top.zbeam = top.zgrid

    # --- zgridprv needs to be updated for the "synchv" step
    # --- Note that zgridprv is also set at the beginning of w3dexe
    top.zgridprv = top.zgrid
    top.time = top.time + top.dt

  #============================================================================
  def generate(self):
    # --- The generate is the normal w3d generate, including resetting
    # --- the E field.
    top.lresetparticlee = true
    w3dgen()
    top.lresetparticlee = false

    # --- Now fetch the initial E field and save it to the old arrays to
    # --- setup for the first part of the first step.
    self.fetche()
    self.saveOldE()
    if w3d.impinterp == 1:
      self.resetparticleb()
      self.addextfields()
      self.savePredEB()

  #==========================================================================
  def printpartextremes(self,leadstring):
    minx = minnd(top.pgroup.xp)
    maxx = maxnd(top.pgroup.xp)
    minz = minnd(top.pgroup.zp)
    maxz = maxnd(top.pgroup.zp)
    print leadstring,minx,maxx,minz,maxz

  #============================================================================
  def step(self):
    """This will be called by the generic step command. It is mostly copied from w3dexe."""
    substarttime = wtime()

    # --- Announce that we're running
    if top.it == 0: remark(" ***  Implicit particle simulation running")

    # --- The reset of the E fields is controlloed by this class
    top.lresetparticlee = false

    # --- Advance the grid frame
    self.advancezgrid()

    # --- Do some diagnostics
    stepid(top.it, top.time+top.dt, top.zgrid)

    # --- set logicals
    top.lfirst = false
    if top.ncall == 1: top.lfirst = true
    top.llast = false
    if top.ncall == top.maxcalls: top.llast = true

    # --- call the routine that does the actual work
    self.dostep()

    # --- Accumulate the run time.
    self.timestep += wtime() - substarttime

  #============================================================================
  def dostep(self):
    #print "STARTING IMPLICIT TIMESTEPPER"
    #print "uzp = ",top.pgroup.uzp[0]
    substarttime = wtime()

    # --- Set the internal lattice variables. This is not generally necessary
    # --- at this point (it is redundant most of the time, the next call to
    # --- setlatt in this subroutine is sufficient). There are cases where
    # --- this is required for consistency. Since it is cheap (time wise),
    # --- it is better to make sure the data is consistent than to save a
    # --- little bit of time. The value of nzl etc must be checked since other
    # --- packages (like WXY or ENV) may have reset it. For example, if the
    # --- ENV package is generated after the W3D package, nzl will be set to
    # --- zero. Switching back to W3D and running step, the internal lattice
    # --- would still be setup for the ENV package and so the step would
    # --- produce erroneaous results.
    if top.nzl == 0:
      top.nzl = top.nzlmax
      top.zlmin = top.zmmin
      top.zlmax = top.zmmax
      top.dzl = (top.zlmax - top.zlmin)/top.nzl
      top.dzli = 1./top.dzl
      top.zlmesh[:] = top.zlmin + dzl*zeros(0,top.nzl+1)

    setlatt()

    # --- This is a special routine needed when there is subcycling.
    #setupevensubcyclingrho(top.it)

    # --- This in effect finishes the implicit advance from the
    # --- previous step, advancing the velocity to n+1/2 and the
    # --- position to n+1 using a(n) = 1/2(a(n-1) + a(n+1))
    # --- Note that there is no gathering of rho since rho(n+1)
    # --- is not used.
    self.restoreOldE()
    top.laccumulate_rho = true # --- XXX
#    execfile("/home/rcohen/warp/scripts/rcdiags.py")
#    print "BEGIN STEP, z,zold,vdzold", printzzoldvzold()
#    print "taking step, Ez = ", top.pgroup.ez[0]
    if w3d.impinterp ==1:
      w3d.ipredcor = 1   # corrector for the blended mover
    if top.lspecial:
      padvnc3d("halfv",top.pgroup)
    else:
      padvnc3d("fullv",top.pgroup)

#    top.pgroup.yp=0.   # to avoid yp out of bounds in solve.  Need y=0 in range.
    top.laccumulate_rho = false # --- XXX
#    print "AFTER HALFV, uzp = ",top.pgroup.uzp[0]

    # --- Inject new particles.
    # --- Note that these particles may have zero for the old E fields
    # --- so the tilde advance won't be quite correct. This also means
    # --- that the average of the old and new E done below will give a
    # --- value that is half too small. This is not easy to fix since
    # --- it would require fetching E for only the new particles.
    inject3d(1,top.pgroup)
    controllers.userinjection()

    # --- Treat particles at boundaries
    particleboundaries3d(top.pgroup)

    # --- Now the position of the grid can be advanced.
    self.advancezbeam()

    # --- The next two variables are the left and right ends of the range
    # --- centered about the end of the current time step plus/minus one half a
    # --- step. The range is used is determining whether diagnostics are done
    # --- which are based on the z location of the beam frame.  The diagnostics
    # --- are done on the time step which ends closest to the value given in
    # --- the controlling arrays.
    # --- The absolute values are taken so that if dt < 0 or vbeamfrm < 0, then
    # --- it will still be true that zbeaml < zbeamr.
    zbeaml = top.zbeam - abs(0.5*top.vbeamfrm*top.dt)
    zbeamr = top.zbeam + abs(0.5*top.vbeamfrm*top.dt)

    # --- Set lattice; this is done just before field solve, and so is
    # --- relative to ZBEAM in the same way that self-fields are.
    setlatt()

    # --- Set logical flags to determine if "always" or "seldom" phase space 
    # --- plots, restart dumps, final timesteps, and moment accumulations
    # --- should be done at the end of this step.
    MACHEPS = 1.0e-14
    NCONTROL = 50
    top.lfinishd = ((top.it >= top.nt) or
                    (top.time >= top.tstop*(1.-MACHEPS)) or
                    (top.zbeam >= top.zstop))
    top.lalways  = (thisstep (top.it       ,top.itplalways,NCONTROL) or
                    thiszbeam(zbeaml,zbeamr,top.zzplalways,NCONTROL) or
                    thisstep (top.it       ,top.itplfreq,  NCONTROL) or
                    thiszbeam(zbeaml,zbeamr,top.zzplfreq,  NCONTROL))
    top.lseldom  = (thisstep (top.it       ,top.itplseldom,NCONTROL) or
                    thiszbeam(zbeaml,zbeamr,top.zzplseldom,NCONTROL) or
                    thisstep (top.it       ,top.itplps,    NCONTROL) or
                    thiszbeam(zbeaml,zbeamr,top.zzplps,    NCONTROL))
    top.lmoments = (thisstep (top.it       ,top.itmomnts,  NCONTROL) or
                    thiszbeam(zbeaml,zbeamr,top.zzmomnts,  NCONTROL))
    if top.nhist != 0:
      top.lhist  = (top.it%top.nhist) == 0
    else:
      top.lhist  = false

    top.ldump    = (top.it%top.itdump) == 0
    top.llabwn   = dolabwn()
    # --- This is not really needed since here, all steps are special, though
    # --- later, that could change.
    top.lspecial = (top.lfinishd or top.lalways or top.lseldom or top.ldump or
                    top.lmoments or top.lhist or top.llabwn or top.llast or
                    (top.it == 0) or top.allspecl)

    #self.printpartextremes("AFTER CORRECTOR, ")
    #print "AFTER CORR STEP, z,zold,vdzold", printzzoldvzold()
    #print "BEFORE SAVEOLDPART, uzp = ",top.pgroup.uzp[0]

    # --- Save the old particle data
    self.saveOldParticleData()

    if self.E0method == 0:
      # --- Zero out the E field, since the first approximation of the
      # --- future E is zero.
      self.zeroE()
    else:
      # --- Fetche the E field from the most recent field solve.
      self.fetche()

    for iter in range(self.niters):
      #print "AFTER SAVEOLD PART, z,zold,vdzold", printzzoldvzold()
      #print "AFTER SAVEOLDPART, uzp = ",top.pgroup.uzp[0]

      # --- Average the old and new E fields to get
      # --- a(n+1) = 1/2(a(n) + a(n+2))
      # --- where a(n+2) is the next estimate of the future fields.
      # --- Note that on the first iteration, a(n+2) = 0.
      # --- The result is stored in top.pgroup.ex etc.
      self.averageOldAndNewEintoE()

      # --- Now do a fullv advance to give the guess at the future x.
      # --- This does not advance any explicit species.
      self.advancezgrid()

      # If we are running the interpolated mover, we are now doing a predictor
      # step.  This should apply the current fields to the corrected x.
      # We can't do it here or it would muck up the Lorentz advance.  So
      # it is done in xpush3d_interp.
      if w3d.impinterp == 1:
        w3d.ipredcor = 0

      #print "BEFORE FULLV, uzp = ",top.pgroup.uzp[0]
      top.pgroup.ldoadvance[:] = top.pgroup.limplicit
      padvnc3d("fullv",top.pgroup)
      #print "AFTER PRED STEP, z,zold,vdzold", printzzoldvzold()
      #print "AFTER FULLV, uzp = ",top.pgroup.uzp[0]
      #self.printpartextremes("AFTER PREDICTOR, ")
      top.pgroup.ldoadvance[:] = true
      self.advancezbeam()

      # --- Treat particles at boundaries
      # --- Note that any particles that have they tilde position lost
      # --- will be permenantly lost, even if the position after the advance
      # --- with the full fields would not be lost.
      # --- This call is needed since parallel and reflecting boundaries must
      # --- be applied, as well as parallel inter-processor boundaries.
      particleboundaries3d(top.pgroup)

      # --- Collect charge density (rho_tilde) and current
      loadrho()
      loadj()

      # --- After the first iteration, or if using the E field from the most
      # --- recent field solve for E0, remove the delsq phi term from
      # --- the right hand side. For E0=0, on the first iteration,
      # --- delsq phi = 0.
      # --- The field solve will calculate del phi.
      # --- Also, save the current phi so the new phi can be constructed
      # --- as the sum of the current phi and del phi.
      if iter > 0 or self.E0method == 1:
        for solver in registeredsolvers():
          solver.removedelsqphi()
          #print "Maximum del rho ",iter,maxnd(abs(solver.rho))
          solver.savepreviouspotential()

      # --- Implicit field-solve for potential.
      # --- This assumes that an appropriate implicit field solver has
      # --- been registered.
      if w3d.lbeforefs: beforefs()
      fieldsol3d(-1)
      bfieldsol3d(-1)
      if w3d.lafterfs: afterfs()

      solver = getregisteredsolver()
      #print "Maximum del phi ",iter,maxnd(abs(solver.phi))

      # --- When the solver is calculating del phi, add together the previous
      # --- phi and del phi to get the new phi
      if iter > 0 or self.E0method == 1:
        for solver in registeredsolvers():
          solver.addinpreviouspotential()
      #plg(solver.rho[:,0,0])

      # --- Fetch the newly calculated implicit field.
      self.fetche()

      # --- Restore the old particle data. This can be done now that the new
      # --- implicit field is fetched at the free streaming positions.
      self.restoreOldParticleData()
      #print "AFTER RESTOREOLDP, uzp = ",top.pgroup.uzp[0]

      # --- If in parallel, the boundary conditions must be reapplied since
      # --- some particles may have crossed processors based on the tilde
      # --- positions - they need to be returned to the processors based
      # --- on the old positions.
      if npes > 0: particleboundaries3d(top.pgroup)

    # --- The newly calculated implicit field is now averaged with
    # --- the old field. This calculates a(n+1) = 1/2(a(n) + a(n+2))
    # --- The result is put into the old field arrays.
    # --- If we are running the interpolated mover we must also save the E,B
    # --- at the particle locations as E_pred,B_pred for use in calculating
    # --- v_drift at current time for next step.
    self.averageOldAndNewEintoOldE()
    if w3d.impinterp == 1:
      self.resetparticleb()
      self.addextfields()
      self.savePredEB()
      # this stores predicted total E, B for use in corrector.

    # --- Set the transverse E fields near any defined apertures.
    #set_aperture_e()

    # --- Is this the correct place?
    # --- Complete constant current and axially directed space-charge limited
    # --- injection with new fields including injected particles.
    inject3d(2,top.pgroup)

    # --- If a flag was set making this a "special" step.
    if top.lspecial:

      # --- Call this here since getzmmnt needs to have lvdts updated.
      # --- Note that it is still called in padvnc3d.
      setuppadvncsubcyclingaveraging(top.it,"synchv",top.pgroup)

      # --- Initialize the moments arrays which are calculated during the
      # --- synchv and gen phases.
      # --- 0. is passed in as a dummy for all of the particles coordinates
      # --- which are not used at this time.
      getzmmnt(1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1,
               top.nplive,0.,0.,0.,1,-1,top.ns,
               top.tempmaxp,top.tempminp,top.tempzmmnts0,top.tempzmmnts)

      # --- Copy the averaged field.
      self.restoreOldE()

      # --- Do a half-advance to bring v to same time level as the particles.
      padvnc3d("synchv",top.pgroup)
      #print "AFTER SYNCHV, uzp = ",top.pgroup.uzp[0]

      # --- Finalize the moments calculation and do other diagnostics.
      getzmmnt(1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,3,
               top.nplive,0.,0.,0.,1,1,top.ns,
               top.tempmaxp,top.tempminp,top.tempzmmnts0,top.tempzmmnts)
      #getlabwn()
      #rhodia()

      # --- Gather moments used in diagnostics at "special" timesteps only. 
      # --- Compute mean beam z velocity from current and line charge density 
      # --- on a 1-d mesh.  Also, calculate the electrostatic energy (getese), 
      # --- electrostatic potential on axis (sphiax), and the axial electric 
      # --- field (sezax).  
      #if w3d.lgetvzofz: getvzofz()
      #gtlchg()
      #srhoax()
      #getese()
      #sphiax()
      #sezax()

    # --- Charge density contour plot diagnostics.
    # --- Note also, that these diagnostics will be made with rho_tilde.
    # --- The actual rho is never calculated, but with multiple iterations,
    # --- rho_tilde will be close to the actual rho.
    if top.lalways or top.lseldom: pltfld3d("rho",always)
    if top.lseldom:                pltfld3d("rho",seldom)

    # --- Electrostatic potential contour plot diagnostics 
    if top.lalways or top.lseldom: pltfld3d("phi",always)
    if top.lseldom:                pltfld3d("phi",seldom)

    # --- 1d array plot diagnostics.
    if top.lalways or top.lseldom: onedplts(always)
    if top.lseldom:                onedplts(seldom)

    # --- Phase space diagnostics
    if top.lalways or top.lseldom: psplots(always)
    if top.lseldom:                psplots(seldom)

    # --- Finally, moment diagnostic printout and history storage
    if top.lspecial: minidiag(top.it,top.time,top.lspecial)

    self.timedostep += wtime() - substarttime

  #============================================================================
  def finish(self):
    gfree("Fields3d")
    gfree("Fields3dParticles")
    gfree("Hist")
    gfree("Win_Moments")
    gfree("Z_Moments")
    gfree("Lab_Moments")
    gfree("Moments")
    gfree("Lattice")
    gfree("LatticeInternal")
    gfree("Z_arrays")

