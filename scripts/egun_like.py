from warp import *
import string
egun_like_version = "$Id: egun_like.py,v 1.11 2002/10/25 23:00:41 dave Exp $"
############################################################################
# EGUN_LIKE algorithm for calculating steady-state behavior in a ion source.
#
# Particles are injected on one time step only and the injection is
# turned off.  Those particles are then tracked through the system until
# there are no particles left.  On each time step, the charge density from
# the particles is accumulated in the array rho.  After all of the
# particles leave the system, the new field is calculated using the
# accumulated charge density.  Also, a selection of particles is saved
# each time step for plotting.
#
# At the end, the charge density is in the array rho, and the particle
# data is in the  particle arrays and ins and nps are set properly.
############################################################################


##############################################################################
print "The command to run is gun(iter,ipsave,save_same_part,maxtime)."
print "For more info type doc(gun)."
print "To do a proper restart, use the command restartgun(file)."
print "To recover from a keyboard interrupt, use the command recovergun()"
print "before doing more iterations."
##############################################################################
# The arguments are preserved by having a shadow with the same name prefixed
# with an underscore. The underscored varaible is the variable used. The one
# one without the underscor is only used as input arguments.
##############################################################################

gun_iter = 0
gun_steps = 0

# --- reduction factor in number of particles saved each timestep
_ipstep = top.npinject

# --- number of particles to save
_ipsave = 0

# --- Sets whether the particles save on each time step are the same particles.
_save_same_part = false

# --- Save value of inject since it is changed during the iterations.
_oinject = top.inject

# --- Save values of nztinjmn and nztinjmx since they are changed also.
if (top.ntinj > 0):
  _onztinjmn = top.nztinjmn
  _onztinjmx = top.nztinjmx

# --- Save value of fstype since it is changed also.
_ofstype = top.fstype

# --- Vz fuzz, particles with uzp bigger are selected.  It should be small
# --- since uzp may be small for newly injected particles.
_vzfuzz = 1.e-20

# --- Read in getzmom script. Only used to make final and initial calls to
# --- getzmmnt routine. Moments are calculated during timesteps and include all
# --- particles.
import getzmom

def gun(iter=1,ipsave=None,save_same_part=None,maxtime=None,
        laccumulate_zmoments=None,rhoparam=None):
  """
Performs steady-state iterations
  - iter=1 number of iterations to perform
  - ipsave=0 number of particles to save from the last iteration
  - save_same_part=0 when true, save same particles each time step instead
    of random particles, i.e. saves particle trajectories
  - maxtime=3*transittime maximum time each iteration will run
  - laccumulate_zmoments=false: When set to true, z-moments are accumulated
    over multiple iterations. Note that getzmom.zmmnt(3) must be called
    by the user to finish the moments calculation.
  - rhoparam=None: Amount of previous rho to mix in with the current rho. This
    can help the relaxation toward a steady state. Caution should be used
    when using this option.
  Note that ipsave and save_same_part are preserved in between calls
  """
  global _oinject,_ofstype,_onztinjmn,_onztinjmx
  global _ipstep,_ipsave,_save_same_part
  global gun_iter,gun_time,gun_steps

  # --- Save general parameters which are modified in this routine
  _oinject = top.inject
  _ofstype = top.fstype
  if (top.ntinj > 0):
    _onztinjmn = top.nztinjmn
    _onztinjmx = top.nztinjmx
  _ifzmmnt = top.ifzmmnt
  _laccumulate_zmoments = top.laccumulate_zmoments
  if laccumulate_zmoments is None:
    laccumulate_zmoments = top.laccumulate_zmoments

  if ipsave: _ipsave = ipsave
  if save_same_part: _save_same_part = save_same_part

  # --- Save current value of top.nhist
  nhist = top.nhist
  top.nhist = 0

  # --- Set injection relaxation parameter if it has not already been set by
  # --- the user.
  if top.inj_param == 1.: top.inj_param = 0.5

  # --- Set logical so that the charge density is accumulated over
  # --- all of the time steps for each iteration.
  top.laccumulate_rho = true

  # --- Turn off rho-diagnostic and calculation of ese
  w3d.lrhodia3d = 0
  w3d.lgetese3d = 0

  # --- Set the attribute 'dump' on the rho array so it is automatically saved
  # --- on a dump. This is done since rho holds the current state of the
  # --- solution.
  if 'dump' not in string.split(w3d.getvarattr('rho')):
    w3d.addvarattr("rho","dump")

  # --- Set verbosity so that the one line diagnostic is not printed out.
  top.verbosity = 1

  # --- Estimate the time that will be required for the particles
  # --- to propagate through the system. It is based off of the Child-Langmuir
  # --- solution for a diode. The diode length is assumed to be (nzfull*dz),
  # --- the diode voltage is assumed to be abs(phi(,,0)-phi(,,nz)).
  transittime = 3.*(w3d.nzfull*w3d.dz)*sqrt(0.5*top.sm[0]/abs(top.sq[0])/ \
                abs(getphi(w3d.ix_axis,w3d.iy_axis,0,bcast=1) - \
                    getphi(w3d.ix_axis,w3d.iy_axis,w3d.nz,bcast=1)))
  # --- Set the default maxtime to 3*transittime. The factor of 3 is a random
  # --- guess at a safety factor. The maxtime is used since in some cases it
  # --- is possible for some particles to get stuck in a low field region,
  # --- requiring a large number of time steps to move out of the system.
  if not maxtime: maxtime = 3*transittime

  # --- make multiple iterations
  for i in xrange(iter):

    # --- set number of particles to save
    # --- assumes the variable 'it' has only been advanced in egun mode and
    # --- that at least one iteration has already been done
    if (_ipsave > 0 and gun_iter > 0):
      if lparallel: npisum = sum(parallelsum(top.npinje_s))
      else: npisum = sum(top.npinje_s)
      _ipstep = gun_steps*npisum/_ipsave
      if (_ipstep < 1): _ipstep = 1

    # --- set number of particles to zero.
    top.nps = 0

    # --- turn on injection (type 1 or 2) for one time step.
    if (top.ntinj > 0):
      top.nztinjmn = _onztinjmn
      top.nztinjmx = _onztinjmx
    top.inject = _oinject

    # --- turn off field solver
    top.fstype = -1

    # --- If rhoparam is not None, then save the previous rho
    if rhoparam is not None: rhoprevious = w3d.rho + 0.

    # --- Zero the charge density array
    w3d.rho = 0.
    if w3d.solvergeom==w3d.RZgeom: reset_rzmgrid_rho()
    top.curr = 0.

    # --- If this is the final iteration and if zmoments are being calculated,
    # --- make the initial call to zero the arrays.
    if ((i == iter-1 or (gun_iter%nhist) == 0) and _ifzmmnt > 0):
      top.ifzmmnt = _ifzmmnt
      getzmom.zmmnt(1)
      # --- Make sure that moments are calculated on each time step. This is
      # --- the only way that the data will make sense.
      top.itmomnts[0:4] = [0,top.nt,1,0]
      # --- Make sure that the moments are accumulated.
      top.laccumulate_zmoments = true
    else:
      # --- Make sure the zmoments are not calculated so the time isn't wasted.
      top.itmomnts[0:4] = [top.nt,top.nt,top.nt,0]
      top.ifzmmnt = 0

    # --- Save current time
    gun_time = top.time

    # --- Make one time step to inject the batch of particles.
    step(1)
    tmp_gun_steps = 1
    print "Number of particles injected = %d"%(top.npinject)

    # --- check if any particles were injected
    if lparallel: npssum = sum(parallelsum(top.nps))
    else: npssum = sum(top.nps)
    if (npssum == 0): raise 'No particles injected'

    # --- only save particles on last iteration
    if (i == iter-1 and _ipsave > 0 and _ipstep > 0):

      # --- Set ins and nps for saved particles
      # --- This only creates the arrays - the actual values don't
      # --- mean anything at this point. They are set later on.
      ins_save = top.ins + top.nps
      nps_save = 0*top.nps

      for js in xrange(top.ns):

        # --- Force particles to the beginning of their block.
        # --- This only needs to be done if there are particles now.
        if (top.nps[js] > 0):
          copypart(top.npmax_s[js]+1,top.nps[js],0,top.ins[js])

        # --- Reset particles counters.
        # --- This always needs to be done since particles may be added
        # --- later.
        top.ins[js] = top.npmax_s[js]+1
        ins_save[js] = top.ins[js] + top.nps[js]
        nps_save[js] = 0

        if (top.nps[js] > 0):
          # --- get indices of live particles.
          if (_save_same_part):
            ip1 = top.ins[js]
            ip2 = top.ins[js]+top.nps[js]-1
            ip3 = _ipstep
            ii = compress(greater(top.uzp[ip1-1:ip2:ip3],_vzfuzz),
                          iota(ip1,ip2,ip3))
          else:
            ip1 = top.ins[js]
            ip2 = top.ins[js]+top.nps[js]-1
            ip3 = 1
            ii = compress(greater(top.uzp[ip1-1:ip2:ip3],_vzfuzz),
                          iota(ip1,ip2,ip3))
            ii = compress(less(ranf(ii),1./_ipstep),ii)

          # --- save data of just injected particles
          if (len(ii) > 0):
            npguess = int(1.5*gun_steps*top.npinje_s[js]/_ipstep)+len(ii),
            nplost = ins_save[js] - top.ins[js] - top.nps[js]
            chckpart(js+1,0,npguess + nplost,false)
            copypart(ins_save[js]+nps_save[js],len(ii),ii,-1)
            nps_save[js] = nps_save[js] + len(ii)

    # --- Turn injection off for remaing time steps. inject is set to a value
    # --- greater than zero so that inject3d subroutine is called so it can
    # --- strip off particles and reduce nps and particles leave the system.
    if (top.ntinj > 0):
      top.nztinjmn = 0
      top.nztinjmx = 0
    top.inject = 100

    # --- Run until all particles are out of the system (no injection, no field
    # --- solves).  Accumulation of charge density and particle moments is done
    # --- automatically in the code.  Save particle data each time step on last
    # --- iteration only.
    maxvz = 2.*_vzfuzz+1.
    while (npssum > 0 and maxvz>_vzfuzz and
           top.time-gun_time < maxtime):
      step()
      tmp_gun_steps = tmp_gun_steps + 1
      # --- only save particles on last iteration
      if (i == iter-1 and _ipsave > 0 and _ipstep > 0):
        for js in xrange(top.ns):
          # --- Make sure that there are actually particles to save
          if (top.nps[js] > 0):
            if (_save_same_part):
              ip1 = top.ins[js]
              ip2 = top.ins[js]+top.nps[js]-1
              ip3 = _ipstep
              ii = compress(greater(top.uzp[ip1-1:ip2:ip3],_vzfuzz),
                            iota(ip1,ip2,ip3))
            else:
              ip1 = top.ins[js]
              ip2 = top.ins[js]+top.nps[js]-1
              ip3 = 1
              ii = compress(greater(top.uzp[ip1-1:ip2:ip3],_vzfuzz),
                            iota(ip1,ip2,ip3))
              ii = compress(less(ranf(ii),1./_ipstep),ii)

            # --- save data of just injected particles
            if (len(ii) > 0):
              npguess = nps_save[js] + len(ii)
              nplost = ins_save[js] - top.ins[js] - top.nps[js]
              chckpart(js+1,0,npguess + nplost,false)
              copypart(ins_save[js]+nps_save[js],len(ii),ii,-1)
              nps_save[js] = nps_save[js] + len(ii)

      if lparallel:
        npssum = sum(parallelsum(top.nps))
        maxvz = parallelmax(top.vzmaxp)
      else:
        npssum = sum(top.nps)
        maxvz = top.vzmaxp

    # --- If rhoparam is not None, mix in the previous rho with the
    # --- new rho
    if rhoparam is not None:
      w3d.rho[:,:,:] = (1.-rhoparam)*w3d.rho + rhoparam*rhoprevious

    # --- Do field solve including newly accumulated charge density.
    # --- The call to perrho3d is primarily needed for the parallel version.
    perrho3d(w3d.rho,w3d.nx,w3d.ny,w3d.nz,top.periinz)
    top.fstype = _ofstype
    fieldsol(-1)
    top.fstype = -1

    # --- Do final work for zmoments calculation
    if ((i == iter-1 or (gun_iter%nhist) == 0) and top.ifzmmnt > 0):
      top.laccumulate_zmoments = laccumulate_zmoments
      getzmom.zmmnt(3)

    # --- Save the history data
    if (gun_iter%nhist) == 0 and top.ifzmmnt > 0:
      top.nhist = gun_iter
      minidiag(gun_iter,gun_time,false)
      top.nhist = 0
      top.hzbeam[top.jhist] = gun_iter

    # --- Print out warning message if needed.
    if top.time-gun_time > maxtime:
      print "Warning: maxtime exceeded - this may be corrected in the next"
      print "iteration. If it is not, increase the value of maxtime argument,"
      print "or look for other problems."
      print "maxtime = ",maxtime

    gun_steps = tmp_gun_steps
    gun_iter = gun_iter + 1
    gun_time = top.time - gun_time
    print "Number of iterations = %d"%gun_iter

  # --- end of multiple iterations

  # --- Set saved particles to be live particles (for diagnostics only).
  if (_ipsave > 0 and _ipstep > 0):
    top.ins[:] = ins_save
    top.nps[:] = nps_save

  # --- Change what is plotted at the bottom of each frame
  stepid(gun_iter,gun_time,top.zbeam)

  # --- Set some additional diagnostic data
  top.vzofz = top.vzbarz

  # --- Return variables to their original values
  top.fstype = _ofstype
  if (top.ntinj > 0):
    top.nztinjmn = _onztinjmn
    top.nztinjmx = _onztinjmx
  top.inject = _oinject
  top.nhist = nhist
  top.ifzmmnt = _ifzmmnt
  top.laccumulate_zmoments = _laccumulate_zmoments


########################################################################
# Restart routine for gun calculations
def restartgun(filename):
  restore(filename)
  fieldsol(0)
  setlatt()

########################################################################
# Recover the gun calculation after a keyboard interrupt
def recovergun():
  # --- Return some variables to the original values
  top.fstype = _ofstype
  if (top.ntinj > 0):
    top.nztinjmn = _onztinjmn
    top.nztinjmx = _onztinjmx
  top.inject = _oinject

