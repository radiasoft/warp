from warp import *
import string
import curses.ascii
import sys
import adjustmesh3d
import __main__
egun_like_version = "$Id: egun_like.py,v 1.35 2005/01/12 17:17:39 dave Exp $"
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
# with an underscore. The underscored variable is the variable used. The one
# without the underscore is only used as input arguments.
##############################################################################

gun_iter = 0
gun_steps = 0

# --- reduction factor in number of particles saved each timestep
_ipstep = 0

# --- number of particles to save
_ipsave = 0

def egungetdata():
  print "_ipstep = ",_ipstep
  print "_ipsave = ",_ipsave

# --- Sets whether the particles save on each time step are the same particles.
_save_same_part = false

# --- Save value of inject since it is changed during the iterations.
_oinject = top.inject

# counter for egundata arrays update
_izdata = 1

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

def plottraces():
  if w3d.solvergeom == w3d.XZgeom:
    warpplp(top.xp[top.ins[0]:top.ins[0]+top.nps[0]],
            top.zp[top.ins[0]:top.ins[0]+top.nps[0]])
  elif w3d.solvergeom == w3d.RZgeom:
    warpplp(sqrt(top.xp[top.ins[0]:top.ins[0]+top.nps[0]]**2+
                 top.yp[top.ins[0]:top.ins[0]+top.nps[0]]**2),
            top.zp[top.ins[0]:top.ins[0]+top.nps[0]])
  else:
    plsys(9)
    warpplp(top.xp[top.ins[0]:top.ins[0]+top.nps[0]],
            top.zp[top.ins[0]:top.ins[0]+top.nps[0]])
    plsys(10)
    warpplp(top.yp[top.ins[0]:top.ins[0]+top.nps[0]],
            top.zp[top.ins[0]:top.ins[0]+top.nps[0]])
  pyg_pending()
  pyg_idler()
  
def gun(iter=1,ipsave=None,save_same_part=None,maxtime=None,
        laccumulate_zmoments=None,rhoparam=None,
        lstatusline=true,insertbeforeiter=None,insertafteriter=None,
        ipstep=None,egundata_window=-1,plottraces_window=-1,
        egundata_nz=None,egundata_zmin=None,egundata_zmax=None):
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
  - lstatusline=1: when try, a line is printed and continuously updated
                   showing the status of the simulation.
  - insertbeforeiter=None: function to be called before each iteration.
  - insertafteriter=None: function to be called after each iteration.
  - egundata_window=2: window in which to display egundata curves for z close 
                       to w3d.zmmax. Set to a negative number to deactivate 
                       plotting.
  - plottraces_window=1: window in which to plot traces.Set to a negative 
                         number to deactivate plotting.
  Note that ipsave and save_same_part are preserved in between calls
  """
  global _oinject,_ofstype,_onztinjmn,_onztinjmx
  global _ipstep,_ipsave,_save_same_part, _izdata
  global gun_iter,gun_time,gun_steps
  global rhoprevious
  global zd, egundata_curr, egundata_xrmsz, egundata_yrmsz
  global egundata_xprmsz, egundata_yprmsz, egundata_epsnxz, egundata_epsnyz

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

  if ipsave is not None: _ipsave = ipsave
  if save_same_part is not None: _save_same_part = save_same_part

  # --- Save current value of top.nhist
  nhist = top.nhist
  top.nhist = 0

  # --- Set injection relaxation parameter if it has not already been set by
  # --- the user.
  if top.inj_param == 1.: top.inj_param = 0.5

  # --- Set logical so that the charge density is accumulated over
  # --- all of the time steps for each iteration.
  top.laccumulate_rho = true

  # --- Set 
  top.clearlostpart = 1

  # --- Turn off rho-diagnostic and calculation of ese
  w3d.lrhodia3d = false
  w3d.lgetese3d = false
  w3d.lgtlchg3d = false
  w3d.lgetvzofz = false
  w3d.lsrhoax3d = false
  w3d.lsphiax3d = false
  w3d.lsezax3d  = false

  # --- setup egundata
  if egundata_nz is not None:
    if egundata_zmin is None:
      egundata_zmin = w3d.zmmin+0.01*(w3d.zmmax-w3d.zmmin)
    if egundata_zmax is None:
      egundata_zmax = w3d.zmmax-0.01*(w3d.zmmax-w3d.zmmin)
    zd = egundata_zmin+arange(egundata_nz)*(egundata_zmax-egundata_zmin)/(egundata_nz-1)
    # --- Store the data in lists. This allows an arbitrary number of
    # --- iterations since the next data is just appended.
    egundata_curr   = []
    egundata_xrmsz  = []
    egundata_yrmsz  = []
    egundata_xprmsz = []
    egundata_yprmsz = []
    egundata_epsnxz = []
    egundata_epsnyz = []

  # --- install plottraces
  if plottraces_window>-1:
    installafterstep(plottraces)


  # --- Set the attribute 'dump' on the rho array so it is automatically saved
  # --- on a dump. This is done since rho holds the current state of the
  # --- solution.
  if 'dump' not in string.split(w3d.getvarattr('rho')):
    w3d.addvarattr("rho","dump")

  # --- Set verbosity so that the one line diagnostic is not printed out.
  top.verbosity = 1

  if not maxtime:
    # --- Estimate the time that will be required for the particles
    # --- to propagate through the system. It is based off of the Child-Langmuir
    # --- solution for a diode. The diode length is assumed to be (nzfull*dz),
    # --- the diode voltage is assumed to be abs(phi(,,0)-phi(,,nz)).
    delphi = abs(getphi(w3d.ix_axis,w3d.iy_axis,0,bcast=1) -
                 getphi(w3d.ix_axis,w3d.iy_axis,w3d.nz,bcast=1))
    if delphi == 0.: delphi = smallpos
    transittime = (3.*(w3d.nzfull*w3d.dz)*sqrt(0.5*top.sm[0]/abs(top.sq[0])/
                                               delphi))
    # --- Set the default maxtime to 3*transittime. The factor of 3 is a random
    # --- guess at a safety factor. The maxtime is used since in some cases it
    # --- is possible for some particles to get stuck in a low field region,
    # --- requiring a large number of time steps to move out of the system.
    maxtime = 3*transittime

  # --- make multiple iterations
  for i in xrange(iter):
    
    # --- plot field and conductors 
    if plottraces_window>-1:
      window(plottraces_window)
      fma()
      if w3d.solvergeom == w3d.XZgeom or w3d.solvergeom == w3d.RZgeom:
        pfzx()
        limits(w3d.zmmin,w3d.zmmax,w3d.xmmin,w3d.xmmax)
      else:
        pfzx(view=9)
        plsys(9);limits(w3d.zmmin,w3d.zmmax,w3d.xmmin,w3d.xmmax)
        pfzy(view=10)
        plsys(10);limits(w3d.zmmin,w3d.zmmax,w3d.ymmin,w3d.ymmax)

    # --- call insertbeforeriter if defined
    if insertbeforeiter is not None:
       insertbeforeiter()

    # --- set number of particles to save
    # --- assumes the variable 'it' has only been advanced in egun mode and
    # --- that at least one iteration has already been done
    if (_ipsave > 0 and gun_iter > 0):
      npisum = sum(parallelsum(top.npinje_s))
      _ipstep = gun_steps*npisum/_ipsave
      if (_ipstep < 1): _ipstep = 1

    # --- If w3d.l_inj_regular is set to true, then always save the
    # --- trajectories of all of the particles.
    if w3d.l_inj_regular and ipsave is None:
      _ipsave = 1
      _ipstep = 1

    if ipstep is not None: _ipstep = ipstep

    # --- set number of particles to zero.
    top.nps = 0
    top.ins[0:]=top.npmax_s[1:]
        
    # --- turn on injection (type 1 or 2) for one time step.
    if (top.ntinj > 0):
      top.nztinjmn = _onztinjmn
      top.nztinjmx = _onztinjmx
    top.inject = _oinject

    # --- turn off field solver
    top.fstype = -1

    # --- If rhoparam is not None, then save the previous rho
    if rhoparam is not None:
     if(w3d.solvergeom<>w3d.RZgeom):
      rhoprevious = w3d.rho + 0.
     else:
       for ig in range(frz.ngrids):
         if(ig==0):
           g = frz.basegrid
           rhoprevious = [g.rho.copy()]
         else:
           try:
             g = g.next
           except:
             g = g.down
           rhoprevious = rhoprevious+[g.rho.copy()]

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
    npssum = sum(parallelsum(top.nps))
    if (npssum == 0): raise 'No particles injected'

    # --- only save particles on last iteration
    if (i == iter-1 and _ipstep > 0):

      # --- Set ins and nps for saved particles
      # --- This only creates the arrays - the actual values don't
      # --- mean anything at this point. They are set later on.
      ins_save = zeros(top.ns)
      nps_save = zeros(top.ns)

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
            ii = iota(ip1,ip2,ip3)
          else:
            ip1 = top.ins[js]
            ip2 = top.ins[js]+top.nps[js]-1
            ip3 = 1
            ii = iota(ip1,ip2,ip3)
            ii = compress(less(ranf(ii),1./_ipstep),ii)

          # --- save data of just injected particles
          if (len(ii) > 0):
            npguess = int(1.5*gun_steps*top.npinje_s[js]/_ipstep)+len(ii)
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
    while (npssum > 0 and top.time-gun_time < maxtime):
      step()
      if lstatusline: statusline()
      tmp_gun_steps = tmp_gun_steps + 1
      # --- only save particles on last iteration
      if (i == iter-1 and _ipstep > 0):
        for js in xrange(top.ns):
          # --- Make sure that there are actually particles to save
          if (top.nps[js] > 0):
            if (_save_same_part):
              ip1 = top.ins[js]
              ip2 = top.ins[js]+top.nps[js]-1
              ip3 = _ipstep
              ii = iota(ip1,ip2,ip3)
            else:
              ip1 = top.ins[js]
              ip2 = top.ins[js]+top.nps[js]-1
              ip3 = 1
              ii = iota(ip1,ip2,ip3)
              ii = compress(less(ranf(ii),1./_ipstep),ii)

            # --- save data of just injected particles
            if (len(ii) > 0):
              if w3d.l_inj_rec_inittime:
                ip1 = top.ins[js]-1
                ip2 = top.ins[js]+top.nps[js]-1
                top.pid[ip1:ip2,top.tpid-1]=top.pid[ip1:ip2,top.tpid-1]-top.dt
              npguess = nps_save[js] + len(ii)
              nplost = ins_save[js] - top.ins[js] - top.nps[js]
              chckpart(js+1,0,npguess + nplost,false)
              copypart(ins_save[js]+nps_save[js],len(ii),ii,-1)
              nps_save[js] = nps_save[js] + len(ii)

      npssum = sum(parallelsum(top.nps))
      maxvz = parallelmax(top.vzmaxp)

    # --- Print a blank line so that the previous status line is not written
    # --- over.
    if lstatusline: print ''

    # --- The call to perrho3d is primarily needed for the parallel version.
    if(w3d.solvergeom==w3d.RZgeom and rhoparam is not None):
      frz.l_distribute = false
    perrho3d(w3d.rho,w3d.nx,w3d.ny,w3d.nz,w3d.bound0,w3d.boundxy)
    if(w3d.solvergeom==w3d.RZgeom and rhoparam is not None):
      frz.l_distribute = true

    # --- If rhoparam is not None, mix in the previous rho with the
    # --- new rho
    if rhoparam is not None:
     if(w3d.solvergeom<>w3d.RZgeom): 
      w3d.rho[:,:,:] = (1.-rhoparam)*w3d.rho + rhoparam*rhoprevious
     else:
      frz.distribute_rho_rz()
      for ig in range(frz.ngrids):
        if(ig==0):
          g = frz.basegrid
        else:
          try:
            g = g.next
          except:
            g = g.down
        g.rho = (1. - rhoparam)*g.rho + rhoparam*rhoprevious[ig]
#       mix_rho_rz(rhoprevious[ig],frz.nrg[ig],frz.nzg[ig],ig+1,rhoparam)  

    # --- Do field solve including newly accumulated charge density.
    top.fstype = _ofstype
    fieldsol(-1)
    top.fstype = -1

    # --- Do final work for zmoments calculation
    if ((i == iter-1 or (gun_iter%nhist) == 0) and top.ifzmmnt > 0):
      top.laccumulate_zmoments = laccumulate_zmoments
      getzmom.zmmnt(3)

    # --- Save the history data
    if nhist>0:
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

    gun.steps = gun_steps
    gun.iter = gun_iter
    gun.time = gun_time

    print "Number of iterations = %d"%gun_iter

    # --- call insertafteriter if defined
    if insertafteriter is not None:
       insertafteriter()

    # --- update egundata arrays
    if egundata_nz is not None:
      zz = zd/w3d.dz
      iz = int(zz)
      dz = zz-iz
      egundata_curr.append((1.-dz)*take(top.curr[...,-1],iz)+dz*take(top.curr[...,-1],iz+1))
      egundata_xrmsz.append((1.-dz)*take(top.xrmsz[...,-1],iz)+dz*take(top.xrmsz[...,-1],iz+1))
      egundata_yrmsz.append((1.-dz)*take(top.yrmsz[...,-1],iz)+dz*take(top.yrmsz[...,-1],iz+1))
      egundata_xprmsz.append((1.-dz)*take(top.xprmsz[...,-1],iz)+dz*take(top.xprmsz[...,-1],iz+1))
      egundata_yprmsz.append((1.-dz)*take(top.yprmsz[...,-1],iz)+dz*take(top.yprmsz[...,-1],iz+1))
      egundata_epsnxz.append((1.-dz)*take(top.epsnxz[...,-1],iz)+dz*take(top.epsnxz[...,-1],iz+1))
      egundata_epsnyz.append((1.-dz)*take(top.epsnyz[...,-1],iz)+dz*take(top.epsnyz[...,-1],iz+1))
      _izdata += 1
  
    # plot egundata 
    if egundata_nz is not None and egundata_window>-1:
      window(egundata_window)
      fma()
      plsys(3)
      # --- The data is plotted this way since it is a list of arrays.
      plg([x[-2] for x in egundata_curr])
      ptitles('Current','Z','',v=3) 
      plsys(4)
      plg([x[-2] for x in egundata_xrmsz])
      plg([x[-2] for x in egundata_yrmsz],color='red')
      ptitles('X, Y RMS','Z','',v=4)
      plsys(5)
      plg([x[-2] for x in egundata_xprmsz])
      plg([x[-2] for x in egundata_yprmsz],color='red')
      ptitles("X', Y' RMS",'Z','',v=5)
      plsys(6)
      plg([x[-2] for x in egundata_epsnxz])
      plg([x[-2] for x in egundata_epsnyz],color='red')
      ptitles('X, Y norm. emittance','Z','',v=6)
      pyg_pending()
      pyg_idler()
      window(0)

  # --- end of multiple iterations

  # --- Set saved particles to be live particles (for diagnostics only).
  if (_ipstep > 0):
    top.ins[:] = ins_save
    top.nps[:] = nps_save

  # --- Change what is plotted at the bottom of each frame
  stepid(gun_iter,gun_time,top.zbeam)

  # --- Set some additional diagnostic data
  if top.nzzarr == top.nzmmnt: top.vzofz[:] = top.vzbarz[:,-1]

  # --- Return variables to their original values
  top.fstype = _ofstype
  if (top.ntinj > 0):
    top.nztinjmn = _onztinjmn
    top.nztinjmx = _onztinjmx
  top.inject = _oinject
  top.nhist = nhist
  top.ifzmmnt = _ifzmmnt
  top.laccumulate_zmoments = _laccumulate_zmoments

  # --- install plottraces
  if plottraces_window>-1:
    uninstallafterstep(plottraces)

  if egundata_nz is not None:
    return  [array(egundata_curr),
             array(egundata_xrmsz), array(egundata_yrmsz), 
             array(egundata_xprmsz), array(egundata_yprmsz),
             array(egundata_epsnxz), array(egundata_epsnyz)]


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

########################################################################
def gunmg(iter=1,itersub=None,ipsave=None,save_same_part=None,maxtime=None,
        laccumulate_zmoments=None,rhoparam=None,
        lstatusline=true,insertbeforeiter=None,insertafteriter=None,
        nmg=0,conductors=None,egundata_window=2,plottraces_window=1,
        egundata_nz=None,egundata_zmin=None,egundata_zmax=None):
  """
Performs steady-state iterations in a cascade using different resolutions.
  - iter=1 number of iterations to perform
  - itersub=1 number of iterations to perform at coarser levels (=iter by default)
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
  - lstatusline=1: when try, a line is printed and continuously updated
                   showing the status of the simulation.
  - insertbeforeiter=None: function to be called before each iteration.
  - insertafteriter=None: function to be called after each iteration.
  - nmg = 0: number of 'multigrid' levels (in addition to main level).
  - conductors: list of conductors.
  - egundata_window=2: window in which to display egundata curves for z close 
                       to w3d.zmmax. Set to a negative number to deactivate 
                       plotting.
  - plottraces_window=1: window in which to plot traces.Set to a negative 
                         number to deactivate plotting.
  Note that ipsave and save_same_part are preserved in between calls
  """
  global rhonext,nrnex,nznext,dxnext,dznext
  global zd, egundata_curr, egundata_xrmsz, egundata_yrmsz
  global egundata_xprmsz, egundata_yprmsz, egundata_epsnxz, egundata_epsnyz
  # if nmg=0, do a normal gun solve
  if(nmg==0):
    return gun(iter,ipsave,save_same_part,maxtime,
        laccumulate_zmoments,rhoparam,
        lstatusline,insertbeforeiter,insertafteriter,
        None,egundata_window,plottraces_window,
        egundata_nz,egundata_zmin,egundata_zmax)
  # initialize itersub
  if itersub is None: itersub = iter
  iterlast = iter
  if(itersub==0):raise('Error in gunmg: called with iter=0 or itersub=0')
  # Define function that makes charge deposition on grid at next level.
  # This is used needed on the last iteration at a given level so that
  # the charge density is transmitted to the next level.
  def setrhonext():
    global rhonext,nrnex,nznext,dxnext,dznext
    frz.dep_rho_rz(1,rhonext,nrnext,nznext,drnext,dznext,0.,w3d.zmmin)
  # declare arrays containing size grids for each mg level
  gunnx = zeros(nmg+1)
  gunny = zeros(nmg+1)
  gunnz = zeros(nmg+1)
  gundt = zeros(nmg+1,'d')
  gunnpinject = zeros(nmg+1)
  # save last level sizes
  i = nmg
  gunnx[nmg] = w3d.nx
  gunny[nmg] = w3d.ny
  gunnz[nmg] = w3d.nz
  gundt[nmg] = top.dt
  gunnpinject[nmg] = top.npinject
  # compute sublevels sizes 
  for i in range(nmg-1,-1,-1):
    gunnx[i] = int(gunnx[i+1]/2)
    gunny[i] = int(gunny[i+1]/2)
    gunnz[i] = int(gunnz[i+1]/2)
    gundt[i] = gundt[i+1]*2
    gunnpinject[i] = int(gunnpinject[i+1]/2)
  if(w3d.solvergeom<>w3d.RZgeom):
    raise('function not yet implemented')
  else:
    # mg loop
    for i in range(0,nmg+1):
      if(i==nmg):
        iter = iterlast
      else:
        iter = itersub
      print 'gunmg: level %g on %g'%(i+1,nmg+1)
      # reset a few variables so that the weights of particles is
      # computed according to the current grid resolution
      top.dt = gundt[i]
      swprev = top.sw*1.
      top.sw[:] = 0.
      top.npinje_s[:] = 0
      top.npinject = gunnpinject[i]
      # resize the mesh and associated arrays
      adjustmesh3d.resizemesh(gunnx[i],0,gunnz[i],0,0,1,1,1,conductors)
      # Except at the coarser level (i=0), the charge density is
      # copied from rhonext where it was stored from at last iteration
      # at the previous level.
      if(i>0):
        frz.basegrid.rho[:,:] = rhonext[:,:]
        rhoprevious = [rhonext]
      # update phi and inj_phi
      fieldsol(-1)
      getinj_phi()
      # Set inj_param=1 (or almost) when starting at a new level since
      # inj_prev has been redimensioned. Do interpolation of old inj_prev
      # to new inj_prev in the future?
      if(i>0):
        top.inj_param = 0.9999999
      else:
        top.inj_param = 0.5
      # reset particle arrays
      top.np_s = 0
      top.npmax_s = 0
      top.npmax = 1
      alotpart()
      # performs all iterations but the last one
      if(iter>1):
        # first iteration is performed with top.inj_param=1 (except i=0)
        gun(1,0,save_same_part,maxtime,
            laccumulate_zmoments,rhoparam,
            lstatusline,insertbeforeiter,insertafteriter,
            None,egundata_window,plottraces_window,
            egundata_nz,egundata_zmin,egundata_zmax)
        # remaining iterations but last one performed with inj_param=0.5
        top.inj_param = 0.5
        if(iter>2):
          gun(iter-2,0,save_same_part,maxtime,
              laccumulate_zmoments,rhoparam,
              lstatusline,insertbeforeiter,insertafteriter,
              None,egundata_window,plottraces_window,
              egundata_nz,egundata_zmin,egundata_zmax)
      # For all sublevels, rhonext is created and setrhonext is installed
      # so that rhonext is eveluated during last iteration at current level.
      if i<nmg:
         nrnext = gunnx[i+1]
         nznext = gunnz[i+1]
         drnext = (w3d.xmmax-w3d.xmmin)/nrnext         
         dznext = (w3d.zmmax-w3d.zmmin)/nznext
         rhonext = fzeros([nrnext+1,nznext+1],'d')
         installafterstep(setrhonext)
      # perform last iteration
      gun(1,ipsave,save_same_part,maxtime,
          laccumulate_zmoments,rhoparam,
          lstatusline,insertbeforeiter,insertafteriter,
          None,egundata_window,plottraces_window,
          egundata_nz,egundata_zmin,egundata_zmax)
      # Uninstall setrhonext if necessary.
      if i<nmg:
         uninstallafterstep(setrhonext)

  if egundata_nz is not None:
    return [array(egundata_curr),array(egundata_xrmsz),array(egundata_yrmsz), 
             array(egundata_xprmsz),array(egundata_yprmsz),
             array(egundata_epsnxz),array(egundata_epsnyz)]

########################################################################
def gunamr(iter=1,itersub=None,ipsave=None,save_same_part=None,maxtime=None,
        nmg=0,AMRlevels=0,
        laccumulate_zmoments=None,rhoparam=None,
        lstatusline=true,insertbeforeiter=None,insertafteriter=None,
        conductors=None,ipstep=None,egundata_window=-1,plottraces_window=-1,
        egundata_nz=None,egundata_zmin=None,egundata_zmax=None):
  """
Performs steady-state iterations in a cascade using different resolutions.
  - iter=1 number of iterations to perform
  - itersub=1 number of iterations to perform at coarser levels (=iter by default)
  - ipsave=5000000 number of particles to save from the last iteration
  - save_same_part=0 when true, save same particles each time step instead
    of random particles, i.e. saves particle trajectories
  - maxtime=3*transittime maximum time each iteration will run
  - laccumulate_zmoments=false: When set to true, z-moments are accumulated
    over multiple iterations. Note that getzmom.zmmnt(3) must be called
    by the user to finish the moments calculation.
  - rhoparam=None: Amount of previous rho to mix in with the current rho. This
    can help the relaxation toward a steady state. Caution should be used
    when using this option.
  - lstatusline=1: when try, a line is printed and continuously updated
                   showing the status of the simulation.
  - insertbeforeiter=None: function to be called before each iteration.
  - insertafteriter=None: function to be called after each iteration.
  - nmg = 0: number of 'multigrid' levels (in addition to main level).
  - conductors: list of conductors.
  - egundata_window=2: window in which to display egundata curves for z close 
                       to w3d.zmmax. Set to a negative number to deactivate 
                       plotting.
  - plottraces_window=1: window in which to plot traces.Set to a negative 
                         number to deactivate plotting.
  Note that ipsave and save_same_part are preserved in between calls
  """
  global zd, egundata_curr, egundata_xrmsz, egundata_yrmsz
  global egundata_xprmsz, egundata_yprmsz, egundata_epsnxz, egundata_epsnyz
  if nmg>0:
   gunmg(itersub,itersub,ipsave,save_same_part,maxtime,
         laccumulate_zmoments,rhoparam,
         lstatusline,insertbeforeiter,insertafteriter,
         nmg,conductors,egundata_window,plottraces_window,
         egundata_nz,egundata_zmin,egundata_zmax)
  else:
   gun(itersub,ipsave,save_same_part,maxtime,
            laccumulate_zmoments,rhoparam,
            lstatusline,insertbeforeiter,insertafteriter,
            None,egundata_window,plottraces_window,
            egundata_nz,egundata_zmin,egundata_zmax)
  if AMRlevels>0:
    w3d.AMRlevels = AMRlevels
    fieldsol()
    tmp = w3d.AMRgenerate_periodicity 
    w3d.AMRgenerate_periodicity = 1
    AMRtree = __main__.__dict__['AMRtree']
    if conductors is not None:
      AMRtree.conductors += conductors
    AMRtree.generate()
    w3d.AMRgenerate_periodicity = 1000000
    fieldsol(-1)
    if rhoparam is not None:
      if iter == 1:
        ipsavetemp = ipsave
        ipsteptemp = ipstep
      else:
        ipsavetemp = None
        ipsteptemp = None
      gun(1,ipsavetemp,save_same_part,maxtime,
          laccumulate_zmoments,None,
          lstatusline,insertbeforeiter,insertafteriter,
          ipsteptemp,egundata_window,plottraces_window,
          egundata_nz,egundata_zmin,egundata_zmax)
      iter = iter - 1
    if iter > 0:
      gun(iter,ipsave,save_same_part,maxtime,
          laccumulate_zmoments,rhoparam,
          lstatusline,insertbeforeiter,insertafteriter,
          None,egundata_window,plottraces_window,
          egundata_nz,egundata_zmin,egundata_zmax)
    w3d.AMRgenerate_periodicity = tmp
  if egundata_nz is not None:
    return [array(egundata_curr),array(egundata_xrmsz),array(egundata_yrmsz), 
             array(egundata_xprmsz),array(egundata_yprmsz),
             array(egundata_epsnxz),array(egundata_epsnyz)]

########################################################################
def statusline():
  """
Prints a running line showing current status of the step.
  """
  if (top.it % 10) == 0:
    CR = curses.ascii.ctrl('m')
    sys.stdout.write("%5d "%top.it)
    nplive = sum(parallelsum(top.nps))
    sys.stdout.write("nplive = %5d "%nplive)
    zz = top.zp[top.ins[0]-1]
    if zz < w3d.zmminglobal: zz = w3d.zmmaxglobal
    sys.stdout.write("zz = %6.4f"%(zz))
    sys.stdout.write(CR)

