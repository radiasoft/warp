from warp import *
from appendablearray import *
singleparticle_version = "$Id: singleparticle.py,v 1.1 2000/10/16 18:34:19 dave Exp $"

# Setups and run the code in single particle mode.

_spissetup = 0
def spsetup(zerophi=1):
  """Sets up for a single particle run: turns diagnostics off, saves some
initial data, and redefines the step command"""
  global _spissetup
  global _itsave,_timesave,_zbeamsave,_zgridsave,_zgridprvsave
  global _vbeamsave,_vbeamfrmsave,_dtsave,_dssave
  _spissetup = 1
  # --- Turn off all of the diagnostics
  wxy.ldiag = 0
  top.ifzmmnt = 0
  top.itmomnts = 0
  top.itplps = 0
  top.itplfreq = 0
  top.zzmomnts = 0
  top.zzplps = 0
  top.zzplfreq = 0
  top.nhist = top.nt
  top.iflabwn = 0
  w3d.lrhodia3d = false
  # --- Turn off the field solver
  top.fstype = -1
  # --- Zero out phi if requested.
  if zerophi:
    w3d.phi = 0.
  # --- Turn off charge deposition
  top.depos = 'none'
  # --- Save initial time-step data
  _itsave = top.it
  _timesave = top.time
  _zbeamsave = top.zbeam
  _zgridsave = top.zgrid
  _zgridprvsave = top.zgridprv
  _vbeamsave = top.vbeam
  _vbeamfrmsave = top.vbeamfrm
  _dtsave = top.dt
  _dssave = wxy.ds
  # --- Add routine after step to save data
  afterstep.append(spsavedata)

# --- Initialize the single particle.
def spinit(x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=None,gi=1.,maxsteps=1000,savedata=1):
  "Initializes one of more particles to run independently"
  global _spx,_spy,_spz,_spvx,_spvy,_spvz,_spgi,_nn
  global _itsave,_timesave,_zbeamsave,_zgridsave,_zgridprvsave
  global _vbeamsave,_vbeamfrmsave,_dtsave,_dssave
  global _savedata
  if not _spissetup: spsetup()
  _savedata = savedata
  # --- Set default value of vz
  if not vz:
    vz = top.vbeam
  if type(x) == type(0.) or type(x) == type(0):
    x = [x]
    y = [y]
    z = [z]
    vx = [vx]
    vy = [vy]
    vz = [vz]
    gi = [gi]
  _nn = len(x)
  # --- make sure there is space
  chckpart(0,0,_nn,false)
  # --- Initialize particle coordinates
  top.ins[0] = 1
  top.nps[0] = _nn
  top.xp[0:_nn] = x
  top.yp[0:_nn] = y
  top.zp[0:_nn] = z
  top.uxp[0:_nn] = vx
  top.uyp[0:_nn] = vy
  top.uzp[0:_nn] = vz
  top.gaminv[0:_nn] = gi
  top.sq[0] = top.zion*top.echarge
  top.sm[0] = top.aion*top.amu
  top.sw[0] = 0.
  # --- Initialize time stepping
  top.it = _itsave
  top.time = _timesave
  top.zbeam = _zbeamsave
  top.zgrid = _zgridsave
  top.zgridprv = _zgridprvsave
  top.vbeam = _vbeamsave
  top.vbeamfrm = _vbeamfrmsave
  top.dt = _dtsave
  wxy.ds = _dssave
  # --- Setup the lattice
  setlatt()
  # --- Create arrays for saving trajectory
  if savedata:
    _spx = []
    _spy = []
    _spz = []
    _spvx = []
    _spvy = []
    _spvz = []
    _spgi = []
    for i in xrange(_nn):
      _spx.append(AppendableArray(maxsteps,'d'))
      _spy.append(AppendableArray(maxsteps,'d'))
      _spz.append(AppendableArray(maxsteps,'d'))
      _spvx.append(AppendableArray(maxsteps,'d'))
      _spvy.append(AppendableArray(maxsteps,'d'))
      _spvz.append(AppendableArray(maxsteps,'d'))
      _spgi.append(AppendableArray(maxsteps,'d'))
    for i in xrange(_nn):
      _spx[i].append(x[i])
      _spy[i].append(y[i])
      _spz[i].append(z[i])
      _spvx[i].append(vx[i])
      _spvy[i].append(vy[i])
      _spvz[i].append(vz[i])
      _spgi[i].append(gi[i])

def spsavedata():
  """Saves data"""
  global _spx,_spy,_spz,_spvx,_spvy,_spvz,_spgi
  if not _savedata: return
  for i in xrange(_nn):
    _spx[i].append(top.xp[i])
    _spy[i].append(top.yp[i])
    _spz[i].append(top.zp[i])
    _spvx[i].append(top.uxp[i])
    _spvy[i].append(top.uyp[i])
    _spvz[i].append(top.uzp[i])
    _spgi[i].append(top.gaminv[i])

def getsavedata():
  "Retrieves the save particle trajectories on a single array"
  global _spx,_spy,_spz,_spvx,_spvy,_spvz,_spgi
  if _nn == 1:
    return array((_spx[0].data(),_spy[0].data(),_spz[0].data(),_spvx[0].data(),_spvy[0].data(),_spvz[0].data(),_spgi[0].data()))
  else:
    r = zeros((_nn,7,_spx[0].len()),'d')
    for i in xrange (_nn):
      r[i,0,:] = _spx[i].data()
      r[i,1,:] = _spy[i].data()
      r[i,2,:] = _spz[i].data()
      r[i,3,:] = _spvx[i].data()
      r[i,4,:] = _spvy[i].data()
      r[i,5,:] = _spvz[i].data()
      r[i,6,:] = _spgi[i].data()
    return r

