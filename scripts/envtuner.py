from warp import *
envtuner_version = "$Id: envtuner.py,v 1.2 2003/08/18 17:40:56 dave Exp $"

def makeenvtuneplot(zl,zu,qz):
  # --- Make the envelope plot, setting plot limits.
  plg(env.aenv,env.zenv,marks=0)
  plg(env.benv,env.zenv,color='red',marks=0)
  envamin = min(env.aenv)
  envbmin = min(env.benv)
  envmin = min(envamin,envbmin)
  envamax = max(env.aenv)
  envbmax = max(env.benv)
  envmax = max(envamax,envbmax)
  plotlimits = (zl,zu,envmin-0.0*envmin,envmax+0.0*envmax,0)
  p=limits(plotlimits)
  p=plotlimits
  # --- Plot vertical guiding lines between quads
  for z in qz:
    plg(array([p[2],p[3]]),array([z,z]),color='white',type='dot',marks=0)
  return p


def envtuner(zl=env.zl,zu=env.zu,scale=1.,mindelta=1.e-3):
  # --- Find the quads within the range.
  iquads = 0
  iquade = 0
  while iquads < top.nquad and top.quadze[iquads]+top.zlatstrt < zl:
    iquads = iquads + 1
  while iquade < top.nquad and top.quadzs[iquade]+top.zlatstrt < zu:
    iquade = iquade + 1
  qz = zeros(iquade - iquads + 1,'d')
  qz[1:] = (top.zlatstrt + 
            0.5*(top.quadze[iquads:iquade] + top.quadzs[iquads+1:iquade+1]))
  if qz[-1] > zu: qz[-1] = zu
  if iquads == 0:
    qz[0] = zl
  else:
    qz[0] = top.zlatstrt + 0.5*(top.quadze[iquads-1] + top.quadzs[iquads])

  # --- Make initial plot
  p = makeenvtuneplot(zl,zu,qz)

  # --- Now, enter the loop.
  done = 0
  while not done:
    m = mouse(1,2,'')
    # --- If null was returned or if right button was clicked, the quit.
    if not m or m[9] == 3:
      done = 1
      continue
    # --- If middle button if pushed, then clear plot.
    if m[9] == 2:
      fma()
      print 'dedx = ' + repr(top.quadde[iquads:iquade+1])
      print 'dbdx = ' + repr(top.quaddb[iquads:iquade+1])
      continue
    # --- Get z position and which quad it is in.
    z = m[0]
    for iq in xrange(len(qz)):
      if qz[iq] < z and z < qz[iq+1]: break
    # --- Get value of delta requested.
    delta = (m[3] - m[1])/(p[3] - p[2])
    if delta == 0:
      if m[1] < 0.5*(p[2]+p[3]):
	delta = -mindelta
      else:
	delta = +mindelta
    # --- Change quad strength appropriately.
    top.quadde[iquads+iq] = top.quadde[iquads+iq]*(1. + delta*scale)
    top.quaddb[iquads+iq] = top.quaddb[iquads+iq]*(1. + delta*scale)
    # --- Recalculate envelope and redo the plot.
    step()
    p = makeenvtuneplot(zl,zu,qz)

