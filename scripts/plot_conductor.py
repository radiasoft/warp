from warp import *
plot_conductor_version = "$Id: plot_conductor.py,v 1.2 2000/11/21 19:54:56 dave Exp $"

######################################################################
# functions to plot the conductor points and subgrid data            #
# in MKS units                                                       #
######################################################################

# --- Convenience function to plot the sub-grid data
def plotsubgrid(izf,nn,ix,iy,iz,delmx,delmy,delpx,delpy,xmin,ymin,dx,dy,color):
  ii = compress(equal(iz[:nn],izf),arange(nn))
  xx = take(ix,ii)*dx+xmin
  yy = take(iy,ii)*dy+ymin
  delmx = take(delmx,ii)*dx
  delmy = take(delmy,ii)*dy
  delpx = take(delpx,ii)*dx
  delpy = take(delpy,ii)*dy
  if lparallel:
    xx = gatherarray(xx)
    yy = gatherarray(yy)
    delmx = gatherarray(delmx)
    delmy = gatherarray(delmy)
    delpx = gatherarray(delpx)
    delpy = gatherarray(delpy)
  for i in xrange(len(xx)):
    if (abs(delmx[i]) < abs(dx)):
      plg([xx[i],xx[i]-delmx[i]],[yy[i],yy[i]],color=color)
    if (abs(delmy[i]) < abs(dy)):
      plg([xx[i],xx[i]],[yy[i],yy[i]-delmy[i]],color=color)
    if (abs(delpx[i]) < abs(dx)):
      plg([xx[i],xx[i]+delpx[i]],[yy[i],yy[i]],color=color)
    if (abs(delpy[i]) < abs(dy)):
      plg([xx[i],xx[i]],[yy[i],yy[i]+delpy[i]],color=color)

# x-y plane
def pfxy(izf=None,contours=None,plotsg=1,scale=1,signx=1,signy=1,
         plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in X-Y plane
  - izf=w3d.iz_axis z index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signx=1 sign of x, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  # --- This logic is needed since in the parallel version, iz_axis already
  # --- has izslave subtracted from it. If the user passes in a value,
  # --- it must be checked for consistency, otherwise coding below could lead
  # --- to a deadlock in the parallel version
  if izf == None:
    izf = w3d.iz_axis
  else:
    if izf < 0 or w3d.nzfull < izf: return
    if lparallel:
      izf = izf - top.izslave[me+1]
  if scale:
    dx = w3d.dx*signx
    dy = w3d.dy*signy
    xmmin = w3d.xmmin
    ymmin = w3d.ymmin
  else:
    dx = 1.*signx
    dy = 1.*signy
    xmmin = 0.
    ymmin = 0.
  if plotphi:
    xx=iota(0,w3d.nx)*dx + xmmin
    yy=iota(0,w3d.ny)*dy + ymmin
    ppp = getphi(iz=izf+top.izslave[me+1])
    ppp = transpose(ppp)
    plotc(ppp,yy,xx,contours=contours,filled=filled,color=blue)
  ii = compress(equal(f3d.izcond[0:f3d.ncond],izf),arange(f3d.ncond))
  yy = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
  xx = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
  if lparallel:
    yy = gatherarray(yy)
    xx = gatherarray(xx)
  if xx:
    plp(yy,xx,color=cyan)
  if (plotsg):
    plotsubgrid(izf,f3d.necndbdy,f3d.iecndy,f3d.iecndx,f3d.iecndz,
                f3d.ecdelmy,f3d.ecdelmx,f3d.ecdelpy,f3d.ecdelpx,
                ymmin,xmmin,dy,dx,green)
    plotsubgrid(izf,f3d.nocndbdy,f3d.iocndy,f3d.iocndx,f3d.iocndz,
                f3d.ocdelmy,f3d.ocdelmx,f3d.ocdelpy,f3d.ocdelpx,
                ymmin,xmmin,dy,dx,red)

# z-x plane
def pfzx(iyf=None,contours=None,plotsg=1,scale=1,signz=1,signx=1,
         plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-X plane
  - iyf=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signx=1 sign of x, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if scale:
    dx = w3d.dx*signx
    dz = w3d.dz*signz
    xmmin = w3d.xmmin
    zmmin = w3d.zmmin
  else:
    dx = 1.*signx
    dz = 1.*signz
    xmmin = 0.
    zmmin = 0.
  if iyf == None: iyf = w3d.iy_axis
  if plotphi:
    xx = iota(0,w3d.nx)*dx + xmmin
    zz = iota(0,w3d.nzfull)*dz + zmmin
    ppp = getphi(iy=iyf)
    plotc(ppp,xx,zz,contours=contours,filled=filled,color=blue)
  ii = compress(equal(f3d.iycond[0:f3d.ncond],iyf),arange(f3d.ncond))
  xx = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
  zz = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  if lparallel:
    xx = gatherarray(xx)
    zz = gatherarray(zz)
  if xx:
    plp(xx,zz,color=cyan)
  if (plotsg):
    plotsubgrid(iyf,f3d.necndbdy,f3d.iecndx,f3d.iecndz,f3d.iecndy,
                f3d.ecdelmx,f3d.ecdelmz,f3d.ecdelpx,f3d.ecdelpz,
                xmmin,zmmin,dx,dz,green)
    plotsubgrid(iyf,f3d.nocndbdy,f3d.iocndx,f3d.iocndz,f3d.iocndy,
                f3d.ocdelmx,f3d.ocdelmz,f3d.ocdelpx,f3d.ocdelpz,
                xmmin,zmmin,dx,dz,red)

# z-y plane
def pfzy(ixf=None,contours=None,plotsg=1,scale=1,signz=1,signy=1,
         plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane
  - ixf=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if scale:
    dy = w3d.dy*signy
    dz = w3d.dz*signz
    ymmin = w3d.ymmin
    zmmin = w3d.zmmin
  else:
    dy = 1.*signy
    dz = 1.*signz
    ymmin = 0.
    zmmin = 0.
  if ixf == None: ixf = w3d.ix_axis
  if plotphi:
    yy = iota(0,w3d.ny)*dy + ymmin
    zz = iota(0,w3d.nzfull)*dz + zmmin
    ppp = getphi(ix=ixf)
    plotc(ppp,yy,zz,contours=contours,filled=filled,color=blue)
  ii = compress(equal(f3d.ixcond[0:f3d.ncond],ixf),arange(f3d.ncond))
  yy = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
  zz = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  if lparallel:
    yy = gatherarray(yy)
    zz = gatherarray(zz)
  if yy:
    plp(yy,zz,color=cyan)
  if (plotsg):
    plotsubgrid(ixf,f3d.necndbdy,f3d.iecndy,f3d.iecndz,f3d.iecndx,
                f3d.ecdelmy,f3d.ecdelmz,f3d.ecdelpy,f3d.ecdelpz,
                ymmin,zmmin,dy,dz,green)
    plotsubgrid(ixf,f3d.nocndbdy,f3d.iocndy,f3d.iocndz,f3d.iocndx,
                f3d.ocdelmy,f3d.ocdelmz,f3d.ocdelpy,f3d.ocdelpz,
                ymmin,zmmin,dy,dz,red)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in grid units                                                      #
######################################################################

# x-y plane
def pfxyg(izf=None,contours=None,plotsg=1,signx=1,signy=1,plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in X-Y plane in grid
frame
  - izf=w3d.iz_axis z index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - signx=1 sign of x, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  pfxy(izf=izf,contours=contours,plotsg=plotsg,scale=0,signx=signx,signy=signy,
       plotphi=plotphi,filled=filled)

# z-x plane
def pfzxg(iyf=None,contours=None,plotsg=1,signz=1,signx=1,plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-X plane in grid
frame
  - iyf=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - signz=1 sign of z, used for plotting symmetry planes
  - signx=1 sign of x, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  pfzx(iyf=iyf,contours=contours,plotsg=plotsg,scale=0,signz=signz,signx=signx,
       plotphi=plotphi,filled=filled)

# z-y plane
def pfzyg(ixf=None,contours=None,plotsg=1,signz=1,signy=1,plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane in grid
frame
  - ixf=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - signz=1 sign of z, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  pfzy(ixf=ixf,contours=contours,plotsg=plotsg,scale=0,signz=signz,signy=signy,
       plotphi=plotphi,filled=filled)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in real units, with inverted x.  Used to make complete plots when  #
# 2-fold symmetry is used.                                           #
######################################################################
 
# z-x plane
def pfzxi(iyf=None,contours=None,plotsg=1,scale=1,signz=1,plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-(-X) plane
  - iyf=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  pfzx(iyf=iyf,contours=contours,plotsg=plotsg,scale=scale,signz=signz,signx=-1,
       plotphi=plotphi,filled=filled)

# z-y plane
def pfzyi(ixf=None,contours=None,plotsg=1,scale=1,signz=1,plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-(-Y) plane
  - ixf=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  pfzy(ixf=ixf,contours=contours,plotsg=plotsg,scale=scale,signz=signz,signy=-1,
       plotphi=plotphi,filled=filled)

# x-y plane
def pfxyi(izf=None,contours=None,plotsg=1,scale=1,plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in X-Y plane,
plotting data in all four quadrants
  - izf=w3d.iz_axis z index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  pfxy(izf=izf,contours=contours,plotsg=plotsg,scale=scale,signx=+1,signy=+1,
       plotphi=plotphi,filled=filled)
  if w3d.l2symtry or w3d.l4symtry:
    pfxy(izf=izf,contours=contours,plotsg=plotsg,scale=scale,signx=+1,signy=-1,
         plotphi=plotphi,filled=filled)
  if w3d.l4symtry:
    pfxy(izf=izf,contours=contours,plotsg=plotsg,scale=scale,signx=-1,signy=+1,
         plotphi=plotphi,filled=filled)
    pfxy(izf=izf,contours=contours,plotsg=plotsg,scale=scale,signx=-1,signy=-1,
         plotphi=plotphi,filled=filled)

############################################################################
# These plot abox at each conductor point
############################################################################

# x-y plane
def pfxybox(izf=None,contours=None,plotsg=1,scale=1,signx=1,signy=1,
            plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in X-Y plane
  - izf=w3d.iz_axis z index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signx=1 sign of x, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if scale:
    dy = w3d.dy*signy
    dx = w3d.dx*signx
    ymmin = w3d.ymmin
    xmmin = w3d.xmmin
  else:
    dy = 1.*signy
    dx = 1.*signx
    ymmin = 0.
    xmmin = 0.
  if not izf: izf = w3d.iz_axis
  if plotphi:
    yy=iota(0,w3d.ny)[:,NewAxis]*ones(w3d.nx+1,'d')*dy + ymmin
    xx=iota(0,w3d.nx)*ones(w3d.ny+1,'d')[:,NewAxis]*dx + xmmin
    ireg=ones((w3d.ny+1,w3d.nx+1))
    ppp = getphi(iz=izf+top.izslave[me+1])
    ppp = transpose(ppp)
    if filled:
      if contours:
        plfc(ppp,yy,xx,ireg,contours=contours)
      else:
        plfc(ppp,yy,xx,ireg)
    else:
      if contours:
        plc(ppp,yy,xx,ireg,color='blue',contours=contours)
      else:
        plc(ppp,yy,xx,ireg,color='blue')
  if (f3d.ncond > 0):
    ii = compress(equal(f3d.izcond[0:f3d.ncond],izf),arange(f3d.ncond))
    if ii:
      y = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
      x = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
      pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
          array([x-dx/2,x+dx/2,x+dx/2,x-dx/2,x-dx/2]),
          color=cyan)

# z-x plane
def pfzxbox(iyf=None,contours=None,plotsg=1,scale=1,signz=1,signx=1,
            plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-X plane
  - iyf=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signx=1 sign of x, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if scale:
    dx = w3d.dx*signx
    dz = w3d.dz*signz
    xmmin = w3d.xmmin
    zmmin = w3d.zmmin
  else:
    dx = 1.*signx
    dz = 1.*signz
    xmmin = 0.
    zmmin = 0.
  if not iyf: iyf = w3d.iy_axis
  if plotphi:
    xx=iota(0,w3d.nx)[:,NewAxis]*ones(w3d.nz+1,'d')*dx + xmmin
    zz=iota(0,w3d.nz)*ones(w3d.nx+1,'d')[:,NewAxis]*dz + zmmin
    ireg=ones((w3d.nx+1,w3d.nz+1))
    ppp = getphi(iy=iyf)
    if filled:
      if contours:
        plfc(ppp,xx,zz,ireg,contours=contours)
      else:
        plfc(ppp,xx,zz,ireg)
    else:
      if contours:
        plc(ppp,xx,zz,ireg,color='blue',contours=contours)
      else:
        plc(ppp,xx,zz,ireg,color='blue')
  if (f3d.ncond > 0):
    ii = compress(equal(f3d.iycond[0:f3d.ncond],iyf),arange(f3d.ncond))
    if ii:
      x = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
      z = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
      pla(array([x-dx/2,x-dx/2,x+dx/2,x+dx/2,x-dx/2]),
          array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
          color=cyan)

# z-y plane
def pfzybox(ixf=None,contours=None,plotsg=1,scale=1,signz=1,signy=1,
            plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-Y plane
  - ixf=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if scale:
    dy = w3d.dy*signy
    dz = w3d.dz*signz
    ymmin = w3d.ymmin
    zmmin = w3d.zmmin
  else:
    dy = 1.*signy
    dz = 1.*signz
    ymmin = 0.
    zmmin = 0.
  if not ixf: ixf = w3d.ix_axis
  if plotphi:
    yy=iota(0,w3d.ny)[:,NewAxis]*ones(w3d.nz+1,'d')*dy + ymmin
    zz=iota(0,w3d.nz)*ones(w3d.ny+1,'d')[:,NewAxis]*dz + zmmin
    ireg=ones((w3d.ny+1,w3d.nz+1))
    ppp = getphi(ix=ixf)
    if filled:
      if contours:
        plfc(ppp,yy,zz,ireg,contours=contours)
      else:
        plfc(ppp,yy,zz,ireg)
    else:
      if contours:
        plc(ppp,yy,zz,ireg,color='blue',contours=contours)
      else:
        plc(ppp,yy,zz,ireg,color='blue')
  if (f3d.ncond > 0):
    ii = compress(equal(f3d.ixcond[0:f3d.ncond],ixf),arange(f3d.ncond))
    if ii:
      y = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
      z = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
      pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
          array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
          color=cyan)

# z-x plane
def pfzxboxi(izf=None,contours=None,plotsg=1,scale=1,signz=1,
             plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-(-X) plane
  - iyf=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  pfzxbox(iyf=iyf,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
          signx=-1,plotphi=plotphi,filled=filled)

# z-y plane
def pfzyboxi(ixf=None,contours=None,plotsg=1,scale=1,signz=1,signy=-1,
             plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-(-Y) plane
  - ixf=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - filled=0 when true, plots filled contours
  """
  pfzybox(ixf=ixf,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
          signy=-1,plotphi=plotphi,filled=filled)






############################################################################
# These plot the conductors in laboratory frame, using the tolabfrm routine
# to convert from code frame to lab frame.  There is also a routine to plot
# the computational x-z grid in lab frame.
############################################################################

# --- plot grid in lab frame (including bends)
def plotgrid(zz=None,ii=2,plotcond=1):
  """
Plots Z-X grid in the lab frame (including bends)
  - zz=top.zbeam is the center position
  - ii=2 is the step size in the grid points plotted
    2 means that every other grid line is plotted
  - plotcond=1 when true, plots conductors
  """
  if not zz: zz=top.zbeam
  # --- declare temporary data space, 2 2-D arrays to hold grid coordinates
  xxx = zeros((w3d.nx/ii+1,w3d.nz/ii+1),'d')
  zzz = zeros((w3d.nx/ii+1,w3d.nz/ii+1),'d')
  for iz in xrange(w3d.nz/ii+1):
    xxx[:,iz] = w3d.xmesh[::ii]
  for ix in xrange(w3d.nx/ii+1):
    zzz[ix,:] = w3d.zmmin + iota(0,w3d.nz,ii)*w3d.dz + w3d.zz

  # --- If in a bend, convert the grid data to the lab frame
  if top.linbend:
    # --- reshape arrays to make 1-D arrays to pass to tolabfrm
    nn = int((w3d.nx/ii+1)*(w3d.nz/ii+1))
    xxx.shape = (nn)
    zzz.shape = (nn)

    # --- Convert data to lab frame
    tolabfrm(zz,nn,xxx,zzz)

    # --- Reshape back into 2-D arrays
    xxx.shape = (w3d.nx/ii+1,w3d.nz/ii+1)
    zzz.shape = (w3d.nx/ii+1,w3d.nz/ii+1)

  # --- Make plots
  pla(xxx,zzz,marks=0)
  pla(transpose(xxx),transpose(zzz),marks=0)

  if plotcond: pfzxlab(zz)

# --- Make pfzx plot in lab frame
def pfzxlab(zz=None):
  """Plots conductors in Z-X lab frame (including bends)
  - zz=top.zbeam is the center position
  """
  if not zz: zz=top.zbeam
  # --- if zz is not equal to zbeam, then calculate conductors for new location
  if (zz != top.zbeam):
    z = top.zbeam
    g = top.zgrid
    top.zbeam = zz
    top.zgrid = zz
    setlatt()
    fieldsol(1)
  # --- gather conductor data
  xxxx=compress(equal(f3d.iycond[0:f3d.ncond],iyf),
                      f3d.ixcond[0:f3d.ncond])*w3d.dx+w3d.xmmin
  zzzz=compress(equal(f3d.iycond[0:f3d.ncond],iyf),
                      f3d.izcond[0:f3d.ncond])*w3d.dz+w3d.zmmin+zz
  # --- convert to lab frame
  tolabfrm(zz,len(xxxx),xxxx,zzzz)   
  # --- make plot
  plg(xxxx,zzzz,marker='\2',color=cyan)
  # --- restore original conductor data at zbeam
  if (zz != top.zbeam):
    top.zbeam = z
    top.zgrid = g
    setlatt()
    fieldsol(1)

