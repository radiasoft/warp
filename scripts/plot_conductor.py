from warp import *
plot_conductor_version = "$Id: plot_conductor.py,v 1.15 2001/06/08 21:11:35 dave Exp $"

def plot_conductordoc():
  print """
The following functions plot contours of the potential in various planes
along with the conductors in that plane. The first three plot with the axis
having units of meters. The suffix 'g' means that it plots with the axis
having units of number of grid cells. The suffix 'i' means that it replicates
the plots in all quadrants (when symmetry is used). The 'box' suffix means
that it plots a box around each grid point inside of a conductor.

pfxy, pfzx, pfzy
pfxyg, pfzxg, pfzyg
pfzxi, pfzyi, pfxyi
pfxybox, pfzxbox, pfzybox, pfzxboxi, pfzyboxi

plotgrid: plots the x-z mesh in the lab frame (including any bends)
pfzxlab: makes the pfzx plot in the lab frame (including any bends)
plotsrfrv: handy command to plot r versus z for a suface of revolution, giving
           the function describing it

plotquadoutline: plots outline of quadrupole structure
  """

######################################################################
# functions to plot the conductor points and subgrid data            #
# in MKS units                                                       #
######################################################################

# --- Convenience function to plot the sub-grid data
def plotsubgrid(iz,nn,ixc,iyc,izc,delmx,delmy,delpx,delpy,xmin,ymin,dx,dy,
                color):
  ii = compress(equal(izc[:nn],iz),arange(nn))
  xx = take(ixc,ii)*dx+xmin
  yy = take(iyc,ii)*dy+ymin
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
def pfxy(iz=None,izf=None,contours=None,plotsg=1,scale=1,signx=1,signy=1,
         plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in X-Y plane
  - iz=w3d.iz_axis z index of plane
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
  if izf != None: iz = izf
  if iz == None:
    iz = w3d.iz_axis
  else:
    if iz < 0 or w3d.nzfull < iz: return
    if lparallel:
      iz = iz - top.izslave[me]
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
    #xx=iota(0,w3d.nx)*dx + xmmin
    #yy=iota(0,w3d.ny)*dy + ymmin
    ppp = getphi(iz=iz+top.izslave[me])
    #ppp = transpose(ppp)
    #plotc(ppp,yy,xx,contours=contours,filled=filled,color=blue)
    ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=blue,
              xmin=xmmin,xmax=xmmin+w3d.nx*dx,
              ymin=ymmin,ymax=ymmin+w3d.ny*dy)
  if f3d.ncond > 0:
    ii = compress(equal(f3d.izcond[0:f3d.ncond],iz),arange(f3d.ncond))
    yy = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
    xx = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    warpplp(yy,xx,color=cyan)
  else:
    warpplp([],[],color=cyan)
  if (plotsg):
    plotsubgrid(iz,f3d.necndbdy,f3d.iecndy,f3d.iecndx,f3d.iecndz,
                f3d.ecdelmy,f3d.ecdelmx,f3d.ecdelpy,f3d.ecdelpx,
                ymmin,xmmin,dy,dx,green)
    plotsubgrid(iz,f3d.nocndbdy,f3d.iocndy,f3d.iocndx,f3d.iocndz,
                f3d.ocdelmy,f3d.ocdelmx,f3d.ocdelpy,f3d.ocdelpx,
                ymmin,xmmin,dy,dx,red)

# z-x plane
def pfzx(iy=None,iyf=None,contours=None,plotsg=1,scale=1,signz=1,signx=1,
         plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-X plane
  - iy=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signx=1 sign of x, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if iyf != None: iy = iyf
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
    if lparallel: zmmin = top.izslave[me]
  if iy == None: iy = w3d.iy_axis
  if plotphi:
    #xx = iota(0,w3d.nx)*dx + xmmin
    #zz = iota(0,w3d.nzfull)*dz + zmmin
    ppp = getphi(iy=iy)
    #plotc(ppp,xx,zz,contours=contours,filled=filled,color=blue)
    ppp = transpose(ppp)
    ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=blue,
              xmin=zmmin,xmax=zmmin+w3d.nzfull*dz,
              ymin=xmmin,ymax=xmmin+w3d.nx*dx)
  if f3d.ncond > 0:
    ii = compress(equal(f3d.iycond[0:f3d.ncond],iy),arange(f3d.ncond))
    xx = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    zz = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
    warpplp(xx,zz,color=cyan)
  else:
    warpplp([],[],color=cyan)
  if (plotsg):
    plotsubgrid(iy,f3d.necndbdy,f3d.iecndx,f3d.iecndz,f3d.iecndy,
                f3d.ecdelmx,f3d.ecdelmz,f3d.ecdelpx,f3d.ecdelpz,
                xmmin,zmmin,dx,dz,green)
    plotsubgrid(iy,f3d.nocndbdy,f3d.iocndx,f3d.iocndz,f3d.iocndy,
                f3d.ocdelmx,f3d.ocdelmz,f3d.ocdelpx,f3d.ocdelpz,
                xmmin,zmmin,dx,dz,red)

# z-y plane
def pfzy(ix=None,ixf=None,contours=None,plotsg=1,scale=1,signz=1,signy=1,
         plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane
  - ix=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if ixf != None: ix = ixf
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
    if lparallel: zmmin = top.izslave[me]
  if ix == None: ix = w3d.ix_axis
  if plotphi:
    #yy = iota(0,w3d.ny)*dy + ymmin
    #zz = iota(0,w3d.nzfull)*dz + zmmin
    ppp = getphi(ix=ix)
    #plotc(ppp,yy,zz,contours=contours,filled=filled,color=blue)
    ppp = transpose(ppp)
    ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=blue,
              xmin=zmmin,xmax=zmmin+w3d.nzfull*dz,
              ymin=ymmin,ymax=ymmin+w3d.ny*dy)
  if f3d.ncond > 0:
    ii = compress(equal(f3d.ixcond[0:f3d.ncond],ix),arange(f3d.ncond))
    yy = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
    zz = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
    warpplp(yy,zz,color=cyan)
  else:
    warpplp([],[],color=cyan)
  if (plotsg):
    plotsubgrid(ix,f3d.necndbdy,f3d.iecndy,f3d.iecndz,f3d.iecndx,
                f3d.ecdelmy,f3d.ecdelmz,f3d.ecdelpy,f3d.ecdelpz,
                ymmin,zmmin,dy,dz,green)
    plotsubgrid(ix,f3d.nocndbdy,f3d.iocndy,f3d.iocndz,f3d.iocndx,
                f3d.ocdelmy,f3d.ocdelmz,f3d.ocdelpy,f3d.ocdelpz,
                ymmin,zmmin,dy,dz,red)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in grid units                                                      #
######################################################################

# x-y plane
def pfxyg(iz=None,izf=None,contours=None,plotsg=1,signx=1,signy=1,plotphi=1,
          filled=0):
  """
Plots conductors and contours of electrostatic potential in X-Y plane in grid
frame
  - iz=w3d.iz_axis z index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - signx=1 sign of x, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if izf != None: iz = izf
  pfxy(iz=iz,contours=contours,plotsg=plotsg,scale=0,signx=signx,signy=signy,
       plotphi=plotphi,filled=filled)

# z-x plane
def pfzxg(iy=None,iyf=None,contours=None,plotsg=1,signz=1,signx=1,plotphi=1,
          filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-X plane in grid
frame
  - iy=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - signz=1 sign of z, used for plotting symmetry planes
  - signx=1 sign of x, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if iyf != None: iy = iyf
  pfzx(iy=iy,contours=contours,plotsg=plotsg,scale=0,signz=signz,signx=signx,
       plotphi=plotphi,filled=filled)

# z-y plane
def pfzyg(ix=None,ixf=None,contours=None,plotsg=1,signz=1,signy=1,plotphi=1,
          filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane in grid
frame
  - ix=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - signz=1 sign of z, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if ixf != None: ix = ixf
  pfzy(ix=ix,contours=contours,plotsg=plotsg,scale=0,signz=signz,signy=signy,
       plotphi=plotphi,filled=filled)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in real units, with inverted x.  Used to make complete plots when  #
# 2-fold symmetry is used.                                           #
######################################################################
 
# z-x plane
def pfzxi(iy=None,iyf=None,contours=None,plotsg=1,scale=1,signz=1,plotphi=1,
          filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-(-X) plane
  - iy=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if iyf != None: iy = iyf
  pfzx(iy=iy,contours=contours,plotsg=plotsg,scale=scale,signz=signz,signx=-1,
       plotphi=plotphi,filled=filled)

# z-y plane
def pfzyi(ix=None,ixf=None,contours=None,plotsg=1,scale=1,signz=1,plotphi=1,
          filled=0):
  """
Plots conductors and contours of electrostatic potential in Z-(-Y) plane
  - ix=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if ixf != None: ix = ixf
  pfzy(ix=ix,contours=contours,plotsg=plotsg,scale=scale,signz=signz,signy=-1,
       plotphi=plotphi,filled=filled)

# x-y plane
def pfxyi(iz=None,izf=None,contours=None,plotsg=1,scale=1,plotphi=1,filled=0):
  """
Plots conductors and contours of electrostatic potential in X-Y plane,
plotting data in all four quadrants
  - iz=w3d.iz_axis z index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if izf != None: iz = izf
  pfxy(iz=iz,contours=contours,plotsg=plotsg,scale=scale,signx=+1,signy=+1,
       plotphi=plotphi,filled=filled)
  if w3d.l2symtry or w3d.l4symtry:
    pfxy(iz=iz,contours=contours,plotsg=plotsg,scale=scale,signx=+1,signy=-1,
         plotphi=plotphi,filled=filled)
  if w3d.l4symtry:
    pfxy(iz=iz,contours=contours,plotsg=plotsg,scale=scale,signx=-1,signy=+1,
         plotphi=plotphi,filled=filled)
    pfxy(iz=iz,contours=contours,plotsg=plotsg,scale=scale,signx=-1,signy=-1,
         plotphi=plotphi,filled=filled)

############################################################################
# These plot abox at each conductor point
############################################################################

# x-y plane
def pfxybox(iz=None,izf=None,contours=None,plotsg=1,scale=1,signx=1,signy=1,
            plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in X-Y plane
  - iz=w3d.iz_axis z index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signx=1 sign of x, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if izf != None: iz = izf
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
  if not iz: iz = w3d.iz_axis
  if plotphi:
    yy=iota(0,w3d.ny)[:,NewAxis]*ones(w3d.nx+1,'d')*dy + ymmin
    xx=iota(0,w3d.nx)*ones(w3d.ny+1,'d')[:,NewAxis]*dx + xmmin
    ireg=ones((w3d.ny+1,w3d.nx+1))
    ppp = getphi(iz=iz+top.izslave[me])
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
  if f3d.ncond > 0:
    ii = compress(equal(f3d.izcond[0:f3d.ncond],iz),arange(f3d.ncond))
    x = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    y = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
  else:
    x = []
    y = []
  if lparallel:
    x = gatherarray(x)
    y = gatherarray(y)
  if len(x) > 0:
    pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
        array([x-dx/2,x+dx/2,x+dx/2,x-dx/2,x-dx/2]),
        color=cyan)

# z-x plane
def pfzxbox(iy=None,iyf=None,contours=None,plotsg=1,scale=1,signz=1,signx=1,
            plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-X plane
  - iy=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signx=1 sign of x, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if iyf != None: iy = iyf
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
    if lparallel: zmmin = top.izslave[me]
  if not iy: iy = w3d.iy_axis
  if plotphi:
    xx=iota(0,w3d.nx)[:,NewAxis]*ones(w3d.nz+1,'d')*dx + xmmin
    zz=iota(0,w3d.nz)*ones(w3d.nx+1,'d')[:,NewAxis]*dz + zmmin
    ireg=ones((w3d.nx+1,w3d.nz+1))
    ppp = getphi(iy=iy)
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
    ii = compress(equal(f3d.iycond[0:f3d.ncond],iy),arange(f3d.ncond))
    x = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    z = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  else:
    x = []
    z = []
  if lparallel:
    x = gatherarray(x)
    z = gatherarray(z)
  if len(x) > 0:
    pla(array([x-dx/2,x-dx/2,x+dx/2,x+dx/2,x-dx/2]),
        array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
        color=cyan)

# z-y plane
def pfzybox(ix=None,ixf=None,contours=None,plotsg=1,scale=1,signz=1,signy=1,
            plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-Y plane
  - ix=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if ixf != None: ix = ixf
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
    if lparallel: zmmin = top.izslave[me]
  if not ix: ix = w3d.ix_axis
  if plotphi:
    yy=iota(0,w3d.ny)[:,NewAxis]*ones(w3d.nz+1,'d')*dy + ymmin
    zz=iota(0,w3d.nz)*ones(w3d.ny+1,'d')[:,NewAxis]*dz + zmmin
    ireg=ones((w3d.ny+1,w3d.nz+1))
    ppp = getphi(ix=ix)
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
    ii = compress(equal(f3d.ixcond[0:f3d.ncond],ix),arange(f3d.ncond))
    y = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
    z = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  else:
    y = []
    z = []
  if lparallel:
    y = gatherarray(y)
    z = gatherarray(z)
  if len(y) > 0:
    pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
        array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
        color=cyan)

# z-x plane
def pfzxboxi(iy=None,iyf=None,contours=None,plotsg=1,scale=1,signz=1,
             plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-(-X) plane
  - iy=w3d.iy_axis y index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  """
  if iyf != None: iy = iyf
  pfzxbox(iy=iy,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
          signx=-1,plotphi=plotphi,filled=filled)

# z-y plane
def pfzyboxi(ix=None,ixf=None,contours=None,plotsg=1,scale=1,signz=1,signy=-1,
             plotphi=1,filled=0):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-(-Y) plane
  - ix=w3d.ix_axis x index of plane
  - contours optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - filled=0 when true, plots filled contours
  """
  if ixf != None: ix = ixf
  pfzybox(ix=ix,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
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
def pfzxlab(zz=None,iy=None):
  """Plots conductors in Z-X lab frame (including bends)
  - zz=top.zbeam is the center position
  """
  if not zz: zz=top.zbeam
  if iy == None: iy = w3d.iy_axis
  # --- if zz is not equal to zbeam, then calculate conductors for new location
  if (zz != top.zbeam):
    z = top.zbeam
    g = top.zgrid
    top.zbeam = zz
    top.zgrid = zz
    setlatt()
    fieldsol(1)
  # --- gather conductor data
  if f3d.ncond > 0:
    xxxx=compress(equal(f3d.iycond[0:f3d.ncond],iy),
                        f3d.ixcond[0:f3d.ncond])*w3d.dx+w3d.xmmin
    zzzz=compress(equal(f3d.iycond[0:f3d.ncond],iy),
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


#####################################################################
def plotsrfrv(srfrv,zmin,zmax,n=1000,color='fg',gridframe=0,rscale=1,zscale=1,
              roff=0,zoff=0,rmin=0.,rmax=top.largepos):
  """Handy function for plotting the r versus z for a surface of revolution
 - srfrv: surface of revolution function to plot
 - zmin,zmax: z range to plot
 - n=1000: number of points to plot
 - color='fg': color of line
 - gridframe=0: when true, plots in grid frame
 - rscale=1: scaling for radius
 - zscale=1: scaling for z
 - roff=0: offset for radius
 - zoff=0: offset for z
 - rmin=0: minimum value of r plotted (before applying rscale and roff)
 - rmax=0: maximum value of r plotted (before applying rscale and roff)
  """
  zz = iota(0,n)*(zmax - zmin)/n + zmin
  rr = ones(n+1,'d')
  for i in range(n+1):
    f3d.srfrv_z = zz[i]
    srfrv()
    rr[i] = f3d.srfrv_r
  if gridframe:
    zz = (zz - w3d.zmmin)/w3d.dz
    rr = (rr)/w3d.dx
  rr = where(less(rr,rmin),rmin,rr)
  rr = where(greater(rr,rmax),rmax,rr)
  plg(rscale*rr+roff,zscale*zz+zoff,color=color)


#####################################################################
def plotelementoutline(color,gridframe,axis,iquad,nquad,
                       ezs,eze,eap,err,erl,egl,egp,eox,eoy,epa,epr,epw,
                       dpal,dpar):
  """Plots the outline of electrostatic elements
  - color: line color
  - gridframe: when true, make plot in grid coordinates
  - axis: selects axis to plot, either 'x' or 'y'
  """
  if axis == 'x': gpsign = 1
  else:           gpsign = -1
  for i in range(iquad,iquad+nquad+1):
    # --- plot rods
    # --- If aperture is zero, then this quad is skipped
    rodap = eap[i]
    if erl[i] > 0.:
      rodlen = erl[i]
      gp = egp[i]*gpsign
      gaplen = egl[i]
    else:
      rodlen = (eze[i] - ezs[i])
      gp = 1*gpsign
      gaplen = 0.
    if err[i] > 0.: rodrr = err[i]
    else:           rodrr = 8./7.*eap[i]
    if axis == 'x': offset = eox[i]
    else:           offset = eoy[i]
    if rodap > 0. and rodlen > 0.:
      rr = rodap + rodrr + rodrr*array([1.,1.,-1.,-1.,1.])
      zz = gp*(-0.5*(rodlen+gaplen) + rodlen*array([0.,1.,1.,0.,0.]))
      rr1 = offset + rr
      rr2 = offset - rr
      zz = 0.5*(eze[i] + ezs[i]) + top.zlatstrt + zz
      if gridframe:
        rr1 = rr1/w3d.dx
        rr2 = rr2/w3d.dx
        zz = (zz - w3d.zmmin)/w3d.dz
      plg(rr1,zz,color=color)
      plg(rr2,zz,color=color)
    # --- Plot end plates
    pw = epw[i]
    if pw > 0.:
      if epa[i] > 0.: pa = epa[i]
      else:           pa = eap[i]
      if epr[i] > 0.: pr = epr[i]
      else:           pr = rodap + 2.*rodrr
      pal = pa + dpal[i]
      par = pa + dpar[i]
      rrl = array([pr,pr,pal,pal,pr])
      rrr = array([pr,pr,par,par,pr])
      zz = pw*array([0.,1.,1.,0.,0.])
      rrl1 = offset + rrl
      rrl2 = offset - rrl
      rrr1 = offset + rrr
      rrr2 = offset - rrr
      zzl = 0.5*(eze[i] + ezs[i]) - 0.5*(rodlen+gaplen) - zz + \
            top.zlatstrt
      zzr = 0.5*(eze[i] + ezs[i]) + 0.5*(rodlen+gaplen) + zz + \
            top.zlatstrt
      if gridframe:
        rrl1 = rrl1/w3d.dx
        rrl2 = rrl2/w3d.dx
        rrr1 = rrr1/w3d.dx
        rrr2 = rrr2/w3d.dx
        zzl = (zzl - w3d.zmmin)/w3d.dz
        zzr = (zzr - w3d.zmmin)/w3d.dz
      plg(rrl1,zzl,color=color)
      plg(rrl2,zzl,color=color)
      plg(rrr1,zzr,color=color)
      plg(rrr2,zzr,color=color)


#####################################################################
def plotquadoutline(iquad=0,nquad=None,color='fg',gridframe=0,axis='x'):
  """Plots the outline of quadrupole elements
  - iquad=0: starting quad to plot
  - nquad=top.nquad: number of quads to plot
  - color='fg': line color
  - gridframe=0: when true, make plot in grid coordinates
  - axis='x': selects axis to plot, either 'x' or 'y'
  """
  if nquad == None: nquad = top.nquad
  plotelementoutline(color,gridframe,axis,iquad,nquad,
                     top.quadzs,top.quadze,top.quadap,top.quadrr,top.quadrl,
                     top.quadgl,top.quadgp,top.qoffx,top.qoffy,
                     top.quadpa,top.quadpr,top.quadpw,
                     top.qdelpal,top.qdelpar)


