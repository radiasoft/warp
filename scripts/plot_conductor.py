from warp import *
plot_conductor_version = "$Id: plot_conductor.py,v 1.28 2002/05/02 17:57:08 dave Exp $"

def plot_conductordoc():
  print """
The following functions plot contours of the potential in various planes
along with the conductors in that plane. The first three plot with the axis
having units of meters. The suffix 'g' means that it plots with the axis
having units of number of grid cells. The 'box' suffix means
that it plots a box around each grid point inside of a conductor.

pfxy, pfzx, pfzy
pfxyg, pfzxg, pfzyg
pfxybox, pfzxbox, pfzybox, pfzxboxi, pfzyboxi

plotgrid: plots the x-z mesh in the lab frame (including any bends)
pfzxlab: makes the pfzx plot in the lab frame (including any bends)
plotsrfrv: handy command to plot r versus z for a suface of revolution, giving
           the function describing it

plotquadoutline: plots outline of quadrupole structure

cleanconductors: not a plot routine, buts removes conductor points not
		 within the the range of the field solve
  """

######################################################################
# functions to plot the conductor points and subgrid data            #
# in MKS units                                                       #
######################################################################

# --- Convenience function to plot the sub-grid data
def plotcond(yy,xx,zz,iz,numb,ymin,xmin,dy,dx,color,lx,ly,lz,signy,signx):
  nn = f3d.ncond
  if nn > 0:
    ixc = eval('f3d.i'+xx+'cond')*signx*lx
    iyc = eval('f3d.i'+yy+'cond')*signy*ly
    izc = eval('f3d.i'+zz+'cond')*lz
    cnumb = f3d.condnumb
    try:
      if xx in ['x','y']: levx = f3d.icondlxy
      else:               levx = f3d.icondlz
      if yy in ['x','y']: levy = f3d.icondlxy
      else:               levy = f3d.icondlz
      if zz in ['x','y']: levz = f3d.icondlxy
      else:               levz = f3d.icondlz
      level = 5*equal(levx,lx) - 2*equal(levy,ly) - 2*equal(levz,lz)
    except:
      level = ones(f3d.ncondmax)
  else:
    ixc = array([])
    iyc = array([])
    izc = array([])
    cnumb = array([])
    level = array([])
  ii = compress(logical_and(equal(izc[:nn],iz),equal(level[:nn],1)),arange(nn))
  xx = take(ixc,ii)*dx+xmin
  yy = take(iyc,ii)*dy+ymin
  if numb is not None: cnumb = take(cnumb,ii)
  warpplp(yy,xx,color=color)

def plotsubgrid(yy,xx,zz,pp,iz,numb,ymin,xmin,dy,dx,color,subgridlen,lx,ly,lz,
                signy,signx):
  assert (pp == 'e' or pp == 'o'),"pp has invalid data"
  nn = eval('f3d.n'+pp+'cndbdy')
  if nn > 0:
    ixc = eval('f3d.i'+pp+'cnd'+xx)*signx
    iyc = eval('f3d.i'+pp+'cnd'+yy)*signy
    izc = eval('f3d.i'+pp+'cnd'+zz)*lz
    delmx = eval('f3d.'+pp+'cdelm'+xx)*signx
    delpx = eval('f3d.'+pp+'cdelp'+xx)*signx
    delmy = eval('f3d.'+pp+'cdelm'+yy)*signy
    delpy = eval('f3d.'+pp+'cdelp'+yy)*signy
    try:
      numbmx = eval('f3d.'+pp+'cnumbm'+xx)
      numbpx = eval('f3d.'+pp+'cnumbp'+xx)
      numbmy = eval('f3d.'+pp+'cnumbm'+yy)
      numbpy = eval('f3d.'+pp+'cnumbp'+yy)
    except:
      numbmx = eval('f3d.'+pp+'cnumb')
      numbpx = eval('f3d.'+pp+'cnumb')
      numbmy = eval('f3d.'+pp+'cnumb')
      numbpy = eval('f3d.'+pp+'cnumb')
    try:
      if xx in ['x','y']: levx = eval('f3d.i'+pp+'cndlxy')
      else:               levx = eval('f3d.i'+pp+'cndlz')
      if yy in ['x','y']: levy = eval('f3d.i'+pp+'cndlxy')
      else:               levy = eval('f3d.i'+pp+'cndlz')
      if zz in ['x','y']: levz = eval('f3d.i'+pp+'cndlxy')
      else:               levz = eval('f3d.i'+pp+'cndlz')
      level = 5*equal(levx,lx) - 2*equal(levy,ly) - 2*equal(levz,lz)
    except:
      level = ones(f3d.ncondmax)
  else:
    ixc = array([])
    iyc = array([])
    izc = array([])
    delmx = array([])
    delpx = array([])
    delmy = array([])
    delpy = array([])
    numbmx = array([])
    numbpx = array([])
    numbmy = array([])
    numbpy = array([])
    level = array([])
  ii = compress(logical_and(equal(izc[:nn],iz),equal(level[:nn],1)),arange(nn))
  dx = dx*lx
  dy = dy*ly
  xx = take(ixc,ii)*dx+xmin
  yy = take(iyc,ii)*dy+ymin
  delmx = take(delmx,ii)*dx
  delpx = take(delpx,ii)*dx
  delmy = take(delmy,ii)*dy
  delpy = take(delpy,ii)*dy
  if numb is not None:
    numbmx = take(numbmx,ii)
    numbpx = take(numbpx,ii)
    numbmy = take(numbmy,ii)
    numbpy = take(numbpy,ii)
  if lparallel:
    xx = gatherarray(xx)
    yy = gatherarray(yy)
    delmx = gatherarray(delmx)
    delpx = gatherarray(delpx)
    delmy = gatherarray(delmy)
    delpy = gatherarray(delpy)
  # --- This code combines all of the individual lines into the list pp.
  # --- This vectorized code avoids slower explicit loops.
  niota = arange(len(xx))
  ii = compress(less(abs(delmy),dy*subgridlen),niota)
  pp = map(lambda x,y,d:[y,y-d,x,x],take(xx,ii),take(yy,ii),take(delmy,ii))
  ii = compress(less(abs(delmx),dx*subgridlen),niota)
  pp = pp + map(lambda x,y,d:[y,y,x,x-d],take(xx,ii),take(yy,ii),take(delmx,ii))
  ii = compress(less(abs(delpy),dy*subgridlen),niota)
  pp = pp + map(lambda x,y,d:[y,y+d,x,x],take(xx,ii),take(yy,ii),take(delpy,ii))
  ii = compress(less(abs(delpx),dx*subgridlen),niota)
  pp = pp + map(lambda x,y,d:[y,y,x,x+d],take(xx,ii),take(yy,ii),take(delpx,ii))
  # --- Convert the list to an array and plot it
  if len(pp) > 0:
    pp = array(pp)
    pldj(pp[:,2],pp[:,0],pp[:,3],pp[:,1],color=color)

# x-y plane
def pfxy(iz=None,izf=None,fullplane=1,plotsg=1,scale=1,
         plotphi=1,subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,lxy=1,lz=1,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in X-Y plane
  - iz=w3d.iz_axis z index of plane
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  - lxy=1,lz=1: level of multigrid to plot data for
  """
  kw.update(kwdict)
  # --- This logic is needed since in the parallel version, iz_axis already
  # --- has izslave subtracted from it. If the user passes in a value,
  # --- it must be checked for consistency, otherwise coding below could lead
  # --- to a deadlock in the parallel version
  if izf is not None: iz = izf
  if iz is None: iz = w3d.iz_axis + top.izslave[me]
  if iz < 0 or w3d.nzfull < iz: return
  izlocal = iz - top.izslave[me]
  if scale:
    dx = w3d.dx
    dy = w3d.dy
    xmmin = w3d.xmmin
    ymmin = w3d.ymmin
    xmmax = w3d.xmmax
    ymmax = w3d.ymmax
  else:
    dx = 1.
    dy = 1.
    xmmin = 0.
    ymmin = 0.
    xmmax = w3d.nx
    ymmax = w3d.ny
  if plotphi:
    if not scale:
      kw['xmin'] = 0
      kw['xmax'] = w3d.nx
      kw['ymin'] = 0
      kw['ymax'] = w3d.ny
    if not kw.has_key('ccolor'): kw['ccolor'] = phicolor
    apply(pcphixy,(iz,fullplane),kw)
  izl = izlocal
  plotcond('y','x','z',izl,numb,ymmin,xmmin,dy,dx,condcolor,lxy,lxy,lz,1,1)
  if fullplane and (w3d.l2symtry or w3d.l4symtry):
    plotcond('y','x','z',izl,numb,ymmin,xmmin,dy,dx,condcolor,lxy,lxy,lz,-1,1)
  if fullplane and w3d.l4symtry:
    plotcond('y','x','z',izl,numb,ymmin,xmmin,dy,dx,condcolor,lxy,lxy,lz,1,-1)
    plotcond('y','x','z',izl,numb,ymmin,xmmin,dy,dx,condcolor,lxy,lxy,lz,-1,-1)
  if (plotsg):
    plotsubgrid('y','x','z','e',izlocal,numb,ymmin,xmmin,dy,dx,evencolor,
                subgridlen,lxy,lxy,lz,1,1)
    plotsubgrid('y','x','z','o',izlocal,numb,ymmin,xmmin,dy,dx,oddcolor,
                subgridlen,lxy,lxy,lz,1,1)
    if fullplane and (w3d.l2symtry or w3d.l4symtry):
      plotsubgrid('y','x','z','e',izlocal,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,lxy,lxy,lz,1,-1)
      plotsubgrid('y','x','z','o',izlocal,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,lxy,lxy,lz,1,-1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('y','x','z','e',izlocal,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,lxy,lxy,lz,-1,1)
      plotsubgrid('y','x','z','o',izlocal,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,lxy,lxy,lz,-1,1)
      plotsubgrid('y','x','z','e',izlocal,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,lxy,lxy,lz,-1,-1)
      plotsubgrid('y','x','z','o',izlocal,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,lxy,lxy,lz,-1,-1)

# z-x plane
def pfzx(iy=None,iyf=None,fullplane=1,lbeamframe=1,plotsg=1,scale=1,
         plotphi=1,subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,lxy=1,lz=1,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in Z-X plane
  - iy=w3d.iy_axis y index of plane
  - fullplane=1: when true, plots all quadrants regardless of symmetries
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  - lxy=1,lz=1: level of multigrid to plot data for
  """
  kw.update(kwdict)
  if iyf is not None: iy = iyf
  if iy is None: iy = w3d.iy_axis
  if iy < 0 or w3d.ny < iy: return
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if scale:
    dx = w3d.dx
    dz = w3d.dz
    xmmin = w3d.xmmin
    zmmin = w3d.zmmin + zbeam
    xmmax = w3d.xmmax
    zmmax = w3d.zmmax + zbeam
  else:
    dx = 1.
    dz = 1.
    xmmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
    xmmax = w3d.nx
    zmmax = w3d.nz
  if plotphi:
    if not scale:
      kw['xmin'] = 0
      kw['xmax'] = w3d.nzfull
      kw['ymin'] = 0.
      kw['ymax'] = w3d.nx
    if not kw.has_key('ccolor'): kw['ccolor'] = phicolor
    apply(pcphizx,(iy,fullplane,lbeamframe),kw)
  plotcond('x','z','y',iy,numb,xmmin,zmmin,dx,dz,condcolor,lxy,lz,lxy,1,1)
  if fullplane and w3d.l4symtry:
    plotcond('x','z','y',iy,numb,xmmin,zmmin,dx,dz,condcolor,lxy,lz,lxy,-1,1)
  if (plotsg):
    plotsubgrid('x','z','y','e',iy,numb,xmmin,zmmin,dx,dz,evencolor,
                subgridlen,lxy,lz,lxy,1,1)
    plotsubgrid('x','z','y','o',iy,numb,xmmin,zmmin,dx,dz,oddcolor,
                subgridlen,lxy,lz,lxy,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('x','z','y','e',iy,numb,xmmin,zmmin,dx,dz,evencolor,
                  subgridlen,lxy,lz,lxy,-1,1)
      plotsubgrid('x','z','y','o',iy,numb,xmmin,zmmin,dx,dz,oddcolor,
                  subgridlen,lxy,lz,lxy,-1,1)

# z-y plane
def pfzy(ix=None,ixf=None,fullplane=1,lbeamframe=1,plotsg=1,scale=1,
         plotphi=1,subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,lxy=1,lz=1,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane
  - ix=w3d.ix_axis x index of plane
  - fullplane=1: when true, plots all quadrants regardless of symmetries
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  - lxy=1,lz=1: level of multigrid to plot data for
  """
  kw.update(kwdict)
  if ixf is not None: ix = ixf
  if ix is None: ix = w3d.ix_axis
  if ix < 0 or w3d.nx < ix: return
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if scale:
    dy = w3d.dy
    dz = w3d.dz
    ymmin = w3d.ymmin
    zmmin = top.zplmin + zbeam
    ymmax = w3d.ymmax
    zmmax = top.zplmax + zbeam
  else:
    dy = 1.
    dz = 1.
    ymmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
    ymmax = w3d.ny
    zmmax = w3d.nz
  if plotphi:
    if not scale:
      kw['xmin'] = 0
      kw['xmax'] = w3d.nzfull
      kw['ymin'] = 0
      kw['ymax'] = w3d.ny
    if not kw.has_key('ccolor'): kw['ccolor'] = phicolor
    apply(pcphizy,(ix,fullplane,lbeamframe),kw)
  plotcond('y','z','x',ix,numb,ymmin,zmmin,dy,dz,condcolor,lxy,lz,lxy,1,1)
  if fullplane and (w3d.l2symtry or w3d.l4symtry):
    plotcond('y','z','x',ix,numb,ymmin,zmmin,dy,dz,condcolor,lxy,lz,lxy,-1,1)
  if (plotsg):
    plotsubgrid('y','z','x','e',ix,numb,ymmin,zmmin,dy,dz,evencolor,
                subgridlen,lxy,lz,lxy,1,1)
    plotsubgrid('y','z','x','o',ix,numb,ymmin,zmmin,dy,dz,oddcolor,
                subgridlen,lxy,lz,lxy,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('y','z','x','e',ix,numb,ymmin,zmmin,dy,dz,evencolor,
                  subgridlen,lxy,lz,lxy,-1,1)
      plotsubgrid('y','z','x','o',ix,numb,ymmin,zmmin,dy,dz,oddcolor,
                  subgridlen,lxy,lz,lxy,-1,1)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in grid units                                                      #
######################################################################

# x-y plane
def pfxyg(iz=None,izf=None,fullplane=1,plotsg=1,plotphi=1,
          phicolor=blue,subgridlen=1.,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,lxy=1,lz=1,**kw):
  """
Plots conductors and contours of electrostatic potential in X-Y plane in grid
frame
Same arguments as pfxy
  """
  if izf is not None: iz = izf
  pfxy(iz=iz,fullplane=fullplane,plotsg=plotsg,scale=0,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,lxy=lxy,lz=lz,kwdict=kw)

# z-x plane
def pfzxg(iy=None,iyf=None,fullplane=1,lbeamframe=1,plotsg=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,lxy=1,lz=1,**kw):
  """
Plots conductors and contours of electrostatic potential in Z-X plane in grid
frame
Same arguments as pfzx
  """
  if iyf is not None: iy = iyf
  pfzx(iy=iy,fullplane=fullplane,lbeamframe=lbeamframe,plotsg=plotsg,scale=0,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,lxy=lxy,lz=lz,kwdict=kw)

# z-y plane
def pfzyg(ix=None,ixf=None,fullplane=1,lbeamframe=1,plotsg=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,lxy=1,lz=1,**kw):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane in grid
frame
Same arguments as pfzy
  """
  if ixf is not None: ix = ixf
  pfzy(ix=ix,fullplane=fullplane,lbeamframe=lbeamframe,plotsg=plotsg,scale=0,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,lxy=lxy,lz=lz,kwdict=kw)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in real units, with inverted x.  Used to make complete plots when  #
# 2-fold symmetry is used.                                           #
######################################################################
 
# x-y plane
def pfxyi(iz=None,izf=None,fullplane=1,plotsg=1,scale=1,plotphi=1,
          phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full X-Y plane,
Same arguments as pfxy
  """
  print "Notice: pfxyi is obsolete is should no longer be used"
  print "        It does the identical thing as pfxy"
  if izf is not None: iz = izf
  pfxy(iz=iz,fullplane=fullplane,plotsg=plotsg,scale=scale,
       plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

# z-x plane
def pfzxi(iy=None,iyf=None,fullplane=1,lbeamframe=1,plotsg=1,scale=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full Z-X plane
Same arguments as pfzx
  """
  print "Notice: pfzxi is obsolete is should no longer be used"
  print "        It does the identical thing as pfzx"
  if iyf is not None: iy = iyf
  pfzx(iy=iy,fullplane=fullplane,lbeamframe=lbeamframe,plotsg=plotsg,
       scale=scale,plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

# z-y plane
def pfzyi(ix=None,ixf=None,fullplane=1,lbeamframe=1,plotsg=1,scale=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full Z-Y plane
Same arguments as pfzy
  """
  print "Notice: pfzyi is obsolete is should no longer be used"
  print "        It does the identical thing as pfzy"
  if ixf is not None: ix = ixf
  pfzy(ix=ix,fullplane=fullplane,lbeamframe=lbeamframe,plotsg=plotsg,
       scale=scale,plotphi=plotphi,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)



############################################################################
# These plot a box at each conductor point
############################################################################

# x-y plane
def pfxybox(iz=None,izf=None,contours=8,plotsg=1,scale=1,signx=1,signy=1,
            plotphi=1,filled=0,phicolor=blue,condcolor=cyan,kwdict={},**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in X-Y plane
  - iz=w3d.iz_axis z index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signx=1 sign of x, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - subgridlen=1 maximum length of subgrid line which are plotted
  """
  kw.update(kwdict)
  if izf is not None: iz = izf
  if not iz: iz = w3d.iz_axis
  if iz < 0 or w3d.nzfull < iz: return
  izlocal = iz - top.izslave[me]
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
  if plotphi:
    ppp = getphi(iz=iz)
    if me == 0:
      if kw.has_key('cellarray') and kw['cellarray']: contours=None
      ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=phicolor,
                xmin=xmmin,xmax=xmmin+w3d.nx*dx,
                ymin=ymmin,ymax=ymmin+w3d.ny*dy,kwdict=kw)
  if f3d.ncond > 0:
    ii = compress(equal(f3d.izcond[0:f3d.ncond],izlocal),arange(f3d.ncond))
    x = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    y = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
  else:
    x = array([])
    y = array([])
  if lparallel:
    x = gatherarray(x)
    y = gatherarray(y)
  if len(x) > 0:
    pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
        array([x-dx/2,x+dx/2,x+dx/2,x-dx/2,x-dx/2]),
        color=condcolor)

# z-x plane
def pfzxbox(iy=None,iyf=None,contours=8,plotsg=1,scale=1,signz=1,signx=1,
            plotphi=1,filled=0,phicolor=blue,condcolor=cyan,kwdict={},**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-X plane
  - iy=w3d.iy_axis y index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signx=1 sign of x, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  """
  kw.update(kwdict)
  if iyf is not None: iy = iyf
  if not iy: iy = w3d.iy_axis
  if iy < 0 or w3d.ny < iy: return
  if scale:
    dx = w3d.dx*signx
    dz = w3d.dz*signz
    xmmin = w3d.xmmin
    zmmin = top.zplmin
  else:
    dx = 1.*signx
    dz = 1.*signz
    xmmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
  if plotphi:
    ppp = getphi(iy=iy)
    ppp = transpose(ppp)
    if me == 0:
      if kw.has_key('cellarray') and kw['cellarray']: contours=None
      ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=phicolor,
                xmin=zmmin,xmax=zmmin+w3d.nzfull*dz,
                ymin=xmmin,ymax=xmmin+w3d.nx*dx,kwdict=kw)
  if (f3d.ncond > 0):
    ii = compress(equal(f3d.iycond[0:f3d.ncond],iy),arange(f3d.ncond))
    x = take(f3d.ixcond[0:f3d.ncond],ii)*dx+xmmin
    z = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  else:
    x = array([])
    z = array([])
  if lparallel:
    x = gatherarray(x)
    z = gatherarray(z)
  if len(x) > 0:
    pla(array([x-dx/2,x-dx/2,x+dx/2,x+dx/2,x-dx/2]),
        array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
        color=condcolor)

# z-y plane
def pfzybox(ix=None,ixf=None,contours=8,plotsg=1,scale=1,signz=1,signy=1,
            plotphi=1,filled=0,phicolor=blue,condcolor=cyan,kwdict={},**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-Y plane
  - ix=w3d.ix_axis x index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - signy=1 sign of y, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  """
  kw.update(kwdict)
  if ixf is not None: ix = ixf
  if not ix: ix = w3d.ix_axis
  if ix < 0 or w3d.nx < ix: return
  if scale:
    dy = w3d.dy*signy
    dz = w3d.dz*signz
    ymmin = w3d.ymmin
    zmmin = top.zplmin
  else:
    dy = 1.*signy
    dz = 1.*signz
    ymmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
  if plotphi:
    ppp = getphi(ix=ix)
    ppp = transpose(ppp)
    if me == 0:
      if kw.has_key('cellarray') and kw['cellarray']: contours=None
      ppgeneric(grid=ppp,contours=contours,filled=filled,ccolor=phicolor,
                xmin=zmmin,xmax=zmmin+w3d.nzfull*dz,
                ymin=ymmin,ymax=ymmin+w3d.ny*dy,kwdict=kw)
  if (f3d.ncond > 0):
    ii = compress(equal(f3d.ixcond[0:f3d.ncond],ix),arange(f3d.ncond))
    y = take(f3d.iycond[0:f3d.ncond],ii)*dy+ymmin
    z = take(f3d.izcond[0:f3d.ncond],ii)*dz+zmmin
  else:
    y = array([])
    z = array([])
  if lparallel:
    y = gatherarray(y)
    z = gatherarray(z)
  if len(y) > 0:
    pla(array([y-dy/2,y-dy/2,y+dy/2,y+dy/2,y-dy/2]),
        array([z-dz/2,z+dz/2,z+dz/2,z-dz/2,z-dz/2]),
        color=condcolor)

# z-x plane
def pfzxboxi(iy=None,iyf=None,contours=8,plotsg=1,scale=1,signz=1,
             plotphi=1,filled=0,phicolor=blue,condcolor=cyan,**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-(-X) plane
  - iy=w3d.iy_axis y index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - plotphi=1 when true, plot contours of potential
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  """
  if iyf is not None: iy = iyf
  pfzxbox(iy=iy,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
          signx=-1,plotphi=plotphi,filled=filled,
          phicolor=phicolor,condcolor=condcolor,kwdict=kw)

# z-y plane
def pfzyboxi(ix=None,ixf=None,contours=8,plotsg=1,scale=1,signz=1,signy=-1,
             plotphi=1,filled=0,phicolor=blue,condcolor=cyan,**kw):
  """
Plots square at conductor points and contours of electrostatic potential
in Z-(-Y) plane
  - ix=w3d.ix_axis x index of plane
  - contours=8 optional number of or list of contours
  - plotsg=1 when true, plots subgrid data
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - signz=1 sign of z, used for plotting symmetry planes
  - filled=0 when true, plots filled contours
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  """
  if ixf is not None: ix = ixf
  pfzybox(ix=ix,contours=contours,plotsg=plotsg,scale=scale,signz=signz,
          signy=-1,plotphi=plotphi,filled=filled,
          phicolor=phicolor,condcolor=condcolor,kwdict=kw)




############################################################################
# These plots plot the conductor points colored based on the conductor
# number. The list of colors is input by the user.

# --- convenience function
def findunique(i):
  ii = sort(i)
  result = list(compress(ii[:-1]!=ii[1:],ii[:-1])) + [ii[-1]]
  return result

def plotcondn(iz,nc,cx,cy,cz,cn,dx,dy,xmmin,ymmin,marker,color):
  ncolor = len(color)
  if f3d.ncond > 0:
    ii = compress(equal(cz[0:nc],iz),arange(nc))
    xx = take(cx[0:nc],ii)*dx+xmmin
    yy = take(cy[0:nc],ii)*dy+ymmin
    nn = take(cn[0:nc],ii)
  else:
    xx = array([])
    yy = array([])
    nn = array([])
  nlist = gatherarray(nn)
  nlist = findunique(nlist)
  nlist = broadcast(nlist)
  for i in range(len(nlist)):
    x = compress(equal(nn,nlist[i]),xx)
    y = compress(equal(nn,nlist[i]),yy)
    warpplp(y,x,color=color[i%ncolor],marker=marker)

def pfzxn(iy=None,numbs=None,colors=None,cmarker=point,smarker=circle,
          scale=1,signz=1,signx=1,subgridlen=1.,fullplane=1):
  if iy is None: iy = w3d.iy_axis
  if iy < 0 or w3d.ny < iy: return
  if colors is None: colors = color
  if scale:
    dx = w3d.dx*signx
    dz = w3d.dz*signz
    xmmin = w3d.xmmin
    zmmin = top.zplmin
  else:
    dx = 1.*signx
    dz = 1.*signz
    xmmin = 0.
    zmmin = 0.
    if lparallel: zmmin = top.izslave[me]
  plotcondn(iy,f3d.ncond,f3d.izcond,f3d.ixcond,f3d.iycond,f3d.condnumb,
            dx,dz,xmmin,zmmin,cmarker,color)
  ncolor = len(colors)
  nlist = gatherarray(f3d.ecnumb[:f3d.necndbdy])
  nlist = findunique(nlist)
  nlist.remove(0)
  nlist = broadcast(nlist)
  for i in range(len(nlist)):
    plotsubgrid('x','z','y','e',iy,nlist[i],xmmin,zmmin,dx,dz,
                colors[i%ncolor],subgridlen,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('x','z','y','e',iy,nlist[i],xmmin,zmmin,dx,dz,
                  colors[i%ncolor],subgridlen,-1,1)
  nlist = gatherarray(f3d.ocnumb[:f3d.nocndbdy])
  nlist = findunique(nlist)
  nlist = broadcast(nlist)
  for i in range(len(nlist)):
    plotsubgrid('x','z','y','o',iy,nlist[i],xmmin,zmmin,dx,dz,
                colors[i%ncolor],subgridlen,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('x','z','y','o',iy,nlist[i],xmmin,zmmin,dx,dz,
                  colors[i%ncolor],subgridlen,-1,1)


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
    zzz[ix,:] = top.zplmin + iota(0,w3d.nz,ii)*w3d.dz + w3d.zz

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
def pfzxlab(zz=None,iy=None,condcolor=cyan):
  """Plots conductors in Z-X lab frame (including bends)
  - zz=top.zbeam is the center position
  - condcolor=cyan color of conductor points inside conductors
  """
  if not zz: zz=top.zbeam
  if iy is None: iy = w3d.iy_axis
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
                        f3d.izcond[0:f3d.ncond])*w3d.dz+top.zplmin+zz
    # --- convert to lab frame
    tolabfrm(zz,len(xxxx),xxxx,zzzz)   
    # --- make plot
    plg(xxxx,zzzz,marker='\2',color=condcolor)
  # --- restore original conductor data at zbeam
  if (zz != top.zbeam):
    top.zbeam = z
    top.zgrid = g
    setlatt()
    fieldsol(1)


#####################################################################
def plotsrfrv(srfrv,zmin,zmax,n=1000,color='fg',gridframe=0,rscale=1,zscale=1,
              roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0):
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
    zz = (zz - top.zplmin)/w3d.dz
    rr = (rr)/w3d.dx + ir_axis
  rr = where(less(rr,rmin),rmin,rr)
  rr = where(greater(rr,rmax),rmax,rr)
  plg(rscale*rr+roff,zscale*zz+zoff,color=color)


#####################################################################
#####################################################################
def plotelementoutline(color,gridframe,axis,ie,ne,
                       ezs,eze,eap,err,erl,egl,egp,eox,eoy,epa,epr,epw,
                       dpal,dpar):
  """Plots the outline of electrostatic elements
  - color: line color
  - gridframe: when true, make plot in grid coordinates
  - axis: selects axis to plot, either 'x' or 'y'
  """
  if axis == 'x': gpsign = 1
  else:           gpsign = -1
  for i in range(ie,ie+ne):
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
        zz = (zz - top.zplmin)/w3d.dz
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
        zzl = (zzl - top.zplmin)/w3d.dz
        zzr = (zzr - top.zplmin)/w3d.dz
      plg(rrl1,zzl,color=color)
      plg(rrl2,zzl,color=color)
      plg(rrr1,zzr,color=color)
      plg(rrr2,zzr,color=color)


#---------------------------------------------------------------------------
def plotquadoutline(iq=0,nq=None,color='fg',gridframe=0,axis='x'):
  """Plots the outline of quadrupole elements
  - iq=0: starting quad to plot
  - nq=top.nquad+1: number of quads to plot
  - color='fg': line color
  - gridframe=0: when true, make plot in grid coordinates
  - axis='x': selects axis to plot, either 'x' or 'y'
  """
  if nq is None: nq = top.nquad + 1
  plotelementoutline(color,gridframe,axis,iq,nq,
                     top.quadzs,top.quadze,top.quadap,top.quadrr,top.quadrl,
                     top.quadgl,top.quadgp,top.qoffx,top.qoffy,
                     top.quadpa,top.quadpr,top.quadpw,
                     top.qdelpal,top.qdelpar)

#---------------------------------------------------------------------------
def plotemltoutline(ie=0,ne=None,color='fg',gridframe=0,axis='x'):
  """Plots the outline of emlt elements
  - ie=0: starting emlt to plot
  - ne=top.nemlt+1: number of emlts to plot
  - color='fg': line color
  - gridframe=0: when true, make plot in grid coordinates
  - axis='x': selects axis to plot, either 'x' or 'y'
  """
  if ne is None: ne = top.nemlt + 1
  plotelementoutline(color,gridframe,axis,ie,ne,
                     top.emltzs,top.emltze,top.emltap,top.emltrr,top.emltrl,
                     top.emltgl,top.emltgp,top.emltox,top.emltoy,
                     top.emltpa,zeros(top.nemlt+1,'d'),top.emltpw,
                     zeros(top.nemlt+1,'d'),zeros(top.nemlt+1,'d'))



#########################################################################
#########################################################################
def cleanconductors():
  """This routine clears out conductor points which are not within the
range of the field solution, w3d.izfsmin and w3d.izfsmax. This is done
for optimization so that time is not wasted on those points.
  """
  if f3d.ncond > 0:
    ii = compress(logical_and(less(w3d.izfsmin-1,f3d.izcond[:f3d.ncond]), \
                              less(f3d.izcond[:f3d.ncond],w3d.izfsmax+1)), \
                  arange(f3d.ncond))
    xx = take(f3d.ixcond,ii)
    yy = take(f3d.iycond,ii)
    zz = take(f3d.izcond,ii)
    vv = take(f3d.condvolt,ii)
    f3d.ncond = len(ii)
    f3d.ncondmax = f3d.ncond
    gchange("PSOR3d")
    if f3d.ncond > 0:
      f3d.ixcond[:] = xx
      f3d.iycond[:] = yy
      f3d.izcond[:] = zz
      f3d.condvolt[:] = vv
  if f3d.necndbdy > 0:
    ii = compress(logical_and(less(w3d.izfsmin-1,f3d.iecndz[:f3d.necndbdy]), \
                              less(f3d.iecndz[:f3d.necndbdy],w3d.izfsmax+1)), \
                  arange(f3d.necndbdy))
    xx = take(f3d.iecndx,ii)
    yy = take(f3d.iecndy,ii)
    zz = take(f3d.iecndz,ii)
    mx = take(f3d.ecdelmx,ii)
    my = take(f3d.ecdelmy,ii)
    mz = take(f3d.ecdelmz,ii)
    px = take(f3d.ecdelpx,ii)
    py = take(f3d.ecdelpy,ii)
    pz = take(f3d.ecdelpz,ii)
    vv = take(f3d.ecvolt,ii)
    f3d.necndbdy = len(ii)
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
    gchange("PSOR3d")
    if f3d.necndbdy > 0:
      f3d.iecndx[:f3d.necndbdy] = xx
      f3d.iecndy[:f3d.necndbdy] = yy
      f3d.iecndz[:f3d.necndbdy] = zz
      f3d.ecdelmx[:f3d.necndbdy] = mx
      f3d.ecdelmy[:f3d.necndbdy] = my
      f3d.ecdelmz[:f3d.necndbdy] = mz
      f3d.ecdelpx[:f3d.necndbdy] = px
      f3d.ecdelpy[:f3d.necndbdy] = py
      f3d.ecdelpz[:f3d.necndbdy] = pz
      f3d.ecvolt[:f3d.necndbdy] = vv
  if f3d.nocndbdy > 0:
    ii = compress(logical_and(less(w3d.izfsmin-1,f3d.iocndz[:f3d.nocndbdy]), \
                              less(f3d.iocndz[:f3d.nocndbdy],w3d.izfsmax+1)), \
                  arange(f3d.nocndbdy))
    xx = take(f3d.iocndx,ii)
    yy = take(f3d.iocndy,ii)
    zz = take(f3d.iocndz,ii)
    mx = take(f3d.ocdelmx,ii)
    my = take(f3d.ocdelmy,ii)
    mz = take(f3d.ocdelmz,ii)
    px = take(f3d.ocdelpx,ii)
    py = take(f3d.ocdelpy,ii)
    pz = take(f3d.ocdelpz,ii)
    vv = take(f3d.ocvolt,ii)
    f3d.nocndbdy = len(ii)
    f3d.ncndmax = max(f3d.necndbdy,f3d.nocndbdy)
    gchange("PSOR3d")
    if f3d.nocndbdy > 0:
      f3d.iocndx[:f3d.nocndbdy] = xx
      f3d.iocndy[:f3d.nocndbdy] = yy
      f3d.iocndz[:f3d.nocndbdy] = zz
      f3d.ocdelmx[:f3d.nocndbdy] = mx
      f3d.ocdelmy[:f3d.nocndbdy] = my
      f3d.ocdelmz[:f3d.nocndbdy] = mz
      f3d.ocdelpx[:f3d.nocndbdy] = px
      f3d.ocdelpy[:f3d.nocndbdy] = py
      f3d.ocdelpz[:f3d.nocndbdy] = pz
      f3d.ocvolt[:f3d.nocndbdy] = vv

#########################################################################
def updatemgconductors():
  """
This routine updates conductors which from old dump files. In older versions
of the code, the conductor coordinates and deltas were stored relative to
the finest grid. In the current version, the data is stored relative the to
grid that data is to be used for.
  """
  # --- Check is one of the MG conductor arrays is allocated. If not, then
  # --- return since the conductors have not been generated for the MG solver.
  # --- This check is used instead of fstype==7 since when the field solve
  # --- is done periodically, fstype may be -1 in the restart dump.
  test = f3d.getpyobject("ecdelmx")
  if test is None: return

  # --- First, copy all of the data
  ixcond = f3d.ixcond[0:f3d.ncond] + 0.
  iycond = f3d.iycond[0:f3d.ncond] + 0.
  izcond = f3d.izcond[0:f3d.ncond] + 0.
  condvolt = f3d.condvolt[0:f3d.ncond] + 0.
  condnumb = f3d.condnumb[0:f3d.ncond] + 0.
  icondlxy = f3d.icondlxy[0:f3d.ncond] + 0
  icondlz  = f3d.icondlz[0:f3d.ncond] + 0
  iecndx = f3d.iecndx[0:f3d.necndbdy] + 0.
  iecndy = f3d.iecndy[0:f3d.necndbdy] + 0.
  iecndz = f3d.iecndz[0:f3d.necndbdy] + 0.
  ecdelmx = f3d.ecdelmx[0:f3d.necndbdy] + 0.
  ecdelpx = f3d.ecdelpx[0:f3d.necndbdy] + 0.
  ecdelmy = f3d.ecdelmy[0:f3d.necndbdy] + 0.
  ecdelpy = f3d.ecdelpy[0:f3d.necndbdy] + 0.
  ecdelmz = f3d.ecdelmz[0:f3d.necndbdy] + 0.
  ecdelpz = f3d.ecdelpz[0:f3d.necndbdy] + 0.
  ecvolt = f3d.ecvolt[0:f3d.necndbdy] + 0.
  ecnumb = f3d.ecnumb[0:f3d.necndbdy] + 0.
  ecvoltmx = f3d.ecvoltmx[0:f3d.necndbdy] + 0.
  ecvoltpx = f3d.ecvoltpx[0:f3d.necndbdy] + 0.
  ecvoltmy = f3d.ecvoltmy[0:f3d.necndbdy] + 0.
  ecvoltpy = f3d.ecvoltpy[0:f3d.necndbdy] + 0.
  ecvoltmz = f3d.ecvoltmz[0:f3d.necndbdy] + 0.
  ecvoltpz = f3d.ecvoltpz[0:f3d.necndbdy] + 0.
  ecnumbmx = f3d.ecnumbmx[0:f3d.necndbdy] + 0.
  ecnumbpx = f3d.ecnumbpx[0:f3d.necndbdy] + 0.
  ecnumbmy = f3d.ecnumbmy[0:f3d.necndbdy] + 0.
  ecnumbpy = f3d.ecnumbpy[0:f3d.necndbdy] + 0.
  ecnumbmz = f3d.ecnumbmz[0:f3d.necndbdy] + 0.
  ecnumbpz = f3d.ecnumbpz[0:f3d.necndbdy] + 0.
  iecndlxy = f3d.iecndlxy[0:f3d.necndbdy] + 0.
  iecndlz = f3d.iecndlz[0:f3d.necndbdy] + 0.
  iocndx = f3d.iocndx[0:f3d.nocndbdy] + 0.
  iocndy = f3d.iocndy[0:f3d.nocndbdy] + 0.
  iocndz = f3d.iocndz[0:f3d.nocndbdy] + 0.
  ocdelmx = f3d.ocdelmx[0:f3d.nocndbdy] + 0.
  ocdelpx = f3d.ocdelpx[0:f3d.nocndbdy] + 0.
  ocdelmy = f3d.ocdelmy[0:f3d.nocndbdy] + 0.
  ocdelpy = f3d.ocdelpy[0:f3d.nocndbdy] + 0.
  ocdelmz = f3d.ocdelmz[0:f3d.nocndbdy] + 0.
  ocdelpz = f3d.ocdelpz[0:f3d.nocndbdy] + 0.
  ocvolt = f3d.ocvolt[0:f3d.nocndbdy] + 0.
  ocnumb = f3d.ocnumb[0:f3d.nocndbdy] + 0.
  ocvoltmx = f3d.ocvoltmx[0:f3d.nocndbdy] + 0.
  ocvoltpx = f3d.ocvoltpx[0:f3d.nocndbdy] + 0.
  ocvoltmy = f3d.ocvoltmy[0:f3d.nocndbdy] + 0.
  ocvoltpy = f3d.ocvoltpy[0:f3d.nocndbdy] + 0.
  ocvoltmz = f3d.ocvoltmz[0:f3d.nocndbdy] + 0.
  ocvoltpz = f3d.ocvoltpz[0:f3d.nocndbdy] + 0.
  ocnumbmx = f3d.ocnumbmx[0:f3d.nocndbdy] + 0.
  ocnumbpx = f3d.ocnumbpx[0:f3d.nocndbdy] + 0.
  ocnumbmy = f3d.ocnumbmy[0:f3d.nocndbdy] + 0.
  ocnumbpy = f3d.ocnumbpy[0:f3d.nocndbdy] + 0.
  ocnumbmz = f3d.ocnumbmz[0:f3d.nocndbdy] + 0.
  ocnumbpz = f3d.ocnumbpz[0:f3d.nocndbdy] + 0.
  iocndlxy = f3d.iocndlxy[0:f3d.nocndbdy] + 0.
  iocndlz = f3d.iocndlz[0:f3d.nocndbdy] + 0.

  # --- Get the list of coarsening levels
  llxy = findunique(icondlxy)
  llz  = findunique(icondlz)

  # --- Create lists to hold the converted data.
  ixcondnew = []
  iycondnew = []
  izcondnew = []
  condvoltnew = []
  iecndxnew = []
  iecndynew = []
  iecndznew = []
  ecdelmxnew = []
  ecdelmynew = []
  ecdelmznew = []
  ecdelpxnew = []
  ecdelpynew = []
  ecdelpznew = []
  ecvoltnew = []
  iocndxnew = []
  iocndynew = []
  iocndznew = []
  ocdelmxnew = []
  ocdelmynew = []
  ocdelmznew = []
  ocdelpxnew = []
  ocdelpynew = []
  ocdelpznew = []
  ocvoltnew = []
  ecvoltmxnew = []
  ecvoltpxnew = []
  ecvoltmynew = []
  ecvoltpynew = []
  ecvoltmznew = []
  ecvoltpznew = []
  ocvoltmxnew = []
  ocvoltpxnew = []
  ocvoltmynew = []
  ocvoltpynew = []
  ocvoltmznew = []
  ocvoltpznew = []
  icondlxynew = []
  icondlznew = []
  iecndlxynew = []
  iecndlznew = []
  iocndlxynew = []
  iocndlznew = []

  # --- Get number of old data points
  ncond = len(ixcond)
  necndbdy = len(iecndx)
  nocndbdy = len(iocndx)

  # --- Loop over coarsening levels, collecting and converting the data
  # --- for each level.
  for j in xrange(len(llxy)):
    ii = compress(logical_and(icondlxy >= llxy[j],icondlz >= llz[j]), \
                  arange(ncond))
    ixcondnew = ixcondnew + list(take(ixcond/llxy[j],ii))
    iycondnew = iycondnew + list(take(iycond/llxy[j],ii))
    izcondnew = izcondnew + list(take(izcond/llz[j],ii))
    condvoltnew = condvoltnew + list(take(condvolt,ii))
    icondlxynew = icondlxynew + list(llxy[j]*ones(len(ii)))
    icondlznew = icondlznew + list(llz[j]*ones(len(ii)))
  
    ii = compress(logical_and(iecndlxy >= llxy[j],iecndlz >= llz[j]), \
                  arange(necndbdy))
    iecndxnew = iecndxnew + list(take(iecndx/llxy[j],ii))
    iecndynew = iecndynew + list(take(iecndy/llxy[j],ii))
    iecndznew = iecndznew + list(take(iecndz/llz[j],ii))
    ecdelmxnew = ecdelmxnew + list(take(ecdelmx/llxy[j],ii))
    ecdelmynew = ecdelmynew + list(take(ecdelmy/llxy[j],ii))
    ecdelmznew = ecdelmznew + list(take(ecdelmz/llz[j],ii))
    ecdelpxnew = ecdelpxnew + list(take(ecdelpx/llxy[j],ii))
    ecdelpynew = ecdelpynew + list(take(ecdelpy/llxy[j],ii))
    ecdelpznew = ecdelpznew + list(take(ecdelpz/llz[j],ii))
    ecvoltnew = ecvoltnew + list(take(ecvolt,ii))
    ecvoltmxnew = ecvoltmxnew + list(take(ecvoltmx,ii))
    ecvoltpxnew = ecvoltpxnew + list(take(ecvoltpx,ii))
    ecvoltmynew = ecvoltmynew + list(take(ecvoltmy,ii))
    ecvoltpynew = ecvoltpynew + list(take(ecvoltpy,ii))
    ecvoltmznew = ecvoltmznew + list(take(ecvoltmz,ii))
    ecvoltpznew = ecvoltpznew + list(take(ecvoltpz,ii))
    iecndlxynew = iecndlxynew + list(llxy[j]*ones(len(ii)))
    iecndlznew = iecndlznew + list(llz[j]*ones(len(ii)))
  
    ii = compress(logical_and(iocndlxy >= llxy[j],iocndlz >= llz[j]), \
                  arange(nocndbdy))
    iocndxnew = iocndxnew + list(take(iocndx/llxy[j],ii))
    iocndynew = iocndynew + list(take(iocndy/llxy[j],ii))
    iocndznew = iocndznew + list(take(iocndz/llz[j],ii))
    ocdelmxnew = ocdelmxnew + list(take(ocdelmx/llxy[j],ii))
    ocdelmynew = ocdelmynew + list(take(ocdelmy/llxy[j],ii))
    ocdelmznew = ocdelmznew + list(take(ocdelmz/llz[j],ii))
    ocdelpxnew = ocdelpxnew + list(take(ocdelpx/llxy[j],ii))
    ocdelpynew = ocdelpynew + list(take(ocdelpy/llxy[j],ii))
    ocdelpznew = ocdelpznew + list(take(ocdelpz/llz[j],ii))
    ocvoltnew = ocvoltnew + list(take(ocvolt,ii))
    ocvoltmxnew = ocvoltmxnew + list(take(ocvoltmx,ii))
    ocvoltpxnew = ocvoltpxnew + list(take(ocvoltpx,ii))
    ocvoltmynew = ocvoltmynew + list(take(ocvoltmy,ii))
    ocvoltpynew = ocvoltpynew + list(take(ocvoltpy,ii))
    ocvoltmznew = ocvoltmznew + list(take(ocvoltmz,ii))
    ocvoltpznew = ocvoltpznew + list(take(ocvoltpz,ii))
    iocndlxynew = iocndlxynew + list(llxy[j]*ones(len(ii)))
    iocndlznew = iocndlznew + list(llz[j]*ones(len(ii)))

  # --- Reset counters and reallocate arrays for the converted data.
  f3d.ncondmax = len(ixcondnew)
  f3d.ncndmax = max(len(iecndxnew),len(iocndxnew))
  gchange("PSOR3d")
  gchange("MultigridConductor3d")
  f3d.ncond = len(ixcondnew)
  f3d.necndbdy = len(iecndxnew)
  f3d.nocndbdy = len(iocndxnew)
  
  # --- Copy converted data into WARP arrays.
  f3d.ixcond[:f3d.ncond] = ixcondnew
  f3d.iycond[:f3d.ncond] = iycondnew
  f3d.izcond[:f3d.ncond] = izcondnew
  f3d.condvolt[:f3d.ncond] = condvoltnew
  f3d.icondlxy[:f3d.ncond]= icondlxynew
  f3d.icondlz[:f3d.ncond]= icondlznew
  f3d.iecndx[:f3d.necndbdy]= iecndxnew
  f3d.iecndy[:f3d.necndbdy]= iecndynew
  f3d.iecndz[:f3d.necndbdy]= iecndznew
  f3d.ecdelmx[:f3d.necndbdy]= ecdelmxnew
  f3d.ecdelmy[:f3d.necndbdy]= ecdelmynew
  f3d.ecdelmz[:f3d.necndbdy]= ecdelmznew
  f3d.ecdelpx[:f3d.necndbdy]= ecdelpxnew
  f3d.ecdelpy[:f3d.necndbdy]= ecdelpynew
  f3d.ecdelpz[:f3d.necndbdy]= ecdelpznew
  f3d.ecvolt[:f3d.necndbdy]= ecvoltnew
  f3d.iocndx[:f3d.nocndbdy]= iocndxnew
  f3d.iocndy[:f3d.nocndbdy]= iocndynew
  f3d.iocndz[:f3d.nocndbdy]= iocndznew
  f3d.ocdelmx[:f3d.nocndbdy]= ocdelmxnew
  f3d.ocdelmy[:f3d.nocndbdy]= ocdelmynew
  f3d.ocdelmz[:f3d.nocndbdy]= ocdelmznew
  f3d.ocdelpx[:f3d.nocndbdy]= ocdelpxnew
  f3d.ocdelpy[:f3d.nocndbdy]= ocdelpynew
  f3d.ocdelpz[:f3d.nocndbdy]= ocdelpznew
  f3d.ocvolt[:f3d.nocndbdy]= ocvoltnew
  f3d.ecvoltmx[:f3d.necndbdy]= ecvoltmxnew
  f3d.ecvoltpx[:f3d.necndbdy]= ecvoltpxnew
  f3d.ecvoltmy[:f3d.necndbdy]= ecvoltmynew
  f3d.ecvoltpy[:f3d.necndbdy]= ecvoltpynew
  f3d.ecvoltmz[:f3d.necndbdy]= ecvoltmznew
  f3d.ecvoltpz[:f3d.necndbdy]= ecvoltpznew
  f3d.ocvoltmx[:f3d.nocndbdy]= ocvoltmxnew
  f3d.ocvoltpx[:f3d.nocndbdy]= ocvoltpxnew
  f3d.ocvoltmy[:f3d.nocndbdy]= ocvoltmynew
  f3d.ocvoltpy[:f3d.nocndbdy]= ocvoltpynew
  f3d.ocvoltmz[:f3d.nocndbdy]= ocvoltmznew
  f3d.ocvoltpz[:f3d.nocndbdy]= ocvoltpznew
  f3d.iecndlxy[:f3d.necndbdy]= iecndlxynew
  f3d.iecndlz[:f3d.necndbdy]= iecndlznew
  f3d.iocndlxy[:f3d.nocndbdy]= iocndlxynew
  f3d.iocndlz[:f3d.nocndbdy]= iocndlznew



#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
# Convenience routines for calling the surface of revolution routines.
# The presence of these allows changing the argument lists of the fortran
# versions without breaking code.
#########################################################################
def pysrfrvout(rofzfunc=" ",volt=0.,zmin=None,zmax=None,xcent=0.,ycent=0.,
             rmax=top.largepos,lfill=false,
             xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
             zmmin=None,zmmax=None,zbeam=None,dx=None,dy=None,dz=None,
             nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
             xmesh=None,ymesh=None,l2symtry=None,l4symtry=None):
  if zmin is None: zmin = w3d.zmmin
  if zmax is None: zmax = w3d.zmmax
  if xmin is None: xmin = w3d.xmmin
  if xmax is None: xmax = w3d.xmmax
  if ymin is None: ymin = w3d.ymmin
  if ymax is None: ymax = w3d.ymmax
  if zmmin is None: zmmin = w3d.zmmin
  if zmmax is None: zmmax = w3d.zmmax
  if zbeam is None: zbeam = top.zbeam
  if dx is None: dx = w3d.dx
  if dy is None: dy = w3d.dy
  if dz is None: dz = w3d.dz
  if nx is None: nx = w3d.nx
  if ny is None: ny = w3d.ny
  if nz is None: nz = w3d.nz
  if ix_axis is None: ix_axis = w3d.ix_axis
  if iy_axis is None: iy_axis = w3d.iy_axis
  if xmesh is None: xmesh = w3d.xmesh
  if ymesh is None: ymesh = w3d.ymesh
  if l2symtry is None: l2symtry = w3d.l2symtry
  if l4symtry is None: l4symtry = w3d.l4symtry

  # --- Now call the fortran version
  f3d.srfrvout(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
               xmin,xmax,ymin,ymax,lshell,zmmin,zmmax,zbeam,dx,dy,dz,
               nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

#---------------------------------------------------------------------------
def pysrfrvin(rofzfunc=" ",volt=0.,zmin=None,zmax=None,xcent=0.,ycent=0.,
            rmin=0.,lfill=false,
            xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
            zmmin=None,zmmax=None,zbeam=None,dx=None,dy=None,dz=None,
            nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
            xmesh=None,ymesh=None,l2symtry=None,l4symtry=None):
  if zmin is None: zmin = w3d.zmmin
  if zmax is None: zmax = w3d.zmmax
  if xmin is None: xmin = w3d.xmmin
  if xmax is None: xmax = w3d.xmmax
  if ymin is None: ymin = w3d.ymmin
  if ymax is None: ymax = w3d.ymmax
  if zmmin is None: zmmin = w3d.zmmin
  if zmmax is None: zmmax = w3d.zmmax
  if zbeam is None: zbeam = top.zbeam
  if dx is None: dx = w3d.dx
  if dy is None: dy = w3d.dy
  if dz is None: dz = w3d.dz
  if nx is None: nx = w3d.nx
  if ny is None: ny = w3d.ny
  if nz is None: nz = w3d.nz
  if ix_axis is None: ix_axis = w3d.ix_axis
  if iy_axis is None: iy_axis = w3d.iy_axis
  if xmesh is None: xmesh = w3d.xmesh
  if ymesh is None: ymesh = w3d.ymesh
  if l2symtry is None: l2symtry = w3d.l2symtry
  if l4symtry is None: l4symtry = w3d.l4symtry

  # --- Now call the fortran version
  f3d.srfrvin(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
              xmin,xmax,ymin,ymax,lshell,zmmin,zmmax,zbeam,dx,dy,dz,
              nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

#---------------------------------------------------------------------------
def pysrfrvinout(rminofz=" ",rmaxofz=" ",volt=0.,zmin=None,zmax=None,
               xcent=0.,ycent=0.,lzend=true,
               xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
               zmmin=None,zmmax=None,zbeam=None,dx=None,dy=None,dz=None,
               nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
               xmesh=None,ymesh=None,l2symtry=None,l4symtry=None):
  if zmin is None: zmin = w3d.zmmin
  if zmax is None: zmax = w3d.zmmax
  if xmin is None: xmin = w3d.xmmin
  if xmax is None: xmax = w3d.xmmax
  if ymin is None: ymin = w3d.ymmin
  if ymax is None: ymax = w3d.ymmax
  if zmmin is None: zmmin = w3d.zmmin
  if zmmax is None: zmmax = w3d.zmmax
  if zbeam is None: zbeam = top.zbeam
  if dx is None: dx = w3d.dx
  if dy is None: dy = w3d.dy
  if dz is None: dz = w3d.dz
  if nx is None: nx = w3d.nx
  if ny is None: ny = w3d.ny
  if nz is None: nz = w3d.nz
  if ix_axis is None: ix_axis = w3d.ix_axis
  if iy_axis is None: iy_axis = w3d.iy_axis
  if xmesh is None: xmesh = w3d.xmesh
  if ymesh is None: ymesh = w3d.ymesh
  if l2symtry is None: l2symtry = w3d.l2symtry
  if l4symtry is None: l4symtry = w3d.l4symtry

  # --- Now call the fortran version
  f3d.srfrvinout(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,lzend,
                 xmin,xmax,ymin,ymax,lshell,zmmin,zmmax,zbeam,dx,dy,dz,
                 nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

#---------------------------------------------------------------------------
