from warp import *
import __main__
plot_conductor_version = "$Id: plot_conductor.py,v 1.56 2003/04/05 01:07:55 dave Exp $"

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

visualizeconductors: create 3-D visualization of conductors based on subgrid
                     data

cleanconductors: not a plot routine, buts removes conductor points not
		 within the the range of the field solve
  """

######################################################################
# functions to plot the conductor points and subgrid data            #
# in MKS units                                                       #
######################################################################

# --- Convenience function to plot the sub-grid data
def plotcond(yy,xx,zz,iz,numb,ymin,xmin,dy,dx,color,mglevel,signy,signx):
  nn = f3d.ncond
  if nn > 0:
    lx = eval('f3d.mglevelsl'+xx)[mglevel]
    ly = eval('f3d.mglevelsl'+yy)[mglevel]
    lz = eval('f3d.mglevelsl'+zz)[mglevel]
    ixc = eval('f3d.i'+xx+'cond')*signx*lx
    iyc = eval('f3d.i'+yy+'cond')*signy*ly
    izc = eval('f3d.i'+zz+'cond')      *lz
    if xx == 'z': ixc[:] = ixc + f3d.mglevelsiz[mglevel]*lx
    if yy == 'z': iyc[:] = iyc + f3d.mglevelsiz[mglevel]*ly
    if zz == 'z': izc[:] = izc + f3d.mglevelsiz[mglevel]*lz
    cnumb = f3d.condnumb
    try:
      level = equal(mglevel,f3d.icondlevel)
    except:
      level = ones(f3d.ncondmax)
  else:
    ixc = array([])
    iyc = array([])
    izc = array([])
    cnumb = array([])
    level = array([])
  ii = compress(logical_and(equal(izc[:nn],iz),level[:nn]),arange(nn))
  if numb is not None:
    cnumb = take(cnumb,ii)
    ii = compress(equal(cnumb,numb),ii)
  xx = take(ixc,ii)*dx+xmin
  yy = take(iyc,ii)*dy+ymin
  warpplp(yy,xx,color=color)

def plotsubgrid(yy,xx,zz,pp,iz,numb,ymin,xmin,dy,dx,color,subgridlen,mglevel,
                signy,signx):
  assert (pp in [0,1]),"pp has invalid data"
  pp = ((pp + f3d.mglevelsiz[mglevel]) % 2)
  pp = ['e','o'][pp]
  nn = eval('f3d.n'+pp+'cndbdy')
  if nn > 0:
    lx = eval('f3d.mglevelsl'+xx)[mglevel]
    ly = eval('f3d.mglevelsl'+yy)[mglevel]
    lz = eval('f3d.mglevelsl'+zz)[mglevel]
    ixc = eval('f3d.i'+pp+'cnd'+xx)*signx
    iyc = eval('f3d.i'+pp+'cnd'+yy)*signy
    izc = eval('f3d.i'+pp+'cnd'+zz)*lz
    if xx == 'z': ixc[:] = ixc + f3d.mglevelsiz[mglevel]
    if yy == 'z': iyc[:] = iyc + f3d.mglevelsiz[mglevel]
    if zz == 'z': izc[:] = izc + f3d.mglevelsiz[mglevel]
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
      lev = eval('f3d.i'+pp+'cndlevel')
      level = equal(mglevel,lev)
    except:
      level = ones(f3d.ncndmax)
  else:
    lx = array([1])
    ly = array([1])
    lz = array([1])
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

def plotcondfill(yy,xx,zz,iz,ymin,xmin,dy,dx,mglevel,signy,signx):
  """
Plots conductors, filling them in with a solid color. The color is given
by the conductor number.
  """
  # --- Get the numbers of conductor points
  nc = f3d.ncond
  ne = f3d.necndbdy
  no = f3d.nocndbdy
  ns = ne + no
  if ns+nc > 0:
    lx = eval('f3d.mglevelsl'+xx)[mglevel]
    ly = eval('f3d.mglevelsl'+yy)[mglevel]
    lz = eval('f3d.mglevelsl'+zz)[mglevel]
    ixc = eval('f3d.i'+xx+'cond')[:nc]
    iyc = eval('f3d.i'+yy+'cond')[:nc]
    izc = eval('f3d.i'+zz+'cond')[:nc]
    numb = f3d.condnumb[:nc]
    # --- The even and odd data are merged into the same list.
    ixs = array(list(eval('f3d.iecnd'+xx)[:ne]) +
                list(eval('f3d.iocnd'+xx)[:no]))
    iys = array(list(eval('f3d.iecnd'+yy)[:ne]) +
                list(eval('f3d.iocnd'+yy)[:no]))
    izs = array(list(eval('f3d.iecnd'+zz)[:ne]) +
                list(eval('f3d.iocnd'+zz)[:no]))*lz
    # --- Add z offset of data. This only applies for the parallel version
    if xx == 'z':
      ixc[:] = ixc + f3d.mglevelsiz[mglevel]
      ixs[:] = ixs + f3d.mglevelsiz[mglevel]
    if yy == 'z':
      iyc[:] = iyc + f3d.mglevelsiz[mglevel]
      iys[:] = iys + f3d.mglevelsiz[mglevel]
    if zz == 'z':
      izc[:] = izc + f3d.mglevelsiz[mglevel]
      izs[:] = izs + f3d.mglevelsiz[mglevel]
    delmx = array(list(eval('f3d.ecdelm'+xx)[:ne]) +
                  list(eval('f3d.ocdelm'+xx)[:no]))
    delpx = array(list(eval('f3d.ecdelp'+xx)[:ne]) +
                  list(eval('f3d.ocdelp'+xx)[:no]))
    delmy = array(list(eval('f3d.ecdelm'+yy)[:ne]) +
                  list(eval('f3d.ocdelm'+yy)[:no]))
    delpy = array(list(eval('f3d.ecdelp'+yy)[:ne]) +
                  list(eval('f3d.ocdelp'+yy)[:no]))
    numbmx = array(list(eval('f3d.ecnumbm'+xx)[:ne]) +
                   list(eval('f3d.ocnumbm'+xx)[:no]))
    numbpx = array(list(eval('f3d.ecnumbp'+xx)[:ne]) +
                   list(eval('f3d.ocnumbp'+xx)[:no]))
    numbmy = array(list(eval('f3d.ecnumbm'+yy)[:ne]) +
                   list(eval('f3d.ocnumbm'+yy)[:no]))
    numbpy = array(list(eval('f3d.ecnumbp'+yy)[:ne]) +
                   list(eval('f3d.ocnumbp'+yy)[:no]))
    try:
      levc = f3d.icondlevel[:nc]
      levelc = equal(mglevel,levc)
      levs = array(list(f3d.iecndlevel[:ne]) + list(f3d.iocndlevel[:no]))
      levels = equal(mglevel,levs)
    except:
      levelc = ones(nc)
      levels = ones(ns)
  else:
    # --- This code block is needed for the parallel version since not all
    # --- processors may have conductors and global communication is needed.
    lx = array([1])
    ly = array([1])
    lz = array([1])
    ixc = array([])
    iyc = array([])
    izc = array([])
    ixs = array([])
    iys = array([])
    izs = array([])
    delmx = array([])
    delpx = array([])
    delmy = array([])
    delpy = array([])
    levelc = array([])
    levels = array([])
  # --- Select out the conductor points in the appropriate slice and in
  # --- the appropriate refinement level.
  iic = compress(logical_and(equal(izc,iz),equal(levelc,1)),arange(nc))
  iis = compress(logical_and(equal(izs,iz),equal(levels,1)),arange(ns))
  dx = dx*lx*signx
  dy = dy*ly*signy
  ixc = take(ixc,iic)
  iyc = take(iyc,iic)
  xxc = ixc*dx+xmin
  yyc = iyc*dy+ymin
  numb = take(numb,iic)
  ixs = take(ixs,iis)
  iys = take(iys,iis)
  xxs = ixs*dx+xmin
  yys = iys*dy+ymin
  delmx = take(delmx,iis)*dx
  delpx = take(delpx,iis)*dx
  delmy = take(delmy,iis)*dy
  delpy = take(delpy,iis)*dy
  numbmx = take(numbmx,iis)
  numbpx = take(numbpx,iis)
  numbmy = take(numbmy,iis)
  numbpy = take(numbpy,iis)
  if lparallel:
    ixc = gatherarray(ixc)
    iyc = gatherarray(iyc)
    xxc = gatherarray(xxc)
    yyc = gatherarray(yyc)
    numb = gatherarray(numb)
    ixs = gatherarray(ixs)
    iys = gatherarray(iys)
    xxs = gatherarray(xxs)
    yys = gatherarray(yys)
    delmx = gatherarray(delmx)
    delpx = gatherarray(delpx)
    delmy = gatherarray(delmy)
    delpy = gatherarray(delpy)
    numbmx = gatherarray(numbmx)
    numbpx = gatherarray(numbpx)
    numbmy = gatherarray(numbmy)
    numbpy = gatherarray(numbpy)
  # --- Now, after the global gather, if there is no data then quit.
  if len(ixc) + len(ixs) == 0: return
  # --- Get max grid point so that an array can be created which covers all
  # --- of the data.
  if len(ixc) == 0:
    maxixc = 0
    maxiyc = 0
  else:
    maxixc = max(ixc)
    maxiyc = max(iyc)
  if len(ixs) == 0:
    maxixs = 0
    maxiys = 0
  else:
    maxixs = max(ixs)
    maxiys = max(iys)
  nx = max(maxixc,maxixs) + 1
  ny = max(maxiyc,maxiys) + 1
  iii = zeros((5,1+nx,1+ny)) - 1
  mx,px,my,py = 1,2,3,4
  # --- Flag grid points where the conductors are.
  for i in xrange(len(ixc)): iii[0,ixc[i],iyc[i]] = i
  for i in xrange(len(ixs)):
    if abs(delmx[i]) < abs(dx): iii[mx,ixs[i],iys[i]] = i
    if abs(delpx[i]) < abs(dx): iii[px,ixs[i],iys[i]] = i
    if abs(delmy[i]) < abs(dy): iii[my,ixs[i],iys[i]] = i
    if abs(delpy[i]) < abs(dy): iii[py,ixs[i],iys[i]] = i
  x = []
  y = []
  z = []
  n = []

  # --- Loop over all conductor points, drawing a fill polygon for each.
  ixall = list(ixc) + list(ixs)
  iyall = list(iyc) + list(iys)
  for ix,iy in map(None,ixall,iyall):
    # --- For each point, deal with the grid cell up and to the right.
    # --- This is done so that no two grid points draw over the same spot.
    i1 = _corner1(ix  ,iy  ,px,py,+delpx,+delpy,numb,numbpx,numbpy,0,
                  iii,xxc,yyc,xxs,yys,x,y,z,n)
    i2 = _corner1(ix  ,iy+1,px,my,+delpx,-delmy,numb,numbpx,numbmy,1,
                  iii,xxc,yyc,xxs,yys,x,y,z,n)
    i3 = _corner1(ix+1,iy+1,mx,my,-delmx,-delmy,numb,numbmx,numbmy,0,
                  iii,xxc,yyc,xxs,yys,x,y,z,n)
    i4 = _corner1(ix+1,iy  ,mx,py,-delmx,+delpy,numb,numbmx,numbpy,1,
                  iii,xxc,yyc,xxs,yys,x,y,z,n)
    if i1[0] + i2[0] + i3[0] + i4[0] > 0:
      n.append(i1[0] + i2[0] + i3[0] + i4[0])
      if   i1[1] is not None: z.append(i1[1])
      elif i2[1] is not None: z.append(i2[1])
      elif i3[1] is not None: z.append(i3[1])
      elif i4[1] is not None: z.append(i4[1])
    # --- One special case needs to be dealt with. If there is no grid
    # --- point in the lower left corner of a grid cell, it would otherwise
    # --- be skipped.
    if (ix > 0 and iii[py,ix,iy] >= 0 and max(iii[:,ix-1,iy]) == -1 and
                   iii[px,ix-1,iy+1] >= 0 and iii[0,ix,iy+1] >= 0):
      i2 = _corner1(ix-1,iy+1,px,my,+delpx,-delmy,numb,numbpx,numbmy,1,
                    iii,xxc,yyc,xxs,yys,x,y,z,n)
      i3 = _corner1(ix  ,iy+1,mx,my,-delmx,-delmy,numb,numbmx,numbmy,0,
                    iii,xxc,yyc,xxs,yys,x,y,z,n)
      i4 = _corner1(ix  ,iy  ,mx,py,-delmx,+delpy,numb,numbmx,numbpy,1,
                    iii,xxc,yyc,xxs,yys,x,y,z,n)
      if i2[0] + i3[0] + i4[0] > 0:
        n.append(i2[0] + i3[0] + i4[0])
        if   i2[1] is not None: z.append(i2[1])
        elif i3[1] is not None: z.append(i3[1])
        elif i4[1] is not None: z.append(i4[1])

  if len(n) > 0:
    # --- Now that the data is gathered, make the plot.
    z = array(z)
    z = z.astype('b')
    plfp(z,y,x,n)

def _corner1(ix,iy,sx,sy,delsx,delsy,nm,nmsx,nmsy,iparity,
             iii,xxc,yyc,xxs,yys,x,y,z,n):
  """
Special routine which gathers the data to draw a polygon in the local grid
cell. It must be called for each of the four corners of a grid cell.
  """
  nn = 0
  numb = None
  if iii[0,ix,iy] >= 0:
    nn = nn + 1
    i0 = iii[0,ix,iy]
    x.append(xxc[i0])
    y.append(yyc[i0])
    numb = nm[i0]
  for i in range(2):
    if i == iparity:
      if iii[sx,ix,iy] >= 0:
        nn = nn + 1
        i0 = iii[sx,ix,iy]
        x.append(xxs[i0]+delsx[i0])
        y.append(yys[i0])
        numb = nmsx[i0]
    else:
      if iii[sy,ix,iy] >= 0:
        nn = nn + 1
        i0 = iii[sy,ix,iy]
        x.append(xxs[i0])
        y.append(yys[i0]+delsy[i0])
        numb = nmsy[i0]
  return nn,numb

######################################################################
######################################################################
# x-y plane
def pfxy(iz=None,izf=None,fullplane=1,
         cond=1,plotsg=1,fill=0,scale=1,plotphi=1,
         subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,mglevel=0,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in X-Y plane
  - iz=w3d.iz_axis z index of plane
  - fullplane=1: when true, plot all quandrants regardless of symmetries
  - cond=1: when true, plot grid point inside of conductors
  - plotsg=1: when true, plots subgrid data
  - fill=0: when true, fills in conductor
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  - mglevel=0: level of multigrid to plot data for
  """
  kw.update(kwdict)
  # --- This logic is needed since in the parallel version, iz_axis already
  # --- has izslave subtracted from it. If the user passes in a value,
  # --- it must be checked for consistency, otherwise coding below could lead
  # --- to a deadlock in the parallel version
  if izf is not None: iz = izf
  if iz is None: iz = w3d.iz_axis + top.izslave[me]
  if iz < 0 or w3d.nzfull < iz: return
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
  if fill:
    plotcondfill('y','x','z',iz,ymmin,xmmin,dy,dx,mglevel,1,1)
    if fullplane and (w3d.l2symtry or w3d.l4symtry):
      plotcondfill('y','x','z',iz,ymmin,xmmin,dy,dx,mglevel,1,-1)
    if fullplane and w3d.l4symtry:
      plotcondfill('y','x','z',iz,ymmin,xmmin,dy,dx,mglevel,-1,1)
      plotcondfill('y','x','z',iz,ymmin,xmmin,dy,dx,mglevel,-1,-1)
  if cond:
    plotcond('y','x','z',iz,numb,ymmin,xmmin,dy,dx,condcolor,mglevel,1,1)
    if fullplane and (w3d.l2symtry or w3d.l4symtry):
      plotcond('y','x','z',iz,numb,ymmin,xmmin,dy,dx,condcolor,mglevel,-1,1)
    if fullplane and w3d.l4symtry:
      plotcond('y','x','z',iz,numb,ymmin,xmmin,dy,dx,condcolor,mglevel,1,-1)
      plotcond('y','x','z',iz,numb,ymmin,xmmin,dy,dx,condcolor,mglevel,-1,-1)
  if plotsg:
    plotsubgrid('y','x','z',0,iz,numb,ymmin,xmmin,dy,dx,evencolor,
                subgridlen,mglevel,1,1)
    plotsubgrid('y','x','z',1,iz,numb,ymmin,xmmin,dy,dx,oddcolor,
                subgridlen,mglevel,1,1)
    if fullplane and (w3d.l2symtry or w3d.l4symtry):
      plotsubgrid('y','x','z',0,iz,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,mglevel,1,-1)
      plotsubgrid('y','x','z',1,iz,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,mglevel,1,-1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('y','x','z',0,iz,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,mglevel,-1,1)
      plotsubgrid('y','x','z',1,iz,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,mglevel,-1,1)
      plotsubgrid('y','x','z',0,iz,numb,ymmin,xmmin,dy,dx,evencolor,
                  subgridlen,mglevel,-1,-1)
      plotsubgrid('y','x','z',1,iz,numb,ymmin,xmmin,dy,dx,oddcolor,
                  subgridlen,mglevel,-1,-1)

# z-x plane
def pfzx(iy=None,iyf=None,fullplane=1,lbeamframe=1,
         cond=1,plotsg=1,fill=0,scale=1,plotphi=1,
         subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,mglevel=0,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in Z-X plane
  - iy=w3d.iy_axis y index of plane
  - fullplane=1: when true, plot all quadrants regardless of symmetries
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  - cond=1: when true, plot grid point inside of conductors
  - plotsg=1: when true, plots subgrid data
  - fill=0: when true, fills in conductor
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  - mglevel=0: level of multigrid to plot data for
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
    zmmin = w3d.zmminglobal + zbeam
  else:
    dx = 1.
    dz = 1.
    xmmin = 0.
    zmmin = 0.
  if plotphi:
    if not scale:
      kw['xmin'] = 0
      kw['xmax'] = w3d.nzfull
      kw['ymin'] = 0.
      kw['ymax'] = w3d.nx
    if not kw.has_key('ccolor'): kw['ccolor'] = phicolor
    apply(pcphizx,(iy,fullplane,lbeamframe),kw)
  if fill:
    plotcondfill('x','z','y',iy,xmmin,zmmin,dx,dz,mglevel,1,1)
    if fullplane and (w3d.l4symtry or w3d.solvergeom == w3d.RZgeom):
      plotcondfill('x','z','y',iy,xmmin,zmmin,dx,dz,mglevel,-1,1)
  if plotsg:
    plotsubgrid('x','z','y',0,iy,numb,xmmin,zmmin,dx,dz,evencolor,
                subgridlen,mglevel,1,1)
    plotsubgrid('x','z','y',1,iy,numb,xmmin,zmmin,dx,dz,oddcolor,
                subgridlen,mglevel,1,1)
    if fullplane and (w3d.l4symtry or w3d.solvergeom == w3d.RZgeom):
      plotsubgrid('x','z','y',0,iy,numb,xmmin,zmmin,dx,dz,evencolor,
                  subgridlen,mglevel,-1,1)
      plotsubgrid('x','z','y',1,iy,numb,xmmin,zmmin,dx,dz,oddcolor,
                  subgridlen,mglevel,-1,1)
  if cond:
    plotcond('x','z','y',iy,numb,xmmin,zmmin,dx,dz,condcolor,mglevel,1,1)
    if fullplane and (w3d.l4symtry or w3d.solvergeom == w3d.RZgeom):
      plotcond('x','z','y',iy,numb,xmmin,zmmin,dx,dz,condcolor,mglevel,-1,1)

# z-y plane
def pfzy(ix=None,ixf=None,fullplane=1,lbeamframe=1,
         cond=1,plotsg=1,fill=0,scale=1,plotphi=1,
         subgridlen=1.,phicolor=blue,condcolor=cyan,
         oddcolor=red,evencolor=green,numb=None,mglevel=0,kwdict={},**kw):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane
  - ix=w3d.ix_axis x index of plane
  - fullplane=1: when true, plot all quadrants regardless of symmetries
  - lbeamframe=1: when true, plot relative to beam frame, otherwise lab frame
  - cond=1: when true, plot grid point inside of conductors
  - plotsg=1: when true, plots subgrid data
  - fill=0: when true, fills in conductor
  - scale=1 when true, plots data in lab frame, otherwise grid frame
  - plotphi=1 when true, plot contours of potential
  - phicolor=blue color of phi contours
  - condcolor=cyan color of conductor points inside conductors
  - oddcolor=red color of odd subgrid points
  - evencolor=green color of even subgrid points
  - subgridlen=1 maximum length of subgrid line which are plotted
  - numb: specify which conductors to plot based on the conductor number
  - mglevel=0: level of multigrid to plot data for
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
    zmmin = w3d.zmminglobal + zbeam
  else:
    dy = 1.
    dz = 1.
    ymmin = 0.
    zmmin = 0.
  if plotphi:
    if not scale:
      kw['xmin'] = 0
      kw['xmax'] = w3d.nzfull
      kw['ymin'] = 0
      kw['ymax'] = w3d.ny
    if not kw.has_key('ccolor'): kw['ccolor'] = phicolor
    apply(pcphizy,(ix,fullplane,lbeamframe),kw)
  if fill:
    plotcondfill('y','z','x',ix,ymmin,zmmin,dy,dz,mglevel,1,1)
    if fullplane and w3d.l4symtry:
      plotcondfill('y','z','x',ix,ymmin,zmmin,dy,dz,mglevel,-1,1)
  if cond:
    plotcond('y','z','x',ix,numb,ymmin,zmmin,dy,dz,condcolor,mglevel,1,1)
    if fullplane and (w3d.l2symtry or w3d.l4symtry):
      plotcond('y','z','x',ix,numb,ymmin,zmmin,dy,dz,condcolor,mglevel,-1,1)
  if plotsg:
    plotsubgrid('y','z','x',0,ix,numb,ymmin,zmmin,dy,dz,evencolor,
                subgridlen,mglevel,1,1)
    plotsubgrid('y','z','x',1,ix,numb,ymmin,zmmin,dy,dz,oddcolor,
                subgridlen,mglevel,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('y','z','x',0,ix,numb,ymmin,zmmin,dy,dz,evencolor,
                  subgridlen,mglevel,-1,1)
      plotsubgrid('y','z','x',1,ix,numb,ymmin,zmmin,dy,dz,oddcolor,
                  subgridlen,mglevel,-1,1)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in grid units                                                      #
######################################################################

# x-y plane
def pfxyg(iz=None,izf=None,fullplane=1,
          cond=1,plotsg=1,fill=0,plotphi=1,
          phicolor=blue,subgridlen=1.,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,mglevel=0,**kw):
  """
Plots conductors and contours of electrostatic potential in X-Y plane in grid
frame
Same arguments as pfxy
  """
  if izf is not None: iz = izf
  pfxy(iz=iz,fullplane=fullplane,scale=0,
       cond=cond,plotsg=plotsg,fill=fill,plotphi=plotphi,
       subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,mglevel=mglevel,kwdict=kw)

# z-x plane
def pfzxg(iy=None,iyf=None,fullplane=1,lbeamframe=1,
          cond=1,plotsg=1,fill=0,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,mglevel=0,**kw):
  """
Plots conductors and contours of electrostatic potential in Z-X plane in grid
frame
Same arguments as pfzx
  """
  if iyf is not None: iy = iyf
  pfzx(iy=iy,fullplane=fullplane,lbeamframe=lbeamframe,scale=0,
       cond=cond,plotsg=plotsg,fill=fill,plotphi=plotphi,
       subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,mglevel=mglevel,kwdict=kw)

# z-y plane
def pfzyg(ix=None,ixf=None,fullplane=1,lbeamframe=1,
          cond=1,plotsg=1,fill=0,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,mglevel=0,**kw):
  """
Plots conductors and contours of electrostatic potential in Z-Y plane in grid
frame
Same arguments as pfzy
  """
  if ixf is not None: ix = ixf
  pfzy(ix=ix,fullplane=fullplane,lbeamframe=lbeamframe,scale=0,
       cond=cond,plotsg=plotsg,fill=fill,plotphi=plotphi,
       subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,mglevel=mglevel,kwdict=kw)

######################################################################
# handy functions to plot the conductor points and subgrid data      #
# in real units, with inverted x.  Used to make complete plots when  #
# 2-fold symmetry is used.                                           #
######################################################################
 
# x-y plane
def pfxyi(iz=None,izf=None,fullplane=1,
          cond=1,plotsg=1,fill=0,scale=1,plotphi=1,
          phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full X-Y plane,
Same arguments as pfxy
  """
  print "Notice: pfxyi is obsolete is should no longer be used"
  print "        It does the identical thing as pfxy"
  if izf is not None: iz = izf
  pfxy(iz=iz,fullplane=fullplane,scale=scale,
       cond=cond,plotsg=plotsg,fill=fill,plotphi=plotphi,
       subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

# z-x plane
def pfzxi(iy=None,iyf=None,fullplane=1,lbeamframe=1,
          cond=1,plotsg=1,fill=0,scale=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full Z-X plane
Same arguments as pfzx
  """
  print "Notice: pfzxi is obsolete is should no longer be used"
  print "        It does the identical thing as pfzx"
  if iyf is not None: iy = iyf
  pfzx(iy=iy,fullplane=fullplane,lbeamframe=lbeamframe,
       cond=cond,plotsg=plotsg,fill=fill,plotphi=plotphi,
       scale=scale,subgridlen=subgridlen,
       phicolor=phicolor,condcolor=condcolor,
       oddcolor=oddcolor,evencolor=evencolor,numb=numb,kwdict=kw)

# z-y plane
def pfzyi(ix=None,ixf=None,fullplane=1,lbeamframe=1,
          cond=1,plotsg=1,fill=0,scale=1,plotphi=1,
          subgridlen=1.,phicolor=blue,condcolor=cyan,
          oddcolor=red,evencolor=green,numb=None,**kw):
  """
Plots conductors and contours of electrostatic potential in full Z-Y plane
Same arguments as pfzy
  """
  print "Notice: pfzyi is obsolete is should no longer be used"
  print "        It does the identical thing as pfzy"
  if ixf is not None: ix = ixf
  pfzy(ix=ix,fullplane=fullplane,lbeamframe=lbeamframe,
       cond=cond,plotsg=plotsg,fill=fill,plotphi=plotphi,
       scale=scale,subgridlen=subgridlen,
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
    izc = f3d.izcond[0:f3d.ncond]+f3d.mglevelsiz[0]
    ii = compress(equal(izc,iz),arange(f3d.ncond))
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
    zmmin = w3d.zmmin
  else:
    dx = 1.*signx
    dz = 1.*signz
    xmmin = 0.
    zmmin = 0.
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
    zmmin = w3d.zmmin
  else:
    dy = 1.*signy
    dz = 1.*signz
    ymmin = 0.
    zmmin = 0.
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

def plotcondn(yy,xx,zz,iz,ymmin,xmmin,dy,dx,mglevel,signy,signx):
  ncolor = len(color)
  if f3d.ncond > 0: nn = f3d.condnumb
  else:             nn = array([])
  nlist = gatherarray(nn)
  nlist = findunique(nlist)
  nlist = broadcast(nlist)
  for i in nlist:
    plotcond(yy,xx,zz,iz,i,ymmin,xmmin,dy,dx,color[i%ncolor],
             mglevel,signy,signx)

def pfzxn(iy=None,numbs=None,colors=None,cmarker=point,smarker=circle,
          scale=1,signz=1,signx=1,subgridlen=1.,fullplane=1,mglevel=0):
  if iy is None: iy = w3d.iy_axis
  if iy < 0 or w3d.ny < iy: return
  if colors is None: colors = color
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
  plotcondn('x','z','y',iy,xmmin,zmmin,dx,dz,mglevel,1,1)
  if fullplane and w3d.l4symtry:
    plotcondn('x','z','y',iy,xmmin,zmmin,dx,dz,mglevel,-1,1)
  ncolor = len(colors)
  nlist = gatherarray(f3d.ecnumb[:f3d.necndbdy])
  nlist = findunique(nlist)
  #nlist.remove(0)
  nlist = broadcast(nlist)
  for i in range(len(nlist)):
    plotsubgrid('x','z','y',0,iy,nlist[i],xmmin,zmmin,dx,dz,
                colors[i%ncolor],subgridlen,mglevel,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('x','z','y',0,iy,nlist[i],xmmin,zmmin,dx,dz,
                  colors[i%ncolor],subgridlen,mglevel,-1,1)
  nlist = gatherarray(f3d.ocnumb[:f3d.nocndbdy])
  nlist = findunique(nlist)
  nlist = broadcast(nlist)
  for i in range(len(nlist)):
    plotsubgrid('x','z','y',1,iy,nlist[i],xmmin,zmmin,dx,dz,
                colors[i%ncolor],subgridlen,mglevel,1,1)
    if fullplane and w3d.l4symtry:
      plotsubgrid('x','z','y',1,iy,nlist[i],xmmin,zmmin,dx,dz,
                  colors[i%ncolor],subgridlen,mglevel,-1,1)


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
                        f3d.izcond[0:f3d.ncond])*w3d.dz+w3d.zmmin+zz
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
def plotsrfrvinout(srfrvin,srfrvout,zmin,zmax,n=1000,color='fg',gridframe=0,
                   rscale=1,zscale=1,
                   roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0,
                   outline=1,fillcolor=None):
  """Handy function for plotting the r versus z for a surface of revolution
 - srfrvin: surface of revolution function to plot
 - srfrvou: surface of revolution function to plot
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
 - outline=1: when true, draw outline
 - fillcolor=None: optionally sets fill color
  """
  zz = iota(0,n)*(zmax - zmin)/n + zmin
  rrin = ones(n+1,'d')
  rrout = ones(n+1,'d')
  for i in range(n+1):
    f3d.srfrv_z = zz[i]
    srfrvin()
    rrin[i] = f3d.srfrv_r
    srfrvout()
    rrout[i] = f3d.srfrv_r
  if gridframe:
    zz = (zz - w3d.zmmin)/w3d.dz
    rrin = (rrin)/w3d.dx + ir_axis
    rrout = (rrout)/w3d.dx + ir_axis
  rr = array(list(rrin) + list(rrout[::-1]))
  zz = array(list(zz) + list(zz[::-1]))
  rr = where(less(rr,rmin),rmin,rr)
  rr = where(greater(rr,rmax),rmax,rr)
  if outline:
    plg(rscale*rr+roff,zscale*zz+zoff,color=color)
  if fillcolor is not None:
    cc = array([fillcolor]).astype('b')
    plfp(cc,rscale*rr+roff,zscale*zz+zoff,[2*n+2])

#####################################################################
def plotsrfrvin(srfrvin,zmin,zmax,n=1000,color='fg',gridframe=0,
                rscale=1,zscale=1,
                roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0,
                outline=1,fillcolor=None):
  """Handy function for plotting the r versus z for a surface of revolution
 - srfrvin: surface of revolution function to plot
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
 - outline=1: when true, draw outline
 - fillcolor=None: optionally sets fill color
  """
  zz = iota(0,n)*(zmax - zmin)/n + zmin
  rr = ones(n+1,'d')
  for i in range(n+1):
    f3d.srfrv_z = zz[i]
    srfrvin()
    rr[i] = f3d.srfrv_r
  if gridframe:
    zz = (zz - w3d.zmmin)/w3d.dz
    rr = (rr)/w3d.dx + ir_axis
  rr = where(less(rr,rmin),rmin,rr)
  rr = where(greater(rr,rmax),rmax,rr)
  rr = rscale*rr+roff
  zz = zscale*zz+zoff
  if outline:
    plg(rr,zz,color=color)
  if fillcolor is not None:
    cc = array([fillcolor]).astype('b')
    rr = array(list(rr) + [rmin,rmin])
    zz = array(list(zz) + [zmax,zmin])
    plfp(cc,rr,zz,[n+3])

#####################################################################
def plotsrfrvout(srfrvin,zmin,zmax,n=1000,color='fg',gridframe=0,
                 rscale=1,zscale=1,
                 roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0,
                 outline=1,fillcolor=None):
  """Handy function for plotting the r versus z for a surface of revolution
 - srfrvout: surface of revolution function to plot
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
 - outline=1: when true, draw outline
 - fillcolor=None: optionally sets fill color
  """
  zz = iota(0,n)*(zmax - zmin)/n + zmin
  rr = ones(n+1,'d')
  for i in range(n+1):
    f3d.srfrv_z = zz[i]
    srfrvout()
    rr[i] = f3d.srfrv_r
  if gridframe:
    zz = (zz - w3d.zmmin)/w3d.dz
    rr = (rr)/w3d.dx + ir_axis
  rr = where(less(rr,rmin),rmin,rr)
  rr = where(greater(rr,rmax),rmax,rr)
  rr = rscale*rr+roff
  zz = zscale*zz+zoff
  if outline:
    plg(rr,zz,color=color)
  if fillcolor is not None:
    cc = array([fillcolor]).astype('b')
    rmax = min(rmax,w3d.xmmax)
    rmax = max(rmax,max(rr))
    rr = array(list(rr) + [rmax,rmax])
    zz = array(list(zz) + [zmax,zmin])
    plfp(cc,rr,zz,[n+3])

#####################################################################
def plotsrfrv(srfrv,zmin,zmax,n=1000,color='fg',gridframe=0,rscale=1,zscale=1,
              roff=0,zoff=0,rmin=0.,rmax=top.largepos,ir_axis=0,
              outline=1,fillcolor=None):
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
 - outline=1: when true, draw outline
 - fillcolor=None: optionally sets fill color
  """
  zz = iota(0,n)*(zmax - zmin)/n + zmin
  rr = ones(n+1,'d')
  for i in range(n+1):
    f3d.srfrv_z = zz[i]
    srfrv()
    rr[i] = f3d.srfrv_r
  if gridframe:
    zz = (zz - w3d.zmmin)/w3d.dz
    rr = (rr)/w3d.dx + ir_axis
  rr = where(less(rr,rmin),rmin,rr)
  rr = where(greater(rr,rmax),rmax,rr)
  if outline:
    plg(rscale*rr+roff,zscale*zz+zoff,color=color)
  if fillcolor is not None:
    cc = array([fillcolor]).astype('b')
    plfp(cc,rscale*rr+roff,zscale*zz+zoff,[n+1])


#####################################################################
#####################################################################
def plotelementoutline(color,gridframe,axis,ie,ne,outline,fillcolor,
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
        zz = (zz - w3d.zmmin)/w3d.dz
      if outline:
        plg(rr1,zz,color=color)
        plg(rr2,zz,color=color)
      if fillcolor is not None:
        cc = array([fillcolor]).astype('b')
        plfp(cc,rr1,zz,[5])
        plfp(cc,rr2,zz,[5])
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
      if outline:
        plg(rrl1,zzl,color=color)
        plg(rrl2,zzl,color=color)
        plg(rrr1,zzr,color=color)
        plg(rrr2,zzr,color=color)
      if fillcolor is not None:
        cc = array([fillcolor]).astype('b')
        plfp(cc,rrl1,zzl,[5])
        plfp(cc,rrl2,zzl,[5])
        plfp(cc,rrr1,zzr,[5])
        plfp(cc,rrr2,zzr,[5])


#---------------------------------------------------------------------------
def plotquadoutline(iq=0,nq=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
  """Plots the outline of quadrupole elements
 - iq=0: starting quad to plot
 - nq=top.nquad+1: number of quads to plot
 - color='fg': line color
 - gridframe=0: when true, make plot in grid coordinates
 - axis='x': selects axis to plot, either 'x' or 'y'
 - outline=1: when true, draw outline
 - fillcolor=None: optionally sets fill color
  """
  if nq is None: nq = top.nquad + 1
  plotelementoutline(color,gridframe,axis,iq,nq,outline,fillcolor,
                     top.quadzs,top.quadze,top.quadap,top.quadrr,top.quadrl,
                     top.quadgl,top.quadgp,top.qoffx,top.qoffy,
                     top.quadpa,top.quadpr,top.quadpw,
                     top.qdelpal,top.qdelpar)

#---------------------------------------------------------------------------
def plotheleoutline(ih=0,nh=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
  """Plots the outline of hele elements
 - ih=0: starting hele to plot
 - nh=top.nhele+1: number of heles to plot
 - color='fg': line color
 - gridframe=0: when true, make plot in grid coordinates
 - axis='x': selects axis to plot, either 'x' or 'y'
 - outline=1: when true, draw outline
 - fillcolor=None: optionally sets fill color
  """
  if nh is None: nh = top.nhele + 1
  plotelementoutline(color,gridframe,axis,ih,nh,outline,fillcolor,
                     top.helezs,top.heleze,top.heleap,top.helerr,top.helerl,
                     top.helegl,top.helegp,top.heleox,top.heleoy,
                     top.helepa,zeros(top.nhele+1,'d'),top.helepw,
                     zeros(top.nhele+1,'d'),zeros(top.nhele+1,'d'))

#---------------------------------------------------------------------------
def plotemltoutline(ie=0,ne=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
  """Plots the outline of emlt elements
 - ie=0: starting emlt to plot
 - ne=top.nemlt+1: number of emlts to plot
 - color='fg': line color
 - gridframe=0: when true, make plot in grid coordinates
 - axis='x': selects axis to plot, either 'x' or 'y'
 - outline=1: when true, draw outline
 - fillcolor=None: optionally sets fill color
  """
  if ne is None: ne = top.nemlt + 1
  plotelementoutline(color,gridframe,axis,ie,ne,outline,fillcolor,
                     top.emltzs,top.emltze,top.emltap,top.emltrr,top.emltrl,
                     top.emltgl,top.emltgp,top.emltox,top.emltoy,
                     top.emltpa,zeros(top.nemlt+1,'d'),top.emltpw,
                     zeros(top.nemlt+1,'d'),zeros(top.nemlt+1,'d'))

#---------------------------------------------------------------------------
def plotpgrdoutline(ie=0,ne=None,color='fg',gridframe=0,axis='x',
                    outline=1,fillcolor=None):
  """Plots the outline of pgrd elements
 - ie=0: starting pgrd to plot
 - ne=top.npgrd+1: number of pgrds to plot
 - color='fg': line color
 - gridframe=0: when true, make plot in grid coordinates
 - axis='x': selects axis to plot, either 'x' or 'y'
 - outline=1: when true, draw outline
 - fillcolor=None: optionally sets fill color
  """
  if ne is None: ne = top.npgrd + 1
  plotelementoutline(color,gridframe,axis,ie,ne,outline,fillcolor,
                     top.pgrdzs,top.pgrdze,top.pgrdap,top.pgrdrr,top.pgrdrl,
                     top.pgrdgl,top.pgrdgp,top.pgrdox,top.pgrdoy,
                     top.pgrdpa,zeros(top.npgrd+1,'d'),top.pgrdpw,
                     zeros(top.npgrd+1,'d'),zeros(top.npgrd+1,'d'))

#########################################################################
#########################################################################
def cleanconductors():
  """This routine clears out conductor points which are not within the
range of the field solution, w3d.izfsmin and w3d.izfsmax. This is done
for optimization so that time is not wasted on those points.
  """
  # --- This routine would corrupt the data for the multigrid fieldsolver.
  if top.fstype == 7: return
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
    gchange("Conductor3d")
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
    gchange("Conductor3d")
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
    gchange("Conductor3d")
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
  # --- Check if one of the MG conductor arrays is allocated. If not, then
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
  icondlxynew = []
  icondlznew = []
  icndxnew = []
  icndynew = []
  icndznew = []
  cdelmxnew = []
  cdelmynew = []
  cdelmznew = []
  cdelpxnew = []
  cdelpynew = []
  cdelpznew = []
  cvoltnew = []
  cvoltmxnew = []
  cvoltpxnew = []
  cvoltmynew = []
  cvoltpynew = []
  cvoltmznew = []
  cvoltpznew = []
  icndlxynew = []
  icndlznew = []

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
    icndxnew = icndxnew + list(take(iecndx/llxy[j],ii))
    icndynew = icndynew + list(take(iecndy/llxy[j],ii))
    icndznew = icndznew + list(take(iecndz/llz[j],ii))
    cdelmxnew = cdelmxnew + list(take(ecdelmx/llxy[j],ii))
    cdelmynew = cdelmynew + list(take(ecdelmy/llxy[j],ii))
    cdelmznew = cdelmznew + list(take(ecdelmz/llz[j],ii))
    cdelpxnew = cdelpxnew + list(take(ecdelpx/llxy[j],ii))
    cdelpynew = cdelpynew + list(take(ecdelpy/llxy[j],ii))
    cdelpznew = cdelpznew + list(take(ecdelpz/llz[j],ii))
    cvoltnew = cvoltnew + list(take(ecvolt,ii))
    cvoltmxnew = cvoltmxnew + list(take(ecvoltmx,ii))
    cvoltpxnew = cvoltpxnew + list(take(ecvoltpx,ii))
    cvoltmynew = cvoltmynew + list(take(ecvoltmy,ii))
    cvoltpynew = cvoltpynew + list(take(ecvoltpy,ii))
    cvoltmznew = cvoltmznew + list(take(ecvoltmz,ii))
    cvoltpznew = cvoltpznew + list(take(ecvoltpz,ii))
    icndlxynew = icndlxynew + list(llxy[j]*ones(len(ii)))
    icndlznew = icndlznew + list(llz[j]*ones(len(ii)))
  
    ii = compress(logical_and(iocndlxy >= llxy[j],iocndlz >= llz[j]), \
                  arange(nocndbdy))
    icndxnew = icndxnew + list(take(iocndx/llxy[j],ii))
    icndynew = icndynew + list(take(iocndy/llxy[j],ii))
    icndznew = icndznew + list(take(iocndz/llz[j],ii))
    cdelmxnew = cdelmxnew + list(take(ocdelmx/llxy[j],ii))
    cdelmynew = cdelmynew + list(take(ocdelmy/llxy[j],ii))
    cdelmznew = cdelmznew + list(take(ocdelmz/llz[j],ii))
    cdelpxnew = cdelpxnew + list(take(ocdelpx/llxy[j],ii))
    cdelpynew = cdelpynew + list(take(ocdelpy/llxy[j],ii))
    cdelpznew = cdelpznew + list(take(ocdelpz/llz[j],ii))
    cvoltnew = cvoltnew + list(take(ocvolt,ii))
    cvoltmxnew = cvoltmxnew + list(take(ocvoltmx,ii))
    cvoltpxnew = cvoltpxnew + list(take(ocvoltpx,ii))
    cvoltmynew = cvoltmynew + list(take(ocvoltmy,ii))
    cvoltpynew = cvoltpynew + list(take(ocvoltpy,ii))
    cvoltmznew = cvoltmznew + list(take(ocvoltmz,ii))
    cvoltpznew = cvoltpznew + list(take(ocvoltpz,ii))
    icndlxynew = icndlxynew + list(llxy[j]*ones(len(ii)))
    icndlznew = icndlznew + list(llz[j]*ones(len(ii)))

  # --- Find odd and even points
  icndxnew = array(icndxnew,typecode='i')
  icndynew = array(icndynew,typecode='i')
  icndznew = array(icndznew,typecode='i')
  ee = compress((icndxnew+icndynew+icndznew)%2==0,arange(len(icndxnew)))
  oo = compress((icndxnew+icndynew+icndznew)%2==1,arange(len(icndxnew)))

  # --- Reset counters and reallocate arrays for the converted data.
  f3d.ncondmax = len(ixcondnew)
  f3d.ncndmax = max(len(ee),len(oo))
  gchange("Conductor3d")
  f3d.ncond = len(ixcondnew)
  f3d.necndbdy = len(ee)
  f3d.nocndbdy = len(oo)
  
  # --- Copy converted data into WARP arrays.
  f3d.ixcond[:f3d.ncond] = ixcondnew
  f3d.iycond[:f3d.ncond] = iycondnew
  f3d.izcond[:f3d.ncond] = izcondnew
  f3d.condvolt[:f3d.ncond] = condvoltnew
  f3d.icondlevel[:f3d.ncond] = nint(log(icondlxynew)/log(2.)+0.5)
  #f3d.icondlxy[:f3d.ncond] = icondlxynew
  #f3d.icondlz[:f3d.ncond] = icondlznew

  f3d.iecndx[:f3d.necndbdy] = take(icndxnew,ee)
  f3d.iecndy[:f3d.necndbdy] = take(icndynew,ee)
  f3d.iecndz[:f3d.necndbdy] = take(icndznew,ee)
  f3d.ecdelmx[:f3d.necndbdy] = take(cdelmxnew,ee)
  f3d.ecdelmy[:f3d.necndbdy] = take(cdelmynew,ee)
  f3d.ecdelmz[:f3d.necndbdy] = take(cdelmznew,ee)
  f3d.ecdelpx[:f3d.necndbdy] = take(cdelpxnew,ee)
  f3d.ecdelpy[:f3d.necndbdy] = take(cdelpynew,ee)
  f3d.ecdelpz[:f3d.necndbdy] = take(cdelpznew,ee)
  f3d.ecvolt[:f3d.necndbdy] = take(cvoltnew,ee)
  f3d.ecvoltmx[:f3d.necndbdy] = take(cvoltmxnew,ee)
  f3d.ecvoltpx[:f3d.necndbdy] = take(cvoltpxnew,ee)
  f3d.ecvoltmy[:f3d.necndbdy] = take(cvoltmynew,ee)
  f3d.ecvoltpy[:f3d.necndbdy] = take(cvoltpynew,ee)
  f3d.ecvoltmz[:f3d.necndbdy] = take(cvoltmznew,ee)
  f3d.ecvoltpz[:f3d.necndbdy] = take(cvoltpznew,ee)
  f3d.iecndlevel[:f3d.ncond] = nint(log(iecndlxynew)/log(2.)+0.5)
  #f3d.iecndlxy[:f3d.necndbdy] = take(icndlxynew,ee)
  #f3d.iecndlz[:f3d.necndbdy] = take(icndlznew,ee)

  f3d.iocndx[:f3d.nocndbdy] = take(icndxnew,oo)
  f3d.iocndy[:f3d.nocndbdy] = take(icndynew,oo)
  f3d.iocndz[:f3d.nocndbdy] = take(icndznew,oo)
  f3d.ocdelmx[:f3d.nocndbdy] = take(cdelmxnew,oo)
  f3d.ocdelmy[:f3d.nocndbdy] = take(cdelmynew,oo)
  f3d.ocdelmz[:f3d.nocndbdy] = take(cdelmznew,oo)
  f3d.ocdelpx[:f3d.nocndbdy] = take(cdelpxnew,oo)
  f3d.ocdelpy[:f3d.nocndbdy] = take(cdelpynew,oo)
  f3d.ocdelpz[:f3d.nocndbdy] = take(cdelpznew,oo)
  f3d.ocvolt[:f3d.nocndbdy] = take(cvoltnew,oo)
  f3d.ocvoltmx[:f3d.nocndbdy] = take(cvoltmxnew,oo)
  f3d.ocvoltpx[:f3d.nocndbdy] = take(cvoltpxnew,oo)
  f3d.ocvoltmy[:f3d.nocndbdy] = take(cvoltmynew,oo)
  f3d.ocvoltpy[:f3d.nocndbdy] = take(cvoltpynew,oo)
  f3d.ocvoltmz[:f3d.nocndbdy] = take(cvoltmznew,oo)
  f3d.ocvoltpz[:f3d.nocndbdy] = take(cvoltpznew,oo)
  f3d.iocndlevel[:f3d.ncond] = nint(log(iocndlxynew)/log(2.)+0.5)
  #f3d.iocndlxy[:f3d.nocndbdy] = take(icndlxynew,oo)
  #f3d.iocndlz[:f3d.nocndbdy] = take(icndlznew,oo)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
# Convenience routines for calling the surface of revolution routines.
# The presence of these allows changing the argument lists of the fortran
# versions without breaking code.
#########################################################################
def srfrvout(rofzfunc=" ",volt=0.,zmin=None,zmax=None,xcent=0.,ycent=0.,
             rmax=top.largepos,lfill=false,
             xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
             zmmin=None,zmmax=None,zbeam=None,dx=None,dy=None,dz=None,
             nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
             xmesh=None,ymesh=None,l2symtry=None,l4symtry=None,
             srfrv_pernz=0,condid=0):
  """
Sets up a conductor represented by the outside of a surface of revolution.
The routine rofzfunc should be of the form
   def rofz():
     f3d.srfrv_r = f(f3d.srfrv_z)
where f() is the radius as a function of z and srfrv_z and srfrv_r are
compiled variables which pass data into and out of rofz.
The begining of the calling sequence would then be
srfrvout(rofz,...)

Input:
  rofzfunc=" ": routine which calculates the radius as function of z.
  volt=0.: voltage on the conductor.
  zmin=w3d.zmmin: minimum z of the conductor.
  zmax=w3d.zmmax: maximum z of the conductor.
  xcent=0.: x center of the conductor.
  ycent=0.: y center of the conductor.
  rmax=LargePos: maximum radius of the conductor.
  lfill=false: logical requesting that the whole conductor be filled
               with points.
  lshell=true: logical requesting that the shell be subgrid resolved
  xmin,xmax,ymin,ymax: min and max transverse extent of conductor.
                       default to w3d.xmmin,xmmax,ymmin,ymmax
  zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
  xmesh,ymesh,l2symtry,l4symtry:
           are all variables describing the grid. Default to variables in w3d
           and top with the same name.
  srfrv_pernz=0: when non-zero, a piece-wise linear approximation is made of
                 the surface (for optimization)
  condid=0: Id number to identify this conductor
Output is put directly into the conductor arrays of Conductor3d.
  """
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
  if srfrv_pernz > 0:
    save_srfrv_pernz = f3d.srfrv_pernz
    f3d.srfrv_pernz = srfrv_pernz

  # --- Make sure the rofzfunc is in main.
  # --- Note that this can only really work if a reference to the function
  # --- is passed in (instead of the name).
  if not f3d.lsrlinr and type(rofzfunc) == FunctionType:
    __main__.__dict__[rofzfunc.__name__] = rofzfunc

  # --- Get the name of the input function if a reference to the function
  # --- was passed in.
  if type(rofzfunc) == FunctionType: rofzfunc = rofzfunc.__name__

  # --- Get srfrv function - in older versions of the code the final 'f'
  # --- was not there
  try:
    srfrvfunc = f3d.srfrvoutf
  except AttributeError:
    srfrvfunc = f3d.srfrvout

  # --- Now call the fortran version
  srfrvfunc(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
                xmin,xmax,ymin,ymax,lshell,zmmin,zmmax,zbeam,dx,dy,dz,
                nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)

  # --- Reset srfrv_pernz if needed
  if srfrv_pernz > 0: f3d.srfrv_pernz = save_srfrv_pernz

#---------------------------------------------------------------------------
def srfrvin(rofzfunc=" ",volt=0.,zmin=None,zmax=None,xcent=0.,ycent=0.,
            rmin=0.,lfill=false,
            xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
            zmmin=None,zmmax=None,zbeam=None,dx=None,dy=None,dz=None,
            nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
            xmesh=None,ymesh=None,l2symtry=None,l4symtry=None,
            srfrv_pernz=0,condid=0):
  """
Sets up a conductor represented by the inside of a surface of revolution.
The routine rofzfunc should be of the form
   def rofz():
     f3d.srfrv_r = f(f3d.srfrv_z)
where f() is the radius as a function of z and srfrv_z and srfrv_r are
compiled variables which pass data into and out of rofz.
The begining of the calling sequence would then be
srfrvout(rofz,...)

Input:
  rofzfunc=" ": routine which calculates the radius as function of z.
  volt=0.: voltage on the conductor.
  zmin=w3d.zmmin: minimum z of the conductor.
  zmax=w3d.zmmax: maximum z of the conductor.
  xcent=0.: x center of the conductor.
  ycent=0.: y center of the conductor.
  rmin=0.: minimum radius of the conductor.
  lfill=false: logical requesting that the whole conductor be filled
               with points.
  lshell=true: logical requesting that the shell be subgrid resolved
  xmin,xmax,ymin,ymax: min and max transverse extent of conductor.
                       default to w3d.xmmin,xmmax,ymmin,ymmax
  zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
  xmesh,ymesh,l2symtry,l4symtry:
           are all variables describing the grid. Default to variables in w3d
           and top with the same name.
  srfrv_pernz=0: when non-zero, a piece-wise linear approximation is made of
                 the surface (for optimization)
  condid=0: Id number to identify this conductor
Output is put directly into the conductor arrays of Conductor3d.
  """
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
  if srfrv_pernz > 0:
    save_srfrv_pernz = f3d.srfrv_pernz
    f3d.srfrv_pernz = srfrv_pernz

  # --- Make sure the rofzfunc is in main.
  # --- Note that this can only really work if a reference to the function
  # --- is passed in (instead of the name).
  if not f3d.lsrlinr and type(rofzfunc) == FunctionType:
    __main__.__dict__[rofzfunc.__name__] = rofzfunc

  # --- Get the name of the input function if a reference to the function
  # --- was passed in.
  if type(rofzfunc) == FunctionType: rofzfunc = rofzfunc.__name__

  # --- Get srfrv function - in older versions of the code the final 'f'
  # --- was not there
  try:
    srfrvfunc = f3d.srfrvinf
  except AttributeError:
    srfrvfunc = f3d.srfrvin

  # --- Now call the fortran version
  srfrvfunc(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
               xmin,xmax,ymin,ymax,lshell,zmmin,zmmax,zbeam,dx,dy,dz,
               nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)

  # --- Reset srfrv_pernz if needed
  if srfrv_pernz > 0: f3d.srfrv_pernz = save_srfrv_pernz

#---------------------------------------------------------------------------
def srfrvinout(rminofz=" ",rmaxofz=" ",volt=0.,zmin=None,zmax=None,
               xcent=0.,ycent=0.,lzend=true,
               xmin=None,xmax=None,ymin=None,ymax=None,lshell=true,
               zmmin=None,zmmax=None,zbeam=None,dx=None,dy=None,dz=None,
               nx=None,ny=None,nz=None,ix_axis=None,iy_axis=None,
               xmesh=None,ymesh=None,l2symtry=None,l4symtry=None,
               srfrv_pernz=0,condid=0):
  """
Sets up a conductor between two surfaces of revolution.
The routines rminofz and rmaxofz should be of the form
   def rofz():
     f3d.srfrv_r = f(f3d.srfrv_z)
where f() is the radius as a function of z and srfrv_z and srfrv_r are
compiled variables which pass data into and out of rofz.

Input:
  rminofz=" ": routine which calculates the inner radius as function of z.
  rmaxofz=" ": routine which calculates the outer radius as function of z.
  volt=0.: voltage on the conductor.
  zmin=w3d.zmmin: minimum z of the conductor.
  zmax=w3d.zmmax: maximum z of the conductor.
  xcent=0.: x center of the conductor.
  ycent=0.: y center of the conductor.
  rmin=0.: minimum radius of the conductor.
  lzend=true: logical requesting that the end of the conductor be included
  lshell=true: logical requesting that the shell be subgrid resolved
  xmin,xmax,ymin,ymax: min and max transverse extent of conductor.
                       default to w3d.xmmin,xmmax,ymmin,ymmax
  zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,ix_axis,iy_axis,
  xmesh,ymesh,l2symtry,l4symtry:
           are all variables describing the grid. Default to variables in w3d
           and top with the same name.
  srfrv_pernz=0: when non-zero, a piece-wise linear approximation is made of
                 the surface (for optimization)
  condid=0: Id number to identify this conductor
Output is put directly into the conductor arrays of Conductor3d.
  """
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
  if srfrv_pernz > 0:
    save_srfrv_pernz = f3d.srfrv_pernz
    f3d.srfrv_pernz = srfrv_pernz

  # --- Make sure the rofzfunc is in main.
  # --- Note that this can only really work if a reference to the function
  # --- is passed in (instead of the name).
  if not f3d.lsrlinr and type(rminofz) == FunctionType:
    __main__.__dict__[rminofz.__name__] = rminofz
  if not f3d.lsrlinr and type(rmaxofz) == FunctionType:
    __main__.__dict__[rmaxofz.__name__] = rmaxofz

  # --- Get the name of the input function if a reference to the function
  # --- was passed in.
  if type(rminofz) == FunctionType: rminofz = rminofz.__name__
  if type(rmaxofz) == FunctionType: rmaxofz = rmaxofz.__name__

  # --- Get srfrv function - in older versions of the code the final 'f'
  # --- was not there
  try:
    srfrvfunc = f3d.srfrvinoutf
  except AttributeError:
    srfrvfunc = f3d.srfrvinout

  # --- Now call the fortran version
  srfrvfunc(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,lzend,
                  xmin,xmax,ymin,ymax,lshell,zmmin,zmmax,zbeam,dx,dy,dz,
                  nx,ny,nz,ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)

  # --- Reset srfrv_pernz if needed
  if srfrv_pernz > 0: f3d.srfrv_pernz = save_srfrv_pernz

#---------------------------------------------------------------------------
def  platepnt(ixmin=None,ixmax=None,iymin=None,iymax=None,
              ix_axis=None,iy_axis=None,dx=None,dy=None,
              aper=None,rmax=None,vvv=None,xoff=None,yoff=None,
              delz_in=None,iz=None,lz_in_plate=None,fuzz=None,condid=0):
  """
Python interface for the platepnt routine. This now just calls
the srfrvout routine. Note that the option lz_in_plate is now ignored.
  ixmin=0: minimum value of ix
  ixmax=w3d.nx: maximum value of ix
  iymin=0: minimum value of iy
  iymax=w3d.ny: maximum value of iy
  ix_axis=w3d.ix_axis: x grid location of beam center
  iy_axis=w3d.iy_axis: y grid location of beam center
  dx=w3d.dx: grid cell size in x
  dy=w3d.dy: grid cell size in y
  aper=0.: inner aperture of the plate
  rmax=LargePos: maximum radius of the plate
  vvv=0.: voltage on the plate
  xoff=0.: x offset of aperture
  yoff=0.: y offset of aperture
  delz_in=0.: fraction of cell, when outside of plate, to edge of plate
  iz=0: axial grid location of the plate
  fuzz=1.e-5*w3d.dz: number smalled compared to grid cell size, used
                     to prevent precision problems
  condid=0: Id number to identify this conductor
  """
  if ixmin is None: ixmin = 0
  if ixmax is None: ixmax = w3d.nx
  if iymin is None: iymin = 0
  if iymax is None: iymax = w3d.ny
  if ix_axis is None: ix_axis = w3d.ix_axis
  if iy_axis is None: iy_axis = w3d.iy_axis
  if dx is None: dx = w3d.dx
  if dy is None: dy = w3d.dy
  if aper is None: aper = 0.
  if rmax is None: rmax = top.largepos
  if vvv is None: vvv = 0.
  if xoff is None: xoff = 0.
  if yoff is None: yoff = 0.
  if delz_in is None: delz_in = 0.
  if iz is None: iz = 0
  if lz_in_plate is None: lz_in_plate = false
  if fuzz is None: fuzz = 1.e-5*w3d.dz
  
  f3d.lsrlinr = true
  f3d.npnts_sr = 2
  gchange("Surface_of_Rev")
  f3d.z_sr[0] = w3d.zmmin + iz*w3d.dz - 1.e-11*w3d.dz
  f3d.z_sr[1] = w3d.zmmin + iz*w3d.dz + 1.e-11*w3d.dz
  f3d.r_sr[:] = aper
  srfrvout(" ",vvv,f3d.z_sr[0],f3d.z_sr[1],xoff,yoff,rmax,true,
           w3d.xmmin+ixmin*w3d.dx,w3d.xmmin+ixmax*w3d.dx,
           w3d.ymmin+iymin*w3d.dy,w3d.ymmin+iymax*w3d.dy,true,
           w3d.zmmin,w3d.zmmax,top.zbeam,dx,dy,w3d.dz,
           w3d.nx,w3d.ny,w3d.nz,ix_axis,iy_axis,condid=condid)

  f3d.lsrlinr = false

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
def setconductorvoltage(voltage,condid=0,discrete=false,setvinject=true):
  """
Sets the voltage on a conductor, given an id.
 - voltage: voltage on conductor
 - condid=0: conductor id number
 - discrete=false: when true, z locations for plus/minus z subgrid
                   points are round up/down.
 - setvinject=false: when true, sets top.vinject
  """

  if w3d.solvergeom == w3d.RZgeom or w3d.solvergeom == w3d.XZgeom:
    if type(voltage) in [ListType,TupleType,ArrayType]:
    # --- Voltage is assumed to be the voltages are the z grid cell locations
    # --- (in the global beam frame).
      setconductorvoltagerz(voltage,w3d.nzfull,top.zmslmin[0],w3d.dz,discrete)
    # --- If setvinject is true, set it to the voltage on the left edge
      if setvinject: top.vinject = voltage[0]
    else:
      setconductorvoltagerz_id(condid,voltage)
    # --- If setvinject is true, set it to the voltage
      if setvinject: top.vinject = voltage
    return

  if f3d.ncond == 0: return

  if type(voltage) in [ListType,TupleType,ArrayType]:
    # --- Voltage is assumed to be the voltages are the z grid cell locations
    # --- (in the global beam frame).

    # --- If setvinject is true, set it to the voltage on the left edge
    if setvinject: top.vinject = voltage[0]

    # --- Get z location of the conductor points (taking into account
    # --- the differing grid cell sizes for coarse levels). Use that to
    # --- gather the value of voltage at those locations.
    icz = f3d.izcond*take(nint(f3d.mglevelslz),f3d.icondlevel) + \
                     take(f3d.mglevelsiz,f3d.icondlevel)
    icz = minimum(icz,len(voltage)-2,icz)
    cv = take(voltage,icz)

    # --- Get z location of even points. For conductors in the transverse
    # --- plane, use the voltage at the grid cell locations.
    iecl = take(nint(f3d.mglevelslz),f3d.iecndlevel)
    iecz = f3d.iecndz*iecl + take(f3d.mglevelsiz,f3d.iecndlevel)
    iecz = minimum(iecz,len(voltage)-2,iecz)
    ecv = take(voltage,iecz)
    ecvmx = ecv
    ecvpx = ecv
    ecvmy = ecv
    ecvpy = ecv

    # --- For conductors to the left, find the location of the conductor
    # --- and linear interpolate from the voltage data. If the discrete flag
    # --- is set, then round down to the nearest grid point.
    ecmz = iecz + where(f3d.ecdelmz < 1.,-f3d.ecdelmz,0)*iecl
    iecmz = ecmz.astype(Int)
    if discrete: wecmz = 0.
    else:        wecmz = ecmz - iecmz
    ecvmz = take(voltage,iecmz)*(1.-wecmz) + take(voltage,iecmz+1)*wecmz

    # --- Same for conductors to the right. If discrete is set, round up.
    ecpz = iecz + where(f3d.ecdelpz < 1.,-f3d.ecdelpz,0)*iecl
    iecpz = ecpz.astype(Int)
    if discrete: wecpz = 1.
    else:        wecpz = ecpz - iecpz
    ecvpz = take(voltage,iecpz)*(1.-wecpz) + take(voltage,iecpz+1)*wecpz

    # --- Repeat for odd conductor points.
    iocl = take(nint(f3d.mglevelslz),f3d.iocndlevel)
    iocz = f3d.iocndz*iocl + take(f3d.mglevelsiz,f3d.iocndlevel)
    iocz = minimum(iocz,len(voltage)-2,iocz)
    ocv = take(voltage,iocz)
    ocvmx = ocv
    ocvpx = ocv
    ocvmy = ocv
    ocvpy = ocv

    ocmz = iocz + where(f3d.ocdelmz < 1.,-1,0)*iocl
    iocmz = ocmz.astype(Int)
    if discrete: wocmz = 0.
    else:        wocmz = ocmz - iocmz
    ocvmz = take(voltage,iocmz)*(1.-wocmz) + take(voltage,iocmz+1)*wocmz

    ocpz = iocz + where(f3d.ocdelpz < 1.,-1,0)*iocl
    iocpz = ocpz.astype(Int)
    if discrete: wocpz = 0.
    else:        wocpz = ocpz - iocpz
    ocvpz = take(voltage,iocpz)*(1.-wocpz) + take(voltage,iocpz+1)*wocpz
  elif type(voltage) == FunctionType:
    # --- Assumes that voltage is a function which takes 3 arguments, the
    # --- coordinates x, y, z, in meters relative to the beam frame.

    # --- If setvinject is true, set it to the voltage at the source center
    if setvinject: top.vinject = voltage(top.xinject,top.yinject,top.zinject)

    zmmin = w3d.zmminglobal
    icx = f3d.ixcond*take(nint(f3d.mglevelslx),f3d.icondlevel)
    icy = f3d.iycond*take(nint(f3d.mglevelsly),f3d.icondlevel)
    icz = f3d.izcond*take(nint(f3d.mglevelslz),f3d.icondlevel) + \
                     take(f3d.mglevelsiz,f3d.icondlevel)
    cx = w3d.xmmin + icx*w3d.dx
    cy = w3d.ymmin + icy*w3d.dy
    cz =     zmmin + icz*w3d.dz
    cv = voltage(cx,cy,cz)

    ieclx = take(nint(f3d.mglevelslx),f3d.iecndlevel)
    iecly = take(nint(f3d.mglevelsly),f3d.iecndlevel)
    ieclz = take(nint(f3d.mglevelslz),f3d.iecndlevel)
    ecx = w3d.xmmin + w3d.dx*f3d.iecndx*ieclx
    ecy = w3d.ymmin + w3d.dy*f3d.iecndy*iecly
    ecz =     zmmin + w3d.dz*(f3d.iecndz*ieclz + \
                              take(f3d.mglevelsiz,f3d.iecndlevel))

    ecxmx = ecx - where(f3d.ecdelmx < 1.,f3d.ecdelmx,0)*ieclx*w3d.dx
    ecxpx = ecx + where(f3d.ecdelpx < 1.,f3d.ecdelpx,0)*ieclx*w3d.dx
    ecymy = ecy - where(f3d.ecdelmy < 1.,f3d.ecdelmy,0)*iecly*w3d.dy
    ecypy = ecy + where(f3d.ecdelpy < 1.,f3d.ecdelpy,0)*iecly*w3d.dy
    eczmz = ecz - where(f3d.ecdelmz < 1.,f3d.ecdelmz,0)*ieclz*w3d.dz
    eczpz = ecz + where(f3d.ecdelpz < 1.,f3d.ecdelpz,0)*ieclz*w3d.dz
    ecvmx = voltage(ecxmx,ecy  ,ecz  )
    ecvpx = voltage(ecxpx,ecy  ,ecz  )
    ecvmy = voltage(ecx  ,ecymy,ecz  )
    ecvpy = voltage(ecx  ,ecypy,ecz  )
    ecvmz = voltage(ecx  ,ecy  ,eczmz)
    ecvpz = voltage(ecx  ,ecy  ,eczpz)

    ioclx = take(nint(f3d.mglevelslx),f3d.iocndlevel)
    iocly = take(nint(f3d.mglevelsly),f3d.iocndlevel)
    ioclz = take(nint(f3d.mglevelslz),f3d.iocndlevel)
    ocx = w3d.xmmin + w3d.dx*f3d.iocndx*ioclx
    ocy = w3d.ymmin + w3d.dy*f3d.iocndy*iocly
    ocz =     zmmin + w3d.dz*(f3d.iocndz*ioclz + \
                              take(f3d.mglevelsiz,f3d.iocndlevel))

    ocxmx = ocx - where(f3d.ocdelmx < 1.,f3d.ocdelmx,0)*ioclx*w3d.dx
    ocxpx = ocx + where(f3d.ocdelpx < 1.,f3d.ocdelpx,0)*ioclx*w3d.dx
    ocymy = ocy - where(f3d.ocdelmy < 1.,f3d.ocdelmy,0)*iocly*w3d.dy
    ocypy = ocy + where(f3d.ocdelpy < 1.,f3d.ocdelpy,0)*iocly*w3d.dy
    oczmz = ocz - where(f3d.ocdelmz < 1.,f3d.ocdelmz,0)*ioclz*w3d.dz
    oczpz = ocz + where(f3d.ocdelpz < 1.,f3d.ocdelpz,0)*ioclz*w3d.dz
    ocvmx = voltage(ocxmx,ocy  ,ocz  )
    ocvpx = voltage(ocxpx,ocy  ,ocz  )
    ocvmy = voltage(ocx  ,ocymy,ocz  )
    ocvpy = voltage(ocx  ,ocypy,ocz  )
    ocvmz = voltage(ocx  ,ocy  ,oczmz)
    ocvpz = voltage(ocx  ,ocy  ,oczpz)

  else:
    # --- If setvinject is true, set it to the voltage
    if setvinject: top.vinject = voltage

    cv = voltage
    ecvmx = ecvpx = ecvmy = ecvpy = ecvmz = ecvpz = voltage
    ocvmx = ocvpx = ocvmy = ocvpy = ocvmz = ocvpz = voltage

  # --- Now, put the voltage data into the fortran arrays.
  if f3d.ncond > 0:
    f3d.condvolt[:] = where(equal(f3d.condnumb,condid),cv,f3d.condvolt)
  if f3d.necndbdy > 0:
    f3d.ecvolt[:]   = where(equal(f3d.ecnumb  ,condid),ecvmx,f3d.ecvolt)
    f3d.ecvoltmx[:] = where(equal(f3d.ecnumbmx,condid),ecvmx,f3d.ecvoltmx)
    f3d.ecvoltmy[:] = where(equal(f3d.ecnumbmy,condid),ecvmy,f3d.ecvoltmy)
    f3d.ecvoltmz[:] = where(equal(f3d.ecnumbmz,condid),ecvmz,f3d.ecvoltmz)
    f3d.ecvoltpx[:] = where(equal(f3d.ecnumbpx,condid),ecvpx,f3d.ecvoltpx)
    f3d.ecvoltpy[:] = where(equal(f3d.ecnumbpy,condid),ecvpy,f3d.ecvoltpy)
    f3d.ecvoltpz[:] = where(equal(f3d.ecnumbpz,condid),ecvpz,f3d.ecvoltpz)
  if f3d.nocndbdy > 0:
    f3d.ocvolt[:]   = where(equal(f3d.ocnumb  ,condid),ocvmx,f3d.ocvolt)
    f3d.ocvoltmx[:] = where(equal(f3d.ocnumbmx,condid),ocvmx,f3d.ocvoltmx)
    f3d.ocvoltmy[:] = where(equal(f3d.ocnumbmy,condid),ocvmy,f3d.ocvoltmy)
    f3d.ocvoltmz[:] = where(equal(f3d.ocnumbmz,condid),ocvmz,f3d.ocvoltmz)
    f3d.ocvoltpx[:] = where(equal(f3d.ocnumbpx,condid),ocvpx,f3d.ocvoltpx)
    f3d.ocvoltpy[:] = where(equal(f3d.ocnumbpy,condid),ocvpy,f3d.ocvoltpy)
    f3d.ocvoltpz[:] = where(equal(f3d.ocnumbpz,condid),ocvpz,f3d.ocvoltpz)

#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
def visualizeconductors(condid=None,color=None,mglevel=0,
                        scene=None,title="Conductor geometry",vrange=None,
                        smooth=0):
  """
Creates 3-D visualization of the conductors based off of the subgrid data.
condid=None: optional conductor ID of object to draw
color=None: color, in the form of a list, [r,g,b], e.g. [1,0,0] to get red
mglevel=0: multigrid level to draw
scene=None: scene to use - by default creates a new one.
            the scene is returned by the function
title="Conductor geometry": window title
vrange=None: range of each dimension - used to scale size of image, in form
             [x,y,z]. e.g. to decrease z by 10, use [1,1,10]

Returns the scene use to draw the image
  """
  if lparallel: return
  try:
    import VPythonobjects
  except ImportError:
    return

  # --- Make sure that the conductor data is properly installed.
  checkconductors(w3d.nx,w3d.ny,w3d.nz,w3d.nzfull,w3d.dx,w3d.dy,w3d.dz,
                  w3d.l2symtry,w3d.l4symtry)

  # --- Save grid size
  nx = w3d.nx
  ny = w3d.ny
  nz = w3d.nz

  # --- Get conductors
  ie = f3d.ecstart[mglevel  ] - 1
  io = f3d.ocstart[mglevel  ] - 1
  ne = f3d.ecstart[mglevel+1] - 1
  no = f3d.ocstart[mglevel+1] - 1
  if condid is None:
    # --- Get all conductors
    iecndx = f3d.iecndx[ie:ne]
    iecndy = f3d.iecndy[ie:ne]
    iecndz = f3d.iecndz[ie:ne]
    ecdelmx = f3d.ecdelmx[ie:ne]
    ecdelpx = f3d.ecdelpx[ie:ne]
    ecdelmy = f3d.ecdelmy[ie:ne]
    ecdelpy = f3d.ecdelpy[ie:ne]
    ecdelmz = f3d.ecdelmz[ie:ne]
    ecdelpz = f3d.ecdelpz[ie:ne]
    iocndx = f3d.iocndx[io:no]
    iocndy = f3d.iocndy[io:no]
    iocndz = f3d.iocndz[io:no]
    ocdelmx = f3d.ocdelmx[io:no]
    ocdelpx = f3d.ocdelpx[io:no]
    ocdelmy = f3d.ocdelmy[io:no]
    ocdelpy = f3d.ocdelpy[io:no]
    ocdelmz = f3d.ocdelmz[io:no]
    ocdelpz = f3d.ocdelpz[io:no]
  else:
    # --- Get only points matching the specified condid
    iecndx = compress(f3d.ecnumb[ie:ne]==condid,f3d.iecndx[ie:ne])
    iecndy = compress(f3d.ecnumb[ie:ne]==condid,f3d.iecndy[ie:ne])
    iecndz = compress(f3d.ecnumb[ie:ne]==condid,f3d.iecndz[ie:ne])
    ecdelmx = compress(f3d.ecnumb[ie:ne]==condid,f3d.ecdelmx[ie:ne])
    ecdelpx = compress(f3d.ecnumb[ie:ne]==condid,f3d.ecdelpx[ie:ne])
    ecdelmy = compress(f3d.ecnumb[ie:ne]==condid,f3d.ecdelmy[ie:ne])
    ecdelpy = compress(f3d.ecnumb[ie:ne]==condid,f3d.ecdelpy[ie:ne])
    ecdelmz = compress(f3d.ecnumb[ie:ne]==condid,f3d.ecdelmz[ie:ne])
    ecdelpz = compress(f3d.ecnumb[ie:ne]==condid,f3d.ecdelpz[ie:ne])
    iocndx = compress(f3d.ocnumb[io:no]==condid,f3d.iocndx[io:no])
    iocndy = compress(f3d.ocnumb[io:no]==condid,f3d.iocndy[io:no])
    iocndz = compress(f3d.ocnumb[io:no]==condid,f3d.iocndz[io:no])
    ocdelmx = compress(f3d.ocnumb[io:no]==condid,f3d.ocdelmx[io:no])
    ocdelpx = compress(f3d.ocnumb[io:no]==condid,f3d.ocdelpx[io:no])
    ocdelmy = compress(f3d.ocnumb[io:no]==condid,f3d.ocdelmy[io:no])
    ocdelpy = compress(f3d.ocnumb[io:no]==condid,f3d.ocdelpy[io:no])
    ocdelmz = compress(f3d.ocnumb[io:no]==condid,f3d.ocdelmz[io:no])
    ocdelpz = compress(f3d.ocnumb[io:no]==condid,f3d.ocdelpz[io:no])

  nn = len(iecndx) + len(iocndx)
  if nn == 0: return
  icndx = concatenate((iecndx,iocndx))
  icndy = concatenate((iecndy,iocndy))
  icndz = concatenate((iecndz,iocndz))
  delmx = concatenate((ecdelmx,ocdelmx))
  delpx = concatenate((ecdelpx,ocdelpx))
  delmy = concatenate((ecdelmy,ocdelmy))
  delpy = concatenate((ecdelpy,ocdelpy))
  delmz = concatenate((ecdelmz,ocdelmz))
  delpz = concatenate((ecdelpz,ocdelpz))
  icnd = array([icndx,icndy,icndz])
  dels = array([delmx,delpx,delmy,delpy,delmz,delpz])

  model = VPythonobjects.VisualModel(twoSided=true,scene=scene,title=title,
                                     vrange=vrange)

  gridmin = array([w3d.xmmin,w3d.ymmin,w3d.zmmin])
  griddd = array([w3d.dx*f3d.mglevelslx[mglevel],
                  w3d.dy*f3d.mglevelsly[mglevel],
                  w3d.dz*f3d.mglevelslz[mglevel]])
  gridnn = array([f3d.mglevelsnx[mglevel],
                  f3d.mglevelsny[mglevel],
                  f3d.mglevelsnz[mglevel]])
    
  # --- This fortran routine generates the triangulated surface. It was
  # --- converted to fortran for speed.
  f3d.ntriangles = 0
  getconductorfacets(nn,icnd,dels,gridnn,griddd,gridmin)
  if smooth: conductorsmoothshading()

  model.triangles = reshape(transpose(f3d.triangles),(3*f3d.ntriangles,3))
  model.normals = reshape(transpose(f3d.normals),(3*f3d.ntriangles,3))
  if color is not None:
    model.colors = (3*f3d.ntriangles)*[color]

  model.Display()
  return model

