from warp import *
expt_diagnostic_version = "$Id: expt_diagnostic.py,v 1.1 2000/11/29 18:38:29 dave Exp $"

#############################################################################
# This script makes plots of x and y phase space like the diagnostic plots
# produced on LBL experiments.
#
# Created by D. P. Grote, April, 1994
# Converted to python, November, 2000
#
#############################################################################

def expt_diagnosticdoc():
  print """
This script makes plots of x and y phase space like the diagnostic plots
produced on LBL experiments.
Available routines
  ppxxp_ticks Makes plot for x-x' phase space
  ppyyp_ticks Makes plot for y-y' phase space
  plot_density_ticks Makes plots for two generic coordinates

Simplest use is ppxxp_ticks(iw=1) where 1 is the desired window number.
  """


#############################################################################
#############################################################################
# Plots density plots using the tick-mark style of the experiments.
def plot_density_ticks(x,y,slope=0.,nx=20,ny=20,
                  xmin=None,xmax=None,ymin=None,ymax=None,
                  line_scale=1.,color='fg',width=1.0,
                  grid=None,contours=None,filled=0,ccolor='fg',
                  titlet=None,titleb=None,titlel=None,titler=None,titles=1):
  """
Makes density tick plots for two generic coordinates
  - x, y: particle data
  - slope=0.: slope to subtract from y coordinate (y-slope*x)
  - nx, ny: grid size, defaults to 20x20
  - xmin, xmax, ymin, ymax: extrema of density grid, defaults to partcle extrema
  - line_scale=1.: scaling factor on line length
  - color='fg': color of ticks
  - width=1.0: width of ticks
  - grid: optional grid to plot (instead of deriving grid from particle data)
  - contours=None: number of countours to plot
  - filled=0: when true, plot filled contours (assumes contours is set)
  - ccolor='fg': contour color (when not filled)
  - titlet, titleb, titlel, titler: plot titles
  - titles=1: when true, plots titles
  """

  # --- If the grid is not passed in, generate it from the inputted
  # --- particle data
  if not grid:
    # --- Create space for data
    if grid == None:
      grid = zeros((1+nx,1+ny),'d')

    # --- Make sure there are particles
    if x == None or y == None or len(x) != len(y):
      raise "both x and y must be specified and of the same length"

    # --- Get slope subtracted value of y
    yms = y - x*slope

    # --- Set grid extrema
    if xmin == None: xmin = min(x)
    if xmax == None: xmax = max(x)
    if ymin == None: ymin = min(yms)
    if ymax == None: ymax = max(yms)

    # --- Get x-y laid down on grid
    setgrid2d(len(x),x,yms,nx,ny,grid,xmin,xmax,ymin,ymax)

  # --- Set plot titles
  if titles: ptitles(titlet,titleb,titlel,titler)

  # --- Grid cell sizes and meshes
  dx = (xmax-xmin)/nx
  dy = (ymax-ymin)/ny
  xmesh = xmin + dx*arange(nx+1)[:,NewAxis]*ones(ny+1,'d')
  ymesh = ymin + dy*arange(ny+1)*ones(nx+1,'d')[:,NewAxis] + slope*xmesh

  # --- make contour plot of grid
  if contours:
    plotc(transpose(grid),transpose(ymesh),transpose(xmesh),
          color=ccolor,contours=contours,filled=filled)

  # --- Set line length
  sss = line_scale*(xmax-xmin)/nx/maxnd(grid)

  # --- Make plot of tick marks
  for ix in range(nx+1):
    for iy in range(ny+1):
      plg(ymesh[ix,iy]+zeros(2),xmesh[ix,iy]+array([0.,sss*grid[ix,iy]]),
          color=color,width=width)



############################################################################
def ppxxp_ticks(iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,
                x=None,xp=None,nx=20,nxp=20,
                xmin=None,xmax=None,xpmin=None,xpmax=None,slope='y',
                line_scale=1.,color='fg',width=1.0,
                grid=None,contours=None,filled=0,ccolor='fg',
                titlet=None,titleb=None,titlel=None,titler=None,titles=1):
  """
Makes density tick plots for x-x' phase space. Particles are selected the
same way as in getx.
Simplest use is ppxxp_ticks(iw=1) where 1 is the desired window number.
Input:
  - iw,js,win,z,iz,wz,zl,zu: all have same meaning as in getx
  - x, xp: optional particle data if it is obtained some other way
  - slope: slope to subtract from x' (x-slope*x'), if not supplied it is
           automatically calculated from the moments data
  - nx, nxp: grid size, defaults to 20x20
  - xmin, xmax, xpmin, xpmax: extrema of density grid, defaults to
                              particle extrema
  - line_scale=1.: scaling factor on line length
  - color='fg': color of ticks
  - width=1.0: width of ticks
  - grid: optional grid to plot (instead of deriving grid from particle data)
  - contours=None: number of countours to plot
  - filled=0: when true, plot filled contours (assumes contours is set)
  - ccolor='fg': contour color (when not filled)
  - titlet, titleb, titlel, titler: plot titles
  - titles=1: when true, plots titles
  """

  # --- Get particles
  if x == None and xp == None:
    x  = getx (iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
    xp = getxp(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  if x == None or xp == None or len(x) != len(xp):
    raise "both x and xp must be specified and of the same length"

  # --- Get the slope of the phase space
  if slope == 'y':
    if 0 <= iz <= w3d.nz:
      slope = (top.xxpbarz[iz]-top.xbarz[iz]*top.xpbarz[iz])/top.xrmsz[iz]**2
      xpo = top.xpbarz[iz]-slope*top.xbarz[iz]
    else:
      iiw = max(0,iw)
      slope = (top.xxpbar[iiw]-top.xbar[iiw]*top.xpbar[iiw])/top.xrms[iiw]**2
      xpo = top.xpbar[iiw]-slope*top.xbar[iiw]

  # --- Set plot titles
  if titles:
    if titlet == None: titlet = "X - X' phase space"
    if titleb == None: titleb = "X (m)"
    if titlel == None: titlel = "X' (rad)"
    if titler == None: titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
    ptitles(titlet,titleb,titlel,titler)

  # --- Call the generic routine to make the plot
  plot_density_ticks(nx=nx,ny=nxp,x=x,y=xp,slope=slope,
                xmin=xmin,xmax=xmax,ymin=xpmin,ymax=xpmax,
                line_scale=line_scale,color=color,width=width,
                grid=grid,contours=contours,filled=filled,
                ccolor=ccolor,
                titlet=titlet,titleb=titleb,titlel=titlel,titler=titler,
                titles=titles)


############################################################################
def ppyyp_ticks(ny=20,nyp=20,
               iw=0,js=1,win=None,z=None,iz=None,wz=1.,zl=None,zu=None,
               y=None,yp=None,
               ymin=None,ymax=None,ypmin=None,ypmax=None,slope='y',
               line_scale=1.,color='fg',width=1.0,
               grid=None,contours=None,filled=0,ccolor='fg',
               titlet=None,titleb=None,titlel=None,titler=None,titles=1):
  """
Makes density tick plots for y-y' phase space. Particles are selected the
same way as in gety.
Simplest use is ppyyp_ticks(iw=1) where 1 is the desired window number.
Input:
  - iw,js,win,z,iz,wz,zl,zu: all have same meaning as in gety
  - y, yp: optional particle data if it is obtained some other way
  - slope: slope to subtract from y' (y-slope*y'), if not supplied it is
           automatically calculated from the moments data
  - ny, nyp: grid size, defaults to 20x20
  - ymin, ymax, ypmin, ypmax: extrema of density grid, defaults to
                              particle extrema
  - line_scale=1.: scaling factor on line length
  - color='fg': color of ticks
  - width=1.0: width of ticks
  - grid: optional grid to plot (instead of deriving grid from particle data)
  - contours=None: number of countours to plot
  - filled=0: when true, plot filled contours (assumes contours is set)
  - ccolor='fg': contour color (when not filled)
  - titlet, titleb, titlel, titler: plot titles
  - titles=1: when true, plots titles
  """

  # --- Get particles
  if y == None and yp == None:
    y  = gety (iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
    yp = getyp(iw=iw,js=js,win=win,z=z,iz=iz,wz=wz,zl=zl,zu=zu)
  if y == None or yp == None or len(y) != len(yp):
    raise "both y and yp must be specified and of the same length"

  # --- Get the slope of the phase space
  if slope == 'y':
    if 0 <= iz <= w3d.nz:
      slope = (top.yypbarz[iz]-top.ybarz[iz]*top.ypbarz[iz])/top.yrmsz[iz]**2
      ypo = top.ypbarz[iz]-slope*top.ybarz[iz]
    else:
      iiw = max(0,iw)
      slope = (top.yypbar[iiw]-top.ybar[iiw]*top.ypbar[iiw])/top.yrms[iiw]**2
      ypo = top.ypbar[iiw]-slope*top.ybar[iiw]

  # --- Set plot titles
  if titles:
    if titlet == None: titlet = "Y - Y' phase space"
    if titleb == None: titleb = "Y (m)"
    if titlel == None: titlel = "Y' (rad)"
    if titler == None: titler = zwintitle(iw=iw,iz=iz,wz=wz,zl=zl,zu=zu)
    ptitles(titlet,titleb,titlel,titler)

  # --- Call the generic routine to make the plot
  plot_density_ticks(nx=ny,ny=nyp,x=y,y=yp,slope=slope,
                xmin=ymin,xmax=ymax,ymin=ypmin,ymax=ypmax,
                line_scale=line_scale,color=color,width=width,
                grid=grid,contours=contours,filled=filled,
                ccolor=ccolor,
                titlet=titlet,titleb=titleb,titlel=titlel,titler=titler,
                titles=titles)

############################################################################
