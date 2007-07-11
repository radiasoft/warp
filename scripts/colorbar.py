from warp import *
colorbar_version = "$Id: "
#############################################################################
#############################################################################
#############################################################################
# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# ed williams' colorbar stuff / modified for Warp by J.-L. Vay on 01/22/2001

def nice_levels(z,n=8) :
  """nice_levels(z,n=8) finds approximately n "nice values"
between min(z) and max(z) for axis labels. n defaults to eight.
  """
  zmax = max(ravel(z))
  zmin = min(ravel(z))
  if zmin == zmax: return array([zmin,zmax])
  finest = abs(zmax - zmin)/float (n)
  # blows up on zmin=zmax
  unit = 10.**floor (log10 (finest))
  finest = finest/unit
  if   finest > 5.0: finest = 10.
  elif finest > 2.:  finest = 5.
  elif finest > 1.:  finest = 2.
  unit = unit*finest
  cmin = unit*ceil(zmin/unit)
  if (abs(cmin - zmin) < 0.01*unit) :
     cmin = cmin
  cmax = unit*floor(zmax/unit)
  if (abs(cmax - zmax) < 0.01*unit) :
     cmax = cmax
  n = int(((cmax - cmin)/unit + 0.5) + 1)
  levs = cmin + arange(n)*unit
  llist = oldnonzero(less(abs(levs),0.1*unit))
  if len(llist) > 0:
     array_set(levs,llist,0.0)
  return levs

#-----------------------------------------------------------------------
def color_bar(zmin,zmax,uselog=0,ncolor=100,view=1):
  """
Plots a color bar to the right of the plot square labelled by the z
values from zmin to zmax.
  - zmin, zmax: lower and upper range for color bar
  - uselog=0: when true, labels are printed in the form 10^x
  - ncolor=100: default number of colors to include
  - view=1: specifies the view that is associated with the color bar
  """
  plsys(0)
  xmin = 0.66
  xmax = 0.68
  ymax = 0.85
  ymin = 0.44
  if type(zmin) == type(zmax) == type(1):
     plotval = reshape(arange(zmin,zmax+1,typecode='b'),(zmax+1-zmin,1))
  else:
     plotval = reshape(arange(ncolor)/(ncolor-1.),(ncolor,1))
  # --- draw the bar
  pli(plotval,xmin,ymin,xmax,ymax)
  # --- Draw a black box around it
  pldj([xmin,xmin,xmin,xmax],[ymin,ymax,ymin,ymin],
       [xmax,xmax,xmin,xmax],[ymin,ymax,ymax,ymax])
  # --- Plot tick marks and labels
  levs = nice_levels(array([zmin,zmax]))
  llev = len(levs)
  scales = []
  if uselog:
    ss = " 10^%.5g"
  else:
    ss = " %.5g"
  for i in xrange(llev):
    scales.append(ss%levs[i])
  if llev==2 and (levs[0] == levs[1]):
    ys = array([ymin,ymax])
  else:
    ys = ymin + (ymax - ymin)*(levs - zmin)/(zmax - zmin)
  for i in xrange(llev):
    plt(scales[i],xmax+0.005,ys[i]-0.005)   # labels
  pldj(llev*[xmin],ys,llev*[xmax+0.005],ys) # ticks
  # --- Return to plot system 1.
  plsys(view)

