from warp import *
from slice3 import *
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
  finest = abs(zmax - zmin)/float (n)
  if zmin == zmax: raise "The min and max of z cannot be the same"
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
  levs = span(cmin,cmax,n)
  llist = nonzero(less(abs(levs),0.1*unit))
  if len(llist) > 0:
     array_set(levs,llist,0.0)
  return levs

#-----------------------------------------------------------------------
def color_bar(zmin,zmax,uselog=0,split=0,ncolor=None):
  """
Plots a color bar to the right of the plot square labelled by the z
values from zmin to zmax.
  - zmin, zmax: lower and upper range for color bar
  - uselog=0: when true, labels are printed in the form 10^x
  - split=0: when non-zero, the color is split into two sections
  - ncolor=100: default number of colors to include
  """
  if ncolor == None: ncolor = 100 + (1 - split)*100
  plsys(0)
  xmin = 0.66
  xmax = 0.68
  ymax = 0.85
  ymin = 0.44
  if type(zmin) == type(zmax) == type(1):
     plotval = reshape(arange(zmin,zmax+1,typecode='b'),(zmax+1-zmin,1))
  elif not split:
     plotval = reshape(span(0,1,ncolor),(ncolor,1))
  else:
     plotval = reshape(split_bytscl(span(0,1,ncolor),0).astype('b'),(ncolor,1))
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
  ys = ymin + (ymax - ymin)*(levs - zmin)/(zmax - zmin)
  for i in xrange(llev):
    plt(scales[i],xmax+0.005,ys[i]-0.005)   # labels
  pldj(llev*[xmin],ys,llev*[xmax+0.005],ys) # ticks
  # --- Return to plot system 1.
  plsys(1)

