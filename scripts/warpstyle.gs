# Gist warpstyle.gs drawing style
# $Id: warpstyle.gs,v 1.2 2000/10/16 19:56:14 dave Exp $

# Multiple viewports
#  1:  Standard, full page, left axis
#  2:  full page, right axis
#  3:  upper left
#  4:  upper right
#  5:  lower left
#  6:  lower right
#  7:  left half
#  8:  right half
#  9:  top half
#  10: bottom half

# See work.gs for description of meanings

# Note:
# The lower x ticks and x labels are drawn by the left axis system,
# but the upper x ticks are drawn by the right axis system.  There is
# no other warning if the x limits of the left system are not the
# same as the x limits of the right system.

landscape= 0

# The default coordinate system template is identical to work.gs
default = {
  legend= 0,

  viewport= { 0.1757, 0.6143, 0.4257, 0.8643 },

  ticks= {

    horiz= {
      nMajor= 5.0,  nMinor= 25.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 3,  gridLevel= 1,  flags= 0x02b,
      tickOff= 0.0007,  labelOff= 0.010,
      tickLen= { 0.009, 0.006, 0.004, 0.0027, 0.0018 },
      tickStyle= { color= -2,  type= 1,  width= 1.0 },
      gridStyle= { color= -2,  type= 3,  width= 1.0 },
      textStyle= { color= -2,  font= 0x08,  height= 0.0182,
        orient= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.395,  yOver= 0.370 },

    vert= {
      nMajor= 5.0,  nMinor= 25.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 4,  gridLevel= 1,  flags= 0x02b,
      tickOff= 0.0007,  labelOff= 0.006,
      tickLen= { 0.009, 0.006, 0.004, 0.0027, 0.0018 },
      tickStyle= { color= -2,  type= 1,  width= 1.0 },
      gridStyle= { color= -2,  type= 3,  width= 1.0 },
      textStyle= { color= -2,  font= 0x08,  height= 0.0182,
        orient= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.150,  yOver= 0.370 },

    frame= 1,
    frameStyle= { color= -2,  type= 1,  width= 1.0 }}}

# The first coordinate system has left ticks and labels
system= {
  #legend= "Left Axis (1)",
  ticks= { horiz= { flags= 0x029 }, vert= { flags= 0x029 } }}

# The second coordinate system has right ticks and labels
system= {
  #legend= "Right Axis (2)",
  ticks= { horiz= { flags= 0x00a }, vert= { flags= 0x04a } }}

# The third coordinate system is only in first quadrant (upper left)
system= {
  #legend= "Upper Left Quadrant(3)",
  viewport= { 0.1757, 0.3511, 0.6889, 0.8643 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

# The fourth coordinate system is only in second quadrant (upper right)
system= {
  #legend= "Upper Right Quadrant(4)",
  viewport= { 0.4389, 0.6143, 0.6889, 0.8643 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

# The fifth coordinate system is only in third quadrant (lower left)
system= {
  #legend= "Lower Left Quadrant(5)",
  viewport= { 0.1757, 0.3511, 0.4257, 0.6011 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

# The sixth coordinate system is only in fourth quadrant (lower right)
system= {
  #legend= "Lower Right Quadrant(6)",
  viewport= { 0.4389, 0.6143, 0.4257, 0.6011 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

# The seventh coordinate system is only in left half
system= {
  #legend= "Left Half(7)",
  viewport= { 0.1757, 0.3511, 0.4257, 0.8643 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

# The eight coordinate system is only in right half
system= {
  #legend= "Right Half(8)",
  viewport= { 0.4389, 0.6143, 0.4257, 0.8643 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

# The nineth coordinate system is only in top half
system= {
  #legend= "Top Half(9)",
  viewport= { 0.1757, 0.6143, 0.6889, 0.8643 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

# The tenth coordinate system is only in botton half
system= {
  #legend= "Bottom Half(10)",
  viewport= { 0.1757, 0.6143, 0.4257, 0.6011 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

# The eleventh coordinate is for full view
system= {
  #legend= "Full View(11)",
  viewport= { 0.1200, 0.6802, 0.3548, 0.9150 },
  ticks= { horiz= { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }},
           vert=  { tickLen= { 0.006, 0.004, 0.002, 0.002, 0.001 }} }}

legends= {
  x= 0.04698,  y= 0.350,  dx= 0.3758,  dy= 0.0,
  textStyle= { color= -2,  font= 0x00,  height= 0.0156,
    orient= 0,  alignH= 1,  alignV= 1,  opaque= 0 },
  nchars= 36,  nlines= 20,  nwrap= 2 }

clegends= { nlines= 0 }
