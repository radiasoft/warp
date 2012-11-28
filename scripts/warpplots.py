"""
| Basic graphics commands
|   setup(): does the work needed to start writing plots to a file automatically
|   winon(): creates X graphic windows and optionally associated cgm or ps files
|   fma(): frame advance, sends plots to the plot file and sets window to be
|        cleared for the next plot
|   hcp(): send current plot to hard-copy file without clearing the display
|   redraw(): redraw the plot display

| All plot commands take the color option. These are possible values:
|   fg, bg, white, black, red, green, blue, cyan, magenta, yellow

| Many commands take the marker option. These are possible values:
|   point, plus, star, circle

| General plotting commands
|   plg(): basic plotting routine, can plot multi-dimensional arrays
|   plp(): plots markers (dots) instead of lines
|   limits(): sets plot limits in order left, right, bottom, top
|   ptitles(): draw plot titles on the current frame
|   zoombox(): when called, use the mouse left button (press and hold) to draw a
|              box around the area to be zoomed to.
|   mouse commmands: left button zoom in, middle shifts, right zoom out

| These return or set a slice out of the rho or phi array.
|   getrho(), getphi(), setrho(), setphi()

| The following plot various particles projections.
|   ppzxy(), ppzx(), ppzy(), ppzr()
|   ppzxp(), ppzvx(), ppzyp(), ppzvy(), ppzvz(), ppzrp(), ppzvr(), ppzvperp()
|   ppxy(), ppxxp(), ppyyp(), ppxpyp(), ppxvx(), ppyvy(), ppxvz(), ppyvz()
|   ppvxvy(), ppvxvz(), ppvyvz(), ppvzvperp()
|   pptrace()
|   pprrp(), pprtp(), pprvr(), pprvz()

| The following plot various particles projections using color.
|   ppzxco(), ppzyco(), ppzxyco(), ppzvzco()

| Plots arbitrary particle projections using color
|   ppco()

| Plots various quantities versus z in beam frame are defined in pzplots.py.
|   (see pzplotsdoc())

| Run histplotsdoc() for a list of history plots.

| Plots solution from envelope code.
|   penv()

| Plots contours of charge density (rho) or electrostatic potential (phi) or self E fields in
| various planes.
|   pcrhozy(), pcrhozx(), pcrhoxy()
|   pcphizy(), pcphizx(), pcphixy()
|   pcselfezy(), pcselfezx(), pcselfexy()
|   pcjzy(), pcjzx(), pcjxy()
|   pcbzy(), pcbzx(), pcbxy()
|   pcazy(), pcazx(), pcaxy()

| Change and query the properties of the plot
|   set_label(): change label properties, such as font and size
|   scale_labels(): scale the label sizes
|   setlinewidth(): set the default line with, including the axis
|                   and tick marks
|   setviewport(): change the plotting viewport
|   wplq(): get information about a plot element in the current view
|   aplq(): get information about all plot elements in the current view

| Dynamically manipulate the color palette using the colorbar
|   changepalette()

| Dynamically view any gist 3-D surface plot
|   viewsurface()

|   pltfld3d(): makes fields plots which have been turned on
|   onedplts(): makes 1-D plots which have been turned on
|   psplots(): makes particle phase space plots which have been turned on

"""
from warp import *

import __main__
__main__.__dict__['with_matplotlib'] = 0
try:
  with_matplotlib = __main__.__dict__['with_matplotlib']
except KeyError:
  with_matplotlib = ('--with-matplotlib' in sys.argv)
  __main__.__dict__['with_matplotlib'] = with_matplotlib

with_gist = not with_matplotlib

if with_matplotlib:
  import pylab
  # --- Set up some defaults to match the basic gist window.
  pylab.rcParams['figure.figsize'] = (8.5,11.)
  pylab.rcParams['figure.subplot.left'] = 0.1757
  pylab.rcParams['figure.subplot.right'] = 0.1757 + 0.4386*11./8.5
  pylab.rcParams['figure.subplot.bottom'] = 0.4257
  pylab.rcParams['figure.subplot.top'] = 0.8643
  pylab.rcParams['font.size'] = 16.0
  maxplotwindows = 64

else:
  try:
    if me == 0 and sys.platform != 'mac':
      import gist
    else:
      import gistdummy as gist
  except ImportError:
    import gistdummy as gist
  try:
    maxplotwindows = gist.GH_NDEVS
  except AttributeError:
    # --- For older versions of pygist
    maxplotwindows = 8

import controllers
import re
import os
import sys
warpplots_version = "$Id: warpplots.py,v 1.278 2012/01/26 23:31:29 grote Exp $"

def warpplotsdoc():
  import warpplots
  print warpplots.__doc__

##########################################################################
top.lpsplots = true
always = top.always
seldom = top.seldom
never = top.never
cgmlogfile = None
numframeslist = ones(maxplotwindows,'l')
_hcp_frame_number = zeros(maxplotwindows,'l')

if with_gist:
  gist.pldefault(marks=0) # --- Set plot defaults, no line marks
  gist.pldefault(legends=0) # --- Turn off the legends in hardcopy output
  #gist.pldefault(width=4.) # --- Set plot defaults, width linewidth

  # --- Set GISTPATH environment variable appropriately if it is not already
  # --- set.
  if "GISTPATH" not in os.environ:
    import warp
    os.environ["GISTPATH"] = os.path.dirname(warp.__file__)

  def active_window(winnum=None):
    if winnum is None:
      return gist.current_window()
    else:
      gist.window(winnum)

if with_matplotlib:
  _matplotwindows = []
  _matplotactivewindow = [0]
  def active_window(winnum=None):
    if winnum is None:
      return _matplotactivewindow[0]
    else:
      _matplotactivewindow[0] = winnum

  def universeaxes():
    # --- Create a new axes which covers the whole plot frame.
    aa = pylab.axes([0.,0.,1.,1.],frameon=False)
    aa.axis([0.,1.,0.,1.])
    aa.set_axis_off()
    return aa

  def palette(cmap):
    '''
Set the default colormap apply to current image if any.
See help(colormaps) for more information and cm.cmapnames() for
a list of colormap names.
    '''
    rc('image', cmap=cmap)
    im = gci()

    if im is not None:
        im.set_cmap(getattr(cm,cmap))
    draw_if_interactive()

# The setup routine does the work needed to start writing plots to a file
# automatically.
def setup(makepsfile=0,prefix=None,cgmlog=1,runcomments='',
          cgmfilesize=100000,pnumb=None,writetodatafile=0,
          lversiontext=1,style='work.gs'):
  """
Does the work needed to start writing plots to a file automatically
  - makepsfile=0: allows the specification of a ps file instead of cgm
  - prefix=None: optional prefix to use for plotfile name instead of runid
  - cgmlog=1: Set to 0 to inhibit cgmlog file creation
  - runcomments='': Additional comments to print on the first plot page
  - cgmfilesize=100000: Max cgmfilesize in units of MBytes.
  - pnumb=None: Optional file name number to be used in place of the
                next available number or the command line argument --pnumb.
                It must be a string.
  - writetodatafile=0: When true, all plot data is written to a data
                       file instead of a gist cgm file.
  - lversiontext=1: When true, write out the version information to the
                    first frame.
  - style='work.gs': Gist style sheet. The style sheet determines
                     the number and location of coordinate systems, tick
                     and label styles, and the like.  Other choices include
                     "axes.gs", "boxed.gs", "work2.gs", and "boxed2.gs", or
                     you can create your own version.
  """
  # --- cgmlogfile is needed elsewhere
  global cgmlogfile

  # --- Setup the plot file name
  if me == 0:
    if prefix is None: prefix = arraytostr(top.runid)
    if writetodatafile:
      suffix = 'pkl'
    else:
      if makepsfile or with_matplotlib:
        suffix = 'ps'
      else:
        suffix = 'cgm'
    if pnumb is None:
      # --- Check command line argument
      if warpoptions.options.pnumb is not None:
        pnumb = warpoptions.options.pnumb
        pname = '%s.%s.%s'%(prefix,pnumb,suffix)
      else:
        # --- Get next available plot file name.
        pname = getnextfilename(prefix,suffix)
        pnumb = pname[-len(suffix)-4:-len(suffix)-1]
    else:
      pname = "%s.%s.%s"%(prefix,pnumb,suffix)
  else:
    # --- This just defines the variables, the values are never used.
    pname = ''
    pnumb = ''

  # --- Save the plotfile name and number, since its not retreivable
  # --- from gist. They are broadcast to the other processors if needed.
  setup.pname = parallel.broadcast(pname)
  setup.pnumb = parallel.broadcast(pnumb)

  # --- Only PE0 (or serial processor) should run the rest of this routine.
  if me > 0: return

  # --- Set cgmfilesize
  if with_gist:
    gist.pldefault(cgmfilesize=cgmfilesize)

  if with_gist:
    # --- Create window(0), but have it only dump to the file pname for now.
    # --- Note that only plots made to window(0) are dumped to the file.
    gist.window(0,display='',hcp=pname,dump=1,style=style)
    # --- Set so all fma's dump plot to file.
    gist.hcpon()
  else:
    # --- Open the file where the plots will be saved.
    _matplotwindows.append(open(setup.pname,'w'))

  if writetodatafile:
    # --- Setup the data file name where plot data is written
    setplotdatafilename(pname)

  print "Plot file name",pname

  if cgmlog:
    # --- Create plot log file and write heading to it.
    #plogname = getnextfilename(prefix,'cgmlog')
    plogname = prefix + '.' + pnumb + '.' + suffix + 'log'
    cgmlogfile = open(plogname,"w")
    cgmlogfile.write("CGMLOG file for "+pname+"\n\n")

  # --- Print the versions to the plot file.
  if with_matplotlib: universeaxes()

  if lversiontext:
    plt(time.ctime(top.starttime)+'\n'+versionstext()+'\n'+runcomments,
        0.15,0.88,justify="LT",local=1)
    fma()

# --- Convenience function to open a window with default value specilized to
# --- WARP. By default, this opens up a window on the current display. If
# --- setup has been called, this just creates a window which is attached to
# --- the already created device. Otherwise, open a window attached to a
# --- new device.
def winon(winnum=0,dpi=100,prefix=None,suffix=None,xon=1,style='work.gs'):
  """
Opens up an X window
  - winnum=0 is the window number
  - dpi=100 is the dots per inch (either 100 or 75)
  - prefix=None: if given, opens a new file with the same suffix number as
                 the one for window 0. Winnum cannot be 0 and setup must have
                 already been called. Warning - this will overwrite a file
                 with the same name.
                 Both prefix and suffix can be specified.
  - suffix=None: if given, opens a new file with the same suffix number as
                 the one for window 0. Winnum cannot be 0 and setup must have
                 already been called. Warning - this will overwrite a file
                 with the same name.
  - xon=1: When true, an X window will be opened.
  - style='work.gs': Gist style sheet. The style sheet determines
                     the number and location of coordinate systems, tick
                     and label styles, and the like.  Other choices include
                     "axes.gs", "boxed.gs", "work2.gs", and "boxed2.gs", or
                     you can create your own version.
  """
  if suffix is None and prefix is None:
    if with_gist:
      if xon and winnum==0 and sys.platform not in ['win32','cygwin']:
        # --- If display isn't set, no X plot window will appear since window0
        # --- is already attached to a device (the plot file).
        # --- The try/except construct takes care of the case where
        # --- the gist package was not compiled with X11.
        try:
          gist.window(winnum,dpi=dpi,display=os.environ['DISPLAY'],style=style)
        except:
          gist.window(winnum,dpi=dpi,style=style)
      else:
        if xon: gist.window(winnum,dpi=dpi,style=style,display=os.environ['DISPLAY'])
        else:   gist.window(winnum,dpi=dpi,style=style,display='')
  else:
    # --- Get the next winnum if it wasn't passed in.
    if winnum == 0:
      try:    winnum = winon.winnum + 1
      except: winnum = 1
    winon.winnum = winnum
    # --- Check input errors
    try: setup.pname
    except AttributeError: raise RuntimeError,'setup has not yet been called'
    assert winnum > 0,'winnum must not be 0'
    # --- Check file name and type from window 0
    pname = setup.pname.split('.')[0]
    filetype = setup.pname.split('.')[-1]
    # --- Create file name
    if prefix is not None: pname = prefix + pname
    if suffix is not None: pname = pname + '_' + suffix
    pname = '%s.%s.%s'%(pname,setup.pnumb,filetype)
    if with_gist:
      # --- Open window
      if xon:
        gist.window(winnum,dpi=dpi,display=os.environ['DISPLAY'],
                    dump=1,hcp=pname,style=style)
      else:
        gist.window(winnum,dpi=dpi,display='',dump=1,hcp=pname,style=style)
    if with_matplotlib:
      # --- Open a new file for the plots and make it active.
      _matplotwindows.append(open(pname,'w'))
      _matplotactivewindow[0] = len(_matplotwindows) - 1

    if _plotdatafilenames[0] is not None:
      # --- If plots for window 0 are bein written to a data file,
      # --- then do the same for this new window.
      setplotdatafilename(pname)

    return winnum

##########################################################################
# Plot run info to the current plot and plot info to the log file.
framet=''
frameb=''
framel=''
framer=''
def plotruninfo():
  "Plot run info to the current plot and plot info to the log file"
  if with_matplotlib:
    # --- Get the current axis.
    ca = pylab.gca()
    # --- Create a new one which covers the whole plot frame.
    aa = pylab.axes([0.,0.,1.,1.],frameon=False)
    aa.axis([0.,1.,0.,1.])
    aa.set_axis_off()
  ss = (arraytostr(top.pline3)+'\n'+
        arraytostr(top.pline2)+'\n'+
        arraytostr(top.pline1))
  if with_gist:
    plt(ss,0.12,0.28,local=1)
  if with_matplotlib:
    aa.text(0.12,0.28,ss)
  runmaker = arraytostr(top.runmaker)
  runtime = arraytostr(top.runtime)
  runid = arraytostr(top.runid)
  try:
    pnumb = '.' + setup.pnumb
  except AttributeError:
    # --- In case setup had not been called yet.
    pnumb = ''
  ss = '%s, %s %s%s'%(runmaker,runtime,runid,pnumb)
  if with_gist:
    # --- Replace _ with !_ so that it prints an underscore instead of making
    # --- a subscript. This is most important for the runid.
    plt(ss.replace('_','!_'),0.12,0.24,local=1)
  if with_matplotlib:
    aa.text(0.12,0.24,ss)
  # --- Increment and print frame number and log
  # --- numframeslist is now incremented in fma
  #numframeslist[active_window()] = numframeslist[active_window()] + 1
  if with_gist:
    if _hcp_frame_number[active_window()] > 0:
      # --- If an hcp has been done, overwrite the frame number using the
      # --- background color. It would better better to delete it but I
      # --- don't think there is a way currently.
      plt(repr(_hcp_frame_number[active_window()]),0.68,0.9,
          justify='RA',local=1,color='bg')
      _hcp_frame_number[active_window()] = 0
    plt(repr(numframeslist[active_window()]),0.68,0.9,justify='RA',local=1)
  if with_matplotlib:
    aa.text(0.68,0.9,repr(numframeslist[active_window()]),
            horizontalalignment='right',
            verticalalignment='top')
  if cgmlogfile:
    cgmlogfile.write('%d %d Step %d %s %s %s %s\n' %
                     (active_window(),numframeslist[active_window()],
                     top.it,framet,frameb,framel,framer))
  if with_matplotlib:
    # --- Restore the previous axis
    pylab.axes(ca)

##########################################################################
if with_gist: _plotpackage = gist
if with_matplotlib: _plotpackage = pylab
def setplotpackage(plotpackage):
  global _plotpackage
  _plotpackage = plotpackage
def getplotpackage():
  return _plotpackage

# --- This is a global switch which toggles between directly calling the gist
# --- plotting routines and accumulating lists of things to plot, which are then
# --- plotted when a final command is given. The functions provide the API for
# --- toggling the switch.
# --- Note that _accumulateplotlists is incremented each time that
# --- accumulateplotlists is called and decrement in makeplotsdirectly.
# --- This allows the switch to be set and unset in multiple places and
# --- ensures that the switch will only be turned off when every call to
# --- accumulateplotlists is matched with a call to makeplotsdirectly.
_accumulateplotlists = 0
def accumulateplotlists():
  global _accumulateplotlists
  _accumulateplotlists += 1
def makeplotsdirectly():
  global _accumulateplotlists
  assert (_accumulateplotlists>0),"makeplotsdirectly should only ever be called after a call to accumulateplotlists"
  _accumulateplotlists -= 1

# --- This allows writing out all plot data into a data file instead of
# --- creating the plot with gist. Each window has its own filename and object.
_plotdatafilenames = maxplotwindows*[None]
_plotdatafileobjects = maxplotwindows*[None]
def setplotdatafilename(filename=None,winnum=None):
  if winnum is None: winnum = active_window()
  _plotdatafilenames[winnum] = filename
  if _plotdatafileobjects[winnum] is not None:
    # --- Some other file was already opened. Make sure it is closed.
    try:
      _plotdatafileobjects[winnum].close()
    except:
      # --- This is probably not a good idea to catch all errors, but this
      # --- if block should never be executed anyway and is only here for
      # --- odd end cases.
      pass
  # --- The file is only opened when something is being written to it.
  _plotdatafileobjects[winnum] = None

# --- This is the global list of the things to be plotted and the function
# --- which actually does the plotting.
_listofthingstoplot = []
def addthingtoplot(pfunc,args,kw):
  """pfunc: name of plotting routine
     args: list of args
     kw: dict of keyword args"""
  _listofthingstoplot.append([pfunc,args,kw])
def plotlistofthings(lturnofflist=0):
  global _listofthingstoplot
  if not _accumulateplotlists: return
  listsofthings = gather(_listofthingstoplot)
  _listofthingstoplot = []
  for things in listsofthings:
    for thing in things:
      handleplotfunctioncall(thing[0],thing[1],thing[2])
  if lturnofflist: makeplotsdirectly()

def callplotfunction(pfunc,args=[],kw={}):
  # --- Note that any None's need to be cleared out since some functions
  # --- don't like them.
  while len(args) > 0 and args[-1] is None: del args[-1]
  if _accumulateplotlists:
    addthingtoplot(pfunc,args,kw)
  else:
    return handleplotfunctioncall(pfunc,args,kw)
def handleplotfunctioncall(pfunc,args=[],kw={}):
  winnum = active_window()
  if _plotdatafilenames[winnum] is None:
    # --- Make direct call to gist or matplotlib routine
    if kw is None:
      return getattr(_plotpackage,pfunc)(*args)
    else:
      return getattr(_plotpackage,pfunc)(*args,**kw)
  else:
    # --- Write data to the file (for later processing)
    # --- Note that for now pickle is used for its simplicity.
    # --- PyPDB is not really usable. It is not robust enough, since the file
    # --- must be closed in order for it to be readable, so a run crash
    # --- would cause all plot data to be lost. But, writing pickled objects
    # --- to the file is not supported in append mode (everything
    # --- written would need to be read in and rewritten out).
    # --- However, pickle is not perfect, since there seems to be an issue
    # --- of moving pickled integers between 32 and 64 bit machines.
    if _plotdatafileobjects[winnum] is None:
      # --- The file is opened here when something is being written to it.
      # --- It is opened in append mode to avoid accidentally deleting data.
      _plotdatafileobjects[winnum] = open(_plotdatafilenames[winnum],'ab')
    try:
      cPickle.dump((numframeslist[active_window()],
                    pfunc,args,kw),_plotdatafileobjects[winnum],protocol=-1)
    except TypeError:
      cPickle.dump((numframeslist[active_window()],
                    pfunc,args,kw),_plotdatafileobjects[winnum])
    _plotdatafileobjects[winnum].flush()

class PlotsFromDataFile(object):
  """This class reads in a plot data file and can be used to make plots
from that data. The available methods are
 - plotframe(framenum=None): plot the frame. framenum defaults to the next frame
 - plotallframes(): plot all frames in the file.
  """
  def __init__(self,filename):
    self.filename = filename
    self.ff = open(filename,'rb')
    self.frames = self._catalogframes()
    self.scaleheight = None
    # --- Current default frame number that plotframe will plot
    self.currentframe = 1
  def _catalogframes(self):
    "Build a catalog of plot frames"
    # --- Build a dictionary with the frame numbers as the keys.
    # --- The value for each key is a list of locations in the file where
    # --- the plot commands that are part of that frame are stored.
    frames = {}
    while 1:
      tell = self.ff.tell()
      try:
        framenum,pfunc,args,kw = cPickle.load(self.ff)
      except EOFError:
        break
      frames.setdefault(framenum,[]).append(tell)
    return frames
  def _getplotdatafromfile(self,framenum,plotcount):
    # --- This reads in the data from the file for the command specified
    # --- by the input arguments.
    # --- For now, no error checking is done. This is expert use only!
    self.ff.seek(self.frames[framenum][plotcount])
    framenum,pfunc,args,kw = cPickle.load(self.ff)
    return pfunc,args,kw
  def plotframe(self,framenum=None):
    """plotframe(framenum=None):
plot the frame.
 - framenum=None: frame number, defaults to the next frame
    """
    if framenum is None: framenum = self.currentframe
    try:
      locs = self.frames[framenum]
    except KeyError:
      # --- The framenum was not found
      raise RuntimeError,"Frame number %d was not found"%framenum
    # --- Set the location in the file to the start of the frame.
    # --- Only this is needed since the plot commands will be contiguous
    # --- in the file.
    self.ff.seek(locs[0])
    for loc in locs:
      framenum,pfunc,args,kw = cPickle.load(self.ff)
      # --- If the fma command is found, that signals the end of this frame.
      # --- Note that this routine does not call fma, allowing further
      # --- changes to the frame.
      # --- Another option would be to skip the last locs element which would
      # --- normally be the 'fma' command. But this wouldn't work in the case
      # --- when the fma command was not given for the last frame at the end
      # --- of a run.
      if pfunc == 'fma': break
      # --- Skip the hcp commands. This is not handled quite right.
      if pfunc == 'hcp': continue
      # --- Scale the height of text if the option is set.
      if pfunc == 'plt' and self.scaleheight is not None:
        try:
          kw['height'] *= self.scaleheight
        except KeyError:
          kw['height'] = 14.*self.scaleheight
      # --- Make the plot, calling the appropriate gist or matplotlib function
      handleplotfunctioncall(pfunc,args,kw)
  def plotallframes(self):
    """plotallframes():
plot all frames stored in the file.
    """
    numframes = max(self.frames.keys())
    # --- All that is needed is to call plotframe the appropriate number
    # --- of times and adding the fma calls.
    for i in range(numframes):
      self.plotframe(i+1)
      # --- The legend is turned off since it was plotted in the
      # --- original run.
      fma(legend=0)
  def gist(self):
    """This offers a rudimentary gist like interface, allowing easy scanning
through the plot file."""
    import fcntl
    import tty
    import termios
    import curses.ascii
    numframes = max(self.frames.keys())
    # --- Get the settings for stdin, so it can be restored later.
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    # --- This gives a carriage return without a linefeed.
    CR = curses.ascii.ctrl('m')
    try:
      # --- Set stdin to be nonblocking.
      rv = fcntl.fcntl(sys.stdin, fcntl.F_SETFL, os.O_NDELAY)
      # --- Set so stdin does not echo the input.
      tty.setraw(sys.stdin.fileno())
      stepstring = ''
      while 1:
        # --- Handle any gist events, such as mouse zooming.
        refresh()
        # --- Check if a character has been entered.
        try:
          c = sys.stdin.read(1)
        except IOError:
          continue
        if c == 'f':
          # --- Advance forward
          self.currentframe += int(stepstring or 1)
          if self.currentframe >= numframes:
            self.currentframe = numframes
            print "At last frame",CR,
          gist.fma()
          oldlimits = limits()
          self.plotframe()
          stepstring = ''
        elif c == 'b':
          # --- Advance backward
          self.currentframe -= int(stepstring or 1)
          if self.currentframe <= 1:
            self.currentframe = 1
            print "At first frame",CR,
          gist.fma()
          oldlimits = limits()
          self.plotframe()
          stepstring = ''
        elif c in ['0','1','2','3','4','5','6','7','8','9']:
          stepstring += c
        elif c == 'r':
          # --- reset
          unzoom()
        elif c == 's':
          # --- Send to cgm file
          hcp()
        elif c == 'q':
          # --- quit
          break
    finally:
      termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)


# Frame advance and redraw routines. The fma routine from gist is replaced
# with one that prints informative text at the bottom of each frame just
# before the normal gist fma is called. Also created are alternate (Basis
# like) names for fma and redraw.
def fma(legend=1):
  """
Frame advance - plots run info on the bottom of the frame, gets graphics window
ready for next plot and sends image to hard copy file if one is opened. Checks
for before and after plot commands.
  - legend=1: when set to 0, the text at the frame bottom is omitted
  """
  plotlistofthings()
  # --- Write the run log info
  if legend: plotruninfo()
  controllers.callafterplotfuncs()
  if with_gist:
    callplotfunction("fma")
    oldlimits = limits()
  if with_matplotlib:
    try:
      pylab.savefig(_matplotwindows[_matplotactivewindow[0]],format='ps')
    except IndexError:
      pass
    pylab.clf()
  # --- Increment frame number
  numframeslist[active_window()] = numframeslist[active_window()] + 1
  controllers.callbeforeplotfuncs()

def hcp(legend=1):
  """
Hardcopy - plots run info on the bottom of the frame and sends image to hard
copy file.
  - legend=1: when set to 0, the text at the frame bottom is omitted
  """
  controllers.callafterplotfuncs()
  # --- Write the run log info
  if legend: plotruninfo()
  if with_gist:
    callplotfunction("hcp")
  if with_matplotlib:
    pylab.savefig(_matplotwindows[_matplotactivewindow[0]],format='ps')
  # --- Save the current frame number so that it can be removed from the
  # --- next frame which will have the incremented frame number.
  _hcp_frame_number[active_window()] = numframeslist[active_window()]
  # --- Increment frame number
  numframeslist[active_window()] = numframeslist[active_window()] + 1
  controllers.callbeforeplotfuncs()

def refresh():
  """
Refresh the current gist windows.
  """
  if with_gist:
    try:
      pyg_pending()
      pyg_idler()
    except:
      ygdispatch()

# --- obsoleted
#nf = fma
#sf = redraw

##########################################################################
# This routine allows plotting of multi-dimensioned arrays.
# It replaces the plg from gist, which can only plot 1-d arrays.
def pla(y,x=None,linetype="solid",local=1,**kw):
  """This comment is replaced with gist.plg.__doc__. The linetype argument is
  only needed for backward compatibility."""
  kw.setdefault('type',linetype)
  if len(shape(y)) == 0: y = [y]
  if x is not None and len(shape(x)) == 0: x = [x]
  y = array(y,copy=False)
  if x is not None:
    x = array(x,copy=False)
    # --- This is the only constraint on the input arrays.
    assert shape(x)[0]==shape(y)[0],\
      'The first dimensions of the two input arrays must be of the same length'
  else:
    # --- If x is not supplied, it is just the integers starting at 0.
    x = arange(y.shape[0],dtype='d')
  if len(shape(x)) > 2:
    # --- Reshape the array, putting all but the 1st dimension into the
    # --- 2nd dimension.
    xx = reshape(x,(x.shape[0],product(x.shape[1:])))
  elif len(shape(x)) == 2:
    # --- The input x is usable as is.
    xx = x
  else:
    # --- Extend xx into a 2-D array, with a second dimension of length 1.
    xx = x[:,newaxis]
  if len(shape(y)) > 2:
    # --- Reshape the array, putting all but the 1st dimension into the
    # --- 2nd dimension.
    yy = reshape(y,(y.shape[0],product(y.shape[1:])))
  elif len(shape(y)) == 2:
    # --- The input y is usable as is.
    yy = y
  else:
    # --- Extend yy into a 2-D array, with a second dimension of length 1.
    yy = y[:,newaxis]
  if not local and lparallel:
    # --- This way is preferred over a gatherarray since, for large data sets,
    # --- it reduces the risk of running out of memory since only part of the
    # --- data is stored on PE0 at a time.
    if me == 0:
      for i in range(0,npes):
        if i > 0:
          yy = mpirecv(i,3)
          xx = mpirecv(i,3)
        if len(xx) > 0 and len(yy)==len(xx):
          plg(yy,xx,local=1,**kw)
    else:
      comm_world.send(yy,0,3)
      comm_world.send(xx,0,3)
  else:
    # --- convert some arguments for pylab
    if with_matplotlib:
      if kw['type'] == 'solid': kw['linestyle'] = '-'
      if kw['type'] == 'none': kw['linestyle'] = 'None'
      del kw['type']
      try:
        if kw['marker'] == "\1": kw['marker'] = ','
      except KeyError:
        pass
      try:
        kw['linewidth'] = kw['msize']
        del kw['msize']
      except KeyError:
        pass
      try:
        kw['linewidth'] = kw['width']
        del kw['width']
      except KeyError:
        pass
      try:
        if kw['color'] == 'fg': kw['color'] = 'k'
      except KeyError:
        pass
    # --- The i%n is used in case the 2nd dimensions are not equal. This
    # --- is most useful if the 2nd dimension of xx is 1, in which case
    # --- all of the plots use that as the abscissa.
    n = shape(xx)[1]
    for i in xrange(yy.shape[1]):
      if len(yy[:,i]) > 0:
        if with_gist:
          callplotfunction("plg",[yy[:,i],xx[:,i%n]],kw)
        if with_matplotlib:
          callplotfunction("plot",[xx[:,i],yy[:,i%n]],kw)

pla.__doc__ = gist.plg.__doc__
plg = pla

# --- This replaces functions from gist, filtering through callplotfunction
def limits(xmin=None,xmax=None,ymin=None,ymax=None,**kw):
  if with_gist:
    if isinstance(xmin,TupleType):
      # --- In this form, gist complains if a keyword dictionary is passed in,
      # --- so pass in None for kw.
      rr = callplotfunction("limits",[xmin],None)
    else:
      rr = callplotfunction("limits",[xmin,xmax,ymin,ymax],kw)
    return rr
  if with_matplotlib:
    callplotfunction("axis",[(xmin,xmax,ymin,ymax)],kw)
limits.__doc__ = gist.limits.__doc__
def pldj(x0,y0,x1,y1,local=1,**kw):
  if not _accumulateplotlists and not local:
    x0 = gatherarray(x0)
    y0 = gatherarray(y0)
    x1 = gatherarray(x1)
    y1 = gatherarray(y1)
  if size(x0) == 0 or size(y0) == 0 or size(x1) == 0 or size(y1) == 0: return
  if with_gist:
    callplotfunction("pldj",[x0,y0,x1,y1],kw)
  if with_matplotlib:
    callplotfunction("plot",[array([x0,x1]),array([y0,y1])],kw)
pldj.__doc__ = gist.pldj.__doc__
def plfp(z,y,x,n,local=1,**kw):
  if not _accumulateplotlists and not local:
    z = gatherarray(z)
    y = gatherarray(y)
    x = gatherarray(x)
    n = gatherarray(n)
  if size(z) == 0 or size(y) == 0 or size(x) == 0 or size(n) == 0: return
  if with_gist:
    callplotfunction("plfp",[z,y,x,n],kw)
  if with_matplotlib:
    i = 0
    for j,c in zip(n,z):
      kw.setdefault('color',(0.,1.-c/200.,c/200.))
      callplotfunction("plot",[x[i:i+j],y[i:i+j]],kw)
      i += j
plfp.__doc__ = gist.plfp.__doc__
def plfc(z,y,x,ireg,local=1,**kw):
  if not _accumulateplotlists and not local:
    z = gatherarray(z)
    y = gatherarray(y)
    x = gatherarray(x)
    ireg = gatherarray(ireg)
  if size(z) == 0 or size(y) == 0 or size(x) == 0 or size(ireg) == 0: return
  if with_gist:
    callplotfunction("plfc",[z,y,x,ireg],kw)
  if with_matplotlib:
    # --- Note: ireg could be handled by making z a masked array
    callplotfunction("contourf",[x,y,z],kw)
plfc.__doc__ = gist.plfc.__doc__
def plc(z,y=None,x=None,ireg=None,local=1,**kw):
  if not _accumulateplotlists and not local:
    z = gatherarray(z)
    if y is not None: y = gatherarray(y)
    if x is not None: x = gatherarray(x)
    if ireg is not None: ireg = gatherarray(ireg)
  if size(z) == 0: return
  if with_gist:
    callplotfunction("plc",[z,y,x,ireg],kw)
  if with_matplotlib:
    # --- Note: ireg could be handled by making z a masked array
    callplotfunction("contour",[x,y,z],kw)
plc.__doc__ = gist.plc.__doc__
def pli(z,x0=None,y0=None,x1=None,y1=None,local=1,**kw):
  if not _accumulateplotlists and not local:
    z = gatherarray(z)
    if x0 is not None: x0 = gatherarray(x0)
    if y0 is not None: y0 = gatherarray(y0)
    if x1 is not None: x1 = gatherarray(x1)
    if y1 is not None: y1 = gatherarray(y1)
  if size(z) == 0: return
  if with_gist:
    callplotfunction("pli",[z,x0,y0,x1,y1],kw)
  if with_matplotlib:
    try:
      kw['vmin'] = kw['cmin']
      del kw['cmin']
    except KeyError:
      pass
    try:
      kw['vmax'] = kw['cmax']
      del kw['cmax']
    except KeyError:
      pass
    try:
      del kw['top']
    except KeyError:
      pass
    ny,nx = z.shape
    x0 = x0 or 0.
    x1 = x1 or (nx - 1)
    y0 = y0 or 0.
    y1 = y1 or (ny - 1)
    dx = (x1 - x0)/(nx-1)
    dy = (y1 - y0)/(ny-1)
    xx = arange(x0-dx/2.,x1+dx/2.,dx)
    yy = arange(y0-dy/2.,y1+dy/2.,dy)
    xg,yg = meshgrid(xx,yy)
    callplotfunction("pcolor",[xg,yg,z],kw)
pli.__doc__ = gist.pli.__doc__
def plf(z,y=None,x=None,ireg=None,local=1,**kw):
  if not _accumulateplotlists and not local:
    z = gatherarray(z)
    if y is not None: y = gatherarray(y)
    if x is not None: x = gatherarray(x)
    if ireg is not None: ireg = gatherarray(ireg)
  if size(z) == 0: return
  if with_gist:
    callplotfunction("plf",[z,y,x,ireg],kw)
  if with_matplotlib:
    # --- ireg not implemented now
    callplotfunction("pcolor",[x,y,z],kw)
plf.__doc__ = gist.plf.__doc__
def plv(vy,vx,y=None,x=None,ireg=None,local=1,**kw):
  if not _accumulateplotlists and not local:
    vy = gatherarray(vy)
    vx = gatherarray(vx)
    if y is not None: y = gatherarray(y)
    if x is not None: x = gatherarray(x)
    if ireg is not None: ireg = gatherarray(ireg)
  if size(vy) == 0 or size(vx) == 0: return
  if with_gist:
    callplotfunction("plv",[vy,vx,y,x,ireg],kw)
  if with_matplotlib:
    callplotfunction("quiver",[x,y,vx,vy],kw)
plv.__doc__ = gist.plv.__doc__
def plt(text,x,y,local=1,**kw):
  if not _accumulateplotlists and not local:
    textlist = gather(text)
    xlist = gather(x)
    ylist = gather(y)
  else:
    textlist = [text]
    xlist = [x]
    ylist = [y]
  if with_matplotlib:
    if 'justify' in kw:
      if kw['justify'] == 'LT':
        kw['horizontalalignment'] = 'left'
        kw['verticalalignment'] = 'top'
        del kw['justify']
  for text,x,y in zip(textlist,xlist,ylist):
    if with_gist:
      callplotfunction("plt",[text,x,y],kw)
    if with_matplotlib:
      callplotfunction("text",[x,y,text],kw)
plt.__doc__ = gist.plt.__doc__
def plsys(n=None,**kw):
  if with_gist:
    if n is None: return getattr(_plotpackage,"plsys")()
    callplotfunction("plsys",[n])
plsys.__doc__ = gist.plsys.__doc__

##########################################################################
# --- Plot particles
circle = '\4'
star = '\3'
plus = '\2'
point = '\1'
def plp(y,x=None,linetype='none',marker="\1",msize=1.0,**kw):
  """Plots particles, same as plg but with different defaults so it plots
markers instead of lines"""
  if len(shape(y)) == 0: y = [y]
  if x is not None and len(shape(x)) == 0: x = [x]
  #if len(y) == 0: return
  kw.setdefault('type',linetype)
  kw['marker'] = marker
  kw['msize'] = msize
  if x is not None:
    plg(y,x,**kw)
  else:
    plg(y,**kw)

# --- Plot history data. Convenience function that is only needed until
# --- the 'limited' capability is implemented.
def plothist(v,iw):
  """
Plots any history versus z
   - v is the history array
   - iw is the window number to plot
  """
  plg(v[iw,0:top.jhist+1],top.hzbeam[0:top.jhist+1])

# --- Simple interface to contour plotting. Only requires the 2-D array
# --- to be plotted.
def plotc(zz,xx=None,yy=None,ireg=None,color='fg',levs=None,contours=8,
          filled=0,width=1.,linetype='solid',cmin=None,cmax=None,
          leveloverlap=None,local=1):
  """
Simple interface to contour plotting, same arguments as plc
  - zz 2-D array to be plotted
  - xx, yy Optional axis. Can either be 1-D or 2-D.
  - ireg Optional region. Must be same shape as zz
  - color='fg'
  - contours=8 Optional number of levels or list of levels
  - filled=0 When 1, draws filled contours
  - cmin, cmax: min and max of contours levels
  - leveloverlap=0. fractional overlap of contour levels
  """
  s = shape(zz)
  if len(s) != 2:
    print 'First argument must be a 2-Dimensional array'
    return
  if xx is None:
    xx = arange(s[0])[:,newaxis]*ones(s[1],'d')
  elif len(shape(xx))==1:
    xx = xx[:,newaxis]*ones(s[1],'d')
  if yy is None:
    yy = arange(s[1])*ones(s[0],'d')[:,newaxis]
  elif len(shape(yy))==1:
    yy = yy*ones(s[0],'d')[:,newaxis]
  if ireg is None:
    ireg = ones(s,'i')
  else:
    assert shape(ireg) == shape(zz),"Shape of ireg must be the same as zz"
  if contours is 0: contours = None
  if levs is not None: contours = levs
  if type(contours) == ListType: contours = array(contours)
  if type(contours) == TupleType: contours = array(contours)
  if type(contours) == type(1):
    # --- cmin and cmax are multiplied by 1. to force them to be standard
    # --- python floats, instead of zero length numpy arrays.
    if cmin is None: cmin = minnd(zz)*1.
    if cmax is None: cmax = maxnd(zz)*1.
    contours = 1.*iota(0,contours)*(cmax-cmin)/contours + cmin
  # --- leveloverlap is handled in a keyword dict so that it won't cause
  # --- problems with older version of pygist that don't support
  # --- the argument.
  if leveloverlap is None:
    kw = {}
  else:
    kw = {'leveloverlap':leveloverlap}
  if filled:
    # --- ireg must be of type integer because some legacy code used
    # --- expects it.
    ireg = ireg.astype('i')
    plfc(zz,xx,yy,ireg,contours=contours,local=local,**kw)
  else:
    plc(zz,xx,yy,ireg,color=color,levs=contours,width=width,type=linetype,
        local=local,**kw)

def plotfc(zz,xx=None,yy=None,ireg=None,contours=8,leveloverlap=None,local=1):
  """
Simple interface to filled contour plotting, same arguments as plfc
  - zz 2-D array to be plotted
  - xx, yy Optional axis. Can either be 1-D or 2-D.
  - ireg Optional region. Must be same shape as zz
  - color='fg'
  - contours Optional number of levels or list of levels
  """
  plotc(zz,xx=xx,yy=yy,ireg=ireg,color=color,contours=contours,filled=1,
        leveloverlap=leveloverlap,local=local)

# --- Define variables names for the allowed colors
fg = 'fg'
bg = 'bg'
white = 'white'
black = 'black'
red = 'red'
green = 'green'
blue = 'blue'
cyan = 'cyan'
magenta = 'magenta'
yellow = 'yellow'

########################################################################
########################################################################
########################################################################
# The next part of this file contains Python routines which emulate compiled
# plot routines.
#
# Here are the plots available so far
#
# ppzx(), ppzy(), ppzr(), ppzxp(), ppzvx(), ppzyp(), ppzvy(), ppzvz()
# ppxy(), ppxxp(), ppyyp(), ppxpyp(), ppxvx(), ppyvy()
# ppxvz(), ppyvz(), pprvz(), ppzxy(), pptrace()
# ppvxvy(), ppvxvz(), ppvyvz()
# ppco(y,x,z;uz,xmin,xmax,ymin,ymax,zmin,zmax)
#
# The following only work properly serially
#
# ppzxco(), ppzyco(), ppzvzco(), ppzxyco()
#
##########################################################################
# List of available named colors.
color = ["red","green","blue","cyan","magenta","yellow"]

########################################################################
# Note: Subtracted off 0.0337 from X position of titlel (10/21/99)
ptitle_placement = [
  [[0.3950, 0.9070], [0.3950, 0.4080], [0.1200, 0.6575], [0.3950, 0.3780]],
  [[0.3950, 0.9070], [0.3950, 0.4080], [0.1200, 0.6575], [0.3950, 0.3780]],
  [[0.2634, 0.8880], [0.2634, 0.6590], [0.1300, 0.7760], [0.3950, 0.3750]],
  [[0.5266, 0.8880], [0.5266, 0.6590], [0.3932, 0.7760], [0.3950, 0.3750]],
  [[0.2634, 0.6230], [0.2634, 0.3970], [0.1300, 0.5160], [0.3950, 0.3750]],
  [[0.5266, 0.6230], [0.5266, 0.3970], [0.3932, 0.5160], [0.3950, 0.3750]],
  [[0.2634, 0.8880], [0.2634, 0.3970], [0.1300, 0.6500], [0.3950, 0.3750]],
  [[0.5266, 0.8880], [0.5266, 0.3970], [0.3932, 0.6500], [0.3950, 0.3750]],
  [[0.3950, 0.8880], [0.3950, 0.6590], [0.1300, 0.7760], [0.3950, 0.3750]],
  [[0.3950, 0.6230], [0.3950, 0.3970], [0.1300, 0.5160], [0.3950, 0.3750]],
  [[0.3950, 0.9070], [0.3950, 0.4080], [0.1200, 0.6575], [0.3950, 0.3780]],
  [[0.3950, 0.9070], [0.3950, 0.4080], [0.1200, 0.6575], [0.3950, 0.3780]],
  [[0.3950, 0.9070], [0.3950, 0.4080], [0.1200, 0.6575], [0.3950, 0.3780]],
  [[0.3950, 0.9053], [0.3950, 0.7486], [0.1300, 0.7975], [0.3950, 0.3780]],
  [[0.3950, 0.7298], [0.3950, 0.5731], [0.1300, 0.6220], [0.3950, 0.3780]],
  [[0.3950, 0.5544], [0.3950, 0.3977], [0.1300, 0.4465], [0.3950, 0.3780]]]
ptitle_heights = [20.,20.,15.,15.,15.,15.,15.,15.,15.,15.,20.,20.,20.,15.,15.,15.]
  
default_titlet=""
default_titleb=""
default_titlel=""
default_titler=""
def settitles(titlet="",titleb="",titlel="",titler=""):
  "Sets titles which are plotted by ptitles"
  global default_titlet,default_titleb,default_titlel,default_titler
  if titlet is not None: default_titlet = titlet
  if titleb is not None: default_titleb = titleb
  if titlel is not None: default_titlel = titlel
  if titler is not None: default_titler = titler
def ptitles(titlet="",titleb="",titlel="",titler="",v=None,height=None):
  "Plots titles, either uses input or titles set by settitles"
  global framet,frameb,framel,framer
  if "ptitles" in _plotpackage.__dict__:
    _plotpackage.ptitles(titlet,titleb,titlel,titler,v)
    return
  if v is None: v = plsys()
  if titlet=="" and default_titlet: titlet = default_titlet
  if titleb=="" and default_titleb: titleb = default_titleb
  if titlel=="" and default_titlel: titlel = default_titlel
  if titler=="" and default_titler: titler = default_titler
  framet=titlet
  frameb=titleb
  framel=titlel
  framer=titler
  if height is None: height = ptitle_heights[v-1]
  if titlet:
    if with_gist:
      plt(titlet,ptitle_placement[v-1][0][0],ptitle_placement[v-1][0][1],
          justify="CC",orient=0,local=1,height=height)
    if with_matplotlib:
      pylab.title(titlet)
  if titleb:
    if with_gist:
      plt(titleb,ptitle_placement[v-1][1][0],ptitle_placement[v-1][1][1],
          justify="CC",orient=0,local=1,height=height)
    if with_matplotlib:
      pylab.xlabel(titleb + '\n' + titler)
  if titlel:
    if with_gist:
      plt(titlel,ptitle_placement[v-1][2][0],ptitle_placement[v-1][2][1],
          justify="CC",orient=1,local=1,height=height)
    if with_matplotlib:
      pylab.ylabel(titlel)
  if titler:
    if with_gist:
      plt(titler,ptitle_placement[v-1][3][0],ptitle_placement[v-1][3][1],
          justify="CC",orient=0,local=1,height=height)
  settitles()
def ptitlebottom(text=""):
  if with_gist:
    plt(text,0.3950,0.37,justify="CC",local=1)

##########################################################################
##########################   UTILITY ROUTINES  ###########################
##########################################################################
##########################################################################
def checkarguments(input,arglist):
  "Compare inputs against and argument list and return list of bad arguments"
  #inputcopy = input.copy()
  #for i in inputcopy:
  #  if i in arglist: del inputcopy[i]
  #return inputcopy
  result = {}
  for k,v in input.iteritems():
    if k not in arglist:
      result[k] = v
  return result
      

##########################################################################
# --- Complete dictionary of possible keywords and their default values
_pptitleright_kwdefaults = {"js":0,"win":None,"z":None,
                "ix":None,"wx":1.,"iy":None,"wy":1.,"iz":None,"wz":1.,
                "zl":None,"zu":None,"zc":None,"slope":0,
                'checkargs':0,'allowbadargs':0}
def pptitleright(iw=0,kwdict={},**kw):
  "Returns right plot title. Takes same arguments as selectparticles"

  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in _pptitleright_kwdefaults.
  kwvalues = _pptitleright_kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)
  #for arg in _pptitleright_kwdefaults: exec(arg+" = kwvalues['"+arg+"']")
  # --- This is MUCH faster than the loop, about 100x, but not as nice looking.
  js = kwvalues['js']
  win = kwvalues['win']
  z = kwvalues['z']
  ix = kwvalues['ix']
  wx = kwvalues['wx']
  iy = kwvalues['iy']
  wy = kwvalues['wy']
  iz = kwvalues['iz']
  wz = kwvalues['wz']
  zl = kwvalues['zl']
  zu = kwvalues['zu']
  zc = kwvalues['zc']
  slope = kwvalues['slope']
  checkargs = kwvalues['checkargs']
  allowbadargs = kwvalues['allowbadargs']

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.
  badargs = checkarguments(kwvalues,_pptitleright_kwdefaults)
  if checkargs: return badargs
  if badargs and not allowbadargs:
    raise TypeError,"bad argument%s"%' '.join(badargs.keys())

  # --- Return appropriate right title
  if zl is not None or zu is not None:
    if z is None: prefix = ""
    else: prefix = "z "
    if zl is None: zl = -top.largepos
    if zu is None: zu = +top.largepos
    result = prefix+"range (%9.4e, %9.4e)"%(zl,zu)
  elif ix is not None:
    xl = w3d.xmmin + ix*w3d.dx - wx*w3d.dx
    xu = w3d.xmmin + ix*w3d.dx + wx*w3d.dx
    result = "ix = %d, x range (%9.4e, %9.4e)"%(ix,xl,xu)
  elif iy is not None:
    yl = w3d.ymmin + iy*w3d.dy - wy*w3d.dy
    yu = w3d.ymmin + iy*w3d.dy + wy*w3d.dy
    result = "iy = %d, y range (%9.4e, %9.4e)"%(iy,yl,yu)
  elif iz is not None:
    zl = w3d.zmmin + iz*w3d.dz - wz*w3d.dz + top.zbeam
    zu = w3d.zmmin + iz*w3d.dz + wz*w3d.dz + top.zbeam
    result = "iz = %d, z range (%9.4e, %9.4e)"%(iz,zl,zu)
  elif zc is not None:
    zl = zc - wz*w3d.dz
    zu = zc + wz*w3d.dz
    result = "zc = %9.4e, z range (%9.4e, %9.4e)"%(zc,zl,zu)
  elif iw < 0:
    result = "subset %d: %d particles"%(-iw,top.npplot[-iw-1])
  else:
    if win is None:
      win = top.zwindows[:,iw] + top.zbeam
      prefix = "z "
    else:
      prefix = ""
    if len(shape(win)) == 2: win = win[:,iw]
    result = prefix+"window%d = %9.4e, %9.4e"%(iw,win[0],win[1])
  if slope != 0:
    result = result + ", slope=%7.4f"%slope
  return result

#-------------------------------------------------------------------------
def ppmoments(text):
  "Plots text in upper right hand corner of the plot"
  plt(text,0.61,.855,justify="RT",height=12,font="courierB",local=1)

#############################################################################
#############################################################################
#############################################################################
#-------------------------------------------------------------------------
# --- Complete dictionary of possible keywords and their default values
_ppgeneric_kwdefaults = {'zz':None,'weights':None,'grid':None,'gridt':None,
                'nx':20,'ny':20,'slope':0.,
                'xoffset':0.,'yoffset':0.,'offset':0.,
                'xscale':1.,'yscale':1.,'titles':1,
                'titlet':'','titleb':'','titlel':'','titler':'',
                'height':20,
                'lframe':0,'xmin':None,'xmax':None,'ymin':None,'ymax':None,
                'pplimits':('e','e','e','e'),
                'particles':0,'uselog':None,'logmin':None,
                'color':'fg','ncolor':None,
                'usepalette':1,'npalette':200,'marker':'\1','msize':1.0,
                'denmin':None,'denmax':None,'chopped':None,
                'hash':0,'line_scale':.9,'hcolor':'fg','width':1.0,
                'contours':None,'filled':0,'ccolor':'fg','leveloverlap':None,
                'cellarray':0,'centering':'node','ctop':199,
                'cmin':None,'cmax':None,'ireg':None,
                'xbound':dirichlet,'ybound':dirichlet,
                'ldensityscale':False,'ldensitycylindrical':False,'gridscale':None,
                'flipxaxis':0,'flipyaxis':0,
                'xcoffset':0.,'ycoffset':0.,
                'view':None,
                'lcolorbar':1,'colbarunitless':0,'colbarlinear':1,'surface':0,
                'xmesh':None,'ymesh':None,
                'returngrid':0,'local':1,
                'checkargs':0,'allowbadargs':0}

def ppgeneric(y=None,x=None,kwdict={},**kw):
  """
Generic particle plotting routine. Allows plotting of particle points, density
contours, and/or density hash marks.
Note that either the x and y coordinates or the grid must be passed in.
  - y, x: optional particle data (instead of using inputted grid)
  - ccolor='fg': contour color (when not filled)
  - cellarray=0: when true, plot grid as cell array, filling each grid cell
                 with a color determined by the value at that cell
  - centering='node': centering of cells with cellarray, other option are
                      'cell' and 'old' (for old incorrect scaling)
  - chopped=None: only particles where r < chopped*maxdensity/density are
                  plotted, where r is a random number between 0 and 1 and
                  density is the density at the particle location
  - cmin=min(grid), cmax=max(grid): min and max of data for coloration for
                                    contours and cellarray. Can be used to
                                    crop the data range for coloration.
  - colbarlinear=1: when true, the color-bar is laid out linearly in density,
                    otherwise each contour level gets an equal sized area.
                    Only in effect when a list of colorbars is specified.
  - colbarunitless=0: when true, color-bar scale is unitless
  - color='fg': color of particles, when color='density', the color of
                particles is determined by the number density. If zz is given
                and color='density', color will be determined by zz.
                Values include red, green, blue, cyan, magenta, yellow.
  - contours=None: number of countour levels to plot
  - ctop=199: max color table index for cellarray plot
  - denmin, denmax: thresholds for removing particles, only particles located
                    where density is between denmin and denmax are plotted
  - filled=0: when true, plot filled contours
  - flipxaxis=0: when true, flips gridded data about the x-axis
  - flipyaxis=0: when true, flips gridded data about the y-axis
  - grid: optional grid to plot (instead of deriving grid from particle data)
  - gridscale=None: scale factor applied to gridded data.
  - gridt: optional grid to plot (instead of deriving grid from particle data)
           The transpose is the grid is plotted.
  - hash=0: flag to turn on or off the hash plot
  - height=20: font size of the axes titles
  - hcolor='fg': color of hash marks for hash plots
  - lcolorbar=1: when plotting colorized data, include a colorbar
  - ldensityscale=False: when true, scale the density by its max.
  - ldensitycylindrical=False: when true, the density is calculated on a
                               cylindrical grid
  - leveloverlap=0.: fractional overlap of filled contour levels
  - lframe=0: when true, the plot limits are set to the plmin and plmax input
              arguments, which default to the plmin and plmax variables from
              the group InDiag
  - line_scale=.9: scaling factor on line length for hash plots
  - local=None: Forces the plotting to be local or non-local (parallel).
                Otherwise, particle plots are non-local and grid plots are
                local. (experts only)
  - logmin=None: when given, and with uselog, values less than logmin are
                 truncated.
  - marker=dot: marker to plot for particles, other options include circle,
                plus, star, or any quoted character, i.e. 'a'.
  - msize=1.: scaled size of marker
  - ncolor=None: when plotting particle color by number density, number of
                 colors to use, defaults to top.ncolor
  - npalette=200: number of colors to use in the palette when plotting
                  particles with colors (passed through to ppco)
  - nx=20, ny=20: grid size when the density is calculated,
  - particles=0: when true, forces plot particles. If a gridded plot (such as
                 contour of cellarray) plot is requested, the particles will
                 also be plotted if this is true.
  - pplimits=None: a tuple of (xplmin, xplmax, yplmin, yplmax), limits of
                   plot range (used when lframe=1)
  - returngrid=0: when true, and when particle data is passed in and a plot
                  which requires a grid is requested (such as a contour
                  plot), no plotting is done and the grid and extrema are
                  returned in a tuple
  - slope=0.: slope to subtract from y coordinate (y-slope*x),
              for example to skew the particles in a phase-space plot.
  - surface=0: when true, a 3-d surface plot is made of the gridded data
               Note: to remove window, use the hidesurfaces() command rather
               than closing the window.
  - titles=1: when true, plot the titles
  - titlet,titleb,titlel,titler='': If specified, added to plot, overriding
                                    other title settings.
  - uselog=None: when given, logarithmic levels of the number density are
                 used.  The value gives the log base, 1 is same as 'e'.
  - usepalette=1: when plotting particles with color, when true, uses palette,
                  otherwise uses colors in the array color (passed through to
                  ppco)
  - view=None: view window to use (experts only)
               When None, uses the current view. Otherwise, the value is passed
               to plsys.
  - width=1.0: width of hash marks for hash plots
  - xbound=dirichlet: sets boundary condition on gridded data for x
  - xcoffset,ycoffset=0: offsets of coordinates in grid plots
  - xmin, xmax, ymin, ymax: extrema of density grid, defaults to particle
                            extrema
  - xoffset=0.: average x of particles
  - xscale=1.: scaling factor applied to x data
  - ybound=dirichlet: sets boundary condition on gridded data for y
  - yoffset=0.: average y of particles
  - yscale=1.: scaling factor applied to y data
  - zz: optional third particle data quantity - when supplied, it is
        deposited on a grid and that is used for contour levels, except if
        color='density' is specified, then zz is used directly to color the
        particles
  """
  # --- Create dictionary of local values and copy it into local dictionary,
  # --- ignoring keywords not listed in _ppgeneric_kwdefaults.
  kwvalues = _ppgeneric_kwdefaults.copy()
  kwvalues.update(kw)
  kwvalues.update(kwdict)

  #for arg in _ppgeneric_kwdefaults: exec(arg+" = kwvalues['"+arg+"']")
  # --- This is MUCH faster than the loop, about 100x, but not as nice looking.
  zz = kwvalues['zz']
  weights = kwvalues['weights']
  grid = kwvalues['grid']
  gridt = kwvalues['gridt']
  nx = kwvalues['nx']
  ny = kwvalues['ny']
  slope = kwvalues['slope']
  xoffset = kwvalues['xoffset']
  yoffset = kwvalues['yoffset']
  offset = kwvalues['offset']
  xscale = kwvalues['xscale']
  yscale = kwvalues['yscale']
  titles = kwvalues['titles']
  titlet = kwvalues['titlet']
  titleb = kwvalues['titleb']
  titlel = kwvalues['titlel']
  titler = kwvalues['titler']
  height = kwvalues['height']
  lframe = kwvalues['lframe']
  xmin = kwvalues['xmin']
  xmax = kwvalues['xmax']
  ymin = kwvalues['ymin']
  ymax = kwvalues['ymax']
  pplimits = kwvalues['pplimits']
  particles = kwvalues['particles']
  uselog = kwvalues['uselog']
  logmin = kwvalues['logmin']
  color = kwvalues['color']
  ncolor = kwvalues['ncolor']
  usepalette = kwvalues['usepalette']
  npalette = kwvalues['npalette']
  marker = kwvalues['marker']
  msize = kwvalues['msize']
  denmin = kwvalues['denmin']
  denmax = kwvalues['denmax']
  chopped = kwvalues['chopped']
  hash = kwvalues['hash']
  line_scale = kwvalues['line_scale']
  hcolor = kwvalues['hcolor']
  width = kwvalues['width']
  contours = kwvalues['contours']
  filled = kwvalues['filled']
  ccolor = kwvalues['ccolor']
  leveloverlap = kwvalues['leveloverlap']
  cellarray = kwvalues['cellarray']
  centering = kwvalues['centering']
  ctop = kwvalues['ctop']
  cmin = kwvalues['cmin']
  cmax = kwvalues['cmax']
  ireg = kwvalues['ireg']
  xbound = kwvalues['xbound']
  ybound = kwvalues['ybound']
  ldensityscale = kwvalues['ldensityscale']
  ldensitycylindrical = kwvalues['ldensitycylindrical']
  gridscale = kwvalues['gridscale']
  flipxaxis = kwvalues['flipxaxis']
  flipyaxis = kwvalues['flipyaxis']
  xcoffset = kwvalues['xcoffset']
  ycoffset = kwvalues['ycoffset']
  view = kwvalues['view']
  lcolorbar = kwvalues['lcolorbar']
  colbarunitless = kwvalues['colbarunitless']
  colbarlinear = kwvalues['colbarlinear']
  surface = kwvalues['surface']
  xmesh = kwvalues['xmesh']
  ymesh = kwvalues['ymesh']
  returngrid = kwvalues['returngrid']
  local = kwvalues['local']
  checkargs = kwvalues['checkargs']
  allowbadargs = kwvalues['allowbadargs']

  # --- Check the argument list for bad arguments.
  # --- 'checkargs' allows this routine to be called only to check the
  # --- input for bad arguments.
  # --- 'allowbadargs' allows this routine to be called with bad arguments.
  # --- These are intentionally undocumented features.
  badargs = checkarguments(kwvalues,_ppgeneric_kwdefaults)
  if checkargs: return badargs
  assert (not badargs or allowbadargs), \
         "bad argument: %s"%' '.join(badargs.keys())

  # --- If gridt is given, take the transpose and put it in grid. Note that
  # --- this will overwrite a grid argument. This is done here to reduce
  # --- the code complexity below. If gridt is specified, it is equivalent
  # --- to specifying grid (except for the transpose).
  if gridt is not None: grid = transpose(gridt)

  # --- If ireg is passed in, get its transpose, since only it will be used.
  if ireg is not None: iregt = transpose(ireg)
  else:                iregt = None

  # --- If y is a 2-d array and x is not input, then assume that the user
  # --- intends to plot gridded data.
  # --- Not sure yet if this is a good idea.
  if y is not None and y.ndim == 2 and x is None:
    grid = y
    y = None

  # --- Do some error checking on the consistency of the input
  assert (type(grid) == ArrayType or \
          (type(x) == ArrayType and type(y) == ArrayType)), \
         "either the grid and/or both x and y must be specified"
  assert (not particles or (type(x) == ArrayType and type(y) == ArrayType)), \
         "both x and y must be specified if particles are to be plotted"
  assert ((type(x) != ArrayType and type(y) != ArrayType) or x.size == y.size),\
         "both x and y must be of the same length"
  assert (zz is None) or (type(zz) == ArrayType and zz.size == x.size),\
         "zz must be the same length as x"
  assert (type(slope) != StringType),"slope must be a number"
  assert (zz is None) or (grid is None),\
         "only one of zz and grid can be specified"
  assert (centering == 'node' or centering == 'cell' or centering == 'old'),\
         "centering must take one of the values 'node', 'cell', or 'old'"
  assert (grid is None or grid.ndim == 2), \
         "the grid specified must be two dimensional"

  # --- If there are no particles and no grid to plot, just return
  if type(x) == ArrayType and type(y) == ArrayType:
    if local:
      np = x.size
    else:
      np = globalsum(x.size)
  else:
    np = 0
  if np == 0 and grid is None: return

  # --- If filled is turned on, but contours is not set, set it to the
  # --- default value of 8.
  if filled and contours is None: contours = 8

  # --- Make sure that contours is not zero, which breaks some code.
  if contours is 0: contours = None

  # --- If particle data was passed in and no specific plots were requested,
  # --- just plot the particles.
  if y is not None and \
     (not hash and contours is None and not surface and not cellarray):
    particles = 1

  # --- If a grid is passed in and no specific plots were requested,
  # --- make a cellarray plot.
  if grid is not None and \
     (not hash and contours is None and not surface and not cellarray
      and not particles):
    cellarray = 1

  # --- Whether a grid plot is parallel or not depends on the input.
  # --- If particle data is input, then the parallel work is handled
  # --- in ppgeneric and the plot is then local. grid_local will be set
  # --- below in that code handling the parallelism. Otherwise, use the
  # --- default value or the one passed in. An input grid can be either
  # --- local or parallel.
  grid_local = local

  # --- Make sure that nothing is not plotted over a surface plot
  if surface:
    particles = 0
    contours = None
    cellarray = 0
    hash = 0
    lframe = 0
    titles = 0

  # -- Set the plotting view window
  if view is not None:
    plsys(view)

  # --- Make sure that the grid size nx and ny are consistent with grid
  # --- is one is input
  if type(grid) == ArrayType:
    nx = shape(grid)[0] - 1
    ny = shape(grid)[1] - 1

  # --- Calculate extrema of the particles
  if type(x) == ArrayType and type(y) == ArrayType:
    # --- Get slope subtracted value of y
    yms = y - x*slope + (xoffset*slope - yoffset - offset)
    # --- For the particles, get mins and maxs that were not supplied by the user.
    if lparallel and not local:
      if xmin is None: xmintemp = globalmin(x)
      if xmax is None: xmaxtemp = globalmax(x)
      if ymin is None: ymintemp = globalmin(yms)
      if ymax is None: ymaxtemp = globalmax(yms)
    else:
      xmintemp = 0.
      xmaxtemp = 0.
      ymintemp = 0.
      ymaxtemp = 0.
      if xmin is None and x.size > 0: xmintemp = x.min()
      if xmax is None and x.size > 0: xmaxtemp = x.max()
      if ymin is None and yms.size > 0: ymintemp = yms.min()
      if ymax is None and yms.size > 0: ymaxtemp = yms.max()
    # --- When neither the min or max are supplied by the user, extend
    # --- extrema by one grid cell so that all particles are within the
    # --- limits of the grid. This is the most common case.
    if xmin is None and xmax is None:
      xmintemp = xmintemp - (xmaxtemp-xmintemp)/(nx-2)
      xmaxtemp = xmaxtemp + (xmaxtemp-xmintemp)/(nx-2)
      if xmintemp == xmaxtemp:
        # --- If there is probably only one particle, leaving
        # --- xmintemp==xmaxtemp, which can cause problems.
        xmintemp -= 1.
        xmaxtemp += 1.
    if ymin is None and ymax is None:
      ymintemp = ymintemp - (ymaxtemp-ymintemp)/(ny-2)
      ymaxtemp = ymaxtemp + (ymaxtemp-ymintemp)/(ny-2)
      if ymintemp == ymaxtemp:
        # --- If there is probably only one particle, leaving
        # --- ymintemp==ymaxtemp, which can cause problems.
        ymintemp -= 1.
        ymaxtemp += 1.
    # --- Now set main versions of min and max
    if xmin is None: xmin = xmintemp
    if xmax is None: xmax = xmaxtemp
    if ymin is None: ymin = ymintemp
    if ymax is None: ymax = ymaxtemp

    # --- Scale the data
    x = x*xscale
    yms = yms*yscale
  else:
    # --- If no particles are inputted and the extrema are not set, then
    # --- can only make a guess.
    if xmin is None: xmin = 0
    if xmax is None: xmax = nx
    if ymin is None: ymin = 0
    if ymax is None: ymax = ny

  # --- Scale the extrema
  xmin = xmin*xscale
  xmax = xmax*xscale
  ymin = ymin*yscale
  ymax = ymax*yscale

  # --- Get grid cell sizes
  if nx != 0: dx = (xmax-xmin)/nx
  else:       dx = 1.
  if ny != 0: dy = (ymax-ymin)/ny
  else:       dy = 1.

  # --- Calculate the density grid. This is needed if a grid or zz quantity
  # --- is not input and a gridded plot is being made. The gridded plots are
  # --- are assumed to be density plots. In this case, grid will be the same
  # --- as densitygrid. The density grid is also needed if chopped or
  # --- denmin or max are specified, which always operate on the density.
  if (((type(grid) != ArrayType and zz is None) and
       (hash or contours is not None or surface or cellarray or
        color=='density'))
      or chopped or denmin or denmax):
    # --- Create space for data
    densitygrid = fzeros((1+nx,1+ny),'d')

    # --- Deposit the density onto the grid.
    if ldensitycylindrical:
      if weights is None:
        setgrid2dcylindrical(x.size,x,yms,nx,ny,densitygrid,xmin,xmax,ymin,ymax)
      else:
        setgrid2dcylindricalw(x.size,x,yms,weights,nx,ny,densitygrid,xmin,xmax,ymin,ymax)
    else:
      if weights is None:
        setgrid2d(x.size,x,yms,nx,ny,densitygrid,xmin,xmax,ymin,ymax)
      else:
        setgrid2dw(x.size,x,yms,weights,nx,ny,densitygrid,xmin,xmax,ymin,ymax)

    # --- If parallel, do a reduction on the grid
    if lparallel and not local:
      try:
        parallelsumrealarray(densitygrid,size(densitygrid))
      except:
        densitygrid = parallelsum(densitygrid)
      # --- Set grid_local so that grid plots will now be done locally
      grid_local = 1

    if (type(grid) != ArrayType and zz is None): grid = densitygrid

  else:
    densitygrid = None

  # --- Calculate a grid based on the input zz quantity when a gridded plot
  # --- is being made. The exception is color=='density', in which case the
  # --- color is taken directly from the zz quantity.
  if ((zz is not None) and
       (hash or contours is not None or surface or cellarray)):

    # --- Create space for data
    grid = fzeros((1+nx,1+ny),'d')
    gridcount = fzeros((1+nx,1+ny),'d')

    # --- Deposit the data onto the grid. itask is 1 so that the parallel
    # --- version can be done properly.
    if(weights is None):
      deposgrid2d(1,x.size,x,yms,zz,nx,ny,grid,gridcount,xmin,xmax,ymin,ymax)
    else:
      deposgrid2dw(1,x.size,x,yms,zz,weights,nx,ny,grid,gridcount,xmin,xmax,ymin,ymax)

    # --- If parallel, do a reduction on the grid
    if lparallel and not local:
      try:
        parallelsumrealarray(grid,size(grid))
        parallelsumrealarray(gridcount,size(gridcount))
      except:
        grid = parallelsum(grid)
        gridcount = parallelsum(gridcount)
      # --- Set grid_local so that grid plots will now be done locally
      grid_local = 1

    # --- Divide out the particle counts by hand.
    grid = grid/where(greater(gridcount,0.),gridcount,1.)

  # --- Enforce boundary conditions on the densitygrid. This operation doesn't
  # --- make sense on anything other than the density grid.
  if densitygrid is not None:
    if xbound == neumann:
      densitygrid[0,:] = 2.*densitygrid[0,:]
      densitygrid[-1,:] = 2.*densitygrid[-1,:]
    elif xbound == periodic:
      densitygrid[0,:] = densitygrid[0,:] + densitygrid[-1,:]
      densitygrid[-1,:] = densitygrid[0,:]
    if ybound == neumann:
      if not ldensitycylindrical:
        densitygrid[:,0] = 2.*densitygrid[:,0]
      densitygrid[:,-1] = 2.*densitygrid[:,-1]
    elif ybound == periodic and not ldensitycylindrical:
      densitygrid[:,0] = densitygrid[:,0] + densitygrid[:,-1]
      densitygrid[:,-1] = densitygrid[:,0]

  # --- If requested, return the grid and extrema, doing no plotting
  if returngrid: return (grid,xmin,xmax,ymin,ymax)

  # --- Scale the grid by its maximum if requested.
  if ldensityscale and densitygrid is not None:
    densitygridmax = maxnd(abs(densitygrid))
    if densitygridmax != 0.:
      densitygrid[:,:] = densitygrid/densitygridmax
  if ldensityscale and grid is not None:
    gridmax = maxnd(abs(grid))
    if gridmax != 0.:
      grid[:,:] = grid/gridmax

  # --- Apply grid scale factor if supplied
  # --- Note that a new array is created so that a grid passed in is not
  # --- modified. Also, any connection between grid and densitygrid is broken.
  if gridscale is not None: grid = grid*gridscale

  # --- If using logarithmic levels, take the log of the grid data.
  if uselog is not None and grid is not None:
    if uselog == 'e' or uselog == 1.: logscale = 1.
    else:                             logscale = log(uselog)
    if grid is densitygrid:
      # --- Take the log, raising all values below 0.1 to 0.1. The
      # --- threshold is used so that none of the elements are zero.
      # --- That value 0.1 is used since values any smaller do not have
      # --- much meaning since a value of 1.0 means that there is already
      # --- only one particle in that cell. The user can reset logmin though.
      if logmin is None: logmin = 0.1
      grid = log(where(less(grid,logmin),logmin,grid))/logscale
    else:
      # --- Before taking the log of the user supplied grid data, make sure
      # --- that there are no negative values. Zero is ok since they will
      # --- be replaced with a minimum value.
      if logmin is None:
        dmax = maxnd(grid)
        logmin = minnd(where(equal(grid,0.),dmax,grid))/10.
        assert logmin>0.,\
               "Can't take log since the grid has negative values"
      grid = log(where(less(grid,logmin),logmin,grid))/logscale

  # --- Flip data and plot limits about axis if requested.
  # --- Note that iregt had already been transposed.
  if flipxaxis:
    xmin = -xmin
    xmax = -xmax
    dx = -dx
    if xmesh is not None: xmesh = -xmesh
  if flipyaxis:
    ymin = -ymin
    ymax = -ymax
    dy = -dy
    if ymesh is not None: ymesh = -ymesh

  # --- Get grid min and max and generate contour levels if needed.
  if grid is not None:
    # --- cmin and cmax are multiplied by 1. to force them to be standard
    # --- python floats, instead of zero length numpy arrays.
    if cmin is None: cmin = minnd(grid)*1.
    if cmax is None: cmax = maxnd(grid)*1.
  elif zz is not None:
    if cmin is None and zz.size > 0: cmin = zz.min()*1.
    if cmax is None and zz.size > 0: cmax = zz.max()*1.
  ppgeneric.cmin = cmin
  ppgeneric.cmax = cmax

  # --- Get grid mesh if it is needed
  if contours is not None or hash or surface or cellarray:
    if xmesh is not None or ymesh is not None: usermesh = 1
    else:                                      usermesh = 0
    # --- The offsets are added in the way they are incase they are arrays.
    # --- Though of course they must be the correct length.
    if xmesh is None:
      xmesh = xmin + dx*arange(nx+1)[:,newaxis]*ones(ny+1,'d') + xcoffset
    else:
      if rank(xmesh) == 1:
        xmesh = xmesh[:,newaxis]*ones(ny+1,'d')
    if ymesh is None:
      ymesh = (ymin + dy*arange(ny+1)*ones(nx+1,'d')[:,newaxis] +
               transpose([ycoffset]))
    else:
      if rank(ymesh) == 1:
        ymesh = ymesh*ones(nx+1,'d')[:,newaxis]

  # --- Make filled contour plot of grid first since it covers everything
  # --- plotted before it.
  if contours is not None and filled and nx > 1 and ny > 1:
    if cmax != cmin:
      plotc(transpose(grid),transpose(ymesh),transpose(xmesh),iregt,
            color=ccolor,contours=contours,filled=filled,cmin=cmin,cmax=cmax,
            leveloverlap=leveloverlap,local=grid_local)

  # --- Make cell-array plot. This also is done early since it covers anything
  # --- done before it. The min and max are adjusted so that the patch for
  # --- each grid cell has the correct centering.
  # --- If the user supplies a mesh, then use plf, a filled mesh plot, since
  # --- the meshes may not be Cartesian.
  if cellarray and nx > 1 and ny > 1:
    if not usermesh:
      if centering == 'node':
        xminc = xmin
        xmaxc = xmax + dx
        yminc = ymin
        ymaxc = ymax + dy
      elif centering == 'cell':
        xminc = xmin - dx/2.
        xmaxc = xmax + dx/2.
        yminc = ymin - dy/2.
        ymaxc = ymax + dy/2.
      elif centering == 'old':
        xminc = xmin
        xmaxc = xmax
        yminc = ymin
        ymaxc = ymax
      pli(transpose(grid),xminc,yminc,xmaxc,ymaxc,top=ctop,cmin=cmin,cmax=cmax,
          local=grid_local)
    else:
      plf(grid,ymesh,xmesh,local=grid_local)

  # --- Plot particles
  if particles:
    if color == 'density':
      if zz is None:
        z1 = zeros(x.size,'d')
        getgrid2d(x.size,x,yms,z1,nx,ny,grid,xmin,xmax,ymin,ymax)
      else:
        z1 = zz
    if chopped or denmin or denmax:
      dd = zeros(x.size,'d')
      getgrid2d(x.size,x,yms,dd,nx,ny,densitygrid,xmin,xmax,ymin,ymax)
      maxdensity = maxnd(densitygrid)
      dd = dd/maxdensity
      ipick = ones(shape(x),'l')
      if chopped:
        ipick[:] = ipick*less(random.random(shape(x)),chopped/dd)
      if denmin:
        ipick[:] = ipick*less(denmin,dd)
      if denmax:
        ipick[:] = ipick*less(dd,denmax)
      x = compress(ipick,x)
      yms = compress(ipick,yms)
      if color == 'density':
        z1 = compress(ipick,z1)
    if color == 'density':
      if ncolor is None: ncolor = top.ncolor
      # --- Plot particles with color based on the density from the grid.
      ppco(yms,x,z1,uz=1.,marker=marker,msize=msize,zmin=cmin,zmax=cmax,
           ncolor=ncolor,usepalette=usepalette,npalette=npalette,local=local)
    else:
      # --- Plot particles as a solid color.
      plp(yms,x,color=color,marker=marker,msize=msize,local=local)

  # --- Now plot unfilled contours, which are easier to see on top of the
  # --- particles
  if contours is not None and not filled and nx > 1 and ny > 1:
    if cmax != cmin:
      plotc(transpose(grid),transpose(ymesh),transpose(xmesh),iregt,
            color=ccolor,contours=contours,filled=filled,cmin=cmin,cmax=cmax,
            leveloverlap=leveloverlap,local=grid_local)

  # --- Plot hash last since it easiest seen on top of everything else.
  if hash:
    # --- Set line length
    if nx != 0 and cmax != cmin:
      sss = line_scale*(xmax-xmin)/nx/(cmax - cmin)
    else:
      sss = 1.
    # --- Make plot of tick marks
    for ix in range(nx+1):
      for iy in range(ny+1):
        plg(ymesh[ix,iy]+zeros(2),xmesh[ix,iy]+array([0.,sss*grid[ix,iy]]),
            color=hcolor,width=width)

  # --- Add colorbar if needed
  if (lcolorbar and
     ((contours is not None and filled==1) or
      (color == 'density') or
      (cellarray))):
    if (contours is not None and filled==1):
      try:
        nc = contours.size + 1
        levs = contours
      except AttributeError:
        nc = contours
        levs = None
    elif (color == 'density'):
      if ncolor is None: ncolor = top.ncolor
      nc = ncolor
      levs = None
    elif (cellarray):
      nc = ctop
      levs = None
    if colbarunitless:
      dmin = 0.
      dmax = 1.0
    elif cellarray:
      dmin = cmin
      dmax = cmax
    else:
      dmin = cmin
      dmax = cmax
    colorbar(dmin,dmax,uselog=uselog,ncolor=nc,view=view,levs=levs,
             colbarlinear=colbarlinear,ctop=ctop)

  # --- Make surface plot
  if surface and me == 0 and nx > 1 and ny > 1:
    try:
      import Opyndx
      if type(color) != ListType: scolor = None
      else:                       scolor = color
      xrange = 1.5*max(abs(xmin),abs(xmax))
      yrange = 1.5*max(abs(ymin),abs(ymax))
      zrange = 1.5*maxnd(abs(grid))
      vo = Opyndx.DXMountainPlot(f=grid,xmin=xmin,ymin=ymin,dx=dx,dy=dy)
    except ImportError:
      from gist import pl3d
      from gist import plwf
      pl3d.orient3()
      pl3d.light3()
      plwf.plwf(grid,xmesh,ymesh,fill=grid,edges=0)
      [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
      #limits(xmin3,xmax3,ymin3,ymax3)

  # --- Finish off the plot, adding titles and setting the frame limits.
  if titles: ptitles(titlet,titleb,titlel,titler,v=view,height=height)
  settitles()
  if (lframe):
    ppp = list(pplimits)
    if ppp[0] != 'e': ppp[0] = ppp[0]*xscale
    if ppp[1] != 'e': ppp[1] = ppp[1]*xscale
    if ppp[2] != 'e': ppp[2] = ppp[2]*yscale
    if ppp[3] != 'e': ppp[3] = ppp[3]*yscale
    limits(ppp[0],ppp[1],ppp[2],ppp[3])


#############################################################################
#############################################################################
def ppvector(gridy=None,gridx=None,ymesh=None,xmesh=None,scale=None,
             xmin=None,xmax=None,ymin=None,ymax=None,**kw):
  """
Generic vector plotting routine.
Note that both the x and y grids must be passed in.
The input can be 1-d or 2-d. If 1-d, ymesh, xmesh and scale must all be
supplied.
  - gridy, gridx: x and y vector comnponents
  - ymesh, xmesh: x and y coordinates
  - scale: scale length of the vectors. This must be supplied when plotting
           1-d data. For 2-d, the default will set the length of the largest
           vector to the grid cell size (scale=(cell size)/(max vector)).
  - xmin,xmax,ymin,ymax: grid extrema; only used if ymesh and xmesh are not
                         specified
  """
  # --- Do some error checking on the consistency of the input
  assert gridx is not None and gridy is not None,"both gridx and gridy must be specified"
  assert len(shape(gridx)) == len(shape(gridy)),"gridx and gridy must have the same number of dimensions"
  assert all(shape(gridx) == shape(gridy)),"gridx and gridy must be the same shape"

  if len(shape(gridx)) == 1:
    # --- handle 1-D arrays
    assert ymesh is not None and xmesh is not None,"if gridx is 1-d, xmesh and ymesh must be supplied"
    assert scale is not None,"if gridx is 1-d, scale must be supplied"
    assert len(gridx) > 0,"at least one point must be given"
    
    # --- This is a bit of hackery. plv is expecting a 2-d array with both
    # --- dimensions at least 2. So, the input arrays need to be reshaped
    # --- appropriately, but this is tricky since the number of points
    # --- must be factorable.
    # --- If the number of points is even, then the shape can be (n/2,2).
    # --- If the number of points is odd, replicate the last point and then
    # --- do shape = ((n+1)/2,2).
    # --- If the number of points is less than 4, replicate as many points
    # --- is needed to get 4 points.
    n = len(gridx)
    while n%2 == 1 or n < 4:
      newgridx = zeros(n+1,dtype=array(gridx,copy=False).dtype)
      newgridy = zeros(n+1,dtype=array(gridy,copy=False).dtype)
      newxmesh = zeros(n+1,dtype=array(xmesh,copy=False).dtype)
      newymesh = zeros(n+1,dtype=array(ymesh,copy=False).dtype)
      newgridx[:n] = gridx
      newgridy[:n] = gridy
      newxmesh[:n] = xmesh
      newymesh[:n] = ymesh
      newgridx[n] = gridx[-1]
      newgridy[n] = gridy[-1]
      newxmesh[n] = xmesh[-1]
      newymesh[n] = ymesh[-1]
      gridx = newgridx
      gridy = newgridy
      xmesh = newxmesh
      ymesh = newymesh
      n += 1

    gridx.shape = (n/2,2)
    gridy.shape = (n/2,2)
    xmesh.shape = (n/2,2)
    ymesh.shape = (n/2,2)
    plv(gridy,gridx,ymesh,xmesh,scale=scale,**kw)

  else:
    # --- handle 2-D arrays
    nx = shape(gridx)[0] - 1
    ny = shape(gridx)[1] - 1

    if xmesh is None:
      if xmin is None: xmin = 0
      if xmax is None: xmax = nx
      dx = (xmax - xmin)/nx
      xmesh = span(xmin,xmax,nx+1)
    else:
      dx = (maxnd(xmesh) - minnd(xmesh))/nx

    if ymesh is None:
      if ymin is None: ymin = 0
      if ymax is None: ymax = ny
      dy = (ymax - ymin)/ny
      ymesh = span(ymin,ymax,ny+1)
    else:
      dy = (maxnd(ymesh) - minnd(ymesh))/ny

    if len(shape(xmesh)) == 1:
      xmesh = xmesh[:,newaxis]*ones(ny+1)[newaxis,:]
    if len(shape(ymesh)) == 1:
      ymesh = ymesh[newaxis,:]*ones(ny+1)[:,newaxis]

    # --- Compute scale
    if scale is None:
      scale = min(dx,dy)/dvnz(max(maxnd(abs(gridx)),maxnd(abs(gridy))))
    #print scale

    # --- Make plot
    plv(gridy,gridx,ymesh,xmesh,scale=scale,**kw)

#############################################################################
#############################################################################
# ed williams' colorbar stuff / modified for Warp by J.-L. Vay on 01/22/2001
def nicelevels(z,n=8) :
  """nicelevels(z,n=8) finds approximately n "nice values"
between min(z) and max(z) for axis labels. n defaults to eight.
  """
  zmin = z.min()
  zmax = z.max()
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
     #array_set(levs,llist,0.0)
     for i in llist: levs[i] = 0.
  return levs

#-----------------------------------------------------------------------
colorbar_placement = [[0.62,0.64,0.4387,0.8713],[0.62,0.64,0.43,0.86],
                      [0.354,0.364,0.692,0.859],[0.617,0.627,0.692,0.859],
                      [0.354,0.364,0.428,0.596],[0.617,0.627,0.428,0.596],
                      [0.354,0.364,0.43,0.86],[0.617,0.627,0.43,0.86],
                      [0.617,0.627,0.692,0.859],[0.617,0.627,0.428,0.596],
                      [],[],[],
                      [0.617,0.627,0.762,0.879],
                      [0.617,0.627,0.587,0.703],
                      [0.617,0.627,0.411,0.528]]

colorbar_fontsize = [14.,14.,11.,11.,11.,11.,11.,11.,11.,11.,14.,14.,14.,8.,8.,8.]

def colorbar(zmin,zmax,uselog=None,ncolor=100,view=None,levs=None,
             colbarlinear=1,ctop=199):
  """
Plots a color bar to the right of the plot square labelled by the z
values from zmin to zmax.
  - zmin, zmax: lower and upper range for color bar
  - uselog=None: when true, labels are printed in the form b^x where b (the
                 base) is given by uselog.
  - ncolor=100: default number of colors to include
  - view=None: specifies the view that is associated with the color bar
               If not given, uses the current view (from plsys).
  - levs: an optional list of color levels
  - ctop=199: number of colors from palette to use
  """
  # --- This is only ever done on processor 0, so otherwise return
  if me > 0: return
  # --- The builtin colorbar is used with pylab
  if with_matplotlib:
    pylab.colorbar(pad=0.02,fraction=0.08)
    return
  if view is None: view = plsys()
  plsys(0)
  xmin,xmax,ymin,ymax = colorbar_placement[view-1]
  fontsize = colorbar_fontsize[view-1]
  # --- draw the bar
  if colbarlinear and levs is not None:
    # --- Use the contour plotting routine plfc for this case. The data
    # --- plotted is uniformly spaced between zmin and zmax. The contour
    # --- levels are those specified. The result is that the colorbar
    # --- shows the contours levels by their values relative to zmin and zmax.
    plotval = span(zmin,zmax,255)[:,newaxis]*ones(2)
    xx = array([xmin,xmax])*ones(255)[:,newaxis]
    yy = span(ymin,ymax,255)[:,newaxis]*ones(2)
    # --- ireg must be of type integer because some legacy code used
    # --- expects it.
    ireg = ones((255,2),'i')
    plfc(plotval,yy,xx,ireg,contours=array(levs),local=1)
  else:
    # --- Use cell array plotting for this case. All of the colors get a block
    # --- of the same size. If levs is not specified, the uniform spacing
    # --- matches the uniform spacing of the contours. If levs is specified,
    # --- each equal sized block represents one contour level, independent of
    # --- the range of the level relative to other levels.
    if type(zmin) == type(zmax) == type(1) and \
       zmin >= 0 and zmax <=199:
       plotval = arange(zmin,zmax+1,typecode=ubyte)[:,newaxis]*ones(2)
    else:
       plotval = (arange(ncolor)/(ncolor-1.))[:,newaxis]*ones(2)
    pli(plotval,xmin,ymin,xmax,ymax,top=ctop,local=1)
  # --- Draw a black box around it
  pldj([xmin,xmin,xmin,xmax],[ymin,ymax,ymin,ymin],
       [xmax,xmax,xmin,xmax],[ymin,ymax,ymax,ymax],local=1)
  # --- Generate nice levels for the labels and tick marks.
  if levs is None:
    # --- Use the nicelevels routine to get evenly spaced labels.
    nicelevs = nicelevels(array([zmin,zmax]))
  else:
    # --- If there are less than 15 specified contour levels, put a label
    # --- at each of the labels. If there are more, pick out roughly 15
    # --- evenly spaced values. Also, if the levels do not extend beyond
    # --- zmin and zmax, add labels at those points too.
    nicelevs = levs
    if zmin < levs[0]:  nicelevs = array([zmin] + list(nicelevs))
    if zmax > levs[-1]: nicelevs = array(list(nicelevs) + [zmax])
  llev = len(nicelevs)
  # --- Create the labels
  labels = []
  # --- Calculate the location of the labels.
  if not colbarlinear and levs is not None:
    # --- The ys are evenly spaced
    ys = ymin + arange(len(nicelevs))/(len(levs)+1.)*(ymax - ymin)
    # --- If the lowest level is less than zmin, then bump up the y's
    # --- by one block size.
    if levs[0] < zmin: ys = ys + 1./(len(levs)+1.)*(ymax - ymin)
  elif llev==2 and (nicelevs[0] == nicelevs[1]):
    ys = array([ymin,ymax])
  else:
    ys = ymin + (ymax - ymin)*(nicelevs - zmin)/(zmax - zmin)
  # --- Plot the labels, skipping ones that are too close together.
  if uselog == 'e' or uselog == 1.: ss = " e^%.5g"
  elif uselog is not None:          ss = " %d^%%.5g"%int(uselog)
  else:                             ss = " %.5g"
  ylast = 0.
  for i in xrange(llev):
    if ys[i] - ylast > (ymax-ymin)/30:
      plt(ss%nicelevs[i],xmax+0.005,ys[i]-0.005,height=fontsize,local=1)
      ylast = ys[i]
  # --- Plot the tick marks
  pldj(llev*[xmax],ys,llev*[xmax+0.005],ys,local=1)
  # --- Return to plot system to its original value.
  plsys(view)

#############################################################################
#############################################################################
def writepalette(filename,r,g,b,comments=None):
  """
Write a palette to the file
  - filename: the file to write the palette to, note that '.gp' is appended
  - r,g,b: lists of colors to write out
           They must be integers between 0 and 255. No checking is done!
  - comments=None: optional comments to write to the file.
  """
  ff = open(filename+'.gp','w')
  ff.write('# gist palette '+filename+'.gp\n')
  if comments is not None: ff.write('# '+comments+'\n')
  ff.write('ncolors = %d\n'%len(r))
  for ri,gi,bi in zip(r,g,b):
    ff.write('%8d%8d%8d\n'%(ri,gi,bi))
  ff.close()

#############################################################################
def changepalette(returnpalette=0,filename='newpalette',help=0,view=None):
  """
Dynamically change the color palette.
  - returnpalette=0: when true, returns tuple of (red, green, blue)
  - filename='newpalette': the palette will be written to the file
                           when requested
  - help=0: when true, prints this message
  """
  print """
Mouse actions:
  Button 1: shifts a point, compressing and stretching the rest of the colors
  Button 2: reset palette to original
  Button 3: shifts a point, sliding the colors up and down
  Control Button 1: add black point
  Control Button 3: add white point
  Shift Button 1: reverse the palette
  Shift Button 2: writes the palette to the file, defaults to newpalette.gp
  Shift Button 3: quits
  """
  # --- Print out help if wanted
  if help: print changepalette.__doc__
  # --- min's and max's are the same as in the colorbar routine
  if view is None: view = plsys()
  xmin,xmax,ymin,ymax = colorbar_placement[view-1]
  # --- Create storate arrays
  # --- rr, gg, bb hold the original palette
  rr = zeros(200,ubyte)
  gg = zeros(200,ubyte)
  bb = zeros(200,ubyte)
  palette(rr,gg,bb,query=1)
  # --- newrr, newgg, newbb hold the new palette
  newrr = zeros(200,ubyte)
  newgg = zeros(200,ubyte)
  newbb = zeros(200,ubyte)
  # --- position relative to the original palette
  cc = arange(0,200)*1.
  newcc = arange(0,200)*1.
  # --- List of black and white points
  blacklist = []
  whitelist = []
  while 1:
    mm = mouse(0,0,"")
    if mm == None or (mm[9] == 3 and mm[10] == 1): break
    # --- Get mouse positions. Skip if outside the colorbar
    (x1, y1, x2, y2) = tuple(mm[:4])
    if x1 < xmin or x1 > xmax or x2 < xmin or x2 > xmax: continue
    if y1 < ymin or y1 > ymax or y2 < ymin or y2 > ymax: continue

    if mm[9] == 1 and mm[10] == 0:
      # --- Button 1, no keys
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      i2 = nint((y2 - ymin)/(ymax - ymin)*200)
      up = (ymax - y1)/(ymax - y2)
      down = (y1 - ymin)/(y2 - ymin)
      for i in xrange(1,i2):
        iold = int(i*down)
        wold =     i*down - iold
        newcc[i] = cc[iold]*(1.-wold) + cc[iold+1]*wold
      for i in xrange(i2,199):
        iold = 199 - int((199-i)*up)
        wold = iold - (199 -    ((199-i)*up))
        newcc[i] = cc[iold]*(1.-wold) + cc[iold-1]*wold

    if mm[9] == 2 and mm[10] == 0:
      # --- Button 2, no keys
      # --- Restore original palette
      newcc = arange(0,200)*1.
      blacklist = []
      whitelist = []

    if mm[9] == 3 and mm[10] == 0:
      # --- Button 3, no keys
      # --- slide whole palette
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      i2 = nint((y2 - ymin)/(ymax - ymin)*200)
      for i in xrange(0,200):
        iold = i - (i2 - i1)
        if iold < 0: newcc[i] = cc[0]
        elif iold > 199: newcc[i] = cc[-1]
        else: newcc[i] = cc[iold]

    if mm[9] == 1 and mm[10] == 1:
      # --- Button 1, shift
      # --- Reverse the palette
      newcc[:] = cc[::-1]

    if mm[9] == 2 and mm[10] == 1:
      # --- Button 2, shift
      print 'Writing palette to '+filename+'.gp'
      writepalette(filename,newrr,newgg,newbb)

    if mm[9] == 1 and mm[10] == 4:
      # --- button 1, control
      # --- Add black point
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      blacklist.append(i1)

    if mm[9] == 3 and mm[10] == 4:
      # --- button 3, control
      # --- Add white point
      i1 = nint((y1 - ymin)/(ymax - ymin)*200)
      whitelist.append(i1)

    # --- Calculate the new palette based on the position relative to the
    # --- original palette.
    for i in xrange(0,200):
      ii = int(newcc[i])
      ww =     newcc[i]  - ii
      iip1 = min(ii+1,199)
      newrr[i] = (nint(rr[ii]*(1.-ww) + rr[iip1]*ww))
      newgg[i] = (nint(gg[ii]*(1.-ww) + gg[iip1]*ww))
      newbb[i] = (nint(bb[ii]*(1.-ww) + bb[iip1]*ww))
    for ii in blacklist: (newrr[ii], newgg[ii], newbb[ii]) = 0,0,0
    for ii in whitelist: (newrr[ii], newgg[ii], newbb[ii]) = 255,255,255
    cc[:] = newcc
    palette(newrr,newgg,newbb)

  if returnpalette: return (newrr,newgg,newbb)

#############################################################################
def makepalette(filename,points,comments=None,ncolor=200):
  """
Creates a palette. A list of rgb points is input and the palette created
varies between the points.
  - filename: the file to write the palette to, note that '.gp' is appended
  - points: list of points, each point must be a list of the rgb values and
            the number of steps to the next point. An optional 5th number
            can be provided which is the exponent of the variation to the next
            point. It defaults to 1, which mean linear.
  - comments: an optional string of comments that is written to the file.
An example:
makepalette('blueorange',[[1,1,1,1],[0,0,.6,79],[0,.6,.6,80],[1,1,0,79],
                          [1,.5,0]])
This makes a palette that varies from blue to orange with points at cyan and
yellow in the middle. Also, note that the vary lowest is point is forced to
white by the first point, [1,1,1,1].
  """
  assert len(points) > 1,'Must specify at least 2 points'
  icolor = 0
  r,g,b = [],[],[]
  for i in range(len(points)-1):
    if len(points[i]) > 3: nc = points[i][3]
    else:                  nc = ncolor - icolor
    if len(points[i]) > 4: power = points[i][4]
    else:                  power = 1.
    s = span(0.,1.,nc+1)**power
    r += list(nint((points[i][0] + (points[i+1][0] - points[i][0])*s[:-1])*255))
    g += list(nint((points[i][1] + (points[i+1][1] - points[i][1])*s[:-1])*255))
    b += list(nint((points[i][2] + (points[i+1][2] - points[i][2])*s[:-1])*255))
    icolor = icolor + nc - 1
  r += [nint(points[-1][0]*255)]
  g += [nint(points[-1][1]*255)]
  b += [nint(points[-1][2]*255)]
  assert len(r) <= 200,'There can be at most 200 colors'
  writepalette(filename,r,g,b,comments)

#############################################################################
#############################################################################
def viewsurface(scale=4.,gnomon=1):
  """
Dynamically view a surface plot. The mouse is used to change to view angle.
With button 1 pushed, the horizontal movement changes the z angle, and
vertical the y angle. With button 2 pressed, horizontal changes the x angle.
When finished, press return in the python window.
  - scale=4.: multiplicative factor to convert mouse movement to angle change
Returns the final values of the parameters that can be passed to pl3d.rot3
to reproduce the same orientation.
  """
  from gist import pl3d
  pl3d.gnomon(gnomon)
  [xmin3min,xmax3max,ymin3min,ymax3max,sys] = limits()
  while 1:
    mm = mouse(0,0,"")
    if mm == None: break
    (xa, ya, za) = (0.,0.,0.)
    if mm[9] == 1:
      ya = - (mm[3] - mm[1])*scale
      za = - (mm[2] - mm[0])*scale
    if mm[9] == 3:
      xa = (mm[2] - mm[0])*scale
    pl3d.rot3(xa,ya,za)
    [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
    xmin3min = min(xmin3min,xmin3)
    xmax3max = max(xmax3max,xmax3)
    ymin3min = min(ymin3min,ymin3)
    ymax3max = max(ymax3max,ymax3)
    limits(xmin3min,xmax3max,ymin3min,ymax3max)
  pl3d.gnomon(gnomon)
  print xa,ya,za

def _viewsurfacetest(scale=4.,gnomon=1):
  """
Dynamically view a surface plot. The mouse is used to change to view angle.
With button 1 pushed, the horizontal movement changes the z angle, and
vertical the y angle. With button 2 pressed, horizontal changes the x angle.
When finished, press return in the python window.
  - scale=4.: multiplicative factor to convert mouse movement to angle change
  """
  from gist import pl3d
  pl3d.gnomon(gnomon)
  pl3d.orient3(phi=0.,theta=0.)
  [xmin3min,xmax3max,ymin3min,ymax3max] = pl3d.draw3(1)
  phi = 0.
  theta = 0.
  (xa, ya, za) = (0.,0.,0.)
  while 1:
    mm = mouse(0,0,"")
    if mm == None: break
    dphi   = (mm[3] - mm[1])*scale
    dtheta = (mm[2] - mm[0])*scale
    print theta,phi
    newxa = xa + dtheta*sin(phi)*cos(theta) + dphi*cos(phi)*cos(theta)
    newya = ya + dtheta*sin(phi)*sin(theta) + dphi*cos(phi)*sin(theta)
    newza = za + dtheta*cos(phi)*cos(theta) + dphi*sin(phi)*sin(theta)
    phi = xa*cos(za) + ya*sin(za)
    theta = za
    pl3d.rot3(newxa-xa,newya-ya,newza-za)
    xa = newxa
    ya = newya
    za = newza
    [xmin3,xmax3,ymin3,ymax3] = pl3d.draw3(1)
    xmin3min = min(xmin3min,xmin3)
    xmax3max = max(xmax3max,xmax3)
    ymin3min = min(ymin3min,ymin3)
    ymax3max = max(ymax3max,ymax3)
    limits(xmin3min,xmax3max,ymin3min,ymax3max)
  pl3d.gnomon(gnomon)

#############################################################################
def zoombox():
  """
When called, use the mouse left button (press and hold) to draw a
box around the area to be zoomed to.
  """
  m1 = mouse(1,1,'')
  xmin = min(m1[0],m1[2])
  xmax = max(m1[0],m1[2])
  ymin = min(m1[1],m1[3])
  ymax = max(m1[1],m1[3])
  limits(xmin,xmax,ymin,ymax)

#############################################################################
#############################################################################
def ppmultispecies(pp,args,kw):
  """checks if js defined and assign it to a list if plotting multispecies.
  Also assign colors accordingly
  """
  if 'js' in kw:
    js = kw['js']
    if js != -1 and type(js) != ListType:
      return false
    else:
      if js == -1: js = range(top.ns)
      ncolor = kw.get('ncolor',240)
      color = kw.get('color',range(0,ncolor,ncolor/len(js)))
      for i in xrange(len(js)):
        kw['js'] = js[i]
        kw['color'] = color[i]
        apply(pp,args,kw)
      return true
  else:
    return false

########################################################################
########################################################################
########################################################################
########################################################################
def checkparticleplotarguments(kw):
  """Convenience routine to check arguments of particle plot routines.
Warning: this has the side affect of adding the arguement allowbadargs to
the kw dictionary. This is done since the calls to these functions here to
make the plots may have unused arguements since the entire kw list passed
into each of the pp plotting routines is passed into each of these
functions.
  """
  badargs = selectparticles(checkargs=1,kwdict=kw)
  badargs = pptitleright(checkargs=1,kwdict=badargs)
  badargs = ppgeneric(checkargs=1,kwdict=badargs)
  badargs = getxxpslope(checkargs=1,kwdict=badargs)
  kw['allowbadargs'] = 1
  if badargs: raise TypeError,"bad arguments%s"%' '.join(badargs.keys())
########################################################################
def ppzxy(iw=0,**kw):
  "Plots Z-X and Z-Y in single page. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzxy,(iw,),kw): return
  kw['view'] = 9
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,win=top.ywindows,z=top.pgroup.yp,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("X vs Z","Z","X",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(getx(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),kwdict=kw)

  kw['view'] = 10
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.yplmin,top.yplmax)
  ii = selectparticles(iw=iw,win=top.xwindows,z=top.pgroup.xp,kwdict=kw)
  settitles("Y vs Z","Z","Y",pptitleright(iw=iw,kwdict=kw))
  ppgeneric(gety(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),kwdict=kw)

##########################################################################
def ppzx(iw=0,**kw):
  "Plots Z-X. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzx,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("X vs Z","Z","X",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getx(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzy(iw=0,**kw):
  "Plots Z-Y. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzy,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.yplmin,top.yplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Y vs Z","Z","Y",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(gety(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzr(iw=0,**kw):
  "Plots Z-R. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzr,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xplmin,top.xplmax)
  kw.setdefault('local',0)
  kw.setdefault('ldensitycylindrical',True)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("R vs Z","Z","R",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getr(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzxp(iw=0,**kw):
  "Plots Z-X'. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzxp,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin,top.xpplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("X' vs Z","Z","X'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getxp(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzvx(iw=0,**kw):
  "Plots Z-Vx. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvx,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vx vs Z","Z","Vx",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvx(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzyp(iw=0,**kw):
  "Plots Z-Y'. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzyp,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.ypplmin,top.ypplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Y' vs Z","Z","Y'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getyp(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzvy(iw=0,**kw):
  "Plots Z-Vy. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvy,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.ypplmin*top.vbeam,top.ypplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vy vs Z","Z","Vy",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvy(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzvz(iw=0,**kw):
  "Plots Z-Vz. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange(kwdict=kw)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,vzmin,vzmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vz vs Z","Z","Vz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvz(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzuz(iw=0,**kw):
  "Plots Z-Uz. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzuz,(iw,),kw): return
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Uz vs Z","Z","Uz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getuz(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzvr(iw=0,**kw):
  "Plots Z-Vr. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvr,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vr vs Z","Z","Vr",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvr(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzvtheta(iw=0,**kw):
  "Plots Z-Vtheta. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvtheta,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vtheta vs Z","Z","Vtheta",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvtheta(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppzvperp(iw=0,**kw):
  "Plots Z-Vperp (sqrt(Vx**2 + Vy**2)). For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvperp,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    vperpmin = min(top.xpplmin*top.vbeam,top.ypplmin*top.vbeam)
    vperpmax = min(top.xpplmax*top.vbeam,top.ypplmax*top.vbeam)
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,vperpmin,vperpmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vperp vs Z","Z","Vperp",pptitleright(iw=iw,kwdict=kw))
  vx = getvx(ii=ii,gather=0,**kw)
  vy = getvy(ii=ii,gather=0,**kw)
  vperp = sqrt(vx**2 + vy**2)
  return ppgeneric(vperp,getz(ii=ii,gather=0,**kw),kwdict=kw)

##########################################################################
def ppzrp(iw=0,**kw):
  "Plots Z-R'. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzrp,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin,top.xpplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("R' vs Z","Z","R'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getrp(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppxex(iw=0,**kw):
  "Plots X-Ex. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxex,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ex vs X","X","Ex",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getex(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppxey(iw=0,**kw):
  "Plots X-Ey. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxey,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ey vs X","X","Ey",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getey(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppxez(iw=0,**kw):
  "Plots X-Ez. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxez,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ez vs X","X","Ez",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getez(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppxbx(iw=0,**kw):
  "Plots X-Bx. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxex,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Bx vs X","X","Bx",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getbx(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppxby(iw=0,**kw):
  "Plots X-By. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxby,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("By vs X","X","By",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getby(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppxbz(iw=0,**kw):
  "Plots X-Bz. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxez,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Bz vs X","X","Bz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getbz(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppyex(iw=0,**kw):
  "Plots Y-Ex. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyex,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ex vs Y","Y","Ex",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getex(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppyey(iw=0,**kw):
  "Plots Y-Ey. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyey,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ey vs Y","Y","Ey",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getey(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppyez(iw=0,**kw):
  "Plots Y-Ez. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyez,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ez vs Y","Y","Ez",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getez(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppybx(iw=0,**kw):
  "Plots Y-Bx. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyex,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Bx vs Y","Y","Bx",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getbx(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppyby(iw=0,**kw):
  "Plots Y-By. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyby,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("By vs Y","Y","By",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getby(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppybz(iw=0,**kw):
  "Plots Y-Bz. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyez,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Bz vs Y","Y","Bz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getbz(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppzex(iw=0,**kw):
  "Plots Z-Ex. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzex,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ex vs Z","Z","Ex",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getex(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppzey(iw=0,**kw):
  "Plots Z-Ey. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzey,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ey vs Z","Z","Ey",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getey(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppzez(iw=0,**kw):
  "Plots Z-Ez. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzez,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Ez vs Z","Z","Ez",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getez(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppzbx(iw=0,**kw):
  "Plots Z-Bx. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzbx,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Bx vs Z","Z","Bx",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getbx(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppzby(iw=0,**kw):
  "Plots Z-By. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzby,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("By vs Z","Z","By",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getby(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppzbz(iw=0,**kw):
  "Plots Z-Bz. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzbz,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.zplmin+top.zbeam,top.zplmax+top.zbeam,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Bz vs Z","Z","Bz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getbz(ii=ii,gather=0,**kw),getz(ii=ii,gather=0,**kw),
                   kwdict=kw)
##########################################################################
def ppexey(iw=0,**kw):
  "Plots Ex-Ey. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppexey,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xpplmin*top.vbeam,top.xpplmax*top.vbeam,
                      top.ypplmin*top.vbeam,top.ypplmax*top.vbeam)
  kw.setdefault('local',0)
  settitles("Ey vs Ex","Ex","Ey",pptitleright(iw=iw,kwdict=kw))
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  return ppgeneric(getey(ii=ii,gather=0,**kw),getex(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppxy(iw=0,**kw):
  "Plots X-Y. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxy,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,top.yplmin,top.yplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Y vs X","X","Y",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(gety(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppxxp(iw=0,**kw):
  "Plots X-X'. If slope='auto', it is calculated from the moments. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxxp,(iw,),kw): return
  if type(kw.get('slope',0.)) == type(''):
    (slope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    kw['slope'] = slope
    kw['yoffset'] = xpoffset
    kw['xoffset'] = xoffset
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("X' vs X","X","X'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getxp(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppyyp(iw=0,**kw):
  "Plots Y-Y'. If slope='auto', it is calculated from the moments. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyyp,(iw,),kw): return
  if type(kw.get('slope',0.)) == type(''):
    (slope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    kw['slope'] = slope
    kw['yoffset'] = ypoffset
    kw['xoffset'] = yoffset
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Y' vs Y","Y","Y'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getyp(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppxpyp(iw=0,**kw):
  "Plots X'-Y'. If slope='auto', it is calculated from the moments. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxpyp,(iw,),kw): return
  slope = kw.get('slope',0.)
  if type(slope) == type(''):
    (xslope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    (yslope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    kw['slope'] = 0.
  else:
    (xslope,xoffset,xpoffset) = (slope,0.,0.)
    (yslope,yoffset,ypoffset) = (slope,0.,0.)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xpplmin,top.xpplmax,top.ypplmin,top.ypplmax)
  kw.setdefault('local',0)
  settitles("Y' vs X'","X'","Y'",pptitleright(iw=iw,kwdict=kw))
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  x = getx(ii=ii,gather=0,**kw)
  y = gety(ii=ii,gather=0,**kw)
  xp = getxp(ii=ii,gather=0,**kw)
  yp = getyp(ii=ii,gather=0,**kw)
  xpms = (xp - xslope*(x-xoffset) - xpoffset)
  ypms = (yp - yslope*(y-yoffset) - ypoffset)
  return ppgeneric(ypms,xpms,kwdict=kw)

##########################################################################
def ppxvx(iw=0,**kw):
  "Plots X-Vx. If slope='auto', it is calculated from the moments. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxvx,(iw,),kw): return
  if type(kw.get('slope',0.)) == type(''):
    (slope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    kw['slope'] = slope*vz
    kw['yoffset'] = xpoffset*vz
    kw['xoffset'] = xoffset
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,
                      top.xpplmin*top.vbeam,top.xpplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vx vs X","X","Vx",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvx(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppyvy(iw=0,**kw):
  "Plots Y-Vy. If slope='auto', it is calculated from the moments. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyvy,(iw,),kw): return
  if type(kw.get('slope',0.)) == type(''):
    (slope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    kw['slope'] = slope*vz
    kw['yoffset'] = ypoffset*vz
    kw['xoffset'] = yoffset
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,
                      top.ypplmin*top.vbeam,top.ypplmax*top.vbeam)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vy vs Y","Y","Vy",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvy(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppxvz(iw=0,**kw):
  "Plots X-Vz. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppxvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange(kwdict=kw)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xplmin,top.xplmax,vzmin,vzmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vz vs X","X","Vz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvz(ii=ii,gather=0,**kw),getx(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppyvz(iw=0,**kw):
  "Plots Y-Vz. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppyvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange(kwdict=kw)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.yplmin,top.yplmax,vzmin,vzmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vz vs Y","Y","Vz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvz(ii=ii,gather=0,**kw),gety(ii=ii,gather=0,**kw),
                   kwdict=kw)

##########################################################################
def ppvxvy(iw=0,**kw):
  "Plots Vx-Vy. If slope='auto', it is calculated from the moments. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppvxvy,(iw,),kw): return
  slope = kw.get('slope',0.)
  kw['slope'] = 0.
  if type(slope) == type(''):
    (xslope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    (yslope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    vxslope = xslope*vz
    vxoffset = xpoffset*vz
    vyslope = yslope*vz
    vyoffset = ypoffset*vz
  else:
    (vxslope,xoffset,vxoffset) = (slope,0.,0.)
    (vyslope,yoffset,vyoffset) = (slope,0.,0.)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xpplmin*top.vbeam,top.xpplmax*top.vbeam,
                      top.ypplmin*top.vbeam,top.ypplmax*top.vbeam)
  kw.setdefault('local',0)
  settitles("Vy vs Vx","Vx","Vy",pptitleright(iw=iw,kwdict=kw))
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  x = getx(ii=ii,gather=0,**kw)
  y = gety(ii=ii,gather=0,**kw)
  vx = getvx(ii=ii,gather=0,**kw)
  vy = getvy(ii=ii,gather=0,**kw)
  vxms = (vx - vxslope*(x-xoffset) - vxoffset)
  vyms = (vy - vyslope*(y-yoffset) - vyoffset)
  return ppgeneric(vyms,vxms,kwdict=kw)

##########################################################################
def ppvxvz(iw=0,**kw):
  "Plots Vx-Vz. If slope='auto', it is calculated from the moments. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppvxvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange(kwdict=kw)
  slope = kw.get('slope',0.)
  kw['slope'] = 0.
  if type(slope) == type(''):
    (xslope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    vxslope = xslope*vz
    vxoffset = xpoffset*vz
  else:
    (vxslope,xoffset,vxoffset) = (slope,0.,0.)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.xpplmin*top.vbeam,top.xpplmax*top.vbeam,vzmin,vzmax)
  kw.setdefault('local',0)
  settitles("Vz vs Vx","Vx","Vz",pptitleright(iw=iw,kwdict=kw))
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  x = getx(ii=ii,gather=0,**kw)
  vx = getvx(ii=ii,gather=0,**kw)
  vxms = (vx - vxslope*(x-xoffset) - vxoffset)
  vz = getvz(ii=ii,gather=0,**kw)
  return ppgeneric(vz,vxms,kwdict=kw)

##########################################################################
def ppvyvz(iw=0,**kw):
  "Plots Vy-Vz. If slope='auto', it is calculated from the moments. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppvyvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange(kwdict=kw)
  slope = kw.get('slope',0.)
  kw['slope'] = 0.
  if type(slope) == type(''):
    (yslope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=kw.get('iz'),kwdict=kw)
    vyslope = yslope*vz
    vyoffset = ypoffset*vz
  else:
    (vyslope,yoffset,vyoffset) = (slope,0.,0.)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (top.ypplmin*top.vbeam,top.ypplmax*top.vbeam,vzmin,vzmax)
  kw.setdefault('local',0)
  settitles("Vz vs Vy","Vy","Vz",pptitleright(iw=iw,kwdict=kw))
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  y = gety(ii=ii,gather=0,**kw)
  vy = getvy(ii=ii,gather=0,**kw)
  vyms = (vy - vyslope*(y-yoffset) - vyoffset)
  vz = getvz(ii=ii,gather=0,**kw)
  return ppgeneric(vz,vyms,kwdict=kw)

##########################################################################
def ppvzvperp(iw=0,**kw):
  "Plots Vz-Vperp (sqrt(Vx**2 + Vy**2)). For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(ppzvperp,(iw,),kw): return
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    vperpmin = min(top.xpplmin*top.vbeam,top.ypplmin*top.vbeam)
    vperpmax = min(top.xpplmax*top.vbeam,top.ypplmax*top.vbeam)
    (vzmin,vzmax) = getvzrange(kwdict=kw)
    kw['pplimits'] = (vzmin,vzmax,vperpmin,vperpmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vperp vs Vz","Vz","Vperp",pptitleright(iw=iw,kwdict=kw))
  vx = getvx(ii=ii,gather=0,**kw)
  vy = getvy(ii=ii,gather=0,**kw)
  vperp = sqrt(vx**2 + vy**2)
  return ppgeneric(vperp,getvz(ii=ii,gather=0,**kw),kwdict=kw)

##########################################################################
def pprrp(iw=0,scale=0,slopejs=-1,**kw):
  """Plots R-R', If slope='auto', it is calculated from the moments.
  - scale=0: when true, scale particle by 2*rms
  - slopejs=-1: Species whose moments are used to calculate the slope
                 -1 means use data combined from all species.
For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`.
  """
  checkparticleplotarguments(kw)
  if ppmultispecies(pprrp,(iw,scale,slopejs),kw): return
  xscale = 1.
  yscale = 1.
  xpscale = 1.
  ypscale = 1.
  if scale:
    iiw = max(0,iw)
    xscale = 2.*top.xrms[iiw,slopejs]
    yscale = 2.*top.yrms[iiw,slopejs]
    xpscale = 2.*top.vxrms[iiw,slopejs]/top.vzbar[iiw,slopejs]
    ypscale = 2.*top.vyrms[iiw,slopejs]/top.vzbar[iiw,slopejs]
  ii = selectparticles(iw=iw,kwdict=kw)
  xx = getx(ii=ii,gather=0,**kw)/xscale
  yy = gety(ii=ii,gather=0,**kw)/yscale
  xp = getxp(ii=ii,gather=0,**kw)/xpscale
  yp = getyp(ii=ii,gather=0,**kw)/ypscale
  rr = sqrt(xx**2 + yy**2)
  tt = arctan2(yy,xx)
  rp = xp*cos(tt) + yp*sin(tt)
  slope = kw.get('slope',0.)
  if type(slope) == type(''):
    aversq = globalave(rr**2)
    averrp = globalave(rr*rp)
    if aversq > 0.:
      slope = averrp/aversq
    else:
      slope = 0.
    kw['slope'] = slope
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                      top.xpplmin/xpscale,top.xpplmax/ypscale)
  kw.setdefault('local',0)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("R' vs R","R","R'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(rp,rr,kwdict=kw)

##########################################################################
def pprtp(iw=0,scale=0,slopejs=-1,**kw):
  """Plots R-Theta', If slope='auto', it is calculated from the moments.
  - scale=0: when true, scale particle by 2*rms
  - slopejs=-1: Species whose moments are used to calculate the slope
                 -1 means use data combined from all species.
For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`.
  """
  checkparticleplotarguments(kw)
  if ppmultispecies(pprtp,(iw,scale,slopejs),kw): return
  xscale = 1.
  yscale = 1.
  xpscale = 1.
  ypscale = 1.
  if scale:
    iiw = max(0,iw)
    xscale = 2.*top.xrms[iiw,slopejs]
    yscale = 2.*top.yrms[iiw,slopejs]
    xpscale = 2.*top.vxrms[iiw,slopejs]/top.vzbar[iiw,slopejs]
    ypscale = 2.*top.vyrms[iiw,slopejs]/top.vzbar[iiw,slopejs]
  ii = selectparticles(iw=iw,kwdict=kw)
  xx = getx(ii=ii,gather=0,**kw)/xscale
  yy = gety(ii=ii,gather=0,**kw)/yscale
  xp = getxp(ii=ii,gather=0,**kw)/xpscale
  yp = getyp(ii=ii,gather=0,**kw)/ypscale
  rr = sqrt(xx**2 + yy**2)
  tt = arctan2(yy,xx)
  tp = -xp*sin(tt) + yp*cos(tt)
  slope = kw.get('slope',0.)
  if type(slope) == type(''):
    aversq = globalave(rr**2)
    avertp = globalave(rr*tp)
    if aversq > 0.:
      slope = avertp/aversq
    else:
      slope = 0.
    kw['slope'] = slope
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                      top.xpplmin/xpscale,top.xpplmax/ypscale)
  kw.setdefault('local',0)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Theta' vs R","R","Theta'",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(tp,rr,kwdict=kw)

##########################################################################
def pprvr(iw=0,scale=0,slopejs=-1,**kw):
  """Plots R-Vr, If slope='auto', it is calculated from the moments.
  - scale=0: when true, scale particle by 2*rms
  - slopejs=-1: Species whose moments are used to calculate the slope
                 -1 means use data combined from all species.
For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`.
  """
  checkparticleplotarguments(kw)
  if ppmultispecies(pprvr,(iw,scale,slopejs),kw): return
  xscale = 1.
  yscale = 1.
  vxscale = 1.
  vyscale = 1.
  if scale:
    iiw = max(0,iw)
    xscale = 2.*top.xrms[iiw,slopejs]
    yscale = 2.*top.yrms[iiw,slopejs]
    vxscale = 2.*top.vxrms[iiw,slopejs]
    vyscale = 2.*top.vyrms[iiw,slopejs]
  ii = selectparticles(iw=iw,kwdict=kw)
  xx = getx(ii=ii,gather=0,**kw)/xscale
  yy = gety(ii=ii,gather=0,**kw)/yscale
  vx = getvx(ii=ii,gather=0,**kw)/vxscale
  vy = getvy(ii=ii,gather=0,**kw)/vyscale
  rr = sqrt(xx**2 + yy**2)
  tt = arctan2(yy,xx)
  vr = vx*cos(tt) + vy*sin(tt)
  slope = kw.get('slope',0.)
  if type(slope) == type(''):
    aversq = globalave(rr**2)
    avervr = globalave(rr*vr)
    if aversq > 0.:
      slope = avervr/aversq
    else:
      slope = 0.
    kw['slope'] = slope
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (0.,max(top.xplmax/xscale,top.yplmax/yscale),
                      top.xpplmin*top.vbeam/vxscale,
                      top.xpplmax*top.vbeam/vyscale)
  kw.setdefault('local',0)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vr vs R","R","Vr",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(vr,rr,kwdict=kw)

##########################################################################
def pprvz(iw=0,**kw):
  "Plots R-Vz. For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`."
  checkparticleplotarguments(kw)
  if ppmultispecies(pprvz,(iw,),kw): return
  (vzmin,vzmax) = getvzrange(kwdict=kw)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (0.,max(top.xplmax,top.yplmax),vzmin,vzmax)
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vz vs R","R","Vz",pptitleright(iw=iw,kwdict=kw))
  return ppgeneric(getvz(ii=ii,gather=0,**kw),getr(ii=ii,gather=0,**kw),kwdict=kw)

##########################################################################
def pptrace(iw=0,normalize=0,**kw):
  """
Plots X-Y, X-X', Y-Y', Y'-X' in single page
If slope='auto', it is calculated from the moments.
pplimits can be a list of up to four tuples, one for each phase space plot.
If any of the tuples are empty, the limits used will be the usual ones for
that plot.
For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`.
  """
  checkparticleplotarguments(kw)
  if ppmultispecies(pptrace,(iw,normalize),kw): return
  kw.setdefault('local',0)
  ii = selectparticles(iw=iw,kwdict=kw)
  x = getx(ii=ii,gather=0,**kw)
  y = gety(ii=ii,gather=0,**kw)
  xp = getxp(ii=ii,gather=0,**kw)
  yp = getyp(ii=ii,gather=0,**kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  slope = kw.get('slope',0.)
  if type(slope)==type(''):
    del kw['slope']
    iz = kw.get('iz',None)
    (xxpslope,xoffset,xpoffset,vz) = getxxpslope(iw=iw,iz=iz,kwdict=kw)
    (yypslope,yoffset,ypoffset,vz) = getyypslope(iw=iw,iz=iz,kwdict=kw)
    xp = xp - xxpslope*(x - xoffset) - xpoffset
    yp = yp - yypslope*(y - yoffset) - ypoffset
  if kw.get('titles',1):
    titler=pptitleright(iw=iw,kwdict=kw)
    ptitles(titler=titler)
  defaultpplimits = [(top.xplmin,top.xplmax,top.yplmin,top.yplmax),
                     (top.yplmin,top.yplmax,top.ypplmin,top.ypplmax),
                     (top.xplmin,top.xplmax,top.xpplmin,top.xpplmax),
                     (top.ypplmin,top.ypplmax,top.xpplmin,top.xpplmax)]
  pplimits = kw.get('pplimits',None)
  if pplimits is None:
    pplimits = defaultpplimits
  else:
    kw['lframe'] = 1
    if type(pplimits[0]) != type(()):
      pplimits = 4*[pplimits]
    else:
      for i in xrange(4):
        if i == len(pplimits): pplimits.append(defaultpplimits[i])
        if not pplimits[i]: pplimits[i] = defaultpplimits[i]

  # --- First make plot with returngrid set to true. If no grids are used,
  # --- then ppgeneric returns None (and makes the plot), otherwise it
  # --- returns the grid. Also, if normalize is false, then set so plots
  # --- are made (and grid is not returned).
  if normalize: rg = 1
  else:         rg = 0
  kw['view'] = 3
  kw['pplimits'] = pplimits[0]
  settitles("Y vs X","X","Y")
  gxy = ppgeneric(y,x,returngrid=rg,kwdict=kw)

  kw['view'] = 4
  kw['pplimits'] = pplimits[1]
  settitles("Y' vs Y","Y","Y'")
  gyyp = ppgeneric(yp,y,returngrid=rg,kwdict=kw)

  kw['view'] = 5
  kw['pplimits'] = pplimits[2]
  settitles("X' vs X","X","X'")
  gxxp = ppgeneric(xp,x,returngrid=rg,kwdict=kw)

  kw['view'] = 6
  kw['pplimits'] = pplimits[3]
  settitles("X' vs Y'","Y'","X'")
  gxpyp = ppgeneric(xp,yp,returngrid=rg,kwdict=kw)

  # --- If the return value is None, then return since plots have already been
  # --- made.
  if gxy is None: return

  # --- If the return value is not None, then call ppgeneric again to
  # --- actually make the plots with the appropriate cmin and cmax
  cmin = kw.get('cmin',None)
  cmax = kw.get('cmax',None)
  if kw.get('cmin',None) is None:
    kw['cmin']=min(minnd(gxy[0]),minnd(gyyp[0]),minnd(gxxp[0]),minnd(gxpyp[0]))
  if kw.get('cmax',None) is None:
    kw['cmax']=max(maxnd(gxy[0]),maxnd(gyyp[0]),maxnd(gxxp[0]),maxnd(gxpyp[0]))

  kw['view'] = 3
  kw['pplimits'] = pplimits[0]
  settitles("Y vs X","X","Y")
  ppgeneric(y,x,kwdict=kw)

  kw['view'] = 4
  kw['pplimits'] = pplimits[1]
  settitles("Y' vs Y","Y","Y'")
  ppgeneric(yp,y,kwdict=kw)

  kw['view'] = 5
  kw['pplimits'] = pplimits[2]
  settitles("X' vs X","X","X'")
  ppgeneric(xp,x,kwdict=kw)

  kw['view'] = 6
  kw['pplimits'] = pplimits[3]
  settitles("X' vs Y'","Y'","X'")
  ppgeneric(xp,yp,kwdict=kw)

##########################################################################
##########################################################################
##########################################################################
def ppzxco(iw=0,ncolor=None,nskipcol=None,nstepcol=None,**kw):
  """Plots Z-X with color based in particle index
 - ncolor=top.ncolor: number of colors to use
 - nskipcol=top.nskipcol:
 - nstepcol=top.nstepcol:
For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`.
  """
  # --- First part copied from ppzx
  if ppmultispecies(ppzxco,(iw,),kw): return
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("X vs Z","Z","X",pptitleright(iw=iw,kwdict=kw))
  x = getx(ii=ii,gather=0,**kw)
  z = getz(ii=ii,gather=0,**kw)
  kw.setdefault('local',0)

  # --- Second part from the original ppzxco
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: nskipcol = top.nskipcol
  if nstepcol is None: nstepcol = top.nstepcol
  inp=1.*len(x)/ncolor
  istep=nskipcol*nstepcol
  istart = 0
  if (inp < istep): istep = 1
  for ij in range(1,istep+1,nskipcol*2):
    for ic in range(1,ncolor+1):
      irs1 = istart+ij+int(inp*(ic-1))
      irs2 = istart+int(inp*ic)
      irs3 = istep
      ii = iota(irs1,irs2,irs3)
      ii = (ii-istart-ij-int(inp*(ic-1)))/istep
      plp(take(x[irs1:irs2:irs3],ii),take(z[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],**kw)
    for ic in range(ncolor,0,-1):
      irs1 = istart+ij+nskipcol+int(inp*(ic-1))
      irs2 = istart+int(inp*ic)
      irs3 = istep
      ii = iota(irs1,irs2,irs3)
      ii = (ii-istart-ij-nskipcol-int(inp*(ic-1)))/istep
      plp(take(x[irs1:irs2:irs3],ii),take(z[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],**kw)

##########################################################################
def ppzyco(iw=0,ncolor=None,nskipcol=None,nstepcol=None,**kw):
  """Plots Z-Y with color based in particle index
 - ncolor=top.ncolor: number of colors to use
 - nskipcol=top.nskipcol:
 - nstepcol=top.nstepcol:
For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`.
  """
  # --- First part copied from ppzy
  if ppmultispecies(ppzyco,(iw,),kw): return
  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Y vs Z","Z","Y",pptitleright(iw=iw,kwdict=kw))
  y = gety(ii=ii,gather=0,**kw)
  z = getz(ii=ii,gather=0,**kw)
  kw.setdefault('local',0)

  # --- Second part from the original ppzyco
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: nskipcol = top.nskipcol
  if nstepcol is None: nstepcol = top.nstepcol
  inp=1.*len(y)/ncolor
  istep=nskipcol*nstepcol
  istart = 0
  if (inp < istep): istep = 1
  for ij in range(1,istep+1,nskipcol*2):
    for ic in range(1,ncolor+1):
      irs1 = istart+ij+int(inp*(ic-1))
      irs2 = istart+int(inp*ic)
      irs3 = istep
      ii = iota(irs1,irs2,irs3)
      ii = (ii-istart-ij-int(inp*(ic-1)))/istep
      plp(take(y[irs1:irs2:irs3],ii),take(z[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],**kw)
    for ic in range(ncolor,0,-1):
      irs1 = istart+ij+nskipcol+int(inp*(ic-1))
      irs2 = istart+int(inp*ic)
      irs3 = istep
      ii = iota(irs1,irs2,irs3)
      ii = (ii-istart-ij-nskipcol-int(inp*(ic-1)))/istep
      plp(take(y[irs1:irs2:irs3],ii),take(z[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],**kw)

##########################################################################
def ppzxyco(iw=0,ncolor=None,nskipcol=None,nstepcol=None,**kw):
  """Plots Z-X and Z-Y in single frame with color based in paricle index
See documentation for ppzxco and ppzyco.
For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`.
  """
  plsys(9)
  ppzxco(ncolor=ncolor,nskipcol=nskipcol,nstepcol=nstepcol,**kw)
  plsys(10)
  ppzyco(ncolor=ncolor,nskipcol=nskipcol,nstepcol=nstepcol,**kw)

##########################################################################
def ppzvzco(iw=0,ncolor=None,nskipcol=None,nstepcol=None,**kw):
  """Plots Z-Vz with color based in particle index
 - ncolor=top.ncolor: number of colors to use
 - nskipcol=top.nskipcol:
 - nstepcol=top.nstepcol:
For particle selection options, see :py:func:`~particles.selectparticles`. For plotting options, see :py:func:`ppgeneric`.
  """
  # --- First part copied from ppzvz
  if ppmultispecies(ppzvzco,(iw,),kw): return
  (vzmin,vzmax) = getvzrange(kwdict=kw)

  ii = selectparticles(iw=iw,kwdict=kw)
  if(top.wpid!=0): kw['weights'] = getpid(id=top.wpid-1,ii=ii,gather=0,**kw)
  settitles("Vz vs Z","Z","Vz",pptitleright(iw=iw,kwdict=kw))
  vz = getvz(ii=ii,gather=0,**kw)
  z = getz(ii=ii,gather=0,**kw)
  kw.setdefault('local',0)

  # --- Second part from the original ppzvzco
  if ncolor is None: ncolor = top.ncolor
  if nskipcol is None: nskipcol = top.nskipcol
  if nstepcol is None: nstepcol = top.nstepcol
  inp=1.*len(vz)/ncolor
  istep=nskipcol*nstepcol
  istart = 0
  if (inp < istep): istep = 1
  for ij in range(1,istep+1,nskipcol*2):
    for ic in range(1,ncolor+1):
      irs1 = int(istart+ij+inp*(ic-1))
      irs2 = int(istart+inp*ic)
      irs3 = istep
      ii = iota(irs1,irs2,irs3)
      ii = (ii-istart-ij-int(inp*(ic-1)))/istep
      plp(take(vz[irs1:irs2:irs3],ii),take(z[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],**kw)
    for ic in range(ncolor,0,-1):
      irs1 = int(istart+ij+nskipcol+inp*(ic-1))
      irs2 = int(istart+inp*ic)
      irs3 = istep
      ii = iota(irs1,irs2,irs3)
      ii = (ii-istart-ij-nskipcol-int(inp*(ic-1)))/istep
      plp(take(vz[irs1:irs2:irs3],ii),take(z[irs1:irs2:irs3],ii),
          color=color[ic%len(color)],**kw)

##########################################################################
def ppco(y,x,z,uz=1.,marker='\1',msize=1.0,zmin=None,zmax=None,
         ncolor=None,usepalette=1,npalette=200,local=1):
  """
Plots y versus x with color based in z
  - y: y coordinate
  - x: x coordinate
  - z: used to calculate the color
  - zmin, zmax: optional bounds on the coloring data
  - ncolor: number of colors to use, defaults to top.ncolor
  - usepalette=1: when true, uses palette, otherwise uses colors in array
                  color
  - npalette=200: number of colors to use in the palette
  """
  # --- Make sure the lengths of the input are the same
  assert (len(y) == len(x) == len(z)),"x, y, and z must all be the same length"

  # --- This routine can be expensive in parallel when there are many
  # --- colors since synchronization is needed for each color.
  # --- So, if there arn't too many particles, transfer everything to PE0
  # --- and let it do the work.
  np = len(y)
  if not local: np = globalsum(np)
  if np < 1000000 and not local:
    local = 1
    y = gatherarray(y)
    x = gatherarray(x)
    z = gatherarray(z)
    if type(uz) == ArrayType: uz = gatherarray(uz)
    if me > 0: return

  # --- Make sure arrays are 1-D
  rx = ravel(x)
  ry = ravel(y)
  rz = ravel(z)

  # --- Find extrema
  if not local:
    if zmin is None: zmin = globalmin(rz)
    if zmax is None: zmax = globalmax(rz)
  else:
    if zmin is None: zmin = rz.min()
    if zmax is None: zmax = rz.max()

  if ncolor is None: ncolor = top.ncolor
  dd = (zmax - zmin)/ncolor
  for ic in xrange(ncolor):
    ii = compress(logical_and(less(zmin+ic*dd,rz),
                              less(rz,zmin+(ic+1)*dd)), iota(0,len(rx)))
    if usepalette:
      if with_matplotlib:
        # --- lut is the number of colors in the current colormap
        lut = pylab.rcParams['image.lut']
        c = nint(lut*ic/(ncolor-1.))
        # --- cm is used to get the current colormap.
        # --- Calling the colormap with an integer, returns an RGB tuple
        # --- for that index.
        from matplotlib import cm
        c = getattr(cm,pylab.rcParams['image.cmap'])(c)
        # --- Note that when plotted this way, the colors of the particle
        # --- is not changed if the colormap is changed. A better way
        # --- would be to use the scatter function, which allows the color
        # --- of each particle to be specified directly, but that is
        # --- really slow.
      else:
        c = nint((npalette-1.)*ic/(ncolor-1.))
    else:
      c = color[ic%len(color)]
    plp(take(y,ii),take(x,ii),color=c,marker=marker,msize=msize,local=local)

##########################################################################
# To be implemented
#ppzx4
#ppzy4
#ppzxp4
#ppzyp4
#ppzvz4
#ppxy4
#ppxxp4
#ppyyp4
#ppxpyp4
#ppxxpco
#ppyypco


##########################################################################
##########################################################################
# history plotting routines have been replaced by those in histplots.py
##########################################################################
##########################################################################

##########################################################################
def penv(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots a and b envelope
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  """
  if not me==0: return
  plg(env.aenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(env.benv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Envelope","Z")
##########################################################################
def penvp(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots a' and b' of envelope
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  """
  if not me==0: return
  plg(env.apenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(env.bpenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Envelope slope","Z")
##########################################################################
def penvaedge(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots a envelope +/- x centroid
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  """
  if not me==0: return
  plg(+env.aenv+env.xenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(-env.aenv+env.xenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("X Envelope edges","Z")
##########################################################################
def penvbedge(color="fg",marks=0,marker=None,msize=1.0,lframe=0,titles=1):
  """
Plots b envelope +/- x centroid
  - color='fg' line color
  - marks=0 turns on identifying marks on the curve
  - marker=None marker type (see gist manual for the list)
  - msize=1.0 marker size
  - lframe=0 specifies whether or not to set plot limits
  - titles=1 specifies whether or not to plot titles
  """
  if not me==0: return
  plg(+env.benv+env.xenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  plg(-env.benv+env.xenv,env.zenv,
      color=color,marks=marks,marker=marker,msize=msize)
  if titles: ptitles("Y Envelope edges","Z")

##########################################################################
def setcmincmaxfromarray(array,kw):
  """Given plot limits, find the min and max of the array within
those limits. Otherwise, the min and max would be over the whole
array, including values outside the plot range.
This also applies the limits. All of the mins and maxes are
extracted from the kw dictionary."""
  xplmin = kw.pop('xplmin',None)
  xplmax = kw.pop('xplmax',None)
  yplmin = kw.pop('yplmin',None)
  yplmax = kw.pop('yplmax',None)

  # --- If no plot limits were specified, then do nothing.
  if xplmin is None and xplmax is None and yplmin is None and yplmax is None: return

  # --- This assumes that the min and max of the data have been specified.
  xmin = kw['xmin']
  xmax = kw['xmax']
  ymin = kw['ymin']
  ymax = kw['ymax']

  # --- Make sure that all plot limits have actual values
  if xplmin == 'e' or xplmin is None: xplmin = xmin
  if xplmax == 'e' or xplmax is None: xplmax = xmax
  if yplmin == 'e' or yplmin is None: yplmin = ymin
  if yplmax == 'e' or yplmax is None: yplmax = ymax

  # --- Apply the limits
  limits(xplmin,xplmax,yplmin,yplmax)

  # --- If cmin and cmax have already been specified, then nothing more needs
  # --- to be done,
  if 'cmin' in kw and 'cmax' in kw: return

  # --- Extract the slice of array that is within the plot limits.
  nxp1,nyp1 = array.shape
  nx,ny = nxp1-1,nyp1-1
  dx,dy = (xmax - xmin)/nx,(ymax - ymin)/ny
  ixmin = max(0,int((xplmin - xmin)/dx))
  ixmax = min(nx,int((xplmax - xmin)/dx) + 1)
  iymin = max(0,int((yplmin - ymin)/dy))
  iymax = min(nx,int((yplmax - ymin)/dy) + 1)

  # --- If the slice is sensible, find the min and max
  # --- setting cmin and cmax (if not already set).
  if (ixmin <= ixmax and iymin <= iymax):
    subarray = array[ixmin:ixmax+1,iymin:iymax+1]
    kw.setdefault('cmin',minnd(subarray))
    kw.setdefault('cmax',maxnd(subarray))

##########################################################################
def pcrhozy(ix=None,fullplane=1,lbeamframe=0,solver=None,local=0,**kw):
  """Plots contours of charge density in the Z-Y plane
  - ix=nint(-xmmin/dx): X index of plane
  - fullplane=1: when true, plots rho in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
For plotting options, see :py:func:`ppgeneric`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if ix is None: ix = nint(-solver.xmmin/solver.dx)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if kw.get('cellarray',1):
    kw.setdefault('contours',20)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (kw['xmin'],kw['xmax'],
                      solver.ymmin,solver.ymmax)
  settitles("Charge density in z-y plane","Z","Y","ix = "+repr(ix))
  rrr = getrho(ix=ix,solver=solver,local=local)
  if me > 0 and not local: rrr = zeros((solver.ny+1,solver.nz+1),'d')
  setcmincmaxfromarray(rrr.T,kw)
  ppgeneric(gridt=rrr,kwdict=kw,local=1)
  if fullplane and (solver.l2symtry or solver.l4symtry):
    ppgeneric(gridt=rrr,kwdict=kw,local=1,flipyaxis=1)
##########################################################################
def pcrhozx(iy=None,fullplane=1,lbeamframe=0,solver=None,local=0,**kw):
  """Plots contours of charge density in the Z-X plane
  - iy=nint(-ymmin/dy): Y index of plane
  - fullplane=1: when true, plots rho in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
For plotting options, see :py:func:`ppgeneric`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iy is None: iy = nint(-solver.ymmin/solver.dy)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.xmminlocal)
    kw.setdefault('ymax',solver.xmmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.xmmin)
    kw.setdefault('ymax',solver.xmmax)
  if kw.get('cellarray',1):
    kw.setdefault('contours',20)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (kw['xmin'],kw['xmax'],
                      solver.xmmin,solver.xmmax)
  settitles("Charge density in z-x plane","Z","X","iy = "+repr(iy))
  rrr = getrho(iy=iy,solver=solver,local=local)
  if me > 0 and not local: rrr = zeros((solver.nx+1,solver.nz+1),'d')
  setcmincmaxfromarray(rrr.T,kw)
  ppgeneric(gridt=rrr,kwdict=kw,local=1)
  if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
    ppgeneric(gridt=rrr,kwdict=kw,local=1,flipyaxis=1)
##########################################################################
def pcrhozr(lbeamframe=0,solver=None,local=0,**kw):
  """Plots contours of charge density in the Z-R plane
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
For plotting options, see :py:func:`ppgeneric`.
  """
  pcrhozx(iy=0,fullplane=0,lbeamframe=lbeamframe,solver=solver,local=local,**kw)
##########################################################################
def pcrhoxy(iz=None,fullplane=1,solver=None,local=0,**kw):
  """Plots contours of charge density in the X-Y plane
  - iz=nint(-zmmin/dz): Z index of plane
  - fullplane=1: when true, plots rho in the symmetric quadrants
For plotting options, see :py:func:`ppgeneric`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iz is None: iz = solver.iz_axis
  if local:
    kw.setdefault('xmin',solver.xmminlocal)
    kw.setdefault('xmax',solver.xmmaxlocal)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.xmmin)
    kw.setdefault('xmax',solver.xmmax)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if kw.get('cellarray',1):
    kw.setdefault('contours',20)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax)
  settitles("Charge density in x-y plane","X","Y","iz = "+repr(iz))
  rrr = getrho(iz=iz,solver=solver,local=local)
  if me > 0 and not local: rrr = zeros((solver.nx+1,solver.ny+1),'d')
  setcmincmaxfromarray(rrr,kw)
  ppgeneric(grid=rrr,kwdict=kw,local=1)
  if fullplane and solver.l4symtry:
    ppgeneric(grid=rrr,kwdict=kw,local=1,flipxaxis=1,flipyaxis=0)
    ppgeneric(grid=rrr,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
    ppgeneric(grid=rrr,kwdict=kw,local=1,flipxaxis=1,flipyaxis=1)
  elif fullplane and solver.l2symtry:
    ppgeneric(grid=rrr,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
##########################################################################
def pcphizy(ix=None,fullplane=1,lbeamframe=0,solver=None,local=0,**kw):
  """Plots contours of electrostatic potential in the Z-Y plane
  - ix=nint(-xmmin/dx): X index of plane
  - fullplane=1: when true, plots phi in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
For plotting options, see :py:func:`ppgeneric`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if ix is None: ix = nint(-solver.xmmin/solver.dx)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if kw.get('cellarray',1):
    kw.setdefault('contours',20)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (kw['xmin'],kw['xmax'],
                      solver.ymmin,solver.ymmax)
  settitles("Electrostatic potential in z-y plane","Z","Y","ix = "+repr(ix))
  ppp = getphi(ix=ix,solver=solver,local=local)
  if me > 0 and not local: ppp = zeros((solver.ny+1,solver.nz+1),'d')
  setcmincmaxfromarray(ppp.T,kw)
  ppgeneric(gridt=ppp,kwdict=kw,local=1)
  if fullplane and (solver.l2symtry or solver.l4symtry):
    ppgeneric(gridt=ppp,kwdict=kw,local=1,flipyaxis=1)
##########################################################################
def pcphizx(iy=None,fullplane=1,lbeamframe=0,solver=None,local=0,**kw):
  """Plots contours of electrostatic potential in the Z-X plane
  - iy=nint(-ymmin/dy): Y index of plane
  - fullplane=1: when true, plots phi in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
For plotting options, see :py:func:`ppgeneric`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iy is None: iy = nint(-solver.ymmin/solver.dy)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.xmminlocal)
    kw.setdefault('ymax',solver.xmmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.xmmin)
    kw.setdefault('ymax',solver.xmmax)
  if kw.get('cellarray',1):
    kw.setdefault('contours',20)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (kw['xmin'],kw['xmax'],
                      solver.xmmin,solver.xmmax)
  settitles("Electrostatic potential in z-x plane","Z","X","iy = "+repr(iy))
  ppp = getphi(iy=iy,solver=solver,local=local)
  if me > 0 and not local: ppp = zeros((solver.nx+1,solver.nz+1),'d')
  setcmincmaxfromarray(ppp.T,kw)
  ppgeneric(gridt=ppp,kwdict=kw,local=1)
  if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
    ppgeneric(gridt=ppp,kwdict=kw,local=1,flipyaxis=1)
##########################################################################
def pcphizr(lbeamframe=0,solver=None,local=0,**kw):
  """Plots contours of electrostatic potential in the Z-R plane
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
For plotting options, see :py:func:`ppgeneric`.
  """
  pcphizx(iy=0,fullplane=0,lbeamframe=lbeamframe,solver=solver,local=local,**kw)
##########################################################################
def pcphixy(iz=None,fullplane=1,solver=None,local=0,**kw):
  """Plots contours of electrostatic potential in the X-Y plane
  - iz=nint(-zmmin/dz): Z index of plane
  - fullplane=1: when true, plots phi in the symmetric quadrants
For plotting options, see :py:func:`ppgeneric`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iz is None: iz = solver.iz_axis
  if local:
    kw.setdefault('xmin',solver.xmminlocal)
    kw.setdefault('xmax',solver.xmmaxlocal)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.xmmin)
    kw.setdefault('xmax',solver.xmmax)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if kw.get('cellarray',1):
    kw.setdefault('contours',20)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax)
  settitles("Electrostatic potential in x-y plane","X","Y","iz = "+repr(iz))
  ppp = getphi(iz=iz,solver=solver,local=local)
  if me > 0 and not local: ppp = zeros((solver.nx+1,solver.ny+1),'d')
  setcmincmaxfromarray(ppp,kw)
  ppgeneric(grid=ppp,kwdict=kw,local=1,flipxaxis=0,flipyaxis=0)
  if fullplane and solver.l4symtry:
    ppgeneric(grid=ppp,kwdict=kw,local=1,flipxaxis=1,flipyaxis=0)
    ppgeneric(grid=ppp,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
    ppgeneric(grid=ppp,kwdict=kw,local=1,flipxaxis=1,flipyaxis=1)
  elif fullplane and solver.l2symtry:
    ppgeneric(grid=ppp,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
##########################################################################
def pcselfezy(comp=None,ix=None,fullplane=1,solver=None,
              lbeamframe=0,vec=0,sz=1,sy=1,local=0,**kw):
  """Plots contours of self E field in the Z-Y plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'E'.
          Use 'E' to get the field magnitude.
  - ix=nint(-xmmin/dx): X index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sy=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if ix is None: ix = nint(-solver.xmmin/solver.dx)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits']=(kw['xmin'],kw['xmax'],
                    solver.ymmin,solver.ymmax)
  if comp == 'E': titlet = 'Emagnitude in z-y plane'
  else:           titlet = 'E%s in z-y plane'%comp
  settitles(titlet,"Z","Y","ix = "+repr(ix))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    eee = getselfe(comp=comp,ix=ix,solver=solver,local=local)
    if me > 0 and not local: eee = zeros((solver.ny+1,solver.nz+1),'d')
    setcmincmaxfromarray(eee.T,kw)
    ppgeneric(gridt=eee,kwdict=kw,local=1)
    if fullplane and (solver.l2symtry or solver.l4symtry):
      ppgeneric(gridt=eee,kwdict=kw,local=1,flipyaxis=1)
  else:
    if fullplane and (solver.l2symtry or solver.l4symtry):
      kw['ymin'] = - kw['ymax']
    ey = getselfe(comp='y',ix=ix,fullplane=fullplane,solver=solver,local=local)
    ez = getselfe(comp='z',ix=ix,fullplane=fullplane,solver=solver,local=local)
    ppvector(transpose(ey[::sy,::sz]),transpose(ez[::sy,::sz]),kwdict=kw,local=1)
##########################################################################
def pcselfezx(comp=None,iy=None,fullplane=1,solver=None,
              lbeamframe=0,vec=0,sz=1,sx=1,local=0,**kw):
  """Plots contours of self E field in the Z-X plane
  - comp: field component to plot, either 'x', 'y', 'z', or 'E'.
          Use 'E' to get the field magnitude.
  - iy=nint(-ymmin/dy): Y index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sx=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iy is None: iy = nint(-solver.ymmin/solver.dy)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.xmminlocal)
    kw.setdefault('ymax',solver.xmmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.xmmin)
    kw.setdefault('ymax',solver.xmmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (kw['xmin'],kw['xmax'],
                      solver.xmmin,solver.xmmax)
  if comp == 'E': titlet = 'Emagnitude in z-x plane'
  else:           titlet = 'E%s in z-x plane'%comp
  settitles(titlet,"Z","X","iy = "+repr(iy))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    eee = getselfe(comp=comp,iy=iy,solver=solver,local=local)
    if me > 0 and not local: eee = zeros((solver.nx+1,solver.nz+1),'d')
    setcmincmaxfromarray(eee.T,kw)
    ppgeneric(gridt=eee,kwdict=kw,local=1)
    if fullplane and solver.l4symtry:
      ppgeneric(gridt=eee,kwdict=kw,local=1,flipyaxis=1)
  else:
    if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
      kw['ymin'] = - kw['ymax']
    ex = getselfe(comp='x',iy=iy,fullplane=fullplane,solver=solver,local=local)
    ez = getselfe(comp='z',iy=iy,fullplane=fullplane,solver=solver,local=local)
    ppvector(transpose(ex[::sx,::sz]),transpose(ez[::sx,::sz]),kwdict=kw,local=1)
##########################################################################
def pcselfezr(comp=None,solver=None,
              lbeamframe=0,vec=0,sz=1,sr=1,local=0,**kw):
  """Plots contours of self E field the Z-R plane
  - comp: field component to plot, either 'r', 'x', 'y', 'z' or 'E'.
          'r' and 'x' are the same thing.
          Use 'E' to get the field magnitude.
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sr=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if comp == 'r': comp = 'x'
  pcselfezx(comp=comp,iy=0,fullplane=0,solver=solver,
            lbeamframe=lbeamframe,vec=vec,sz=sz,sx=sr,local=local,**kw)
##########################################################################
def pcselfexy(comp=None,iz=None,fullplane=1,solver=None,vec=0,sx=1,sy=1,
              local=0,**kw):
  """Plots contours of self E field in the X-Y plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'E'.
          Use 'E' to get the field magnitude.
  - iz=nint(-zmmin/dz): Z index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - vec=0: when true, plots E field vectors
  - sx,sy=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iz is None: iz = solver.iz_axis
  if local:
    kw.setdefault('xmin',solver.xmminlocal)
    kw.setdefault('xmax',solver.xmmaxlocal)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.xmmin)
    kw.setdefault('xmax',solver.xmmax)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax)
  if comp == 'E': titlet = 'Emagnitude in x-y plane'
  else:           titlet = 'E%s in x-y plane'%comp
  settitles(titlet,"X","Y","iz = "+repr(iz))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    eee = getselfe(comp=comp,iz=iz,solver=solver,local=local)
    if me > 0 and not local: eee = zeros((solver.nx+1,solver.ny+1),'d')
    setcmincmaxfromarray(eee,kw)
    ppgeneric(grid=eee,kwdict=kw,local=1)
    if fullplane and solver.l4symtry:
      ppgeneric(grid=eee,kwdict=kw,local=1,flipxaxis=1,flipyaxis=0)
      ppgeneric(grid=eee,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
      ppgeneric(grid=eee,kwdict=kw,local=1,flipxaxis=1,flipyaxis=1)
    elif fullplane and solver.l2symtry:
      ppgeneric(grid=eee,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
  else:
    if fullplane and solver.l4symtry:
      kw['ymin'] = - kw['ymax']
      kw['xmin'] = - kw['xmax']
    elif fullplane and solver.l2symtry:
      kw['ymin'] = - kw['ymax']
    ex = getselfe(comp='x',iz=iz,fullplane=fullplane,solver=solver,local=local)
    ey = getselfe(comp='y',iz=iz,fullplane=fullplane,solver=solver,local=local)
    ppvector(ey[::sx,::sy],ex[::sx,::sy],kwdict=kw,local=1)
##########################################################################
def pcjzy(comp=None,ix=None,fullplane=1,solver=None,
          lbeamframe=0,vec=0,sz=1,sy=1,local=0,**kw):
  """Plots contours of current density in the Z-Y plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'J'.
          Use 'J' to get the field magnitude.
  - ix=nint(-xmmin/dx): X index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sy=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if ix is None: ix = nint(-solver.xmmin/solver.dx)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits']=(kw['xmin'],kw['xmax'],
                    solver.ymmin,solver.ymmax)
  if comp == 'J': titlet = 'Current Density Jmagnitude in z-y plane'
  else:           titlet = 'Current Density J%s in z-y plane'%comp
  settitles(titlet,"Z","Y","ix = "+repr(ix))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    j = getj(comp=comp,ix=ix,solver=solver,local=local)
    if me > 0 and not local: j = zeros((solver.ny+1,solver.nz+1),'d')
    setcmincmaxfromarray(j.T,kw)
    ppgeneric(gridt=j,kwdict=kw,local=1)
    if fullplane and (solver.l2symtry or solver.l4symtry):
      ppgeneric(gridt=j,kwdict=kw,local=1,flipyaxis=1)
  else:
    if fullplane and (solver.l2symtry or solver.l4symtry):
      kw['ymin'] = - kw['ymax']
    jy = getj(comp='y',ix=ix,fullplane=fullplane,solver=solver,local=local)
    jz = getj(comp='z',ix=ix,fullplane=fullplane,solver=solver,local=local)
    ppvector(transpose(jy[::sy,::sz]),transpose(jz[::sy,::sz]),kwdict=kw,local=1)
##########################################################################
def pcjzx(comp=None,iy=None,fullplane=1,solver=None,
              lbeamframe=0,vec=0,sz=1,sx=1,local=0,**kw):
  """Plots contours of current density in the Z-X plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'J'.
          Use 'J' to get the field magnitude.
  - iy=nint(-ymmin/dy): Y index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sx=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iy is None: iy = nint(-solver.ymmin/solver.dy)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.xmminlocal)
    kw.setdefault('ymax',solver.xmmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.xmmin)
    kw.setdefault('ymax',solver.xmmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (kw['xmin'],kw['xmax'],
                      solver.xmmin,solver.xmmax)
  if comp == 'J': titlet = 'Current Density Jmagnitude in z-x plane'
  else:           titlet = 'Current Density J%s in z-x plane'%comp
  settitles(titlet,"Z","X","iy = "+repr(iy))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    j = getj(comp=comp,iy=iy,solver=solver,local=local)
    if me > 0 and not local: j = zeros((solver.nx+1,solver.nz+1),'d')
    setcmincmaxfromarray(j.T,kw)
    ppgeneric(gridt=j,kwdict=kw,local=1)
    if fullplane and solver.l4symtry:
      ppgeneric(gridt=j,kwdict=kw,local=1,flipyaxis=1)
  else:
    if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
      kw['ymin'] = - kw['ymax']
    jx = getj(comp='x',iy=iy,fullplane=fullplane,solver=solver,local=local)
    jz = getj(comp='z',iy=iy,fullplane=fullplane,solver=solver,local=local)
    ppvector(transpose(jx[::sx,::sz]),transpose(jz[::sx,::sz]),kwdict=kw,local=1)
##########################################################################
def pcjxy(comp=None,iz=None,fullplane=1,solver=None,vec=0,sx=1,sy=1,
          local=0,**kw):
  """Plots contours of current density in the X-Y plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'J'.
          Use 'J' to get the field magnitude.
  - iz=nint(-zmmin/dz): Z index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - vec=0: when true, plots E field vectors
  - sx,sy=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iz is None: iz = solver.iz_axis
  if local:
    kw.setdefault('xmin',solver.xmminlocal)
    kw.setdefault('xmax',solver.xmmaxlocal)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.xmmin)
    kw.setdefault('xmax',solver.xmmax)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax)
  if comp == 'J': titlet = 'Current Density Jmagnitude in x-y plane'
  else:           titlet = 'Current Density J%s in x-y plane'%comp
  settitles(titlet,"X","Y","iz = "+repr(iz))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    j = getj(comp=comp,iz=iz,solver=solver,local=local)
    if me > 0 and not local: j = zeros((solver.nx+1,solver.ny+1),'d')
    setcmincmaxfromarray(j,kw)
    ppgeneric(grid=j,kwdict=kw,local=1)
    if fullplane and solver.l4symtry:
      ppgeneric(grid=j,kwdict=kw,local=1,flipxaxis=1,flipyaxis=0)
      ppgeneric(grid=j,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
      ppgeneric(grid=j,kwdict=kw,local=1,flipxaxis=1,flipyaxis=1)
    elif fullplane and solver.l2symtry:
      ppgeneric(grid=j,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
  else:
    if fullplane and solver.l4symtry:
      kw['ymin'] = - kw['ymax']
      kw['xmin'] = - kw['xmax']
    elif fullplane and solver.l2symtry:
      kw['ymin'] = - kw['ymax']
    jx = getj(comp='x',iz=iz,fullplane=fullplane,solver=solver,local=local)
    jy = getj(comp='y',iz=iz,fullplane=fullplane,solver=solver,local=local)
    ppvector(jy[::sx,::sy],jx[::sx,::sy],kwdict=kw,local=1)
##########################################################################
def pcbzy(comp=None,ix=None,fullplane=1,solver=None,
          lbeamframe=0,vec=0,sz=1,sy=1,local=0,**kw):
  """Plots contours of the magnetic field in the Z-Y plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'B'.
          Use 'B' to get the field magnitude.
  - ix=nint(-xmmin/dx): X index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sy=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if ix is None: ix = nint(-solver.xmmin/solver.dx)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits']=(kw['xmin'],kw['xmax'],
                    solver.ymmin,solver.ymmax)
  if comp == 'B': titlet = 'Magnetic Field Bmagnitude in z-y plane'
  else:           titlet = 'Magnetic Field B%s in z-y plane'%comp
  settitles(titlet,"Z","Y","ix = "+repr(ix))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    b = getb(comp=comp,ix=ix,solver=solver,local=local)
    if me > 0 and not local: b = zeros((solver.ny+1,solver.nz+1),'d')
    setcmincmaxfromarray(b.T,kw)
    ppgeneric(gridt=b,kwdict=kw,local=1)
    if fullplane and (solver.l2symtry or solver.l4symtry):
      ppgeneric(gridt=b,kwdict=kw,local=1,flipyaxis=1)
  else:
    if fullplane and (solver.l2symtry or solver.l4symtry):
      kw['ymin'] = - kw['ymax']
    by = getb(comp='y',ix=ix,fullplane=fullplane,solver=solver,local=local)
    bz = getb(comp='z',ix=ix,fullplane=fullplane,solver=solver,local=local)
    ppvector(transpose(by[::sy,::sz]),transpose(bz[::sy,::sz]),kwdict=kw,local=1)
##########################################################################
def pcbzx(comp=None,iy=None,fullplane=1,solver=None,
          lbeamframe=0,vec=0,sz=1,sx=1,local=0,**kw):
  """Plots contours of the magnetic field in the Z-X plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'B'.
          Use 'B' to get the field magnitude.
  - iy=nint(-ymmin/dy): Y index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sx=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iy is None: iy = nint(-solver.ymmin/solver.dy)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.xmminlocal)
    kw.setdefault('ymax',solver.xmmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.xmmin)
    kw.setdefault('ymax',solver.xmmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (kw['xmin'],kw['xmax'],
                      solver.xmmin,solver.xmmax)
  if comp == 'B': titlet = 'Magnetic Field Bmagnitude in z-x plane'
  else:           titlet = 'Magnetic Field B%s in z-x plane'%comp
  settitles(titlet,"Z","X","iy = "+repr(iy))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    b = getb(comp=comp,iy=iy,solver=solver,local=local)
    if me > 0 and not local: b = zeros((solver.nx+1,solver.nz+1),'d')
    setcmincmaxfromarray(b.T,kw)
    ppgeneric(gridt=b,kwdict=kw,local=1)
    if fullplane and solver.l4symtry:
      ppgeneric(gridt=b,kwdict=kw,local=1,flipyaxis=1)
  else:
    if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
      kw['ymin'] = - kw['ymax']
    bx = getb(comp='x',iy=iy,fullplane=fullplane,solver=solver,local=local)
    bz = getb(comp='z',iy=iy,fullplane=fullplane,solver=solver,local=local)
    ppvector(transpose(bx[::sx,::sz]),transpose(bz[::sx,::sz]),kwdict=kw,local=1)
##########################################################################
def pcbxy(comp=None,iz=None,fullplane=1,solver=None,vec=0,sx=1,sy=1,
          local=0,**kw):
  """Plots contours of the magnetic field in the X-Y plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'B'.
          Use 'B' to get the field magnitude.
  - iz=nint(-zmmin/dz): Z index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - vec=0: when true, plots E field vectors
  - sx,sy=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iz is None: iz = solver.iz_axis
  if local:
    kw.setdefault('xmin',solver.xmminlocal)
    kw.setdefault('xmax',solver.xmmaxlocal)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.xmmin)
    kw.setdefault('xmax',solver.xmmax)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax)
  if comp == 'B': titlet = 'Magnetic Field Bmagnitude in x-y plane'
  else:           titlet = 'Magnetic Field B%s in x-y plane'%comp
  settitles(titlet,"X","Y","iz = "+repr(iz))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    b = getb(comp=comp,iz=iz,solver=solver,local=local)
    if me > 0 and not local: b = zeros((solver.nx+1,solver.ny+1),'d')
    setcmincmaxfromarray(b,kw)
    ppgeneric(grid=b,kwdict=kw,local=1)
    if fullplane and solver.l4symtry:
      ppgeneric(grid=b,kwdict=kw,local=1,flipxaxis=1,flipyaxis=0)
      ppgeneric(grid=b,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
      ppgeneric(grid=b,kwdict=kw,local=1,flipxaxis=1,flipyaxis=1)
    elif fullplane and solver.l2symtry:
      ppgeneric(grid=b,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
  else:
    if fullplane and solver.l4symtry:
      kw['ymin'] = - kw['ymax']
      kw['xmin'] = - kw['xmax']
    elif fullplane and solver.l2symtry:
      kw['ymin'] = - kw['ymax']
    bx = getb(comp='x',iz=iz,fullplane=fullplane,solver=solver,local=local)
    by = getb(comp='y',iz=iz,fullplane=fullplane,solver=solver,local=local)
    ppvector(by[::sx,::sy],bx[::sx,::sy],kwdict=kw,local=1)
##########################################################################
def pcazy(comp=None,ix=None,fullplane=1,solver=None,
          lbeamframe=0,vec=0,sz=1,sy=1,local=0,**kw):
  """Plots contours of the magnetic vector potential in the Z-Y plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'A'.
          Use 'A' to get the field magnitude.
  - ix=nint(-xmmin/dx): X index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sy=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if ix is None: ix = nint(-solver.xmmin/solver.dx)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits']=(kw['xmin'],kw['xmax'],
                    solver.ymmin,solver.ymmax)
  if comp == 'A': titlet = 'Magnetic Vector Potential Amagnitude in z-y plane'
  else:           titlet = 'Magnetic Vector Potential A%s in z-y plane'%comp
  settitles(titlet,"Z","Y","ix = "+repr(ix))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    a = geta(comp=comp,ix=ix,solver=solver,local=local)
    if me > 0 and not local: a = zeros((solver.ny+1,solver.nz+1),'d')
    setcmincmaxfromarray(a.T,kw)
    ppgeneric(gridt=a,kwdict=kw,local=1)
    if fullplane and (solver.l2symtry or solver.l4symtry):
      ppgeneric(gridt=a,kwdict=kw,local=1,flipyaxis=1)
  else:
    if fullplane and (solver.l2symtry or solver.l4symtry):
      kw['ymin'] = - kw['ymax']
    ay = geta(comp='y',ix=ix,fullplane=fullplane,solver=solver,local=local)
    az = geta(comp='z',ix=ix,fullplane=fullplane,solver=solver,local=local)
    ppvector(transpose(ay[::sy,::sz]),transpose(az[::sy,::sz]),kwdict=kw,local=1)
##########################################################################
def pcazx(comp=None,iy=None,fullplane=1,solver=None,
          lbeamframe=0,vec=0,sz=1,sx=1,local=0,**kw):
  """Plots contours of the magnetic vector potential in the Z-X plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'A'.
          Use 'A' to get the field magnitude.
  - iy=nint(-ymmin/dy): Y index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - lbeamframe=0: when true, plot relative to beam frame, otherwise lab frame
  - xplmin,xplmax,yplmin,yplmax: optional plot limits.
  - vec=0: when true, plots E field vectors
  - sz,sx=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iy is None: iy = nint(-solver.ymmin/solver.dy)
  if lbeamframe: zbeam = 0.
  else:          zbeam = top.zbeam
  if local:
    kw.setdefault('xmin',solver.zmminlocal + zbeam)
    kw.setdefault('xmax',solver.zmmaxlocal + zbeam)
    kw.setdefault('ymin',solver.xmminlocal)
    kw.setdefault('ymax',solver.xmmaxlocal)
  else:
    kw.setdefault('xmin',solver.zmmin + zbeam)
    kw.setdefault('xmax',solver.zmmax + zbeam)
    kw.setdefault('ymin',solver.xmmin)
    kw.setdefault('ymax',solver.xmmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (kw['xmin'],kw['xmax'],
                      solver.xmmin,solver.xmmax)
  if comp == 'A': titlet = 'Magnetic Vector Potential Amagnitude in z-x plane'
  else:           titlet = 'Magnetic Vector Potential A%s in z-x plane'%comp
  settitles(titlet,"Z","X","iy = "+repr(iy))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    a = geta(comp=comp,iy=iy,solver=solver,local=local)
    if me > 0 and not local: a = zeros((solver.nx+1,solver.nz+1),'d')
    setcmincmaxfromarray(a.T,kw)
    ppgeneric(gridt=a,kwdict=kw,local=1)
    if fullplane and solver.l4symtry:
      ppgeneric(gridt=a,kwdict=kw,local=1,flipyaxis=1)
  else:
    if fullplane and (solver.l4symtry or solver.solvergeom == w3d.RZgeom):
      kw['ymin'] = - kw['ymax']
    ax = geta(comp='x',iy=iy,fullplane=fullplane,solver=solver,local=local)
    az = geta(comp='z',iy=iy,fullplane=fullplane,solver=solver,local=local)
    ppvector(transpose(ax[::sx,::sz]),transpose(az[::sx,::sz]),kwdict=kw,local=1)
##########################################################################
def pcaxy(comp=None,iz=None,fullplane=1,solver=None,vec=0,sx=1,sy=1,
          local=0,**kw):
  """Plots contours of the magnetic vector potential in the X-Y plane
  - comp: field component to plot, either 'x', 'y', 'z' or 'A'.
          Use 'A' to get the field magnitude.
  - iz=nint(-zmmin/dz): Z index of plane
  - fullplane=1: when true, plots E in the symmetric quadrants
  - vec=0: when true, plots E field vectors
  - sx,sy=1: step size in grid for plotting fewer points
For plotting options, see :py:func:`ppgeneric` or :py:func:`ppvector`.
  """
  if solver is None: solver = (getregisteredsolver() or w3d)
  if iz is None: iz = solver.iz_axis
  if local:
    kw.setdefault('xmin',solver.xmminlocal)
    kw.setdefault('xmax',solver.xmmaxlocal)
    kw.setdefault('ymin',solver.ymminlocal)
    kw.setdefault('ymax',solver.ymmaxlocal)
  else:
    kw.setdefault('xmin',solver.xmmin)
    kw.setdefault('xmax',solver.xmmax)
    kw.setdefault('ymin',solver.ymmin)
    kw.setdefault('ymax',solver.ymmax)
  if 'pplimits' in kw:
    kw['lframe'] = 1
  else:
    kw['pplimits'] = (solver.xmmin,solver.xmmax,solver.ymmin,solver.ymmax)
  if comp == 'A': titlet = 'Magnetic Vector Potential Amagnitude in x-y plane'
  else:           titlet = 'Magnetic Vector Potential A%s in x-y plane'%comp
  settitles(titlet,"X","Y","iz = "+repr(iz))
  if not vec:
    if kw.get('cellarray',1):
      kw.setdefault('contours',20)
    a = geta(comp=comp,iz=iz,solver=solver,local=local)
    if me > 0 and not local: a = zeros((solver.nx+1,solver.ny+1),'d')
    setcmincmaxfromarray(a,kw)
    ppgeneric(grid=a,kwdict=kw,local=1)
    if fullplane and solver.l4symtry:
      ppgeneric(grid=a,kwdict=kw,local=1,flipxaxis=1,flipyaxis=0)
      ppgeneric(grid=a,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
      ppgeneric(grid=a,kwdict=kw,local=1,flipxaxis=1,flipyaxis=1)
    elif fullplane and solver.l2symtry:
      ppgeneric(grid=a,kwdict=kw,local=1,flipxaxis=0,flipyaxis=1)
  else:
    if fullplane and solver.l4symtry:
      kw['ymin'] = - kw['ymax']
      kw['xmin'] = - kw['xmax']
    elif fullplane and solver.l2symtry:
      kw['ymin'] = - kw['ymax']
    ax = geta(comp='x',iz=iz,fullplane=fullplane,solver=solver,local=local)
    ay = geta(comp='y',iz=iz,fullplane=fullplane,solver=solver,local=local)
    ppvector(ay[::sx,::sy],ax[::sx,::sy],kwdict=kw,local=1)
##########################################################################
##########################################################################
def ppdecompositionz(scale=1.,minscale=0.,gap=0.2):
  """Shows the domain decomposition in a graphical way. For each
processor, the total mesh extent is plotted as a filled rectangle
covering the z-length and with height determined by 'scale' and the
number of processors. Another filled rectangle is plotted in the top
half showing the particle domains, and one on the lower half shows the
field domain.
  - scale=1.: the maximum vertical extent of the graph
  - minscale=0.: the minimum vertical extent of the graph
  - gap=0.2: fractional vertical gap between rectangles
  """
  z = []
  x = []
  y = []
  dd = 1.*scale/top.nprocs
  mm = 1. - gap
  for i in xrange(top.nprocs):
    z = z + [1.]
    zmin = top.izfsslave[i]*w3d.dz + w3d.zmmin
    zmax = (top.izfsslave[i] + top.nzfsslave[i])*w3d.dz + w3d.zmmin
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([-mm,-mm,mm,mm,-mm]))
  for i in xrange(top.nprocs):
    z = z + [2.]
    zmin = top.zpslmin[i]
    zmax = top.zpslmax[i]
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([0,0,mm,mm,0]))
  for i in xrange(top.nprocs):
    z = z + [3.]
    zmin = top.izfsslave[i]*w3d.dz
    zmax = top.izfsslave[i]*w3d.dz + top.nzfsslave[i]*w3d.dz
    x = x + [zmin,zmax,zmax,zmin,zmin]
    y = y + list(i*dd + 0.5*dd*array([-mm,-mm,0,0,-mm]))
  plfp(array(z),y,x,5*ones(len(z),'l'),cmin=0,cmax=4,local=1)
  for i in xrange(len(z)):
    pldj(x[i*5:i*5+4],y[i*5:i*5+4],x[i*5+1:i*5+5],y[i*5+1:i*5+5],local=1)

def _ppdecomposition_work(ix,nx,iz,nz):
  for izproc in xrange(len(iz)):
    for ixproc in xrange(len(ix)):
      ix1 = ix[ixproc]
      ix2 = ix[ixproc]+nx[ixproc]
      iz1 = iz[izproc]
      iz2 = iz[izproc]+nz[izproc]
      plg([ix1,ix1,ix2,ix2,ix1],[iz1,iz2,iz2,iz1,iz1],
          color=color[(ixproc+izproc*len(ix))%len(color)])

def ppdecompzx(decomp=None,scale=1):
  """Plots the decomposition in the z-x plane.
  - decomp=top.fsdecomp: a Decomposition object. Another example is
    top.ppdecomp, the particle decomposition
  - scale=true: When true, plot in in meters, otherwise grid cells
  """
  if decomp is None: decomp=top.fsdecomp
  ix = decomp.ix
  nx = decomp.nx
  iz = decomp.iz
  nz = decomp.nz
  if scale:
    ix = w3d.xmmin + ix*w3d.dx
    nx = nx*w3d.dx
    iz = w3d.zmmin + iz*w3d.dz
    nz = nz*w3d.dz
  _ppdecomposition_work(ix,nx,iz,nz)

def ppdecompzy(decomp=None,scale=1):
  """Plots the decomposition in the z-y plane.
  - decomp=top.fsdecomp: a Decomposition object. Another example is
    top.ppdecomp, the particle decomposition
  - scale=true: When true, plot in in meters, otherwise grid cells
  """
  if decomp is None: decomp=top.fsdecomp
  iy = decomp.iy
  ny = decomp.ny
  iz = decomp.iz
  nz = decomp.nz
  if scale:
    iy = w3d.ymmin + iy*w3d.dy
    ny = ny*w3d.dy
    iz = w3d.zmmin + iz*w3d.dz
    nz = nz*w3d.dz
  _ppdecomposition_work(iy,ny,iz,nz)

def ppdecompxy(decomp=None,scale=1):
  """Plots the decomposition in the x-y plane.
  - decomp=top.fsdecomp: a Decomposition object. Another example is
    top.ppdecomp, the particle decomposition
  - scale=true: When true, plot in in meters, otherwise grid cells
  """
  if decomp is None: decomp=top.fsdecomp
  ix = decomp.ix
  nx = decomp.nx
  iy = decomp.iy
  ny = decomp.ny
  if scale:
    ix = w3d.xmmin + ix*w3d.dx
    nx = nx*w3d.dx
    iy = w3d.ymmin + iy*w3d.dy
    ny = ny*w3d.dy
  _ppdecomposition_work(iy,ny,ix,nx)

##########################################################################
##########################################################################
def pltfld3d(fld='phi',freqflag=always):
  """Makes fields plots which have been turned on
  - fld='phi' quantity to plot, either 'phi' or 'rho'
  - freqflag=always frequency flag, either always, seldom, or never"""
  currentwindow = active_window()
  active_window(0)
  nwindows = 9
  for i in xrange(nwindows):
    if (w3d.icrhoxy[i] == freqflag and fld == "rho"): pcrhoxy[i]
    if (w3d.icrhozx[i] == freqflag and fld == "rho"): pcrhozx[i]
    if (w3d.icrhozy[i] == freqflag and fld == "rho"): pcrhozy[i]
    if (w3d.icphixy[i] == freqflag and fld == "phi"): pcphixy[i]
    if (w3d.icphizx[i] == freqflag and fld == "phi"): pcphizx[i]
    if (w3d.icphizy[i] == freqflag and fld == "phi"): pcphizy[i]
  #if (top.icrhoxy4 == freqflag and fld == "rho"): pcrhoxy4
  #if (top.icrhozx4 == freqflag and fld == "rho"): pcrhozx4
  #if (top.icrhozy4 == freqflag and fld == "rho"): pcrhozy4
  #if (top.icphixy4 == freqflag and fld == "phi"): pcphixy4
  #if (top.icphizx4 == freqflag and fld == "phi"): pcphizx4
  #if (top.icphizy4 == freqflag and fld == "phi"): pcphizy4
  oldlimits = limits()
  active_window(currentwindow)

##########################################################################
def onedplts(freqflag=always):
  """Makes 1-D plots which have been turned on
  - freqflag=always frequency flag, either always, seldom, or never"""
  currentwindow = active_window()
  active_window(0)
  if freqflag == top.ipcurr: pzcurr()
  if freqflag == top.ipegap: pzegap()
  if freqflag == top.iplchg: pzlchg()
  if freqflag == top.ipvzofz: pzvzofz()
  if freqflag == top.iprhoax: pzrhoax()
  if freqflag == top.ipphiax: pzphiax()
  if freqflag == top.ipezax: pzezax()
  oldlimits = limits()
  active_window(currentwindow)

# --- Thses are defined for the fortran interface. If WARP is not imported
# --- main, then the functions and the always and seldom parameters will
# --- not be accessible from the fortran call. This way avoids that by
# --- declaring parameterless functions and explicitly adding them to main.
def onedpltsalways():
  onedplts(always)
def onedpltsseldom():
  onedplts(seldom)

##########################################################################
def psplots(freqflag=always,js=0):
  """Makes particle phase space plots which have been turned on
  - freqflag=always frequency flag, either always, seldom, or never
  - js=0 specifies the species of particles to plot"""
  # --- Phase space plots, both "frequent" ones and others
  # --- Do z-x,y 2-to-a-page subset and all-particle plots
  bb = wtime()
  # --- Save current device and set active device to window(0). This
  # --- ensures that plots created by this routine will be dumped to
  # --- the appropriate plot file.
  currentwindow = active_window()
  active_window(0)

  nsubsets = 3
  nwindows = 9

  for i in xrange(-nsubsets,1):
    if (top.ipzxy[i] == freqflag):
      ppzxy(i,lframe=true)
      fma()

  # --- Do z-x,y 2-to-a-page in color, skipping NSKIPCOL particles
  if (top.ipzxyco == freqflag):
    ppzxyco(js,lframe=true)
    fma()

  # --- Do z-vz in color, skipping NSKIPCOL particles
  if (top.ipzvzco == freqflag):
    ppzvzco(js,lframe=true)
    fma()

  # --- Do x-xp in color, skipping NSKIPCOL particles
  for i in xrange(nwindows+1):
   if (top.ipxxpco[i] == freqflag):
     ppxxpco(i,lframe=true)
     fma()

  # --- Do y-yp in color, skipping NSKIPCOL particles
  for i in xrange(nwindows+1):
   if (top.ipyypco[i] == freqflag):
     ppyypco(i,lframe=true)
     fma()

  # --- Do z-x and z-xp subset and y-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzx[i] == freqflag):
      ppzx(i,lframe=true)
      fma()
  #if (top.ipzx4 == freqflag):
    #ppzx4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzxp[i] == freqflag):
      ppzxp(i,lframe=true)
      fma()
  #if (top.ipzxp4 == freqflag):
    #ppzxp4

  # --- Do z-y and z-yp subset and x-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzy[i] == freqflag):
      ppzy(i,lframe=true)
      fma()
  #if (top.ipzy4 == freqflag):
    #ppzy4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzyp[i] == freqflag):
      ppzyp(i,lframe=true)
      fma()
  #if (top.ipzyp4 == freqflag):
    #ppzyp4
    #fma()

  # --- Do z-vz subset and r-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipzvz[i] == freqflag):
      ppzvz(i,lframe=true)
      fma()
  #if (top.ipzvz4 == freqflag):
    #ppzvz4
    #fma()

  # --- Do transverse phase-space subset and z-window plots
  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxy[i] == freqflag):
      ppxy(i,lframe=true)
      fma()
  #if (top.ipxy4 == freqflag):
    #ppxy4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxxp[i] == freqflag):
      ppxxp(i,lframe=true)
      fma()
  #if (top.ipxxp4 == freqflag):
    #ppxxp4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipyyp[i] == freqflag):
      ppyyp(i,lframe=true)
      fma()
  #if (top.ipyyp4 == freqflag):
    #ppyyp4
    #fma()

  for i in xrange(-nsubsets,nwindows+1):
    if (top.ipxpyp[i] == freqflag):
      ppxpyp(i,lframe=true)
      fma()
  #if (top.ipxpyp4 == freqflag):
    #ppxpyp4
    #fma()

  # --- Do trace-space z-window plots
  for i in xrange(nwindows+1):
    if (top.iptrace[i] == freqflag and i >= 0):
      pptrace(i,lframe=true)
      fma()

  # --- Do the user defined plots
  if freqflag == always: controllers.callplalwaysfuncs()
  if freqflag == seldom: controllers.callplseldomfuncs()

  # --- Reset the current window to it previous value.
  active_window(currentwindow)

  # --- Accumulate time
  aa = wtime()
  try: psplots.time = psplots.time + (aa - bb)
  except: psplots.time = 0.

# --- Thses are defined for the fortran interface. If WARP is not imported
# --- main, then the functions and the always and seldom parameters will
# --- not be accessible from the fortran call. This way avoids that by
# --- declaring parameterless functions and explicitly adding them to main.
def psplotsalways():
  psplots(always)
def psplotsseldom():
  psplots(seldom)

##########################################################################
def gstyle():
    global gist_style
    gist_style = gist.get_style()
    for i in range(0,len(gist_style['systems'])):
      gist_style['systems'][i]['legend']=''

##########################################################################
def set_label(height=None,font=None,bold=0,italic=0,axis='all',system=None):
    """change plots label attributes
    - height=None,
    - scale=1.,
    - font=None ('Courier'=0,'Times'=1,'Helvetica'=2,'Symbol'=3,'New Century'=4),
    - bold=0
    - italic=0
    - axis='all'
    - system='all'"""
    try:
      gist_style
    except:
      gstyle()

    if font is not None:
      if(type(font)==type('string')):
        if font == 'Courier':     font = 0
        if font == 'Times':       font = 1
        if font == 'Helvetica':   font = 2
        if font == 'Symbol':      font = 3
        if font == 'New Century': font = 4
      font=4*font+bold+2*italic
    if system is None:
      # --- Not sure why plsys is called twice, but the rewrite just below
      # --- is needed to be consistent with the wrapped version of plsys
      # --- defined above.
      #systems = [plsys(plsys())-1]
      view = gist.plsys()
      gist.plsys(view)
      view = gist.plsys()
      systems = [view-1]
    else:
      if(system=='all'):
        systems = range(0,len(gist_style['systems']))
      else:
        systems = [system-1]
    for i in systems:
      if(axis=='x' or axis=='all'):
        if height is not None:
          gist_style ['systems'][i]['ticks']['horizontal']['textStyle']['height']= height
        if font is not None:
          gist_style ['systems'][i]['ticks']['horizontal']['textStyle']['font']=font
      if(axis=='y' or axis=='all'):
        if height is not None:
          gist_style ['systems'][i]['ticks']['vertical']['textStyle']['height']= height
        if font is not None:
          gist_style ['systems'][i]['ticks']['vertical']['textStyle']['font']=font
    set_style(gist_style)

##########################################################################
def scale_labels(scale):
  """Scale the label sizes. This changes the sizes of the lables and shifts them
to make space for the changed font size.
  - scale: scale factor.
  """
  set_label(height=0.0182*scale)

  colorbar_fontsize = array([14.,14.,8.,8.,8.,8.,8.,8.,8.,8.])*scale

  for i in range(1):
    # --- This numbers only work for the systems with plots covering the full
    # --- page.
    ptitle_placement[i][0][1] += (0.9 - 0.895)*scale
    ptitle_placement[i][1][1] -= (0.4 - 0.3927)*scale
    ptitle_placement[i][2][0] -= (0.14 - 0.12)*scale
    ptitle_placement[i][3][1] -= (0.403 - 0.3927)*scale

  # --- Change the default value of the height argument.
  d = list(ptitles.func_defaults)
  d[-1] = 20.*scale
  ptitles.func_defaults = tuple(d)

##########################################################################
def setlinewidth(width=1.):
  """Set the line width for the axis and the data.
Note that this does not affect contour lines.
 - width=1.: desired line width
  """
  if with_matplotlib:
    print "Not yet implemented for matplotlib"
  else:
    gist.pldefault(width=width)
    style = gist.get_style()
    for s in style['systems']:
      s['ticks']['horizontal']['tickStyle']['width'] = width
      s['ticks']['horizontal']['gridStyle']['width'] = width
      s['ticks']['vertical']['tickStyle']['width'] = width
      s['ticks']['vertical']['gridStyle']['width'] = width
      s['ticks']['frameStyle']['width'] = width
      # --- This is needed to get around a bug.
      s['legend']=''
    gist.set_style(style)

##########################################################################
def setviewport(x0=None,x1=None,y0=None,y1=None,view=None,
                movecurrenttitles=False):
    """Set the viewport for a view. This will automatically move the
titles and colorbars. (Though, the change to the color bar location
will only take affect on the next plot.)
 - x0,x1,y0,y1: The viewport in universal coordinates (ranging between
                0 and 1). If any are not given, the current value is used.
 - view=current view: View to change.
 - movecurrenttitles=False: When true, and if view is the current view,
                            then actively move any titles currently
                            being displayed.
    """
    if view is None: view = plsys()
    print plsys()
    try:
        gist_style
    except:
        gstyle()

    # --- Get the any unspecified view port dimensions. These are
    # --- needed for the title placements.
    if x0 is None: x0 = gist_style['systems'][view-1]['viewport'][0]
    if x1 is None: x1 = gist_style['systems'][view-1]['viewport'][1]
    if y0 is None: y0 = gist_style['systems'][view-1]['viewport'][2]
    if y1 is None: y1 = gist_style['systems'][view-1]['viewport'][3]

    # --- Set the new viewport
    gist_style['systems'][view-1]['viewport'][0] = x0
    gist_style['systems'][view-1]['viewport'][1] = x1
    gist_style['systems'][view-1]['viewport'][2] = y0
    gist_style['systems'][view-1]['viewport'][3] = y1
    currentview = plsys()
    set_style(gist_style)
    plsys(currentview)

    # --- Save the current location of the titles.
    t = copy.copy(ptitle_placement[view-1][0])
    b = copy.copy(ptitle_placement[view-1][1])
    l = copy.copy(ptitle_placement[view-1][2])
    r = copy.copy(ptitle_placement[view-1][3])

    # --- Set the new location of the titles.
    ptitle_placement[view-1][0] = [(x0+x1)/2.,y1 + 0.0327]
    ptitle_placement[view-1][1] = [(x0+x1)/2.,y0 - 0.0277]
    ptitle_placement[view-1][2] = [x0 - 0.0557,(y0+y1)/2.]
    ptitle_placement[view-1][3] = [(x0+x1)/2.,y0 - 0.0577]

    if movecurrenttitles and view == plsys():
        # --- If changing the current view, then move any existing
        # --- plot titles.
        currentview = plsys()
        plsys(0)
        try:
            # --- This assumes that the titles are the first four elements
            # --- in system zero. dx and dy specify the change in location.
            pledit(1,dx=ptitle_placement[view-1][0][0] - t[0],
                     dy=ptitle_placement[view-1][0][1] - t[1])
            pledit(2,dx=ptitle_placement[view-1][1][0] - b[0],
                     dy=ptitle_placement[view-1][1][1] - b[1])
            pledit(3,dx=ptitle_placement[view-1][2][0] - l[0],
                     dy=ptitle_placement[view-1][2][1] - l[1])
            pledit(4,dx=ptitle_placement[view-1][3][0] - r[0],
                     dy=ptitle_placement[view-1][3][1] - r[1])
        except:
            pass
        plsys(currentview)

    # --- Update the location of the colorbar. Note that the colorbar
    # --- is not redrawn.
    colorbar_placement[view-1] = [x1+0.0057,x1+0.0257,y0+0.003,y1-0.003]

    # --- Return the new viewport
    return (x0,x1,y0,y1)

##########################################################################
class _getstdout:
    def __init__(self):
        self.out = []
    def write(self,s):
        self.out.append(s)
    def clear(self):
        self.out = []
    def flush(self):
        pass

##########################################################################
def wplq(i):
    """return dictionary of plot options"""
    stdout = _getstdout()
    sys.stdout = stdout
    try:
      plq(i)
    finally:
      # --- Make sure that this reassignment happens no matter what.
      sys.stdout = sys.__stdout__
    result = {}
    for line in stdout.out:
      words = line.split()
      k = 0
      while k < len(words):
        if words[k].find('=') > 0:
          arg = words[k].replace('=','')
          val = words[k+1].replace(',','')
          try:
            val=float(val)
          except ValueError:
            try:
              val=int(val)
            except ValueError:
              val = val.replace('"','')
              val = val.replace("'",'')
              pass
          result[arg] = val
          k += 2
        else:
          k += 1
    return result

##########################################################################
def aplq():
    """return list of dictionaries for all elements in active window"""
    result = []
    i = 1
    while True:
      # --- wplq generally will not raise an error (only if i > maxint)
      d = wplq(i)
      if len(d) == 0: break
      result += [d]
      i += 1
    return result

##########################################################################
def plellipse(l,h,np=100,thetamin=0.,thetamax=2.*pi,xcent=0.,ycent=0.,**kw):
  """Plot ellipse
       - l,               : length
       - h,               : height
       - np=100,          : nb points
       - thetamin = 0.,   : min theta
       - thetamax = 2.*pi : max theta
  """
  dtheta = (thetamax-thetamin)/(np-1)
  theta = arange(thetamin,thetamax+dtheta/2,dtheta)
  x = 0.5*l*cos(theta) + xcent
  y = 0.5*h*sin(theta) + ycent
  pla(y,x,**kw)

