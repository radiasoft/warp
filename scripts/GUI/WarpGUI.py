#!/usr/bin/env python
#Boa:App:BoaApp

from wxPython.wx import *
from warp import *
import WarpRun
import __main__

modules ={'ConsoleClass':     [0, '', 'ConsoleClass.py'],
          'wxDialog_proto':   [0, '', 'wxDialog_proto.py'],
          'EnvelopeGUI':      [0, '', 'EnvelopeGUI.py'],
          'ParticlePlotsGUI': [0, '', 'ParticlePlotsGUI.py'],
          'PzplotsGUI':       [0, '', 'PzplotsGUI.py'],
          'WarpGUIInfo':      [0, '', 'WarpGUIInfo.py'], 
          'pygistDialog':     [0, '', 'pygistDialog.py'],
          'MatchingGUI':      [0, '', 'MatchingGUI.py'],
          'WarpRun':          [1, 'Main frame of Application', 'WarpRun.py']}

class BoaApp(wxApp):
    def OnInit(self):
        self.main = WarpRun.create(None)
        return true

panels=[]
wgui=None
wgui=BoaApp(0)
wgui.initialized=False
wgui.closed=False

def add_panel(panel,name):
    global panels
    panels+=[[panel,name]]
    
def process_gui_events():
    while(wgui.Pending()):
        wgui.Dispatch()
__main__.process_gui_events = process_gui_events
__main__.wgui = wgui

def gui():
  if wgui.initialized:
    print 'The GUI is already running.'
    return
  if wgui.closed:
    print 'The GUI has already been opened and closed once and cannot be reopened in this session. Sorry.'
    return
  wgui.main.init()
  wgui.main.Show()
  if not wgui.initialized:
    wgui.SetTopWindow(wgui.main)
    for i in panels:
      wgui.main.add_panel(i[0],i[1])
    installafterstep(process_gui_events)
    wgui.initialized=True
  wgui.MainLoop()
  

if __name__ == '__main__':
    gui()
