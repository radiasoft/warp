#!/usr/bin/env python
#Boa:App:BoaApp

from wxPython.wx import *
from warp import *

import WarpRun

modules ={'ConsoleClass': [0, '', 'ConsoleClass.py'],
 'EnvelopeGUI': [0, '', 'EnvelopeGUI.py'],
 'ParticlePlotsGUI': [0, '', 'ParticlePlotsGUI.py'],
 'PzplotsGUI': [0, '', '../../boa/PzplotsGUI.py'],
 'WarpGUIInfo': [0, '', 'WarpGUIInfo.py'],
 'WarpRun': [1, 'Main frame of Application', 'WarpRun.py']}

class BoaApp(wxApp):
    def OnInit(self):
        self.main = WarpRun.create(None)
        self.main.Show(true)
        self.SetTopWindow(self.main)
        return true

def main():
    application = BoaApp(0)
    application.MainLoop()

if __name__ == '__main__':
    main()
